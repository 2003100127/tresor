import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from scipy.stats import nbinom, chi2

TRUTH_CMI_30 = "GGGAAACCCTTTGGGCCCTTTAAACCCTTT"  # known homotrimer CMI (30bp)

# -------------------------------
# 1) Per-read: error positions & counts
# -------------------------------

def compute_cmi_errors(
    df: pd.DataFrame,
    cmi_col: str,
    truth_30: str = TRUTH_CMI_30,
    drop_bad_len: bool = True
) -> pd.DataFrame:
    """
    For each read, compare observed 30-mer CMI to the known truth and compute:
      - error_positions: 1-based indices where obs != truth
      - error_count: number of mismatches (0..30)
    If drop_bad_len=True, rows whose CMI length != 30 are dropped.
    Returns a copy of df with two new columns.
    """
    out = df.copy()
    pos_list = []
    cnt_list = []

    for s in out[cmi_col].astype(str):
        if len(s) != len(truth_30):
            if drop_bad_len:
                pos_list.append(None); cnt_list.append(None)
                continue
            else:
                # pad/trim (not recommended); here we just count as mismatched
                s = s[:30].ljust(30, "N")

        errs = [i+1 for i, (a, b) in enumerate(zip(s, truth_30)) if a != b]  # 1-based
        pos_list.append(errs)
        cnt_list.append(len(errs))

    out["error_positions"] = pos_list
    out["error_count"] = cnt_list
    if drop_bad_len:
        out = out.dropna(subset=["error_count"]).reset_index(drop=True)
        out["error_count"] = out["error_count"].astype(int)
    return out

# -------------------------------
# 2) NB fitting (truncated to 0..K=30), GOF, plotting
# -------------------------------

def _pool_bins(observed: np.ndarray, expected: np.ndarray, min_exp: float = 5.0):
    """Pool adjacent bins so every expected >= min_exp (needed for Pearson's chi-square)."""
    obs = observed.astype(float).copy()
    exp = expected.astype(float).copy()

    # merge from left
    while len(exp) > 1 and exp[0] < min_exp:
        exp[1] += exp[0]; obs[1] += obs[0]
        exp = exp[1:];    obs = obs[1:]
    # merge from right
    while len(exp) > 1 and exp[-1] < min_exp:
        exp[-2] += exp[-1]; obs[-2] += obs[-1]
        exp = exp[:-1];     obs = obs[:-1]
    # inner small bins
    i = 1
    while i < len(exp)-1:
        if exp[i] < min_exp:
            # merge into the larger neighbor
            if exp[i-1] >= exp[i+1]:
                exp[i-1] += exp[i]; obs[i-1] += obs[i]
                exp = np.delete(exp, i); obs = np.delete(obs, i)
            else:
                exp[i+1] += exp[i]; obs[i+1] += obs[i]
                exp = np.delete(exp, i); obs = np.delete(obs, i)
        else:
            i += 1
    return obs, exp

def _nb_from_moments(counts: np.ndarray) -> Tuple[float, float]:
    """
    Method-of-moments estimate for NB with support {0,1,2,...}:
      r = mu^2 / (var - mu),  p = r / (r + mu)
    """
    mu = counts.mean()
    var = counts.var(ddof=0)
    if var <= mu + 1e-12:
        # nearly Poisson/Binomial; push a tiny overdispersion to avoid degenerate r
        r = (mu**2) / max(var - mu, 1e-6)
    else:
        r = (mu**2) / (var - mu)
    r = max(r, 1e-6)
    p = r / (r + mu) if (r + mu) > 0 else 0.999999
    p = min(max(p, 1e-9), 1 - 1e-9)
    return float(r), float(p)

def _truncated_nb_expected_hist(
    r: float, p: float, N: int, K: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Expected histogram under NB(r,p) truncated to {0,...,K}.
    Returns (pmf_trunc, expected_counts).
    """
    k = np.arange(K+1)
    pmf_full = nbinom.pmf(k, r, p)  # SciPy: pmf(k) = C(k+r-1,k) * (1-p)^k * p^r
    Z = pmf_full.sum()  # sum over 0..K
    pmf_trunc = pmf_full / max(Z, 1e-300)
    expected = pmf_trunc * N
    return pmf_trunc, expected

def _chisq_gof_truncated_nb(
    obs_hist: np.ndarray,
    r: float, p: float, K: int,
    n_params: int = 2,
    min_exp: float = 5.0
) -> Dict:
    N = int(obs_hist.sum())
    pmf_trunc, expected = _truncated_nb_expected_hist(r, p, N, K)
    obs_p, exp_p = _pool_bins(obs_hist, expected, min_exp=min_exp)
    mask = exp_p > 0
    obs_p, exp_p = obs_p[mask], exp_p[mask]
    chi2_stat = float(np.sum((obs_p - exp_p)**2 / exp_p))
    df = max(len(exp_p) - 1 - n_params, 1)
    pval = float(1 - chi2.cdf(chi2_stat, df))
    return {"chi2": chi2_stat, "df": df, "p": pval, "pmf_trunc": pmf_trunc, "expected": expected}

def validate_nbinom_and_plot(
    df: pd.DataFrame,
    cmi_col: str,
    truth_30: str = TRUTH_CMI_30,
    use_moments: bool = True,
    bootstrap_B: int = 0,
    ax: Optional[plt.Axes] = None
) -> Dict:
    """
    End-to-end:
      1) compute per-read error_count vs truth
      2) fit NB(r,p) by moments (default) or MLE via scipy.fit (optional extension)
      3) Pearson chi-square GOF with truncated NB to 0..30
      4) (optional) parametric bootstrap p-value
      5) plot Observed vs Expected histogram

    Returns a dict with r, p, mu, var, chi2/df/p, (optional) p_boot, and the per-read table.
    """
    # per-read errors
    df2 = compute_cmi_errors(df, cmi_col=cmi_col, truth_30=truth_30)
    counts = df2["error_count"].values.astype(int)
    N = len(counts)
    K = 30
    obs_hist = np.bincount(counts, minlength=K+1)

    # fit NB
    if use_moments:
        r, p = _nb_from_moments(counts)
    else:
        # SciPy fit with loc fixed at 0; then convert to r,p (scipy returns n,p,loc)
        n_hat, p_hat, loc_hat = nbinom.fit(counts, floc=0)
        r, p = float(n_hat), float(p_hat)

    mu = r * (1 - p) / p
    var = mu + (mu**2) / r
    phi = 1.0 / r  # dispersion (common in GLM)

    # GOF (truncated to 0..30)
    gof = _chisq_gof_truncated_nb(obs_hist, r, p, K=K, n_params=2, min_exp=5.0)

    # optional parametric bootstrap under truncated NB
    p_boot = None
    if bootstrap_B and bootstrap_B > 0:
        rng = np.random.default_rng(1)
        stat_obs = gof["chi2"]
        ge = 0
        # simulate from full NB then truncate by resampling until within 0..K
        # (efficient when tail mass > 0 is tiny; otherwise, sample via rejection with cap)
        pmf = gof["pmf_trunc"]
        cdf = np.cumsum(pmf)
        for _ in range(bootstrap_B):
            u = rng.random(N)
            sim = np.searchsorted(cdf, u).astype(int)  # inverse-CDF sampling on 0..K
            sim_hist = np.bincount(sim, minlength=K+1)
            stat_sim = _chisq_gof_truncated_nb(sim_hist, r, p, K=K, n_params=2)["chi2"]
            ge += (stat_sim >= stat_obs)
        p_boot = ge / bootstrap_B

    # plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    xs = np.arange(K+1)
    ax.bar(xs - 0.2, obs_hist, width=0.4, label="Observed")
    ax.bar(xs + 0.2, gof["expected"], width=0.4, label="Expected")
    ax.set_xlabel("Errors per read (0..30)")
    ax.set_ylabel("Frequency")
    ax.set_title(f"NB fit (truncated 0..30): r={r:.3g}, p={p:.3g} | "
                 f"chi2={gof['chi2']:.1f}, df={gof['df']}, p={gof['p']:.3g}")
    ax.legend()
    plt.tight_layout()

    return {
        "per_read": df2[["error_positions", "error_count"]],
        "hist_observed": obs_hist,
        "nb_r": r, "nb_p": p,
        "mu": mu, "var": var, "phi": phi,  # phi=1/r
        "chi2": gof["chi2"], "df": gof["df"], "p_value": gof["p"],
        "p_bootstrap": p_boot,
    }

# -------------------------------
# Minimal example
# -------------------------------
if __name__ == "__main__":
    # Example: fake data (replace with your real dataframe `df` and column name)
    # df = pd.DataFrame({
    #     "CMI": [
    #         TRUTH_CMI_30,  # 0 error
    #         TRUTH_CMI_30[:10] + "G" + TRUTH_CMI_30[11:],  # 1 error
    #         TRUTH_CMI_30[:5] + "T" + TRUTH_CMI_30[6:15] + "A" + TRUTH_CMI_30[16:]  # 2 errors
    #     ] * 1000
    # })
    bam_fpn = "/mnt/d/Document/Programming/Python/umiche/umiche/data/r1/cmi/Aligned_GACTGCTACT_gene_sorted_primary_tagged.bam"
    from umiche.bam.Reader import ReaderChunk

    df_bam = ReaderChunk(
        bam_fpn=bam_fpn,
        bam_fields=['contig', 'pos', 'CIGAR', 'seq', 'read'],
        tag_whitelist=['CB', 'MB', 'XT', 'XS'],
        verbose=True,
    ).todf(chunk_size=1_000_000)
    from tresor.util.Console import Console
    Console(True).df_column_summary(df=df_bam)
    # df_bam = df_bam.reset_index(drop=True)
    df_bam = df_bam.dropna(subset=["XT"]).reset_index(drop=True)
    Console(True).df_column_summary(df=df_bam)
    df_bam = df_bam[df_bam["MB"] != ""]
    Console(True).df_column_summary(df=df_bam)
    df_bam = df_bam[df_bam["MB"].str.len() == 30]
    Console(True).df_column_summary(df=df_bam)

    print(df_bam)
    # res = validate_nbinom_and_plot(df, cmi_col="CMI", bootstrap_B=0)  # set B=1000 for bootstrap p
    # print(res["nb_r"], res["nb_p"], res["chi2"], res["df"], res["p_value"])
    # plt.show()
