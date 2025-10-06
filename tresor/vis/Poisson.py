import math
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


@dataclass
class GOFResult:
    model: str
    params: Dict[str, float]
    loglik: float
    aic: float
    bic: float
    chi2: float
    df: int
    p_value: float
    n_bins_used: int
    notes: str = ""


def _merge_small_expected(obs: np.ndarray, exp: np.ndarray, min_exp: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
    """Merge adjacent bins until each merged bin has expected >= min_exp."""
    merged_obs, merged_exp = [], []
    run_o = run_e = 0.0
    for o, e in zip(obs, exp):
        run_o += o; run_e += e
        if run_e >= min_exp:
            merged_obs.append(run_o); merged_exp.append(run_e)
            run_o = run_e = 0.0
    if run_e > 0:
        if merged_obs:
            merged_obs[-1] += run_o; merged_exp[-1] += run_e
        else:
            merged_obs.append(run_o); merged_exp.append(run_e)
    return np.asarray(merged_obs), np.asarray(merged_exp)


def _poisson_expected_arrays(counts: np.ndarray, lam: float, cover_q: float = 0.9995):
    """Observed/expected arrays over k=0..K where K covers most mass and observed max."""
    N = len(counts)
    if lam <= 1e-12:   # all mass at k=0
        ks = np.array([0], dtype=int)
        pmf = np.array([1.0])
        exp_full = pmf * N
        obs_full = np.bincount(counts, minlength=1).astype(float)
        return ks, obs_full, exp_full
    K_q = int(stats.poisson.ppf(cover_q, lam))
    K = max(int(np.max(counts)), K_q)
    ks = np.arange(K + 1)
    pmf = stats.poisson.pmf(ks, lam)
    exp_full = pmf * N
    obs_full = np.bincount(counts, minlength=K + 1).astype(float)
    return ks, obs_full, exp_full


def _poisson_gof(counts: np.ndarray, lam: float, min_exp: float = 5.0):
    N = len(counts)
    # Special-case lam≈0
    if lam <= 1e-12:
        obs0 = np.array([np.sum(counts == 0)], dtype=float)
        exp0 = np.array([float(N)])
        chi2 = 0.0
        df = 1
        pval = 1.0
        return chi2, df, pval, 1, np.array([0]), obs0, exp0

    ks, obs_full, exp_full = _poisson_expected_arrays(counts, lam)
    obs_m, exp_m = _merge_small_expected(obs_full, exp_full, min_exp=min_exp)
    chi2 = float(((obs_m - exp_m) ** 2 / np.clip(exp_m, 1e-12, None)).sum())
    # fitted one parameter (lambda)
    df = max(int(len(obs_m) - 1 - 1), 1)
    pval = float(1 - stats.chi2.cdf(chi2, df))
    return chi2, df, pval, len(obs_m), ks, obs_full, exp_full


def validate_poisson_and_plot(
    s: pd.Series,
    title: str = "Per-read error counts (Poisson validation)",
    min_exp: float = 5.0,
    show: bool = True,
    make_plots: bool = True,
) -> Dict[str, GOFResult]:
    """
    Validate whether counts follow a Poisson distribution and plot.

    Args:
        s: Pandas Series of non-negative integers (errors per read).
        min_exp: minimum expected count per merged bin for chi-square.
        make_plots/show: control plotting behavior.

    Returns:
        Dict with key 'poisson' (GOFResult).
    """
    import seaborn as sns

    sns.set(font="Helvetica")
    sns.set_style("ticks")

    counts = s.to_numpy(dtype=int)
    N = len(counts)

    # MLE for lambda is the sample mean
    lam = float(np.mean(counts))
    # Log-likelihood, AIC, BIC
    ll = float(np.sum(stats.poisson.logpmf(counts, lam)))
    k = 1
    aic = 2 * k - 2 * ll
    bic = k * math.log(N) - 2 * ll

    chi2, df, pval, nbins, ks, obs_full, exp_full = _poisson_gof(counts, lam, min_exp=min_exp)

    out = {
        "poisson": GOFResult(
            model="Poisson",
            params={"lambda_hat": lam},
            loglik=ll, aic=aic, bic=bic,
            chi2=chi2, df=df, p_value=pval, n_bins_used=nbins,
            notes=""
        )
    }

    if make_plots:
        rel_obs = obs_full / max(N, 1)
        # Histogram + PMF
        plt.figure(figsize=(8, 5))
        plt.bar(ks, rel_obs, width=0.9, alpha=0.6, label="Observed (relative freq)")
        pmf = stats.poisson.pmf(ks, lam)
        plt.plot(ks, pmf, marker="o", linestyle="-", label=f"Poisson PMF (λ={lam:.3g})")
        plt.xlabel("Errors per read")
        plt.ylabel("Probability")
        plt.title(f"{title} – Histogram vs PMF")
        plt.legend()
        plt.tight_layout()
        ax = plt.gca()
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")

        # CDF comparison
        plt.figure(figsize=(8, 5))
        emp_cdf = np.cumsum(rel_obs)
        plt.step(ks, emp_cdf, where="post", label="Empirical CDF")
        cdf_p = stats.poisson.cdf(ks, lam)
        plt.plot(ks, cdf_p, marker="o", linestyle="-", label="Poisson CDF")
        plt.xlabel("Errors per read")
        plt.ylabel("CDF")
        plt.title(f"{title} – CDF comparison")
        plt.legend()
        plt.tight_layout()
        ax = plt.gca()
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")

        if show:
            plt.show()

    return out


def plot_poisson_panels(
    data_by_cycle: Dict[int, pd.Series],
    cols: int = 4,
    min_exp: float = 5.0,
    read_len: int = 80,
    fs_sub_row: float = 4.2,
    fs_sub_col: float = 3.8,
    suptitle: str = "Poisson — per PCR cycle",
):
    """
    Draw per-cycle panels: (Histogram+PMF) and (CDF), mirroring your NB/BB panels.
    """
    keys = sorted(data_by_cycle.keys())
    m = len(keys)
    rows = math.ceil(m / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(cols * fs_sub_row, rows * fs_sub_col), squeeze=False)
    for i, k in enumerate(keys):
        r, c = divmod(i, cols); ax = axes[r][c]
        s = data_by_cycle[k].astype(int)
        res = validate_poisson_and_plot(s, title="", min_exp=min_exp, show=False, make_plots=False)
        pois = res["poisson"]
        lam = float(pois.params["lambda_hat"])
        N = len(s)
        N_nt = N * read_len
        n_err_nt = int(s.sum())
        pct_err_nt = n_err_nt / max(N_nt, 1)
        n_err_read = int((s != 0).sum())

        ks, obs_full, exp_full = _poisson_expected_arrays(s.to_numpy(dtype=int), lam)
        rel_obs = obs_full / max(N, 1)
        pmf = stats.poisson.pmf(ks, lam)

        ax.bar(ks, rel_obs, width=0.9, alpha=0.55, label="Observed")
        ax.plot(ks, pmf, marker="o", linestyle="-", label=f"Poisson PMF (λ={lam:.3g})")

        ax.set_xlabel("Errors per read", fontsize=12)
        ax.set_ylabel("Prob.", fontsize=12)
        ax.set_title(
            f"PCR cycle {k}: N-read={N} | N-err-nt={n_err_nt}\n"
            f"N-err-read={n_err_read} | %-err-nt={pct_err_nt:.2e}\n"
            f"χ² p={pois.p_value:.3g}",
            pad=6
        )
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")
        if i == 0:
            ax.legend(loc="upper right", fontsize=10, frameon=False)

    for j in range(m, rows * cols):
        r, c = divmod(j, cols); axes[r][c].axis("off")

    fig.tight_layout(rect=[0, 0.03, 1, 0.90])
    fig.suptitle(suptitle + " (Histogram vs PMF)", y=0.98, fontsize=14)

    fig2, axes2 = plt.subplots(rows, cols, figsize=(cols * fs_sub_row, rows * fs_sub_col), squeeze=False)
    for i, k in enumerate(keys):
        r, c = divmod(i, cols); ax = axes2[r][c]
        s = data_by_cycle[k].astype(int)
        res = validate_poisson_and_plot(s, title="", min_exp=min_exp, show=False, make_plots=False)
        lam = float(res["poisson"].params["lambda_hat"])
        N = len(s)

        ks, obs_full, _ = _poisson_expected_arrays(s.to_numpy(dtype=int), lam)
        emp_cdf = np.cumsum(obs_full / max(N, 1))
        cdf_p = stats.poisson.cdf(ks, lam)

        ax.step(ks, emp_cdf, where="post", label="Empirical CDF")
        ax.plot(ks, cdf_p, marker="o", linestyle="-", label="Poisson CDF")

        ax.set_xlabel("Errors per read", fontsize=12)
        ax.set_ylabel("CDF", fontsize=12)
        ax.set_title(f"PCR cycle {k}", pad=6)
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")
        if i == 0:
            ax.legend(loc="lower right", fontsize=10, frameon=False)

    for j in range(m, rows * cols):
        r, c = divmod(j, cols); axes2[r][c].axis("off")

    fig2.tight_layout(rect=[0, 0.03, 1, 0.90])
    fig2.suptitle(suptitle + " (CDF comparison)", y=0.98, fontsize=14)


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    lam_demo = 0.8
    demo = rng.poisson(lam=lam_demo, size=10000)
    s = pd.Series(demo, name="num_err_per_read")

    res = validate_poisson_and_plot(s, title="Poisson demo")
    pois = res["poisson"]
    print(f"Poisson: λ̂={pois.params['lambda_hat']:.4f}, "
          f"loglik={pois.loglik:.2f}, AIC={pois.aic:.2f}, BIC={pois.bic:.2f}, "
          f"χ²={pois.chi2:.2f}, df={pois.df}, p={pois.p_value:.3g}, bins={pois.n_bins_used}")
