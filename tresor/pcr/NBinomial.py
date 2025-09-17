import math
from dataclasses import dataclass
from typing import Optional, Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize


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


def _nbinom_loglik(counts: np.ndarray, r: float, mu: float) -> float:
    """Negative Binomial fitting / GOF"""
    # Convert (mu, r) -> (n=r, p=r/(r+mu))
    p = r / (r + mu)
    p = float(np.clip(p, 1e-12, 1 - 1e-12))
    return float(np.sum(stats.nbinom.logpmf(counts, r, p)))


def _nbinom_mle(counts: np.ndarray) -> Tuple[float, float, str]:
    """
    MLE for NB using parameterization (mu>0, r>0). Optimizes on (log_mu, log_r).
    Returns (mu_hat, r_hat, note).
    """
    m = counts.mean()
    v = counts.var(ddof=1) if len(counts) > 1 else counts.var()
    note = "MLE with log-params."
    # MoM init
    if v <= m + 1e-9:
        # No over-dispersion vs Poisson; NB collapses to Poisson.
        mu0 = max(m, 1e-9)
        r0 = 1e6
        note += " Var<=Mean → NB≈Poisson; initialized r large."
    else:
        r0 = max((m**2) / (v - m), 1e-6)
        mu0 = max(m, 1e-9)

    def nll(xy):
        log_mu, log_r = xy
        mu = np.exp(log_mu)
        r = np.exp(log_r)
        return -_nbinom_loglik(counts, r, mu)

    res = minimize(
        nll,
        x0=np.array([np.log(mu0), np.log(r0)], dtype=float),
        method="L-BFGS-B",
        bounds=[(np.log(1e-9), np.log(1e12)), (np.log(1e-9), np.log(1e12))],
    )
    if not res.success:
        note += f" | Optimizer: {res.message}"
    mu_hat = float(np.exp(res.x[0]))
    r_hat  = float(np.exp(res.x[1]))
    return mu_hat, r_hat, note


def _nbinom_expected_arrays(counts: np.ndarray, r: float, mu: float, cover_q: float = 0.9995):
    """Build full observed/expected arrays over k=0..K where K covers cover_q of mass."""
    N = len(counts)
    p = r / (r + mu)
    # Choose K to cover most mass & include observed max
    K_q = int(stats.nbinom.ppf(cover_q, r, p))
    K = max(int(np.max(counts)), K_q)
    ks = np.arange(K + 1)
    pmf = stats.nbinom.pmf(ks, r, p)
    exp_full = pmf * N
    obs_full = np.bincount(counts, minlength=K + 1).astype(float)
    return ks, obs_full, exp_full


def _nbinom_gof(counts: np.ndarray, r: float, mu: float, d_params: int = 2, min_exp: float = 5.0):
    ks, obs_full, exp_full = _nbinom_expected_arrays(counts, r, mu)
    obs_m, exp_m = _merge_small_expected(obs_full, exp_full, min_exp=min_exp)
    chi2 = float(((obs_m - exp_m) ** 2 / np.clip(exp_m, 1e-12, None)).sum())
    df = max(int(len(obs_m) - 1 - d_params), 1)
    pval = float(1 - stats.chi2.cdf(chi2, df))
    return chi2, df, pval, len(obs_m), ks, obs_full, exp_full


def _poisson_fit_and_gof(counts: np.ndarray, min_exp: float = 5.0):
    N = len(counts)
    lam = float(np.mean(counts))
    K_q = int(stats.poisson.ppf(0.9995, lam))
    K = max(int(np.max(counts)), K_q)
    ks = np.arange(K + 1)
    pmf = stats.poisson.pmf(ks, lam)
    exp_full = pmf * N
    obs_full = np.bincount(counts, minlength=K + 1).astype(float)
    obs_m, exp_m = _merge_small_expected(obs_full, exp_full, min_exp=min_exp)
    chi2 = float(((obs_m - exp_m) ** 2 / np.clip(exp_m, 1e-12, None)).sum())
    df = max(int(len(obs_m) - 1 - 1), 1)  # fitted λ → 1 parameter
    pval = float(1 - stats.chi2.cdf(chi2, df))
    ll = float(np.sum(stats.poisson.logpmf(counts, lam)))
    k = 1
    aic = 2 * k - 2 * ll
    bic = k * math.log(N) - 2 * ll
    res = GOFResult(
        model="Poisson",
        params={"lambda_hat": lam},
        loglik=ll, aic=aic, bic=bic, chi2=chi2, df=df, p_value=pval, n_bins_used=len(obs_m),
    )
    return res, ks, obs_full, exp_full

def validate_nbinom_and_plot(
        s: pd.Series,
        fit_poisson: bool = True,
        title: str = "Per-read error counts (NB validation)",
        min_exp: float = 5.0,
        show: bool = True,
        make_plots: bool = True,
) -> Dict[str, GOFResult]:
    """
    Validate whether counts follow a Negative Binomial (over-dispersed) distribution and plot.

    Args:
        s: Pandas Series of non-negative integers (error counts per read).
        fit_poisson: Also fit Poisson as a baseline for comparison.
        title: Plot title prefix.
        min_exp: Minimum expected count per merged bin for chi-square test.
        show: Whether to plt.show().

    Returns:
        Dict of GOFResult: {'nbinom': ..., 'poisson': ... (optional)}.
    """
    import seaborn as sns

    sns.set(font="Helvetica")
    sns.set_style("ticks")

    counts = s.to_numpy(dtype=int)
    N = len(counts)
    out: Dict[str, GOFResult] = {}

    mu_hat, r_hat, note = _nbinom_mle(counts)
    ll_nb = _nbinom_loglik(counts, r_hat, mu_hat)
    k_nb = 2
    aic_nb = 2 * k_nb - 2 * ll_nb
    bic_nb = k_nb * math.log(N) - 2 * ll_nb
    chi2, df, pval, nbins, ks, obs_full, exp_full = _nbinom_gof(counts, r_hat, mu_hat, d_params=2, min_exp=min_exp)
    out["nbinom"] = GOFResult(
        model="Negative Binomial",
        params={"mu_hat": float(mu_hat), "r_hat": float(r_hat), "p_hat": float(r_hat / (r_hat + mu_hat))},
        loglik=float(ll_nb),
        aic=float(aic_nb),
        bic=float(bic_nb),
        chi2=float(chi2),
        df=int(df),
        p_value=float(pval),
        n_bins_used=int(nbins),
        notes=note
    )

    if fit_poisson:
        pois_res, ks_p, obs_p, exp_p = _poisson_fit_and_gof(counts, min_exp=min_exp)
        out["poisson"] = pois_res

    if make_plots:
        # Histogram + PMFs
        plt.figure(figsize=(8, 5))
        rel_obs = obs_full / max(N, 1)
        plt.bar(ks, rel_obs, width=0.9, alpha=0.6, label="Observed (relative freq)")
        pmf_nb = stats.nbinom.pmf(ks, r_hat, r_hat / (r_hat + mu_hat))
        plt.plot(ks, pmf_nb, marker="o", linestyle="-", label=f"NB PMF (mu={mu_hat:.3g}, r={r_hat:.3g})")
        if fit_poisson:
            lam = out["poisson"].params["lambda_hat"]
            pmf_p = stats.poisson.pmf(ks, lam)
            plt.plot(ks, pmf_p, marker="x", linestyle="--", label=f"Poisson PMF (λ={lam:.3g})")
        plt.xlabel("Errors per read")
        plt.ylabel("Probability")
        plt.title(f"{title} – Histogram vs PMF")
        plt.legend()
        plt.tight_layout()

        ax = plt.gca()
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")

        # CDF comparison
        plt.figure(figsize=(8, 5))
        emp_cdf = np.cumsum(rel_obs)
        plt.step(ks, emp_cdf, where="post", label="Empirical CDF")
        cdf_nb = stats.nbinom.cdf(ks, r_hat, r_hat / (r_hat + mu_hat))
        plt.plot(ks, cdf_nb, linestyle="-", marker="o", label="NB CDF")
        if fit_poisson:
            lam = out["poisson"].params["lambda_hat"]
            cdf_p = stats.poisson.cdf(ks, lam)
            plt.plot(ks, cdf_p, linestyle="--", marker="x", label="Poisson CDF")
        plt.xlabel("Errors per read")
        plt.ylabel("CDF")
        plt.title(f"{title} – CDF comparison")
        plt.legend()
        plt.tight_layout()
        ax = plt.gca()
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")

        if show:
            plt.show()

    return out


def plot_nbinom_panels(
        data_by_cycle: dict,
        cols: int = 4,
        min_exp: float = 5.0,
        read_len=80,
        fs_sub_row=4.2,
        fs_sub_col=3.8,
        suptitle: str = "Negative Binomial — per PCR cycle",
        show_poisson: bool = False
):
    """
    """
    keys = sorted(data_by_cycle.keys())
    m = len(keys)
    rows = math.ceil(m / cols)

    # Histogram + PMF
    fig, axes = plt.subplots(rows, cols, figsize=(cols * fs_sub_row, rows * fs_sub_col), squeeze=False)
    for i, k in enumerate(keys):
        r, c = divmod(i, cols)
        ax = axes[r][c]
        s = data_by_cycle[k].astype(int)
        res = validate_nbinom_and_plot(
            s,
            fit_poisson=show_poisson,
            title="",
            min_exp=min_exp,
            show=False,
            make_plots=False,
        )
        nb = res["nbinom"]
        mu = float(nb.params["mu_hat"]); rsize = float(nb.params["r_hat"])
        p = rsize / (rsize + mu)
        N = len(s)
        N_nt = len(s) * read_len
        n_err_nt = s.sum()
        pct_err_nt = n_err_nt / N_nt
        n_err_read = (s != 0).sum()

        # Upper limit support: It also includes the maximum observation value and covers 99.95% of theoretical quality
        K = int(max(s.max(), stats.nbinom.ppf(0.9995, rsize, p)))
        ks = np.arange(K + 1)
        pmf_nb = stats.nbinom.pmf(ks, rsize, p)
        obs = np.bincount(s, minlength=K + 1).astype(float) / max(N, 1)

        ax.bar(ks, obs, width=0.9, alpha=0.55, label="Observed")
        ax.plot(ks, pmf_nb, marker="o", linestyle="-", label=f"NB PMF (μ={mu:.3g}, r={rsize:.3g})")
        if show_poisson and "poisson" in res:
            lam = float(res["poisson"].params["lambda_hat"])
            pmf_p = stats.poisson.pmf(ks, lam)
            ax.plot(ks, pmf_p, marker="x", linestyle="--", label=f"Poisson PMF (λ={lam:.3g})")

        ax.set_xlabel("Errors per read", fontsize=12)
        ax.set_ylabel("Prob.", fontsize=12)
        ax.set_title(f"PCR cycle {k}: N-read={N} | N-err-nt={n_err_nt} \n N-err-read={n_err_read} | %-err-nt={pct_err_nt:.2e} \n χ² p={nb.p_value:.3g}")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")
        # if i == 0:
        ax.legend(loc="upper right", fontsize=12, frameon=False)

    for j in range(m, rows * cols):
        r, c = divmod(j, cols)
        axes[r][c].axis("off")

    fig.suptitle(suptitle + " (Histogram vs PMF)", y=0.98, fontsize=14)
    fig.tight_layout()

    # CDF comparison
    fig2, axes2 = plt.subplots(rows, cols, figsize=(cols * fs_sub_row, rows * fs_sub_col), squeeze=False)
    for i, k in enumerate(keys):
        r, c = divmod(i, cols)
        ax = axes2[r][c]
        s = data_by_cycle[k].astype(int)
        res = validate_nbinom_and_plot(
            s,
            fit_poisson=show_poisson,
            title="",
            min_exp=min_exp,
            show=False,
            make_plots=False,
        )
        nb = res["nbinom"]
        mu = float(nb.params["mu_hat"]); rsize = float(nb.params["r_hat"])
        p = rsize / (rsize + mu)
        N = len(s)

        K = int(max(s.max(), stats.nbinom.ppf(0.9995, rsize, p)))
        ks = np.arange(K + 1)
        obs_full = np.bincount(s, minlength=K + 1).astype(float) / max(N, 1)
        emp_cdf = np.cumsum(obs_full)
        cdf_nb = stats.nbinom.cdf(ks, rsize, p)

        ax.step(ks, emp_cdf, where="post", label="Empirical CDF")
        ax.plot(ks, cdf_nb, marker="o", linestyle="-", label="NB CDF")
        if show_poisson and "poisson" in res:
            lam = float(res["poisson"].params["lambda_hat"])
            cdf_p = stats.poisson.cdf(ks, lam)
            ax.plot(ks, cdf_p, marker="x", linestyle="--", label="Poisson CDF")

        ax.set_xlabel("Errors per read", fontsize=12)
        ax.set_ylabel("CDF", fontsize=12)
        ax.set_title(f"PCR cycle {k}")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")
        # if i == 0:
        ax.legend(loc="lower right", fontsize=12, frameon=False)

    for j in range(m, rows * cols):
        r, c = divmod(j, cols)
        axes2[r][c].axis("off")

    fig2.suptitle(suptitle + " (CDF comparison)", y=0.98, fontsize=14)
    fig2.tight_layout()


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    mu_demo, r_demo, N_demo = 0.6, 3.0, 10000
    p_demo = r_demo / (r_demo + mu_demo)
    demo = rng.negative_binomial(n=r_demo, p=p_demo, size=N_demo)
    s = pd.Series(demo, name="num_err_per_read")

    res = validate_nbinom_and_plot(s, fit_poisson=True, title="NB demo")
    # print(res['aic'])
    for k, v in res.items():
        print(f"{v.model}: {v.params}, loglik={v.loglik:.2f}, AIC={v.aic:.2f}, BIC={v.bic:.2f}, "
              f"chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, bins={v.n_bins_used}. {v.notes}")
