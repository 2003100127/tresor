import math
from dataclasses import dataclass
from typing import Optional, Dict, Any, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize


@dataclass
class GOFResult:
    model: str
    n: int
    params: Dict[str, float]
    loglik: float
    aic: float
    bic: float
    chi2: float
    df: int
    p_value: float
    n_bins_used: int
    notes: str = ""


def _method_of_moments_binomial(counts: np.ndarray) -> Tuple[int, float, str]:
    """Estimate n and p for Binomial using method-of-moments."""
    m = counts.mean()
    v = counts.var(ddof=1)
    if v > m:
        # Over-dispersion relative to Binomial; MoM would give negative p.
        return int(counts.max()), max(1e-12, m / max(int(counts.max()), 1)), \
               "Var > Mean ⇒ over-dispersion; MoM(n,p) unreliable. Falling back to n=max(k), p=m/n."
    # MoM formulas: p = 1 - v/m; n = m / p
    p_hat = max(1e-12, min(1 - (v / m if m > 0 else 0), 1 - 1e-12)) if m > 0 else 1e-12
    n_hat = max(int(round(m / p_hat)) if p_hat > 0 else int(counts.max()), int(counts.max()))
    # Recompute p with integer n
    p_hat = np.clip(m / max(n_hat, 1), 1e-12, 1 - 1e-12) if n_hat > 0 else 1e-12
    return n_hat, p_hat, "MoM estimation."


def _merge_small_expected(obs: np.ndarray, exp: np.ndarray, min_exp: float = 5.0) -> Tuple[np.ndarray, np.ndarray]:
    """Merge adjacent bins so that each merged bin has expected count >= min_exp."""
    merged_obs, merged_exp = [], []
    run_o, run_e = 0.0, 0.0
    for o, e in zip(obs, exp):
        run_o += o
        run_e += e
        if run_e >= min_exp:
            merged_obs.append(run_o); merged_exp.append(run_e)
            run_o, run_e = 0.0, 0.0
    if run_e > 0:  # tail
        # Merge tail into the previous bin if exists; otherwise keep as is.
        if merged_obs:
            merged_obs[-1] += run_o
            merged_exp[-1] += run_e
        else:
            merged_obs.append(run_o); merged_exp.append(run_e)
    return np.asarray(merged_obs), np.asarray(merged_exp)


def _binomial_gof(counts: np.ndarray, n: int, p: float, d_params: int = 1, min_exp: float = 5.0) -> Tuple[float, int, float, int]:
    N = len(counts)
    obs_full = np.bincount(counts, minlength=n + 1).astype(float)
    pmf = stats.binom.pmf(np.arange(n + 1), n, p)
    exp_full = pmf * N
    # Merge tails for chi-square validity
    obs, exp = _merge_small_expected(obs_full, exp_full, min_exp=min_exp)
    chi2 = ((obs - exp) ** 2 / np.clip(exp, 1e-12, None)).sum()
    # df = bins - 1 - number_of_fitted_params
    df = max(int(len(obs) - 1 - d_params), 1)
    pval = 1 - stats.chi2.cdf(chi2, df)
    return chi2, df, pval, len(obs)


def _loglik_binomial(counts: np.ndarray, n: int, p: float) -> float:
    k = counts
    # Use logpmf; guard p
    p = float(np.clip(p, 1e-12, 1 - 1e-12))
    return np.sum(stats.binom.logpmf(k, n, p))


def _mle_betabinom(counts: np.ndarray, n: int) -> Tuple[float, float, str]:
    """MLE for Beta-Binomial(alpha, beta) with fixed n. Returns (alpha, beta, note)."""
    # Method-of-moments init
    m = counts.mean()
    v = counts.var(ddof=1)
    p = m / max(n, 1)
    p = float(np.clip(p, 1e-9, 1 - 1e-9))
    base_var = n * p * (1 - p)
    note = "MoM init."
    if v <= base_var + 1e-9:
        # No over-dispersion vs Binomial → alpha,beta → very large (approx Binomial). Start with big τ.
        tau = 1e6
    else:
        A = v / max(base_var, 1e-12)
        # A = (n + τ)/(1 + τ) → τ = (n - A)/(A - 1)
        tau = (n - A) / max(A - 1, 1e-9)
        if tau <= 0:
            tau = 10.0
            note = "MoM produced non-positive τ; reset to 10."
    alpha0 = max(p * tau, 1e-6)
    beta0 = max((1 - p) * tau, 1e-6)

    def nll(xy):
        a, b = xy
        if a <= 0 or b <= 0:
            return 1e50
        return -np.sum(stats.betabinom.logpmf(counts, n, a, b))

    res = minimize(nll, x0=np.array([alpha0, beta0]),
                   bounds=[(1e-9, 1e9), (1e-9, 1e9)],
                   method="L-BFGS-B")
    if not res.success:
        note += f" | Optimizer: {res.message}"
    a_hat, b_hat = float(res.x[0]), float(res.x[1])
    return a_hat, b_hat, note


def validate_binomial_and_plot(
        s: pd.Series,
        n_trials: Optional[int] = None,
        p_true: Optional[float] = None,
        fit_betabinomial: bool = True,
        title: str = "Per-read error counts",
        min_exp: float = 5.0,
        show: bool = True,
        make_plots: bool = True,
) -> Dict[str, GOFResult]:
    """
    Validate whether the counts follow a Binomial distribution and plot distributions.

    Args:
        s: Pandas Series of non-negative integers, each = number of errors per read.
        n_trials: If known, the number of Bernoulli trials per read (e.g., read length L).
        p_true: If known, the per-base error rate used in simulation; used for reporting.
        fit_betabinomial: Also fit Beta-Binomial(n, alpha, beta) for comparison (AIC/BIC).
        title: Title prefix for plots.
        min_exp: Minimum expected count per merged bin for chi-square test.
        show: Whether to call plt.show().

    Returns:
        Dict of GOFResult for 'binom' (and 'betabinom' if requested).
    """
    import seaborn as sns

    sns.set(font="Helvetica")
    sns.set_style("ticks")

    counts = s.to_numpy(dtype=int)
    N = len(counts)
    results: Dict[str, GOFResult] = {}

    # Estimate (n, p)
    note = ""
    if n_trials is None:
        n_hat, p_hat, note = _method_of_moments_binomial(counts)
        d_params = 2
    else:
        n_hat = int(n_trials)
        # MLE for p with fixed n
        p_hat = np.clip(counts.mean() / max(n_hat, 1), 1e-12, 1 - 1e-12)
        d_params = 1

    # Binomial GOF
    ll = _loglik_binomial(counts, n_hat, p_hat)
    k_params = d_params
    aic = 2 * k_params - 2 * ll
    bic = k_params * math.log(N) - 2 * ll
    chi2, df, pval, n_bins_used = _binomial_gof(counts, n_hat, p_hat, d_params=d_params, min_exp=min_exp)
    results["binom"] = GOFResult(
        model="Binomial",
        n=n_hat,
        params={"p_hat": float(p_hat), "p_true": float(p_true) if p_true is not None else np.nan},
        loglik=float(ll),
        aic=float(aic),
        bic=float(bic),
        chi2=float(chi2),
        df=int(df),
        p_value=float(pval),
        n_bins_used=int(n_bins_used),
        notes=note.strip()
    )

    # Optionally fit Beta-Binomial for comparison
    if fit_betabinomial:
        a_hat, b_hat, bb_note = _mle_betabinom(counts, n_hat)
        ll_bb = np.sum(stats.betabinom.logpmf(counts, n_hat, a_hat, b_hat))
        k_bb = 2  # alpha, beta (n is fixed)
        aic_bb = 2 * k_bb - 2 * ll_bb
        bic_bb = k_bb * math.log(N) - 2 * ll_bb
        # Chi-square vs expected from betabinom
        pmf_bb = stats.betabinom.pmf(np.arange(n_hat + 1), n_hat, a_hat, b_hat)
        exp_full = pmf_bb * N
        obs_full = np.bincount(counts, minlength=n_hat + 1).astype(float)
        obs_m, exp_m = _merge_small_expected(obs_full, exp_full, min_exp=min_exp)
        chi2_bb = ((obs_m - exp_m) ** 2 / np.clip(exp_m, 1e-12, None)).sum()
        df_bb = max(int(len(obs_m) - 1 - 2), 1)
        pval_bb = 1 - stats.chi2.cdf(chi2_bb, df_bb)
        results["betabinom"] = GOFResult(
            model="Beta-Binomial",
            n=n_hat,
            params={"alpha_hat": float(a_hat), "beta_hat": float(b_hat)},
            loglik=float(ll_bb),
            aic=float(aic_bb),
            bic=float(bic_bb),
            chi2=float(chi2_bb),
            df=int(df_bb),
            p_value=float(pval_bb),
            n_bins_used=int(len(obs_m)),
            notes=bb_note
        )

    if make_plots:
        # Histogram (normalized) + PMF overlay
        ks = np.arange(0, n_hat + 1)
        pmf_bin = stats.binom.pmf(ks, n_hat, p_hat)
        # Observed relative frequencies over full support
        obs_full = np.bincount(counts, minlength=n_hat + 1).astype(float)
        rel_obs = obs_full / max(N, 1)

        plt.figure(figsize=(8, 5))
        plt.bar(ks, rel_obs, width=0.9, alpha=0.6, label="Observed (relative freq)")
        plt.plot(ks, pmf_bin, marker="o", linestyle="-", label=f"Binomial PMF (n={n_hat}, p={p_hat:.4g})")
        if "betabinom" in results:
            pmf_bb = stats.betabinom.pmf(ks, n_hat, results["betabinom"].params["alpha_hat"], results["betabinom"].params["beta_hat"])
            plt.plot(ks, pmf_bb, marker="x", linestyle="--", label="Beta-Binomial PMF (fitted)")
        plt.xlabel("Errors per read", fontsize=12)
        plt.ylabel("Probability", fontsize=12)
        plt.title(f"{title} – Histogram vs PMF")
        plt.legend()
        plt.tight_layout()
        ax = plt.gca()
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")

        # CDF (empirical vs theoretical)
        plt.figure(figsize=(8, 5))
        emp_cdf = np.cumsum(rel_obs)
        plt.step(ks, emp_cdf, where="post", label="Empirical CDF")
        cdf_bin = stats.binom.cdf(ks, n_hat, p_hat)
        plt.plot(ks, cdf_bin, linestyle="-", marker="o", label="Binomial CDF")
        if "betabinom" in results:
            cdf_bb = stats.betabinom.cdf(ks, n_hat, results["betabinom"].params["alpha_hat"], results["betabinom"].params["beta_hat"])
            plt.plot(ks, cdf_bb, linestyle="--", marker="x", label="Beta-Binomial CDF")
        plt.xlabel("Errors per read", fontsize=12)
        plt.ylabel("CDF", fontsize=12)
        plt.title(f"{title} – CDF comparison")
        plt.legend()
        plt.tight_layout()
        ax = plt.gca()
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")

        if show:
            plt.show()

    return results


def plot_binomial_panels(
        data_by_cycle: dict,
        n_trials: int | None = None,
        cols: int = 4,
        read_len=80,
        fs_sub_row=4.2,
        fs_sub_col=3.8,
        min_exp: float = 5.0,
        suptitle: str = "Binomial — per PCR cycle"
):
    """
    """
    keys = sorted(data_by_cycle.keys())
    m = len(keys)
    rows = math.ceil(m / cols)

    # ---------- Page 1: Histogram + PMF ----------
    fig, axes = plt.subplots(rows, cols, figsize=(cols * fs_sub_row, rows * fs_sub_col), squeeze=False)
    for i, k in enumerate(keys):
        r, c = divmod(i, cols)
        ax = axes[r][c]
        s = data_by_cycle[k].astype(int)
        res = validate_binomial_and_plot(
            s,
            n_trials=n_trials,
            fit_betabinomial=False,
            min_exp=min_exp,
            title="",
            show=False,
            make_plots=False,
        )
        binres = res["binom"]
        n_hat = int(binres.n)
        p_hat = float(binres.params["p_hat"])
        N = len(s)

        ks = np.arange(n_hat + 1)
        pmf = stats.binom.pmf(ks, n_hat, p_hat)
        obs = np.bincount(s, minlength=n_hat + 1).astype(float) / max(N, 1)

        ax.bar(ks, obs, width=0.9, alpha=0.55, label="Observed")
        ax.plot(ks, pmf, marker="o", linestyle="-", label=f"Binom PMF (n={n_hat}, p={p_hat:.3g})")

        ax.set_xlabel("Errors per read", fontsize=12)
        ax.set_ylabel("Prob.", fontsize=12)
        ax.set_title(f"PCR cycle {k}  N={N}\nχ² p={binres.p_value:.3g}")
        # 轴样式
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")

        if i == 0:
            ax.legend(loc="upper right", fontsize=10, frameon=False)

    # hide empty subplots
    for j in range(m, rows * cols):
        r, c = divmod(j, cols)
        axes[r][c].axis("off")

    fig.suptitle(suptitle + " (Histogram vs PMF)", y=0.98, fontsize=14)
    fig.tight_layout()

    fig2, axes2 = plt.subplots(rows, cols, figsize=(cols * fs_sub_row, rows * fs_sub_col), squeeze=False)
    for i, k in enumerate(keys):
        r, c = divmod(i, cols)
        ax = axes2[r][c]
        s = data_by_cycle[k].astype(int)
        res = validate_binomial_and_plot(
            s,
            n_trials=n_trials,
            fit_betabinomial=False,
            min_exp=min_exp,
            title="",
            show=False,
            make_plots=False,
        )
        binres = res["binom"]
        n_hat = int(binres.n); p_hat = float(binres.params["p_hat"])
        N = len(s)

        ks = np.arange(n_hat + 1)
        obs_full = np.bincount(s, minlength=n_hat + 1).astype(float) / max(N, 1)
        emp_cdf = np.cumsum(obs_full)
        cdf_th = stats.binom.cdf(ks, n_hat, p_hat)

        ax.step(ks, emp_cdf, where="post", label="Empirical CDF")
        ax.plot(ks, cdf_th, marker="o", linestyle="-", label="Binomial CDF")
        ax.set_xlabel("Errors per read", fontsize=12)
        ax.set_ylabel("CDF", fontsize=12)
        ax.set_title(f"PCR cycle {k}")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")
        if i == 0:
            ax.legend(loc="lower right", fontsize=10, frameon=False)

    for j in range(m, rows * cols):
        r, c = divmod(j, cols)
        axes2[r][c].axis("off")

    fig2.suptitle(suptitle + " (CDF comparison)", y=0.98, fontsize=14)
    fig2.tight_layout()


if __name__ == "__main__":
    demo = pd.Series([0,0,0,2,1] + [0]*12081 + [1,0,0,0], name="num_err_per_read")
    print(demo)
    results = validate_binomial_and_plot(demo, n_trials=None, p_true=None, title="PCR per-read error counts")
    for k, v in results.items():
        print(f"{v.model}: n={v.n}, params={v.params}, loglik={v.loglik:.2f}, "
              f"AIC={v.aic:.2f}, BIC={v.bic:.2f}, chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, bins={v.n_bins_used}. {v.notes}")
