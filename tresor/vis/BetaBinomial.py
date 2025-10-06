import math
from dataclasses import dataclass
from typing import Optional, Dict, Tuple, Union

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


def _estimate_n_mom_binomial(counts: np.ndarray) -> Tuple[int, float, str]:
    """
    Safe MoM guess for 'effective n' when n_trials is unknown.
    Returns (n_hat, p_hat, note). Never explodes n.
    """
    m = float(np.mean(counts))
    v = float(np.var(counts, ddof=1)) if len(counts) > 1 else float(np.var(counts))
    max_k = int(np.max(counts))

    if m <= 0:
        n0 = max(max_k, 1)
        p0 = 1e-12
        return n0, p0, "Zero mean; set n=max(k)."

    # 过度离散或小样本抖动时 v ≥ m：不要用 m/(1 - v/m)（会爆 n）
    if v >= m:
        # 用一个保守的上界：覆盖观测最大值，并留出少量余量
        margin = int(math.ceil(4.0 * math.sqrt(max(v, 1e-9))))  # ~2σ 的安全余量
        n0 = max(max_k, int(math.ceil(m)) + margin)
        n0 = int(min(n0, max_k + 50))   # 再加一道硬上限，绝不超过 max(k)+50
        p0 = np.clip(m / max(n0, 1), 1e-9, 1 - 1e-9)
        return n0, float(p0), "Over-dispersion: use n≈max(k)+margin."

    # v < m：按二项式矩法近似
    p0 = max(1e-9, min(1 - v / m, 1 - 1e-9))
    n0 = int(max(max_k, round(m / p0)))
    p0 = np.clip(m / max(n0, 1), 1e-9, 1 - 1e-9)
    return n0, float(p0), "MoM estimation (binomial-like)."


def _betabinom_loglik(counts: np.ndarray, n: int, a: float, b: float) -> float:
    """Sum log-likelihood under Beta–Binomial(n, a, b)."""
    if a <= 0 or b <= 0:
        return -np.inf
    return float(np.sum(stats.betabinom.logpmf(counts, n, a, b)))


def _betabinom_mle(counts: np.ndarray, n: int) -> Tuple[float, float, str]:
    """
    MLE for Beta–Binomial with fixed n.
    Optimize on log(alpha), log(beta) for stability; MoM for initialization.
    """
    m = float(np.mean(counts))
    v = float(np.var(counts, ddof=1)) if len(counts) > 1 else float(np.var(counts))
    p = np.clip(m / max(n, 1), 1e-9, 1 - 1e-9)

    # MoM init for tau = alpha+beta from: v = n p(1-p) * (n + tau)/(1 + tau)
    base = n * p * (1 - p)
    note = "MLE with log-params."
    if v <= base + 1e-12:  # no over-dispersion vs binomial
        tau = 1e6
        note += " Var≤n·p·(1-p) → approx Binomial; initialized tau large."
    else:
        A = v / max(base, 1e-12)
        tau = (n - A) / max(A - 1, 1e-9)
        if tau <= 0:
            tau = 10.0
            note += " MoM produced non-positive tau; reset to 10."

    a0 = max(p * tau, 1e-6)
    b0 = max((1 - p) * tau, 1e-6)

    def nll(xy):
        loga, logb = xy
        a, b = float(np.exp(loga)), float(np.exp(logb))
        return -_betabinom_loglik(counts, n, a, b)

    res = minimize(nll, x0=np.array([np.log(a0), np.log(b0)], dtype=float),
                   bounds=[(np.log(1e-9), np.log(1e12)), (np.log(1e-9), np.log(1e12))],
                   method="L-BFGS-B")
    if not res.success:
        note += f" | Optimizer: {res.message}"
    a_hat, b_hat = float(np.exp(res.x[0])), float(np.exp(res.x[1]))
    return a_hat, b_hat, note


def _betabinom_expected_arrays(counts: np.ndarray, n: int, a: float, b: float):
    """Observed/expected arrays over full support k=0..n."""
    N = len(counts)
    ks = np.arange(n + 1)
    pmf = stats.betabinom.pmf(ks, n, a, b)
    exp_full = pmf * N
    obs_full = np.bincount(counts, minlength=n + 1).astype(float)
    return ks, obs_full, exp_full


def _betabinom_gof(counts: np.ndarray, n: int, a: float, b: float, d_params: int, min_exp: float = 5.0):
    ks, obs_full, exp_full = _betabinom_expected_arrays(counts, n, a, b)
    obs_m, exp_m = _merge_small_expected(obs_full, exp_full, min_exp=min_exp)
    chi2 = float(((obs_m - exp_m) ** 2 / np.clip(exp_m, 1e-12, None)).sum())
    df = max(int(len(obs_m) - 1 - d_params), 1)
    pval = float(1 - stats.chi2.cdf(chi2, df))
    return chi2, df, pval, len(obs_m), ks, obs_full, exp_full


def _binom_fit_and_gof(counts: np.ndarray, n: int, min_exp: float = 5.0):
    N = len(counts)
    p_hat = float(np.clip(np.mean(counts) / max(n, 1), 1e-12, 1 - 1e-12))
    ks = np.arange(n + 1)
    pmf = stats.binom.pmf(ks, n, p_hat)
    exp_full = pmf * N
    obs_full = np.bincount(counts, minlength=n + 1).astype(float)
    obs_m, exp_m = _merge_small_expected(obs_full, exp_full, min_exp=min_exp)
    chi2 = float(((obs_m - exp_m) ** 2 / np.clip(exp_m, 1e-12, None)).sum())
    df = max(int(len(obs_m) - 1 - 1), 1)  # fitted p only
    pval = float(1 - stats.chi2.cdf(chi2, df))
    ll = float(np.sum(stats.binom.logpmf(counts, n, p_hat)))
    aic = 2 * 1 - 2 * ll
    bic = 1 * math.log(N) - 2 * ll
    res = GOFResult(
        model="Binomial",
        params={"n": float(n), "p_hat": p_hat},
        loglik=ll, aic=aic, bic=bic, chi2=chi2, df=df, p_value=pval, n_bins_used=len(obs_m),
    )
    return res, ks, obs_full, exp_full


def validate_betabinom_and_plot(
    s: pd.Series,
    n_trials: Optional[int] = None,
    title: str = "Per-read error counts (Beta–Binomial validation)",
    min_exp: float = 5.0,
    show: bool = True,
    make_plots: bool = True,
    fit_binomial: bool = True,
) -> Dict[str, GOFResult]:
    """
    Validate whether counts follow a Beta–Binomial distribution and plot.

    Args:
        s: Series of non-negative integers (errors per read).
        n_trials: number of Bernoulli trials per read (e.g., read length L).
                  If None, an 'effective n' is estimated via method-of-moments.
        min_exp: minimum expected count per merged bin for chi-square.
        make_plots/show: control plotting behavior.
        fit_binomial: also fit a Binomial baseline with the same n.

    Returns:
        Dict with keys: 'betabinom' (and optionally 'binom').
    """
    # optional style
    import seaborn as sns

    sns.set(font="Helvetica")
    sns.set_style("ticks")

    counts = s.to_numpy(dtype=int)
    N = len(counts)

    note_n = ""
    if n_trials is None:
        n_hat0, p0, note_n = _estimate_n_mom_binomial(counts)
        n = int(max(n_hat0, int(np.max(counts))))
        d_params = 3  # (n, alpha, beta) all estimated from data
    else:
        n = int(n_trials)
        d_params = 2  # (alpha, beta) only

    max_k = int(np.max(counts))
    # 缺省把 n 限制在 [max_k, max_k+50] 区间（你可以把 +50 调大/调小）
    if n < max_k:
        n = max_k
    if n > max_k + 50:
        notes_extra = f" | n capped from {n} to {max_k + 50}"
        n = max_k + 50
    else:
        notes_extra = ""

    # MLE for (alpha, beta) with fixed n
    a_hat, b_hat, note_ab = _betabinom_mle(counts, n)
    ll_bb = _betabinom_loglik(counts, n, a_hat, b_hat)
    k_bb = 2 if d_params == 2 else 3
    aic_bb = 2 * k_bb - 2 * ll_bb
    bic_bb = k_bb * math.log(N) - 2 * ll_bb
    chi2, df, pval, nbins, ks, obs_full, exp_full = _betabinom_gof(
        counts, n, a_hat, b_hat, d_params=d_params, min_exp=min_exp
    )

    out: Dict[str, GOFResult] = {}
    out["betabinom"] = GOFResult(
        model="Beta–Binomial",
        params={
            "n": float(n),
            "alpha_hat": float(a_hat),
            "beta_hat": float(b_hat),
            "p_bar_hat": float(a_hat / (a_hat + b_hat)),
            "rho_hat": float(1.0 / (a_hat + b_hat + 1.0)),
        },
        loglik=float(ll_bb),
        aic=float(aic_bb),
        bic=float(bic_bb),
        chi2=float(chi2),
        df=int(df),
        p_value=float(pval),
        n_bins_used=int(nbins),
        notes=(note_n + " | " + note_ab + notes_extra).strip(" |")
    )

    # Binomial baseline with same n
    if fit_binomial:
        bin_res, ks_b, obs_b, exp_b = _binom_fit_and_gof(counts, n, min_exp=min_exp)
        out["binom"] = bin_res

    if make_plots:
        # Histogram + PMF
        plt.figure(figsize=(8, 5))
        rel_obs = obs_full / max(N, 1)
        plt.bar(ks, rel_obs, width=0.9, alpha=0.6, label="Observed (relative freq)")
        pmf_bb = stats.betabinom.pmf(ks, n, a_hat, b_hat)
        plt.plot(ks, pmf_bb, marker="o", linestyle="-",
                 label=f"Beta–Binomial PMF (n={n}, α={a_hat:.3g}, β={b_hat:.3g})")
        if fit_binomial and "binom" in out:
            p_hat = out["binom"].params["p_hat"]
            pmf_bin = stats.binom.pmf(ks, n, p_hat)
            plt.plot(ks, pmf_bin, marker="x", linestyle="--",
                     label=f"Binomial PMF (p={p_hat:.3g})")
        plt.xlabel("Errors per read")
        plt.ylabel("Probability")
        plt.title(f"{title} – Histogram vs PMF")
        plt.legend()
        plt.tight_layout()
        ax = plt.gca()
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")

        # CDF
        plt.figure(figsize=(8, 5))
        emp_cdf = np.cumsum(rel_obs)
        plt.step(ks, emp_cdf, where="post", label="Empirical CDF")
        cdf_bb = stats.betabinom.cdf(ks, n, a_hat, b_hat)
        plt.plot(ks, cdf_bb, linestyle="-", marker="o", label="Beta–Binomial CDF")
        if fit_binomial and "binom" in out:
            p_hat = out["binom"].params["p_hat"]
            cdf_bin = stats.binom.cdf(ks, n, p_hat)
            plt.plot(ks, cdf_bin, linestyle="--", marker="x", label="Binomial CDF")
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


def plot_betabinom_panels(
        data_by_cycle: Dict[int, pd.Series],
        n_trials: Optional[Union[int, Dict[int, int]]] = None,   # int for all cycles or dict{cycle:n}
        cols: int = 4,
        min_exp: float = 5.0,
        read_len: int = 80,
        fs_sub_row: float = 4.2,
        fs_sub_col: float = 3.8,
        suptitle: str = "Beta–Binomial — per PCR cycle",
        show_binomial: bool = True,
):
    """
    Draw per-cycle panels (Histogram+PMF) and (CDF) just like your NB panels.
    """
    keys = sorted(data_by_cycle.keys())
    m = len(keys)
    rows = math.ceil(m / cols)

    # Histogram + PMF
    fig, axes = plt.subplots(rows, cols, figsize=(cols * fs_sub_row, rows * fs_sub_col), squeeze=False)
    for i, k in enumerate(keys):
        r, c = divmod(i, cols); ax = axes[r][c]
        s = data_by_cycle[k].astype(int)
        n_k = n_trials[k] if isinstance(n_trials, dict) else n_trials
        res = validate_betabinom_and_plot(
            s, n_trials=n_k, title="", min_exp=min_exp,
            show=False, make_plots=False, fit_binomial=show_binomial
        )
        bb = res["betabinom"]
        n = int(bb.params["n"])
        a = float(bb.params["alpha_hat"]); b = float(bb.params["beta_hat"])
        pbar = float(bb.params["p_bar_hat"])
        rho = float(bb.params["rho_hat"])
        N = len(s)
        N_nt = N * read_len
        n_err_nt = int(s.sum())
        pct_err_nt = n_err_nt / max(N_nt, 1)
        n_err_read = int((s != 0).sum())

        ks = np.arange(n + 1)
        pmf_bb = stats.betabinom.pmf(ks, n, a, b)
        obs = np.bincount(s, minlength=n + 1).astype(float) / max(N, 1)

        ax.bar(ks, obs, width=0.9, alpha=0.55, label="Observed")
        ax.plot(ks, pmf_bb, marker="o", linestyle="-",
                label=f"BB PMF (n={n}, α={a:.3g}, β={b:.3g})")
        if show_binomial and "binom" in res:
            p_hat = float(res["binom"].params["p_hat"])
            pmf_bin = stats.binom.pmf(ks, n, p_hat)
            ax.plot(ks, pmf_bin, marker="x", linestyle="--",
                    label=f"Binom PMF (p={p_hat:.3g})")

        ax.set_xlabel("Errors per read", fontsize=12)
        ax.set_ylabel("Prob.", fontsize=12)
        ax.set_title(
            f"PCR cycle {k}: N-read={N} | N-err-nt={n_err_nt}\n"
            f"N-err-read={n_err_read} | %-err-nt={pct_err_nt:.2e}\n"
            f"χ² p={bb.p_value:.3g} | p̄̂={pbar:.3g}, ρ̂={rho:.3g}",
            pad=6
        )
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")
        ax.legend(loc="upper right", fontsize=10, frameon=False)

    for j in range(m, rows * cols):
        r, c = divmod(j, cols); axes[r][c].axis("off")

    fig.tight_layout(rect=[0, 0.03, 1, 0.90])
    fig.suptitle(suptitle + " (Histogram vs PMF)", y=0.98, fontsize=14)

    # CDF
    fig2, axes2 = plt.subplots(rows, cols, figsize=(cols * fs_sub_row, rows * fs_sub_col), squeeze=False)
    for i, k in enumerate(keys):
        r, c = divmod(i, cols); ax = axes2[r][c]
        s = data_by_cycle[k].astype(int)
        n_k = n_trials[k] if isinstance(n_trials, dict) else n_trials
        res = validate_betabinom_and_plot(
            s, n_trials=n_k, title="", min_exp=min_exp,
            show=False, make_plots=False, fit_binomial=show_binomial
        )
        bb = res["betabinom"]
        n = int(bb.params["n"])
        a = float(bb.params["alpha_hat"]); b = float(bb.params["beta_hat"])
        N = len(s)

        ks = np.arange(n + 1)
        obs_full = np.bincount(s, minlength=n + 1).astype(float) / max(N, 1)
        emp_cdf = np.cumsum(obs_full)
        cdf_bb = stats.betabinom.cdf(ks, n, a, b)

        ax.step(ks, emp_cdf, where="post", label="Empirical CDF")
        ax.plot(ks, cdf_bb, marker="o", linestyle="-", label="Beta–Binomial CDF")
        if show_binomial and "binom" in res:
            p_hat = float(res["binom"].params["p_hat"])
            cdf_bin = stats.binom.cdf(ks, n, p_hat)
            ax.plot(ks, cdf_bin, marker="x", linestyle="--", label="Binomial CDF")

        ax.set_xlabel("Errors per read", fontsize=12)
        ax.set_ylabel("CDF", fontsize=12)
        ax.set_title(f"PCR cycle {k}", pad=6)
        for side in ("top", "right"): ax.spines[side].set_visible(False)
        ax.tick_params(axis="both", direction="out")
        ax.legend(loc="lower right", fontsize=10, frameon=False)

    for j in range(m, rows * cols):
        r, c = divmod(j, cols); axes2[r][c].axis("off")

    fig2.tight_layout(rect=[0, 0.03, 1, 0.90])
    fig2.suptitle(suptitle + " (CDF comparison)", y=0.98, fontsize=14)


if __name__ == "__main__":
    # Demo with synthetic Beta–Binomial data
    rng = np.random.default_rng(0)
    n_demo = 12
    a_demo, b_demo = 50.0, 1500.0
    demo = stats.betabinom.rvs(n_demo, a_demo, b_demo, size=10000, random_state=rng)
    s = pd.Series(demo, name="num_err_per_read")

    res = validate_betabinom_and_plot(s, n_trials=n_demo, title="Beta–Binomial")
    for k, v in res.items():
        print(f"{v.model}: {v.params}, loglik={v.loglik:.2f}, AIC={v.aic:.2f}, "
              f"BIC={v.bic:.2f}, chi2={v.chi2:.2f}, df={v.df}, p={v.p_value:.3g}, "
              f"bins={v.n_bins_used}. {v.notes}")
