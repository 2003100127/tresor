from pathlib import Path
from typing import Dict, Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# ---------- utilities ----------
def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def _beautify(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def _gammaln(x):
    """Vector-safe log-gamma: SciPy if available else vectorized math.lgamma."""
    try:
        from scipy.special import gammaln  # type: ignore
        return gammaln(x)
    except Exception:
        import math
        vlgamma = np.vectorize(math.lgamma, otypes=[float])
        return vlgamma(x)

def _binom_pmf(n: int, p: float, k: np.ndarray) -> np.ndarray:
    from math import comb
    k = np.asarray(k, dtype=int)
    pmf = np.zeros_like(k, dtype=float)
    valid = (k >= 0) & (k <= n)
    kv = k[valid]
    pmf[valid] = [comb(n, int(kk)) * (p ** int(kk)) * ((1.0 - p) ** (n - int(kk))) for kk in kv]
    return pmf

def _nbinom_pmf_mu_theta(k: np.ndarray, mu: float, theta: float) -> np.ndarray:
    """
    NB2: Var = mu + mu^2/theta; r=theta, p=r/(r+mu), pmf = C(k+r-1,k)*p^r*(1-p)^k
    """
    if theta <= 0:
        theta = 1e12
    r = float(theta)
    p = r / (r + mu)
    q = 1.0 - p
    k = np.asarray(k, dtype=float)
    logpmf = (_gammaln(k + r) - _gammaln(r) - _gammaln(k + 1.0)) + r * np.log(p) + k * np.log(q)
    pmf = np.exp(logpmf)
    pmf[~np.isfinite(pmf)] = 0.0
    return pmf

def _chi_square_gof(obs_counts: np.ndarray, exp_probs: np.ndarray, min_exp: float = 5.0) -> Tuple[float, int]:
    """
    Pearson χ² 粗略拟合优度；把期望<min_exp 的相邻桶合并。
    返回 (chi2, dof)。p值如需可自己用 SciPy 计算。
    """
    n = obs_counts.sum()
    exp_counts = exp_probs * n
    k = len(obs_counts)
    merged_obs, merged_exp = [], []
    i = 0
    while i < k:
        o = float(obs_counts[i]); e = float(exp_counts[i]); j = i
        while e < min_exp and j + 1 < k:
            j += 1
            o += float(obs_counts[j]); e += float(exp_counts[j])
        merged_obs.append(o); merged_exp.append(e)
        i = j + 1
    mo = np.array(merged_obs); me = np.array(merged_exp)
    mask = me > 0
    chi2 = float(np.sum((mo[mask] - me[mask]) ** 2 / me[mask]))
    dof = max(int(mask.sum() - 1 - 1), 1)  # -1 sum-to-N，再 -1 近似扣掉1个已估参数
    return chi2, dof


# ---------- main ----------
def analyze_tresor_pcr_counts(
    num_err_per_read_dict: Dict[int, pd.Series],
    # out_png: str,
    readlen: Optional[int] = None,  # 若给读长，则启用二项拟合与φ、χ²(bin)
    p_theory: Optional[float] = None,  # 给理论 p；否则用合并均值/读长估计
    title: str = "Tresor PCR: pooled per-read error counts"
):
    """
    num_err_per_read_dict: {cycle_idx: pandas.Series(num_errors_per_read)}, 各循环独立的 per-read 错误个数。
    本函数会将所有循环的 per-read 错误“合并”为一个总体样本，并与 Binomial/NB 对比；
    另外在一个子图里展示“每循环”的 mean-vs-variance 散点及参考曲线。

    生成一张 2x3 的子图总览（去掉替换矩阵）。
    """
    # ---- 1) 合并样本 & 基本统计 ----
    # 按循环顺序拼接（忽略 index），代表“PCR 全过程中的总体 per-read 错误分布”
    pooled = pd.concat([num_err_per_read_dict[k].astype(int).reset_index(drop=True)
                        for k in sorted(num_err_per_read_dict.keys())], ignore_index=True)
    counts = pooled.values
    n_reads = len(counts)
    max_k = int(counts.max()) if n_reads else 0
    hist_obs = np.bincount(counts, minlength=max_k + 1) if n_reads else np.zeros(1, dtype=int)

    mean_obs = float(np.mean(counts)) if n_reads else 0.0
    var_obs  = float(np.var(counts, ddof=0)) if n_reads else 0.0

    # ---- 2) 每循环 mean/var（用于散点）----
    cyc = []
    for c, s in sorted(num_err_per_read_dict.items()):
        v = s.astype(int).values
        m = float(np.mean(v)) if len(v) else 0.0
        vv = float(np.var(v, ddof=0)) if len(v) else 0.0
        cyc.append((c, m, vv, len(v)))
    df_cyc = pd.DataFrame(cyc, columns=["cycle", "mean", "var", "n"])

    # ---- 3) 二项理论（若有读长） & NB 拟合 ----
    have_binom = readlen is not None and readlen > 0
    if have_binom:
        n = int(readlen)
        p_est = (mean_obs / n) if p_theory is None else float(p_theory)
        mu_bin = n * p_est
        var_bin = n * p_est * (1 - p_est)
        phi = (var_obs / var_bin) if var_bin > 0 else np.nan

        # Binomial pmf（用于叠加/GOF）
        k_emp = np.arange(0, len(hist_obs))
        pmf_bin = _binom_pmf(n, p_est, k_emp)
        chi2_bin, dof_bin = _chi_square_gof(hist_obs.astype(float), pmf_bin)
    else:
        n = None
        p_est = None
        mu_bin = np.nan
        var_bin = np.nan
        phi = np.nan
        k_emp = np.arange(0, len(hist_obs))

    # NB：矩估计 θ̂（以“合并样本”的 mean/var 为准）
    mu_nb = mean_obs if not have_binom else mu_bin
    if var_obs > mu_nb:
        theta_hat = (mu_nb ** 2) / (var_obs - mu_nb)
        theta_hat = float(max(theta_hat, 1e-6))
    else:
        theta_hat = float(1e12)  # 无过散，NB≈二项/泊松

    pmf_nb = _nbinom_pmf_mu_theta(k_emp, mu=mu_nb, theta=theta_hat)
    # NB 在理论上支持 k>max_k；图内只画到 max_k；χ² 里我们用相同长度做合并，
    # 由于我们的观测没有 k>max_k 的桶，这里不追加尾部质量
    chi2_nb, dof_nb = _chi_square_gof(hist_obs.astype(float), pmf_nb)

    # ---- 4) 绘图（2x3，总览）----
    # out_path = Path(out_png).resolve()
    # _ensure_dir(out_path.parent)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    ax1, ax2, ax3, ax4, ax5, ax6 = axes.ravel()

    # (1) 直方图（灰色）
    bins = np.arange(0, max_k + 2) - 0.5
    ax1.hist(counts, bins=bins, color="gray", edgecolor="none")
    ax1.set_xlabel("Errors per read")
    ax1.set_ylabel("Reads")
    ax1.set_title("Pooled histogram")
    _beautify(ax1)

    # (2) 经验 CDF（灰色）
    if n_reads > 0:
        xs = np.sort(counts)
        ys = np.arange(1, n_reads + 1) / n_reads
        ax2.step(xs, ys, where="post", color="gray")
    ax2.set_xlabel("Errors per read")
    ax2.set_ylabel("Empirical CDF")
    ax2.set_ylim(0, 1.0)
    ax2.set_title("Pooled empirical CDF")
    _beautify(ax2)

    # (3) 每循环 mean-vs-var 散点 + 参考曲线
    ax3.scatter(df_cyc["mean"], df_cyc["var"], s=30, color="gray")
    mx = float(max(df_cyc["mean"].max() if len(df_cyc) else 1.0, mean_obs, 1.0))
    xline = np.linspace(0, max(mx, mean_obs) * 1.1, 200)
    # 泊松参考线 var=mean
    ax3.plot(xline, xline, linestyle="--", linewidth=1.2, label="Poisson: var=mean")
    # 二项上界（若有读长）：var = mean * (1 - mean/n)
    if have_binom:
        ax3.plot(xline, xline * (1 - xline / n), linestyle="-.", linewidth=1.2,
                 label=f"Binomial bound (n={n}): var=mean*(1-mean/n)")
    ax3.set_xlabel("Mean errors/read (per cycle)")
    ax3.set_ylabel("Variance (per cycle)")
    ax3.set_title("Per-cycle mean–variance")
    ax3.legend(frameon=False)
    _beautify(ax3)

    # (4) 空面板（你的需求：去掉替换矩阵）
    # ax4.axis("off")
    # ax4.text(0.5, 0.5, "Substitution matrix: N/A", ha="center", va="center")

    # (5) 合并样本：经验 vs Binomial & NB
    rel_emp = hist_obs / max(1, n_reads)
    ax4.bar(k_emp, rel_emp, width=0.9, alpha=0.5, color="gray", label="Empirical")
    if have_binom:
        ax4.plot(k_emp, pmf_bin, marker="o", linewidth=1.5,
                 label=f"Binomial(n={n}, p={p_est:.4g})")
    ax4.plot(k_emp, pmf_nb, marker="s", linestyle="--", linewidth=1.5,
             label=f"NB(μ={mu_nb:.3f}, θ={theta_hat:.3g})")
    ax4.set_xlabel("Errors per read (k)")
    ax4.set_ylabel("Probability / Relative freq.")
    ax4.set_title("Pooled: empirical vs Binomial vs NB")
    ax4.legend(frameon=False)
    _beautify(ax4)

    # (6) 指标面板
    ax5.axis("off")
    lines = [
        title,
        f"reads = {n_reads}",
        f"mean(obs) = {mean_obs:.3f}, var(obs) = {var_obs:.3f}",
        f"NB θ̂ = {theta_hat:.3g}, χ²(NB) = {chi2_nb:.2f} (dof {dof_nb})",
    ]
    if have_binom:
        lines += [
            f"readlen (n) = {n}",
            f"p = {p_est:.6f}" if p_theory is None else f"p (theory) = {p_est:.6f}",
            f"mean(bin) = {mu_bin:.3f}, var(bin) = {var_bin:.3f}",
            f"φ = Var(obs)/Var(bin) = {(var_obs/var_bin):.3f}" if var_bin>0 else "φ = N/A",
            f"χ²(bin) = {chi2_bin:.2f} (dof {dof_bin})",
        ]
    ax5.text(0.02, 0.98, "\n".join(lines), va="top", ha="left", family="monospace")

    fig = ax5.get_figure()
    fig.suptitle("PCR error counts: pooled over cycles", y=1.02, fontsize=14)
    # fig.savefig(out_path, dpi=150)
    # plt.close(fig)
    # print(f"[OK] figure -> {out_path}")

    ax6.axis("off")
    # ax46.text(0.5, 0.5, "Substitution matrix: N/A", ha="center", va="center")

    plt.show()


# --------- example usage ----------
if __name__ == "__main__":
    # 假设你已经有 num_err_per_read_dict（见你的粘贴）
    # 示例：把它命名为 pcr_counts_dict 后直接调用：
    #
    # analyze_tresor_pcr_counts(
    #     num_err_per_read_dict=pcr_counts_dict,
    #     out_png="tresor_pcr_pooled.png",
    #     readlen=100,        # 若已知读长；不清楚就设 None
    #     p_theory=None,      # 有理论 p 就传入，否则自动用 mean/n
    #     title="Tresor PCR (cycles 1..8)"
    # )
    pass
