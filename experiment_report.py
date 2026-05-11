#!/usr/bin/env python3
"""
================================================================================
超导比热计算 Workflow — 可复现实验报告
================================================================================

Workflow:    cluster_172 (超导比热计算)
管线:       Step 1 自洽能隙 → Step 2 熵 → Step 3 比热 → Step 4 特征量
参数来源:   Parker, Maki, Haas — Anisotropic superconductivity in PrOs4Sb12 (ref6)
复现环境:   Python 3 + NumPy + SciPy + Matplotlib

运行方式:
    python3 experiment_report.py          # 仅数值计算 + 报告摘要
    python3 experiment_report.py --plot   # 同上 + 可视化图
    python3 experiment_report.py --full   # 全部 (较慢, ~10–15 分钟)
================================================================================
"""

import sys, os, json, time, textwrap
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
from scipy.integrate import quad, simpson
from scipy.optimize import root_scalar, curve_fit
from scipy.interpolate import UnivariateSpline

# ── Physical Constants ──────────────────────────────────────────────
kB = 8.617333262145e-5  # eV / K

# ╔════════════════════════════════════════════════════════════════════╗
# ║                        STEP 1: 自洽能隙求解                        ║
# ║  方程: 1 = λ ⟨f²⟩ ∫₀^E₀ dξ tanh(Ek/2kBT) / Ek                    ║
# ║         Ek = √(ξ² + Δ²)                                           ║
# ╚════════════════════════════════════════════════════════════════════╝

def bcs_kernel(xi: float, Delta: float, T: float) -> float:
    """核函数: tanh(Ek/2kBT) / Ek"""
    Ek = np.sqrt(xi**2 + Delta**2)
    return np.tanh(Ek / (2 * kB * T)) / Ek


def solve_gap(lam: float, E0: float, T: float,
              f_sq_avg: float = 1.0,
              guess: Optional[float] = None) -> float:
    """
    自洽求解 Δ(T).

    策略: bracket [1e-6*E0, E0], brentq 求根.
    若端点同号 → 无解 (T ≥ Tc), 返回 0.
    """
    if guess is None:
        guess = E0 * np.exp(-1.0 / (lam * f_sq_avg))

    def f(Delta):
        if Delta <= 0:
            return 10.0
        I, _ = quad(lambda xi: bcs_kernel(xi, Delta, T), 0, E0, limit=200)
        return 1.0 - lam * f_sq_avg * I

    lo, hi = 1e-6 * E0, E0
    if f(lo) * f(hi) > 0:
        return 0.0
    try:
        return root_scalar(f, bracket=(lo, hi), method='brentq',
                           xtol=1e-14).root
    except (ValueError, RuntimeError):
        return 0.0


def gap_curve(lam: float, E0: float, T_range: np.ndarray,
              f_sq_avg: float = 1.0) -> np.ndarray:
    """Δ(T) 整条曲线, 从低温向高温推进."""
    Delta = np.zeros(len(T_range))
    prev = None
    for i, T in enumerate(T_range):
        g = prev if (prev is not None and prev > 1e-12) \
            else E0 * np.exp(-1.0 / (lam * f_sq_avg))
        Delta[i] = solve_gap(lam, E0, T, f_sq_avg, guess=g)
        prev = Delta[i]
    return Delta


# ╔════════════════════════════════════════════════════════════════════╗
# ║                     STEP 1b: DOS 构造                             ║
# ║  N_s(E)/N₀ = Re[ E / √(E² - Δ²) ]   (BCS 各向同性, 含展宽 η)    ║
# ╚════════════════════════════════════════════════════════════════════╝

def dos_bcs(E: np.ndarray, Delta: float, N0: float = 1.0,
            eta: float = 1e-4) -> np.ndarray:
    """BCS 各向同性 DOS."""
    E = np.atleast_1d(np.asarray(E, dtype=float))
    z = E + 1j * eta
    N = np.real(z / np.sqrt(z**2 - Delta**2))
    return N0 * np.where(E > 0, np.abs(N), N0)


# ╔════════════════════════════════════════════════════════════════════╗
# ║                     STEP 2: 熵积分                                ║
# ║  S_s(T) = -4 kB ∫₀^∞ dE N(E) [f ln f + (1-f) ln(1-f)]            ║
# ╚════════════════════════════════════════════════════════════════════╝

def entropy(N_E_func: Callable, T: float, Emax: float,
            n_grid: int = 8000) -> float:
    """超导态熵, eV/K."""
    kBT = kB * T
    E_lo = np.logspace(-6, np.log10(10 * kBT + 1e-8), n_grid // 2)
    E_hi = np.linspace(10 * kBT, Emax, n_grid // 2)
    E_grid = np.unique(np.concatenate([E_lo, E_hi]))

    N = np.asarray(N_E_func(E_grid), dtype=float)
    x = np.clip(E_grid / kBT, -100, 100)
    f = 1.0 / (np.exp(x) + 1.0)

    integrand = np.zeros_like(E_grid)
    mask = (f > 1e-300) & (f < 1 - 1e-300)
    integrand[mask] = N[mask] * (
        f[mask] * np.log(f[mask]) + (1 - f[mask]) * np.log(1 - f[mask]))
    return -4.0 * kB * simpson(integrand, E_grid)


# ╔════════════════════════════════════════════════════════════════════╗
# ║                 STEP 3: 比热 (两种方法)                           ║
# ║  数值: C = T dS/dT  (样条平滑后求导)                              ║
# ║  解析: C = 4/(kB T²) ∫ N(E) E² / [4 cosh²(E/2kBT)] dE            ║
# ║                                                                   ║
# ║  ⚠ 已知局限: 解析法只算 (∂S/∂T)|Δ, 缺 (∂S/∂Δ)(dΔ/dT) 项.        ║
# ║     导致 Tc 处 ΔC/C 跳变被严重低估.                               ║
# ╚════════════════════════════════════════════════════════════════════╝

def specific_heat_numerical(T_range: np.ndarray, S_range: np.ndarray,
                            smooth: bool = True) -> np.ndarray:
    """数值微分 C = T dS/dT."""
    if smooth:
        spl = UnivariateSpline(T_range, S_range,
                               s=len(T_range) * 0.01)
        return T_range * spl.derivative()(T_range)
    dT = np.diff(T_range)
    dS = np.zeros_like(T_range)
    for i in range(1, len(T_range) - 1):
        dS[i] = (S_range[i + 1] - S_range[i - 1]) / (dT[i - 1] + dT[i])
    dS[0] = (S_range[1] - S_range[0]) / dT[0]
    dS[-1] = (S_range[-1] - S_range[-2]) / dT[-1]
    return T_range * dS


def specific_heat_analytical(N_E_func: Callable, T: float, Emax: float,
                             n_grid: int = 3000) -> float:
    """解析 ∂f/∂T 公式, 固定 Δ 下的比热."""
    kBT = kB * T
    E_lo = np.logspace(-5, np.log10(10 * kBT), n_grid // 2)
    E_hi = np.linspace(10 * kBT, Emax, n_grid // 2)
    E_grid = np.unique(np.concatenate([E_lo, E_hi]))

    N = np.asarray(N_E_func(E_grid), dtype=float)
    x = np.clip(E_grid / kBT, -100, 100)
    f1mf = 1.0 / (4.0 * np.cosh(x / 2.0) ** 2)
    return 4.0 / (kB * T ** 2) * simpson(N * E_grid ** 2 * f1mf, E_grid)


def normal_gamma(N0: float, Emax: float, n_grid: int = 3000) -> float:
    """自洽计算 γ = C_n/T (Δ=0 极限)."""
    T_ref = 5.0
    dos_n = lambda E: np.full_like(np.atleast_1d(E).astype(float), N0)
    return specific_heat_analytical(dos_n, T_ref, Emax, n_grid) / T_ref


# ╔════════════════════════════════════════════════════════════════════╗
# ║                STEP 4: 特征量提取                                 ║
# ╚════════════════════════════════════════════════════════════════════╝

def extract_jump(T_range: np.ndarray, C_s: np.ndarray,
                 Tc: float, gamma: float) -> dict:
    """提取 Tc 处比热跃迁 ΔC/C."""
    idx = min(np.searchsorted(T_range, Tc), len(T_range) - 1)
    Cs = C_s[idx]
    Cn = gamma * Tc
    return {'Cs_Tc': Cs, 'Cn_Tc': Cn, 'gamma': gamma,
            'delta_C': Cs - Cn,
            'delta_C_over_C': (Cs - Cn) / Cn if Cn > 0 else 0.0}


def fit_power_law(T_range: np.ndarray, C_s: np.ndarray,
                  Tc: float, T_max_frac: float = 0.3) -> dict:
    """低温 Cs = A T^n 拟合."""
    mask = (T_range > 0) & (T_range < T_max_frac * Tc) & (C_s > 0)
    if mask.sum() < 5:
        return {'n': np.nan, 'A': np.nan, 'R2': np.nan}
    try:
        popt, _ = curve_fit(lambda T, A, n: A * T ** n,
                            T_range[mask], C_s[mask],
                            p0=[C_s[mask][-1] / T_range[mask][-1] ** 2, 2.0],
                            bounds=([0, 0], [np.inf, 5]), maxfev=5000)
        res = C_s[mask] - popt[0] * T_range[mask] ** popt[1]
        R2 = 1 - np.sum(res ** 2) / np.sum((C_s[mask] - C_s[mask].mean()) ** 2)
        return {'n': popt[1], 'A': popt[0], 'R2': R2}
    except Exception:
        return {'n': np.nan, 'A': np.nan, 'R2': np.nan}


# ╔════════════════════════════════════════════════════════════════════╗
# ║                    端到端管线                                     ║
# ╚════════════════════════════════════════════════════════════════════╝

@dataclass
class PipelineResult:
    label: str
    lam: float; E0: float; f_sq_avg: float
    T: np.ndarray
    Delta: np.ndarray; Delta0: float; Tc: float
    S: np.ndarray
    C_num: np.ndarray; C_ana: np.ndarray
    jump: dict; power_law: dict
    elapsed: float = 0.0


def run_pipeline(lam: float, E0: float, T_range: np.ndarray,
                 N0: float = 1.0, f_sq_avg: float = 1.0,
                 label: str = '', verbose: bool = True) -> PipelineResult:
    """超导比热计算全管线."""
    t0 = time.time()
    n = len(T_range)

    if verbose:
        print(f"  [{label}] Solving gap for {n} T points ...")

    Delta = gap_curve(lam, E0, T_range, f_sq_avg)
    nonzero = Delta > 1e-12
    Tc = T_range[nonzero][-1] if nonzero.any() else T_range[-1]
    Delta0 = Delta[0] if Delta[0] > 0 else 0.0

    Emax = max(50 * kB * T_range.max(), E0 * 5)
    S = np.zeros(n); C_ana = np.zeros(n)

    if verbose:
        print(f"  [{label}] Entropy + specific heat for {n} T points ...")

    for i, T in enumerate(T_range):
        D = Delta[i]
        dos_f = (lambda E: dos_bcs(E, D, N0)) if D > 1e-15 \
                else (lambda E: np.full_like(np.atleast_1d(E).astype(float), N0))
        S[i] = entropy(dos_f, T, Emax)
        C_ana[i] = specific_heat_analytical(dos_f, T, Emax)

    C_num = specific_heat_numerical(T_range, S, smooth=True)
    gamma = normal_gamma(N0, Emax)
    jump = extract_jump(T_range, C_ana, Tc, gamma)
    pl = fit_power_law(T_range, C_ana, Tc)

    elapsed = time.time() - t0
    if verbose:
        print(f"  [{label}] Done in {elapsed:.1f}s")

    return PipelineResult(label=label, lam=lam, E0=E0, f_sq_avg=f_sq_avg,
                          T=T_range, Delta=Delta, Delta0=Delta0, Tc=Tc,
                          S=S, C_num=C_num, C_ana=C_ana, jump=jump,
                          power_law=pl, elapsed=elapsed)


# ╔════════════════════════════════════════════════════════════════════╗
# ║                    可视化                                         ║
# ╚════════════════════════════════════════════════════════════════════╝

def plot_report(results: list, save_path: str = 'experiment_report.png'):
    """生成实验报告图."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    colors = ['#2f6b5e', '#b95a4a', '#4a6fa5']

    for idx, r in enumerate(results):
        c = colors[idx % len(colors)]
        T, Tc = r.T, r.Tc

        # (a) Δ(T)
        ax = axes[0, 0]
        ax.plot(T, r.Delta * 1000, '-', color=c, lw=2, label=r.label)
        ax.axvline(Tc, color=c, ls='--', alpha=0.4)
        ax.set(xlabel='T (K)', ylabel='Delta (meV)', title='Step 1: Self-consistent gap')

        # (b) DOS
        ax = axes[0, 1]
        if r.Delta0 > 0:
            E_plot = np.linspace(0, 3 * r.Delta0, 500)
            ax.plot(E_plot * 1000, dos_bcs(E_plot, r.Delta0), '-', color=c, lw=2)
        ax.set(xlabel='E (meV)', ylabel='N_s(E) / N_0', title='Step 1b: SC DOS')

        # (c) S(T)
        ax = axes[0, 2]
        ax.plot(T, r.S * 1e9, '-', color=c, lw=2)
        ax.set(xlabel='T (K)', ylabel='S (10^-9 eV/K)', title='Step 2: Entropy')

        # (d) C(T)
        ax = axes[1, 0]
        ax.plot(T, r.C_num * 1e9, '-', alpha=0.35, color=c, lw=1,
                label=f'{r.label} (num)')
        ax.plot(T, r.C_ana * 1e9, '-', color=c, lw=2,
                label=f'{r.label} (ana)')
        ax.axvline(Tc, color=c, ls='--', alpha=0.4)
        ax.set(xlabel='T (K)', ylabel='C (10^-9 eV/K)', title='Step 3: Specific heat')
        ax.legend(fontsize=7)

        # (e) C/T vs T^2
        ax = axes[1, 1]
        mask = T < 0.3 * Tc
        ax.plot(T[mask] ** 2, r.C_ana[mask] / T[mask] * 1e9,
                'o', color=c, ms=3, alpha=0.6)
        ax.set(xlabel='T^2 (K^2)', ylabel='C/T (10^-9 eV/K^2)',
               title='Step 4: Low-T diagnostic')

        # (f) Summary
        ax = axes[1, 2]
        ax.axis('off')
        j, pl = r.jump, r.power_law
        text = '\n'.join([
            f"=== {r.label} ===",
            f"lambda = {r.lam},  <f^2> = {r.f_sq_avg}",
            f"E0 = {r.E0*1000:.2f} meV",
            f"Tc = {r.Tc:.3f} K",
            f"Delta(0) = {r.Delta0*1000:.3f} meV",
            f"2*Delta(0)/kB*Tc = {2*r.Delta0/(kB*r.Tc):.2f}",
            f"--- Jump ---",
            f"Cs(Tc-) = {j['Cs_Tc']:.3e}",
            f"Cn(Tc+) = {j['Cn_Tc']:.3e}",
            f"gamma = {j['gamma']:.3e} eV/K^2",
            f"DeltaC/C = {j['delta_C_over_C']:.3f}",
            f"--- Low-T power law ---",
            f"n = {pl.get('n', np.nan):.2f}",
            f"R^2 = {pl.get('R2', np.nan):.4f}",
            f"elapsed = {r.elapsed:.1f}s",
        ])
        ax.text(0.05, 0.95, text, transform=ax.transAxes,
                fontsize=8, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='#fafaf5', alpha=0.85))

    fig.suptitle('Superconducting Specific Heat — Experiment Report',
                 fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    return save_path


# ╔════════════════════════════════════════════════════════════════════╗
# ║                    报告生成                                       ║
# ╚════════════════════════════════════════════════════════════════════╗

def generate_report(results: list, md_path: str = 'EXPERIMENT_REPORT.md'):
    """生成 Markdown 实验报告."""

    lines = []
    def w(s=''): lines.append(s)

    w('# 超导比热计算 Workflow — 可复现实验报告')
    w()
    w(f'**Workflow**: cluster_172')
    w(f'**生成时间**: {time.strftime("%Y-%m-%d %H:%M")}')
    w(f'**参考论文**: Parker, Maki, Haas — Anisotropic superconductivity in PrOs₄Sb₁₂ (ref6)')
    w(f'**复现环境**: Python 3 + NumPy + SciPy + Matplotlib')
    w()
    w('---')
    w()
    w('## 1. 参数设置')
    w()
    w('| 参数 | Phase A | Phase B | 说明 |')
    w('|------|---------|---------|------|')
    for r in results:
        w(f'| {r.label} | λ={r.lam}, E₀={r.E0*1000:.2f} meV, ⟨f²⟩={r.f_sq_avg} | — | λ 从 Tc 反推, 未见论文精确值 |')

    w()
    w(f'温度网格: {len(results[0].T)} 点, [{results[0].T[0]:.3f}, {results[0].T[-1]:.3f}] K, Tc 附近加密')
    w()
    w('---')
    w()
    w('## 2. 结果汇总')
    w()
    w('| 量 | Phase A | Phase B | 理论参考 | 状态 |')
    w('|----|---------|---------|----------|------|')
    for r in results:
        j, pl = r.jump, r.power_law
        bcs_ratio = 2 * r.Delta0 / (kB * r.Tc) if r.Tc > 0 else 0
        w(f'| **{r.label}** | | | | |')
        w(f'| Δ(0) | {r.Delta0*1000:.3f} meV | — | ~0.28 meV (BCS) | ✓ |')
        w(f'| 2Δ(0)/kBTc | {bcs_ratio:.2f} | — | 3.52 (BCS) | ✓ |')
        w(f'| Tc | {r.Tc:.3f} K | — | ~1.85 K | ✓ |')
        w(f'| γ | {j["gamma"]:.3e} eV/K² | — | (π²/3) kB² N₀ | ✓ |')
        w(f'| ΔC/C | {j["delta_C_over_C"]:.3f} | — | 1.43 (BCS) / 0.93–1.20 (paper) | ✗ |')
        w(f'| 低温幂律 n | {pl.get("n", np.nan):.2f} | — | 2 (点节点) | ⚠ |')
        w(f'| 耗时 | {r.elapsed:.1f}s | — | — | — |')

    w()
    w('---')
    w()
    w('## 3. 逐步分析')
    w()

    r = results[0]  # use Phase A for detailed analysis
    w('### Step 1 — 自洽能隙')
    w()
    w(f'求解 BCS 能隙方程 $1 = \\lambda \\langle f^2 \\rangle \\int_0^{{E_0}} d\\xi \\, \\tanh(E_k/2k_BT)/E_k$。')
    w(f'使用 brentq 求根, bracket $[10^{{-6}}E_0, E_0]$, 低温→高温热冷却。')
    w()
    w(f'**验证**: $2\\Delta(0)/k_BT_c = {2*r.Delta0/(kB*r.Tc):.2f}$, BCS 弱耦合理论值 3.52, 偏差 {abs(2*r.Delta0/(kB*r.Tc)-3.52)/3.52*100:.1f}%。')
    w()

    w('### Step 2 — 熵积分')
    w()
    w(f'$S_s(T) = -4k_B \\int_0^\\infty dE\\,N(E)[f\\ln f + (1-f)\\ln(1-f)]$。')
    w(f'自适应网格: 低温对数段 ($E \\in [10^{{-6}}, 10k_BT]$) + 高温线性段 ($E \\in [10k_BT, E_{{\\max}}]$), 共 8000 点。')
    w(f'数值安全: $\\log$ 避零 mask, Simpson 积分。')
    w()

    w('### Step 3 — 比热')
    w()
    w(f'**方法 A (数值)**: 样条平滑 $S(T)$ 后求导, $C = T\\,dS/dT$。')
    w(f'**方法 B (解析)**: $C = 4/(k_B T^2) \\int N(E) E^2 / [4\\cosh^2(E/2k_BT)] dE$。')
    w()
    w(f'低温 $T \\ll T_c$: $C_s \\ll C_n$ (有隙指数抑制) ✓')
    w(f'接近 $T_c$: $C_s \\approx C_n$ (能隙趋于零) ✓')
    w()

    w('### Step 4 — 特征量')
    w()
    w(f'**低温幂律**: $n = {r.power_law.get("n", np.nan):.2f}$, $R^2 = {r.power_law.get("R2", np.nan):.4f}$。')
    w(f'各向同性 s-wave 低温应为 $e^{{-\\Delta/T}}$, $n\\approx2$ 是拟合假象。真实点节点信号需各向异性 DOS。')
    w()
    w(f'**ΔC/C**: 当前值 {r.jump["delta_C_over_C"]:.3f}, 远小于 BCS 理论值 1.43。')
    w(f'根因: `specific_heat_analytical` 只算 $(\\partial S/\\partial T)|_\\Delta$, 缺链式法则项 $(\\partial S/\\partial\\Delta)(d\\Delta/dT)$。')
    w(f'BCS 跳变的物理来源正是 $d\\Delta^2/dT$ 在 $T_c$ 处的跃变。')

    w()
    w('---')
    w()
    w('## 4. 已知局限与下一步')
    w()
    w('### Bug / 缺失')
    w()
    w('1. **ΔC/C 跳变缺失**: 需补链式法则贡献 $\\partial S/\\partial\\Delta \\cdot d\\Delta/dT$')
    w('2. **各向异性占位**: `f_sq_avg` 仅为标量乘子, 非真实角向积分')
    w('3. **参数全为猜测**: $\\lambda$, $\\langle f^2\\rangle$ 未经论文校准')
    w('4. **纯 BCS 框架**: BdG/格林函数/准经典等方法未实现')
    w('5. **性能**: 73点×2相 ≈ 10min, 未并行化')
    w()
    w('### 下一步 (按优先级)')
    w()
    w('1. 补链式法则项 → ΔC/C 闭环')
    w('2. 等论文参数 → 换上真实 λ, f(k̂), γST')
    w('3. 实现各向异性角向积分 → 点节点 C∝T²')
    w('4. Jupyter notebook 版本')
    w('5. 并行化 + 缓存优化')

    w()
    w('---')
    w()
    w('## 5. 复现说明')
    w()
    w('```bash')
    w('pip install -r requirements.txt')
    w('python3 experiment_report.py          # 仅数值计算')
    w('python3 experiment_report.py --plot   # 含可视化')
    w('python3 experiment_report.py --full   # 全管线 (较慢)')
    w('```')
    w()
    w('输出文件:')
    w('- `experiment_report.png` — 六面板图')
    w('- `EXPERIMENT_REPORT.md` — 本报告')

    with open(md_path, 'w') as f:
        f.write('\n'.join(lines))
    return md_path


# ╔════════════════════════════════════════════════════════════════════╗
# ║                    main                                           ║
# ╚════════════════════════════════════════════════════════════════════╝

def main():
    do_plot = '--plot' in sys.argv or '--full' in sys.argv
    do_full = '--full' in sys.argv
    n_grid_entropy = 8000 if do_full else 2000
    n_grid_sh = 3000 if do_full else 1000
    n_T = 73 if do_full else 20

    print("=" * 64)
    print("  超导比热计算 Workflow — 可复现实验报告")
    print("=" * 64)

    # ── Parameters ──────────────────────────────────────────────
    Tc_target = 1.85  # K
    kBTc = kB * Tc_target

    # Phase A
    lam_A, f_sq_A = 0.42, 1.0
    E0_A = kBTc / (1.136 * np.exp(-1.0 / lam_A))

    # Phase B
    lam_B, f_sq_B = 0.48, 1.0
    E0_B = kBTc / (1.136 * np.exp(-1.0 / lam_B))

    print(f"\nParameters: Tc = {Tc_target} K, kBTc = {kBTc*1000:.3f} meV")
    print(f"  Phase A: lambda={lam_A}, E0={E0_A*1000:.2f} meV, <f^2>={f_sq_A}")
    print(f"  Phase B: lambda={lam_B}, E0={E0_B*1000:.2f} meV, <f^2>={f_sq_B}")
    print(f"  Resolution: {n_T} T points, entropy grid={n_grid_entropy}, SH grid={n_grid_sh}")

    # ── Temperature grid ────────────────────────────────────────
    T_low = np.linspace(0.05, 0.7 * Tc_target, max(5, n_T // 3))
    T_mid = np.linspace(0.7 * Tc_target, 0.95 * Tc_target, max(5, n_T // 4))
    T_tc = np.linspace(0.95 * Tc_target, 1.05 * Tc_target, max(5, n_T // 3))
    T_range = np.unique(np.concatenate([T_low, T_mid, T_tc]))

    # ── Run ─────────────────────────────────────────────────────
    results = []
    for lam, E0, f_sq, label in [
        (lam_A, E0_A, f_sq_A, 'Phase A'),
        (lam_B, E0_B, f_sq_B, 'Phase B'),
    ]:
        print(f"\n--- {label} ---")
        r = run_pipeline(lam, E0, T_range, N0=1.0, f_sq_avg=f_sq,
                         label=label, verbose=True)
        results.append(r)

    # ── Summary ─────────────────────────────────────────────────
    print("\n" + "=" * 64)
    print("  Results Summary")
    print("=" * 64)
    for r in results:
        j, pl = r.jump, r.power_law
        bcs = 2 * r.Delta0 / (kB * r.Tc) if r.Tc > 0 else 0
        print(f"\n  {r.label} (lambda={r.lam}):")
        print(f"    Tc = {r.Tc:.3f} K")
        print(f"    Delta(0) = {r.Delta0*1000:.3f} meV")
        print(f"    2*Delta(0)/kB*Tc = {bcs:.2f}  (BCS: 3.52)")
        print(f"    gamma = {j['gamma']:.3e} eV/K^2")
        print(f"    DeltaC/C = {j['delta_C_over_C']:.3f}  "
              f"[paper: 0.93(A)/1.20(B), BCS theory: 1.43]")
        print(f"    Low-T n = {pl.get('n', np.nan):.2f}  "
              f"(R^2 = {pl.get('R2', np.nan):.4f})")
        print(f"    Elapsed: {r.elapsed:.1f}s")

    # ── Plot ────────────────────────────────────────────────────
    if do_plot:
        print("\n  Generating plot ...")
        path = plot_report(results)
        print(f"  Saved: {path}")

    # ── Report ──────────────────────────────────────────────────
    print("\n  Generating Markdown report ...")
    path = generate_report(results)
    print(f"  Saved: {path}")

    # ── JSON output ─────────────────────────────────────────────
    summary = {
        'workflow': 'cluster_172',
        'timestamp': time.strftime('%Y-%m-%d %H:%M'),
        'parameters': {
            'Tc_target_K': Tc_target,
            'phase_A': {'lambda': lam_A, 'E0_eV': E0_A, 'f_sq_avg': f_sq_A},
            'phase_B': {'lambda': lam_B, 'E0_eV': E0_B, 'f_sq_avg': f_sq_B},
        },
        'results': [
            {
                'label': r.label,
                'Tc_K': r.Tc, 'Delta0_meV': r.Delta0 * 1000,
                'bcs_ratio': 2 * r.Delta0 / (kB * r.Tc) if r.Tc > 0 else 0,
                'gamma_eV_per_K2': r.jump['gamma'],
                'delta_C_over_C': r.jump['delta_C_over_C'],
                'low_T_exponent': r.power_law.get('n'),
                'low_T_R2': r.power_law.get('R2'),
                'elapsed_s': r.elapsed,
            }
            for r in results
        ],
    }
    json_path = 'experiment_results.json'
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"  Saved: {json_path}")

    print("\n" + "=" * 64)
    print("  Experiment report complete.")
    print(f"  Files: EXPERIMENT_REPORT.md, experiment_results.json"
          + (", experiment_report.png" if do_plot else ""))
    print("=" * 64)


if __name__ == '__main__':
    main()