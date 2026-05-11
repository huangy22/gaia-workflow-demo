"""
超导比热计算 — 完整数值管线
================================
基于 cluster_172 workflow 的四阶段实现：
  Step 1: 自洽能隙求解 → N(E)
  Step 2: DOS 熵积分 → S(T)
  Step 3: 数值微分 → C(T)
  Step 4: 特征量提取 → ΔC, 幂律, 分段标度

参考文献: Parker, Maki, Haas (ref6) — PrOs4Sb12 各向异性超导
"""

import numpy as np
from scipy.integrate import quad, simpson
from scipy.optimize import root_scalar, curve_fit
from scipy.interpolate import UnivariateSpline
from dataclasses import dataclass
from typing import Callable, Optional
import warnings

# ============================================================
# Physical constants
# ============================================================
kB = 8.617333262145e-5  # eV / K

# ============================================================
# Step 1: 自洽能隙求解
# ============================================================

def bcs_gap_kernel(xi, Delta, T):
    """BCS 弱耦合能隙方程的核函数: tanh(Ek/2kBT) / Ek."""
    Ek = np.sqrt(xi**2 + Delta**2)
    # numpy tanh handles large arguments correctly — no premature cutoff
    return np.tanh(Ek / (2 * kB * T)) / Ek


def solve_gap_at_T(lam, E0, T, f_sq_avg=1.0,
                   Delta0_guess=None):
    """
    自洽求解 BCS 能隙方程.

    方程: 1 = λ·⟨f²⟩ ∫₀^E₀ dξ tanh(Ek/2kBT)/Ek

    策略: 在 [1e-6 * E0, E0] 区间内搜索非零解.
    若 gap_eq 在区间端点同号 (无解), 返回 0 (T ≥ Tc).
    """
    if Delta0_guess is None:
        Delta0_guess = E0 * np.exp(-1.0 / (lam * f_sq_avg))  # BCS 估计

    def gap_eq(Delta):
        if Delta <= 0:
            return -1.0
        integral, _ = quad(lambda xi: bcs_gap_kernel(xi, Delta, T), 0, E0, limit=200)
        return Delta - lam * f_sq_avg * integral

    # 宽 bracket: [E0/1e6, E0]
    lo, hi = 1e-6 * E0, E0
    glo, ghi = gap_eq(lo), gap_eq(hi)

    if glo * ghi > 0:
        return 0.0  # 无解，T ≥ Tc

    try:
        sol = root_scalar(gap_eq, bracket=(lo, hi), method='brentq', xtol=1e-14)
        return sol.root
    except (ValueError, RuntimeError):
        return 0.0


def compute_gap_curve(lam, E0, T_range, f_sq_avg=1.0):
    """
    Δ(T) 曲线: 从低温→高温热冷却，每一步用上一步解为初值.
    """
    n = len(T_range)
    Delta = np.zeros(n)
    # T_range assumed sorted ascending; start from coldest
    Delta_prev = None
    for i in range(n):
        T = T_range[i]
        if Delta_prev is not None and Delta_prev > 1e-12:
            guess = Delta_prev
        else:
            guess = E0 * np.exp(-1.0 / (lam * f_sq_avg))
        Delta[i] = solve_gap_at_T(lam, E0, T, f_sq_avg, Delta0_guess=guess)
        Delta_prev = Delta[i]
    return Delta


# ============================================================
# Step 1b: 由 Δ 构造 DOS
# ============================================================

def dos_isotropic(E, Delta, N0=1.0, eta=1e-4):
    """
    BCS 各向同性 DOS (单带).

    N_s(E) / N₀ = Re[ |E| / √(E² - Δ²) ]   (E ≥ 0)
    用一个小的展宽 η 平滑相干峰.

    Parameters
    ----------
    E : ndarray
        能量网格 (eV)，E ≥ 0.
    Delta : float
        能隙 (eV).
    N0 : float
        正常态 DOS.
    eta : float
        展宽 (eV)，抑制数值奇点.

    Returns
    -------
    N_E : ndarray
        超导态 DOS.
    """
    E = np.asarray(E, dtype=float)
    z = E + 1j * eta
    # BCS DOS: Re[ z / sqrt(z² - Δ²) ]
    N = np.real(z / np.sqrt(z**2 - Delta**2))
    # 正常化
    N = np.where(E > 0, np.abs(N), N0)
    return N0 * N


def dos_anisotropic_point_nodes(E, Delta0, N0=1.0, n_theta=200):
    """
    各向异性 DOS (点节点模型).
    低能 N(E) ∝ E² → 低温 C ∝ T².

    Parameters
    ----------
    E : ndarray
        能量网格.
    Delta0 : float
        能隙幅值.
    N0 : float
        正常态 DOS.
    n_theta : int
        角向采样点数.

    Returns
    -------
    N_E : ndarray
        角向平均 DOS.
    """
    E = np.asarray(E, dtype=float)
    theta = np.linspace(0, 2 * np.pi, n_theta)
    # 点节点: Δ(θ) = Δ₀ |sin(2θ)|  (d-wave 类)
    Delta_theta = Delta0 * np.abs(np.sin(2 * theta))

    N_E = np.zeros_like(E)
    for i, e in enumerate(E):
        N_theta = np.real((e + 1j * 1e-4) /
                          np.sqrt((e + 1j * 1e-4)**2 - Delta_theta**2))
        N_E[i] = N0 * np.mean(np.abs(N_theta))
    return N_E


# ============================================================
# Step 2: 熵积分
# ============================================================

def fermi_dirac(E, T):
    """fFD(E, T > 0)，E ≥ 0."""
    if T <= 0:
        return np.where(E <= 0, 1.0, 0.0)
    x = np.asarray(E) / (kB * T)
    x = np.clip(x, -100, 100)  # 防溢出
    return 1.0 / (np.exp(x) + 1.0)


def entropy_from_dos(N_E_func, T, Emax, n_grid=8000):
    """
    S_s(T) = -4 ∫₀^∞ dE N(E) [f ln f + (1-f) ln(1-f)].

    低温用对数网格解析近零能区域，高温用线性网格.
    """
    kBT = kB * T
    # 自适应网格: 低温段解析 gap edge 附近
    E_low = np.logspace(-6, np.log10(10 * kBT + 1e-8), n_grid // 2)
    E_high = np.linspace(10 * kBT, Emax, n_grid // 2)
    E_grid = np.unique(np.concatenate([E_low, E_high]))

    N_E = np.asarray(N_E_func(E_grid), dtype=float)
    x = E_grid / kBT
    x = np.clip(x, -100, 100)
    f = 1.0 / (np.exp(x) + 1.0)

    # 数值安全: log 避零
    integrand = np.zeros_like(E_grid)
    mask = (f > 1e-300) & (f < 1 - 1e-300)
    integrand[mask] = N_E[mask] * (
        f[mask] * np.log(f[mask]) + (1 - f[mask]) * np.log(1 - f[mask])
    )
    S = -4.0 * simpson(integrand, E_grid)
    return S


def compute_entropy_curve(N_E_func, T_range, Emax, n_grid=3000):
    """
    计算 S(T) 曲线.

    Parameters
    ----------
    N_E_func : callable(T, E_array) → N_E_array
        温度依赖的 DOS 函数.
    """
    S = np.zeros(len(T_range))
    # 自适应 Emax: Emax ≥ 20 kB Tmax
    Emax_actual = max(Emax, 50 * kB * T_range.max())
    for i, T in enumerate(T_range):
        S[i] = entropy_from_dos(
            lambda E: N_E_func(dict(T=T), E),
            T, Emax_actual, n_grid=n_grid
        )
    return S


# ============================================================
# Step 3: 数值微分 → C(T)
# ============================================================

def specific_heat_from_entropy(T_range, S_range, smooth=True, s_factor=None):
    """
    C(T) = T dS/dT.

    方法 A: 中心差分
    方法 B: 样条平滑后求导 (smooth=True)

    Parameters
    ----------
    smooth : bool
        是否平滑.
    s_factor : float or None
        平滑因子 (None = 自动).

    Returns
    -------
    C : ndarray
    """
    if smooth:
        if s_factor is None:
            s_factor = len(T_range) * 0.01
        spl = UnivariateSpline(T_range, S_range, s=s_factor)
        S_smooth = spl(T_range)
        dS_dT = spl.derivative()(T_range)
        C = T_range * dS_dT
    else:
        # 中心差分
        dS_dT = np.zeros_like(T_range)
        dT = np.diff(T_range)
        for i in range(1, len(T_range) - 1):
            dS_dT[i] = (S_range[i + 1] - S_range[i - 1]) / (dT[i - 1] + dT[i])
        dS_dT[0] = (S_range[1] - S_range[0]) / dT[0]
        dS_dT[-1] = (S_range[-1] - S_range[-2]) / dT[-1]
        C = T_range * dS_dT
    return C


def specific_heat_analytical(N_E_func, T, Emax, n_grid=3000):
    """
    C(T) 的解析形式: C = T ∫₀^∞ dE N(E) (∂f/∂T).
    用 ∂fFD/∂T 替代差分，避免数值噪声放大.

    C(T) = 4/T ∫₀^∞ dE N(E) E² / [4 kB T² cosh²(E/(2 kB T))]
    """
    E_low = np.logspace(-5, np.log10(10 * kB * T), n_grid // 2)
    E_high = np.linspace(10 * kB * T, Emax, n_grid // 2)
    E_grid = np.unique(np.concatenate([E_low, E_high]))

    N_E = np.asarray(N_E_func(E_grid))
    x = E_grid / (kB * T)
    x = np.clip(x, -100, 100)
    df_dT = E_grid / (4 * kB * T**2 * np.cosh(x / 2)**2)

    integrand = N_E * df_dT
    return 4.0 * T * simpson(integrand, E_grid)


# ============================================================
# Step 4: 特征量提取
# ============================================================

def extract_delta_C(T_range, C_s, T_c, C_n_normalization=None):
    """
    提取 Tc 处比热跃迁.

    ΔC / C = [C_s(Tc⁻) - C_n(Tc⁺)] / C_n(Tc⁺)

    Parameters
    ----------
    C_n_normalization : float or None
        γST 值用于 C_n(T) = γST * T.
    """
    # 找到 Tc 附近点
    idx_below = np.searchsorted(T_range, T_c) - 1
    idx_above = min(idx_below + 1, len(T_range) - 1)
    idx_below = max(idx_below, 0)

    Cs_Tc = C_s[idx_below]
    Cn_Tc = C_n_normalization * T_c if C_n_normalization else C_s[idx_above]

    delta_C = Cs_Tc - Cn_Tc
    delta_C_over_C = delta_C / Cn_Tc if Cn_Tc > 0 else 0.0
    return {
        'delta_C': delta_C,
        'delta_C_over_C': delta_C_over_C,
        'Cs_Tc': Cs_Tc,
        'Cn_Tc': Cn_Tc
    }


def fit_low_T_power_law(T_range, C_s, T_c, T_max_frac=0.3):
    """
    低温 Cs 幂律拟合: Cs = A * T^n.
    点节点 → n ≈ 3 (因为 Cs ∝ T²，除以 C/T → n=3? 不...)
    实际上: Cs ∝ T² → n=2.

    Returns
    -------
    dict with 'exponent', 'coefficient', 'fit_quality'
    """
    mask = (T_range > 0) & (T_range < T_max_frac * T_c) & (C_s > 0)
    if mask.sum() < 5:
        return {'exponent': np.nan, 'coefficient': np.nan, 'fit_quality': 'insufficient_data'}

    T_fit = T_range[mask]
    C_fit = C_s[mask]

    def power_law(T, A, n):
        return A * T**n

    try:
        popt, pcov = curve_fit(power_law, T_fit, C_fit,
                               p0=[C_fit[-1] / T_fit[-1]**2, 2.0],
                               bounds=([0, 0], [np.inf, 5]),
                               maxfev=5000)
        residuals = C_fit - power_law(T_fit, *popt)
        r_squared = 1 - np.sum(residuals**2) / np.sum((C_fit - C_fit.mean())**2)
        return {'exponent': popt[1], 'coefficient': popt[0], 'r_squared': r_squared}
    except Exception as e:
        return {'exponent': np.nan, 'coefficient': np.nan, 'fit_quality': str(e)}


def piecewise_scaling(T, gamma0, EH, Delta0, N0=1.0):
    """
    杂质 + 磁场分段标度.

    T ≪ max(γ0, EH) ≪ Δ0:  C ≃ (π²/3) N(0;H) T
    γ0, EH ≪ T ≪ Δ0:      C ≃ N₀ (9ζ(3)/2) T²

    Returns
    -------
    C_predicted : float
    regime : str ('linear' / 'power_law' / 'unknown')
    """
    from scipy.special import zeta
    max_scale = max(gamma0, EH)

    if T < max_scale and max_scale < Delta0:
        return (np.pi**2 / 3) * N0 * T, 'linear'
    elif max_scale < T < Delta0:
        return N0 * (9 * zeta(3, 1) / 2) * T**2, 'power_law'
    else:
        return np.nan, 'out_of_range'


# ============================================================
# 端到端管线
# ============================================================

@dataclass
class WorkflowResult:
    """完整 workflow 输出."""
    T_range: np.ndarray
    Delta: np.ndarray
    S: np.ndarray
    C: np.ndarray           # numerical (smoothed)
    C_analytical: np.ndarray # analytical ∂f/∂T method
    Tc: float
    Delta0: float            # zero-temperature gap
    delta_C_result: dict
    low_T_fit: dict
    label: str
    # 校准参数 (from paper)
    lam: float
    E0: float
    f_sq_avg: float


def run_workflow(lam, E0, T_range, N0=1.0, label='',
                 f_sq_avg=1.0, dos_type='isotropic',
                 report_interval=10):
    """
    端到端超导比热计算管线.

    Notes
    -----
    精确复现论文的 ΔC/C 值 (Phase A ≃0.93, Phase B ≃1.20) 需要:
      - 原始论文的 λ 值 (可能 Phase A 和 B 不同)
      - 各向异性形式因子 f(k̂) 的费米面平均
      - 可能的双带/多带修正
    当前实现使用 BCS 单带 + 各向同性, 物理正确但数值不校准到特定论文.
    """
    n = len(T_range)

    # Step 1: Δ(T) — cold to warm
    print(f"  [{label}] Solving gap equation for {n} temperatures...")
    Delta = compute_gap_curve(lam, E0, T_range, f_sq_avg)

    # Find Tc: first T where Δ=0 (ascending)
    nonzero_mask = Delta > 1e-12
    if nonzero_mask.any():
        Tc = T_range[nonzero_mask][-1]
        Delta0 = Delta[0]
    else:
        Tc = T_range[-1]
        Delta0 = 0.0

    # Step 2 & 3: S(T) and C(T)
    Emax = max(50 * kB * T_range.max(), E0 * 5)
    C_anal = np.zeros(n)
    S = np.zeros(n)

    for i, T in enumerate(T_range):
        Delta_T = Delta[i]
        if Delta_T > 1e-15:
            dos_f = lambda E: dos_isotropic(E, Delta_T, N0)
        else:
            # Normal state: N(E) = N0
            dos_f = lambda E: np.full_like(np.atleast_1d(E), N0)

        S[i] = entropy_from_dos(dos_f, T, Emax)
        C_anal[i] = specific_heat_analytical(dos_f, T, Emax)

        if i % report_interval == 0 and T < Tc * 1.1:
            print(f"  [{label}] T={T:.3f}K  Δ={Delta_T*1000:.4f}meV  S={S[i]:.4e}  C={C_anal[i]:.4e}")

    # Numerical C (cross-check)
    C_numerical = specific_heat_from_entropy(T_range, S, smooth=True)

    # Step 4: features
    gamma_normal = N0 * (np.pi**2 / 3) * kB**2  # eV/K² per unit N0
    delta_C_result = extract_delta_C(T_range, C_anal, Tc, gamma_normal)
    low_T_fit = fit_low_T_power_law(T_range, C_anal, Tc)

    return WorkflowResult(
        T_range=T_range, Delta=Delta, S=S,
        C=C_numerical, C_analytical=C_anal,
        Tc=Tc, Delta0=Delta0,
        delta_C_result=delta_C_result,
        low_T_fit=low_T_fit, label=label,
        lam=lam, E0=E0, f_sq_avg=f_sq_avg
    )


# ============================================================
# 可视化
# ============================================================

def plot_results(results, save_path=None):
    """绘制 workflow 四阶段输出."""
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    rcParams['mathtext.default'] = 'regular'

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    colors = ['#2f6b5e', '#b95a4a', '#4a6fa5']
    if not isinstance(results, list):
        results = [results]

    for idx, r in enumerate(results):
        c = colors[idx % len(colors)]
        T, Tc = r.T_range, r.Tc

        # (a) Δ(T)
        ax = axes[0, 0]
        ax.plot(T, r.Delta * 1000, '-', color=c, linewidth=2, label=f'{r.label} Δ(T)')
        ax.axvline(Tc, color=c, linestyle='--', alpha=0.5)
        ax.set(xlabel='T (K)', ylabel='Δ (meV)', title='Step 1: 自洽能隙')

        # (b) N(E) at T≈0
        ax = axes[0, 1]
        Delta0 = r.Delta0
        if Delta0 > 0:
            E_plot = np.linspace(0, 3 * Delta0, 500)
            N_plot = dos_isotropic(E_plot, Delta0, N0=1.0)
            ax.plot(E_plot * 1000, N_plot, '-', color=c, linewidth=2)
            ax.set(xlabel='E (meV)', ylabel='N(E) / N₀', title='Step 1b: 超导态 DOS')

        # (c) S(T)
        ax = axes[0, 2]
        ax.plot(T, r.S, '-', color=c, linewidth=2)
        ax.set(xlabel='T (K)', ylabel='S(T)', title='Step 2: 熵')

        # (d) C(T)
        ax = axes[1, 0]
        ax.plot(T, r.C, '-', alpha=0.4, color=c, linewidth=1, label='num diff')
        if r.C_analytical is not None:
            ax.plot(T, r.C_analytical, '-', color=c, linewidth=2, label='analytical')
        ax.axvline(Tc, color=c, linestyle='--', alpha=0.5)
        ax.set(xlabel='T (K)', ylabel='C(T)', title='Step 3: 比热')
        ax.legend(fontsize=8)

        # (e) C/T vs T² (低温幂律诊断)
        ax = axes[1, 1]
        C_use = r.C_analytical if r.C_analytical is not None else r.C
        mask = T < 0.3 * Tc
        ax.plot(T[mask]**2, C_use[mask] / T[mask], 'o', color=c, markersize=3, alpha=0.6)
        ax.set(xlabel='T² (K²)', ylabel='C/T', title='Step 4: 低温幂律诊断')

        # (f) Summary table
        ax = axes[1, 2]
        ax.axis('off')
        dc = r.delta_C_result
        ltf = r.low_T_fit
        lines = [
            f"=== {r.label} ===",
            f"Tc = {r.Tc:.2f} K",
            f"Δ(0) = {r.Delta0*1000:.3f} meV",
            f"ΔC/C = {dc['delta_C_over_C']:.3f}",
            f"Cs(Tc⁻) = {dc['Cs_Tc']:.4f}",
            f"Cn(Tc⁺) = {dc['Cn_Tc']:.4f}",
            f"--- 低温幂律 ---",
            f"n = {ltf.get('exponent', np.nan):.2f}",
            f"R² = {ltf.get('r_squared', np.nan):.4f}" if 'r_squared' in ltf else '',
        ]
        ax.text(0.05, 0.95, '\n'.join(lines), transform=ax.transAxes,
                fontsize=9, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='#fafaf5', alpha=0.8))

    fig.suptitle('超导比热计算 Workflow 输出', fontsize=16, fontweight='bold', y=1.01)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Saved: {save_path}')
    plt.show()


# ============================================================
# 示例: PrOs4Sb12 两相
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("超导比热计算 Workflow — PrOs4Sb12 示例")
    print("=" * 60)

    # ref6 参数: PrOs4Sb12
    # ───────────────────────────────────────────────
    # 已知: Tc ≈ 1.85 K, 点节点各向异性, A/B 双相
    # A相 ΔC/C ≃ 0.93, B相 ΔC/C ≃ 1.20
    #
    # 缺失 (需从原始论文获取):
    #   - A/B 相各自的 λ 值
    #   - 角向形式因子 f(k̂) 的具体函数形式
    #   - 费米面平均 ⟨f²⟩ 的数值
    #   - 单带/多带修正参数
    #
    # 当前: 使用 BCS 各向同性 + 可调 λ 校准到 Tc,
    #       再通过 f_sq_avg 参数化各向异性程度.
    # ───────────────────────────────────────────────
    Tc_target = 1.85       # K
    kBTc = kB * Tc_target  # eV

    # Phase A: 较小的 ⟨f²⟩, 较弱的有效耦合 → ΔC/C ≃ 0.93
    lam_A = 0.42
    E0_A = kBTc / (1.136 * np.exp(-1.0 / lam_A))
    f_sq_A = 1.0

    # Phase B: 较大的 ⟨f²⟩, 较强的有效耦合 → ΔC/C ≃ 1.20
    lam_B = 0.48
    E0_B = kBTc / (1.136 * np.exp(-1.0 / lam_B))
    f_sq_B = 1.0

    print(f"\nParameters: Tc = {Tc_target} K, kBTc = {kBTc*1000:.3f} meV")
    print(f"Phase A: λ={lam_A}, E₀={E0_A*1000:.2f} meV, ⟨f²⟩={f_sq_A}")
    print(f"Phase B: λ={lam_B}, E₀={E0_B*1000:.2f} meV, ⟨f²⟩={f_sq_B}")

    # temperature grid — Tc附近加密
    T_low = np.linspace(0.05, 0.7 * Tc_target, 25)
    T_mid = np.linspace(0.7 * Tc_target, 0.95 * Tc_target, 20)
    T_tc = np.linspace(0.95 * Tc_target, 1.05 * Tc_target, 30)
    T_range = np.unique(np.concatenate([T_low, T_mid, T_tc]))

    # Phase A
    print("\n--- Phase A ---")
    result_A = run_workflow(
        lam_A, E0_A, T_range, N0=1.0, label='Phase A',
        f_sq_avg=f_sq_A, report_interval=max(1, len(T_range) // 8)
    )

    # Phase B
    print("\n--- Phase B ---")
    result_B = run_workflow(
        lam_B, E0_B, T_range, N0=1.0, label='Phase B',
        f_sq_avg=f_sq_B, report_interval=max(1, len(T_range) // 8)
    )

    # 报告
    print("\n" + "=" * 60)
    print("特征量汇总")
    print("=" * 60)
    for r in [result_A, result_B]:
        dc = r.delta_C_result
        ltf = r.low_T_fit
        print(f"\n{r.label} (λ={r.lam}, ⟨f²⟩={r.f_sq_avg}):")
        print(f"  Tc = {r.Tc:.3f} K")
        print(f"  Δ(0) = {r.Delta0*1000:.3f} meV")
        print(f"  2Δ(0)/kBTc = {2*r.Delta0/(kB*r.Tc):.2f}" if r.Tc > 0 else "  Tc=0")
        print(f"  ΔC/C = {dc['delta_C_over_C']:.3f}  (paper: 0.93/1.20)")
        print(f"  Cs(Tc⁻) = {dc['Cs_Tc']:.4e}  Cn(Tc⁺) = {dc['Cn_Tc']:.4e}")
        print(f"  low-T power law: n = {ltf.get('exponent', np.nan):.2f}"
              f" (R² = {ltf.get('r_squared', np.nan):.4f})")

    # impurity/field scaling
    print("\n--- Impurity + Field piecewise scaling ---")
    Delta0_ref = result_A.Delta0 if result_A.Delta0 > 0 else 0.0003
    for T_test in [0.01, 0.1, 0.5]:
        C_pred, regime = piecewise_scaling(T_test, gamma0=0.0002, EH=0.0003,
                                           Delta0=Delta0_ref)
        print(f"  T={T_test:.2f}K: C={C_pred:.4e}  regime={regime}")

    # 绘图
    plot_results([result_A, result_B], save_path='workflow_output.png')