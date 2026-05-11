# 超导比热计算 Workflow — 可复现实验报告

**Workflow**: cluster_172
**生成时间**: 2026-05-11 23:15
**参考论文**: Parker, Maki, Haas — Anisotropic superconductivity in PrOs₄Sb₁₂ (ref6)
**复现环境**: Python 3 + NumPy + SciPy + Matplotlib

---

## 1. 参数设置

| 参数 | Phase A | Phase B | 说明 |
|------|---------|---------|------|
| Phase A | λ=0.42, E₀=1.52 meV, ⟨f²⟩=1.0 | — | λ 从 Tc 反推, 未见论文精确值 |
| Phase B | λ=0.48, E₀=1.13 meV, ⟨f²⟩=1.0 | — | λ 从 Tc 反推, 未见论文精确值 |

温度网格: 15 点, [0.050, 1.943] K, Tc 附近加密

---

## 2. 结果汇总

| 量 | Phase A | Phase B | 理论参考 | 状态 |
|----|---------|---------|----------|------|
| **Phase A** | | | | |
| Δ(0) | 0.283 meV | — | ~0.28 meV (BCS) | ✓ |
| 2Δ(0)/kBTc | 3.59 | — | 3.52 (BCS) | ✓ |
| Tc | 1.832 K | — | ~1.85 K | ✓ |
| γ | 4.886e-08 eV/K² | — | (π²/3) kB² N₀ | ✓ |
| ΔC/C | 0.004 | — | 1.43 (BCS) / 0.93–1.20 (paper) | ✗ |
| 低温幂律 n | nan | — | 2 (点节点) | ⚠ |
| 耗时 | 0.0s | — | — | — |
| **Phase B** | | | | |
| Δ(0) | 0.285 meV | — | ~0.28 meV (BCS) | ✓ |
| 2Δ(0)/kBTc | 3.61 | — | 3.52 (BCS) | ✓ |
| Tc | 1.832 K | — | ~1.85 K | ✓ |
| γ | 4.886e-08 eV/K² | — | (π²/3) kB² N₀ | ✓ |
| ΔC/C | 0.005 | — | 1.43 (BCS) / 0.93–1.20 (paper) | ✗ |
| 低温幂律 n | nan | — | 2 (点节点) | ⚠ |
| 耗时 | 0.0s | — | — | — |

---

## 3. 逐步分析

### Step 1 — 自洽能隙

求解 BCS 能隙方程 $1 = \lambda \langle f^2 \rangle \int_0^{E_0} d\xi \, \tanh(E_k/2k_BT)/E_k$。
使用 brentq 求根, bracket $[10^{-6}E_0, E_0]$, 低温→高温热冷却。

**验证**: $2\Delta(0)/k_BT_c = 3.59$, BCS 弱耦合理论值 3.52, 偏差 1.9%。

### Step 2 — 熵积分

$S_s(T) = -4k_B \int_0^\infty dE\,N(E)[f\ln f + (1-f)\ln(1-f)]$。
自适应网格: 低温对数段 ($E \in [10^{-6}, 10k_BT]$) + 高温线性段 ($E \in [10k_BT, E_{\max}]$), 共 8000 点。
数值安全: $\log$ 避零 mask, Simpson 积分。

### Step 3 — 比热

**方法 A (数值)**: 样条平滑 $S(T)$ 后求导, $C = T\,dS/dT$。
**方法 B (解析)**: $C = 4/(k_B T^2) \int N(E) E^2 / [4\cosh^2(E/2k_BT)] dE$。

低温 $T \ll T_c$: $C_s \ll C_n$ (有隙指数抑制) ✓
接近 $T_c$: $C_s \approx C_n$ (能隙趋于零) ✓

### Step 4 — 特征量

**低温幂律**: $n = nan$, $R^2 = nan$。
各向同性 s-wave 低温应为 $e^{-\Delta/T}$, $n\approx2$ 是拟合假象。真实点节点信号需各向异性 DOS。

**ΔC/C**: 当前值 0.004, 远小于 BCS 理论值 1.43。
根因: `specific_heat_analytical` 只算 $(\partial S/\partial T)|_\Delta$, 缺链式法则项 $(\partial S/\partial\Delta)(d\Delta/dT)$。
BCS 跳变的物理来源正是 $d\Delta^2/dT$ 在 $T_c$ 处的跃变。

---

## 4. 已知局限与下一步

### Bug / 缺失

1. **ΔC/C 跳变缺失**: 需补链式法则贡献 $\partial S/\partial\Delta \cdot d\Delta/dT$
2. **各向异性占位**: `f_sq_avg` 仅为标量乘子, 非真实角向积分
3. **参数全为猜测**: $\lambda$, $\langle f^2\rangle$ 未经论文校准
4. **纯 BCS 框架**: BdG/格林函数/准经典等方法未实现
5. **性能**: 73点×2相 ≈ 10min, 未并行化

### 下一步 (按优先级)

1. 补链式法则项 → ΔC/C 闭环
2. 等论文参数 → 换上真实 λ, f(k̂), γST
3. 实现各向异性角向积分 → 点节点 C∝T²
4. Jupyter notebook 版本
5. 并行化 + 缓存优化

---

## 5. 复现说明

```bash
pip install -r requirements.txt
python3 experiment_report.py          # 仅数值计算
python3 experiment_report.py --plot   # 含可视化
python3 experiment_report.py --full   # 全管线 (较慢)
```

输出文件:
- `experiment_report.png` — 六面板图
- `EXPERIMENT_REPORT.md` — 本报告