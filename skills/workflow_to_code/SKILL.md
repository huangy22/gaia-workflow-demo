# SKILL.md — Workflow JSON → Python 代码管线

## 概述

将一个"工作流聚类综述"的结构化 JSON（符合 `schema.json`）转化为可执行的 Python 数值管线。

## 输入

- `review_cluster_*.json`：结构化 workflow 数据
- 原始 Markdown 文件（可选，用于补全公式细节）

## 输出

单个 `.py` 文件，包含：
1. 物理常数
2. 每个 step 的模块化函数
3. 端到端 `run_workflow()` 管线
4. 可视化模块 `plot_results()`
5. `__main__` 示例入口

## 流程

### Phase 1：识别可代码化的 Step

遍历 `steps[]`，判断哪些 step 是纯计算逻辑（方程求解、积分、求导、拟合），哪些是非代码步骤（数据采集、人工决策）。

**可代码化**的标志：
- description 含数学公式（LaTeX）
- substeps 描述具体的数值操作
- key_parameters 给出可取值
- 输出是数值量

**不可直接代码化**（改为注释/接口桩）：
- 依赖外部实验数据输入
- 人为判断节点
- 未给出具体算法的通用描述

### Phase 2：抽取公式

从 step/substep 的 `description` 字段中提取 LaTeX 公式，作为代码的核心数学逻辑。

**要点**：
- 每个 `$...$` 和 `$$...$$` 块对应一个数值操作
- 公式中的符号映射到 Python 变量名（`\Delta` → `Delta`，`\int` → `quad/simpson`）
- 保留公式编号作为代码注释引用

### Phase 3：确定数值方法

根据公式类型选择实现策略：

| 数学操作 | Python 实现 |
|----------|------------|
| 定积分 $\int_a^b f(x)dx$ | `scipy.integrate.quad` (光滑函数) 或 `scipy.integrate.simpson` (离散网格) |
| 非线性求根 $f(x)=0$ | `scipy.optimize.root_scalar` (brentq / secant) |
| 导数 $df/dx$ | 中心差分 或 `scipy.interpolate.UnivariateSpline.derivative()` |
| 曲线拟合 | `scipy.optimize.curve_fit` |
| 角向平均 $\langle\cdot\rangle_{\hat{k}}$ | 数值角向积分, 离散求和 |

### Phase 4：处理数值陷阱

从 `common_pitfalls` 字段提取信息，在代码中加入防御：

| 陷阱 | 代码对策 |
|------|---------|
| 低温 $f_{FD}$ 数值病态 | $\log f$ 加 mask 防 $\log(0)$；低温段对数网格加密 |
| 积分截断不足 | $E_{\max} \ge 20$–$50\,k_B T_{\max}$，做收敛扫描 |
| 数值微分放大噪声 | 样条平滑后求导；或使用 $\partial f_{FD}/\partial T$ 解析代入 |
| 自洽方程收敛到平凡解 $\Delta=0$ | 提供非零初值猜测；检查 bracket 端点符号 |
| 简并因子/前因子约定不一致 | 明确 $N(E)$ 是否含自旋简并，在注释中标明 |

### Phase 5：单位一致性

这是最容易出 bug 的环节。必须明确全局单位体系：

1. **统一能量单位**：所有公式使用一致的单位（推荐 eV）
2. **温度→能量转换**：$k_B T$ 处处一致，$k_B = 8.617\times10^{-5}$ eV/K
3. **耦合常数 $\lambda$ 无量纲**
4. **$E_0$ 截断能**：从 $k_B T_c = 1.136\,E_0\,e^{-1/\lambda}$ 反推时注意单位（$E_0$ 是能量，不是温度）
5. **检查比值**：如 $2\Delta(0)/k_B T_c$，BCS 弱耦合 ≈ 3.52

### Phase 6：代码结构

```python
"""
{workflow.name} — 完整数值管线
参考: {paper_refs[foundational]}
"""

# === 物理常数 ===
kB = 8.617333262145e-5  # eV/K

# === Step 1: {name} ===
def step1_core(...):
    """公式来源: step_*.substeps[*].description"""
    ...

# === Step 2: {name} ===
def step2_core(...):
    ...

# === 端到端管线 ===
@dataclass
class WorkflowResult:
    T_range: ndarray
    Delta: ndarray
    S: ndarray
    C: ndarray
    Tc: float
    ...

def run_workflow(lam, E0, T_range, ...):
    ...

# === 可视化 ===
def plot_results(results):
    ...

# === 示例 ===
if __name__ == '__main__':
    # 参数来源: workflow.claims + 文献
    ...
    result = run_workflow(...)
    print(result)
    plot_results(result)
```

### Phase 7：校准与验证

1. **检查已知结果**：如 $\Delta C/C$、$T_c$ 公式
2. **反推缺失参数**：若论文报告了可观测值但没有直接给 $\lambda$，从 BCS 公式反推
3. **明确标注参数来源**：
   - 直接来自论文 → 写论文出处
   - 反推估计 → 标注"fit to paper result"
   - 缺失 → 标注"需要原始论文获取"，给出合理默认值

## 常见错误

1. **单位混用** — 最常见：$E_0$ 从 $T_c$ 公式得到的是温度还是能量？必须显式乘 $k_B$
2. **数值溢出** — $e^{E/k_B T} \to \infty$ 在低温下，用 `np.clip(x, -100, 100)` 保护
3. **平凡解** — 自洽方程 $\Delta=0$ 总是数学解，需要用 bracket 排除
4. **过早截断** — 对核函数加 `Ek/kBT > threshold` 的 return 0 可能切掉 gap edge 贡献
5. **熵在低温趋于 0 是物理的** — 对有隙超导体，$T \ll T_c$ 时 $S \sim e^{-\Delta/T}$，不是 bug