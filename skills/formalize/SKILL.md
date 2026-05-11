# SKILL.md — Workflow Review 结构化抽取

## 概述

本 skill 用于将一篇"工作流聚类综述"的 Markdown 文件抽取为结构化 JSON，遵循 `schema.json` 定义的契约。

## 输入

一篇 Markdown 文件，结构特征：
- 开头有引言，描述工作流的科学问题与整体框架
- 主体按阶段（step）分章，每章内含若干子节（将成为 substep）
- 每章末尾有：推荐工具 / 关键参数 / 常见陷阱 / 案例研究
- 末尾有结论与参考文献列表

## 输出

一个 JSON 文件，包含以下顶层键：

| 键 | 说明 | 何时为空 |
|---|---|---|
| `workflow` | 工作流总览 + 结论 + claims 数组 | 永不为空 |
| `steps` | 阶段节点，含 substeps + step 级字段 | 永不为空 |
| `tools` | 具体软件工具 | **只放原文明确报告的具名软件**，否则留空 |
| `methods` | 理论/算法/协议 | 永不为空（至少一条） |
| `data` | 公开数据集 | **只放原文明确引用的具名数据集**，否则留空 |
| `paper_refs` | 引用文献池 | 永不为空 |

## 抽取规则

### 1. Workflow

- `description`：引言中描述工作流整体流程的段落，保留原文丰富度
- `step_ids`：按章节顺序排列的 step id 列表
- `member_paper_ids`：所有聚类成员文献的 id
- `conclusions`：结论章节的自然语言总结（可保留为单段 markdown）
- `open_issues`：局限性 / 未来工作中与 workflow 直接相关的部分，去掉套话
- `claims`：将结论拆成原子声明数组，每条有独立 `id`（c1, c2...）和 `text`

### 2. Steps

- 每个 step 对应原文一个主章节
- `id` 命名：`step_<英文关键词>`，如 `step_input`、`step_thermo`
- `description`：该章开头概述段落，1–2 句话
- `substeps`：该章的子节，每个只有 `name`（子节标题）和 `description`（子节正文，保留 LaTeX 公式）
- `method_ids`：该 step 涉及的方法 id 列表
- `paper_ref_ids`：该 step 引用的文献 id 列表
- **Step 级字段**（从原文每章末尾的"推荐工具/关键参数/常见陷阱/案例研究"块提取）：
  - `recommended_tools`：保留原文 bullet points，不压缩
  - `key_parameters`：保留原文参数名 + 建议值
  - `common_pitfalls`：保留原文错误 + 解决方案对
  - `case_study`：该章的复现案例描述

### 3. Tools

**关键原则：只收录"可以点名的具体软件/工具"**，不收录通用类别。

- 正确：「VASP 6.3」「Quantum ESPRESSO」「SciPy quad」
- 错误：「数值积分器」「非线性方程求解器」
- 如果原文只说"使用数值积分器（未报告具体软件）"，则 `tools` 数组留空，工具类别信息放入对应 step 的 `recommended_tools` 字段

### 4. Methods

每个 method 是一个独立对象：
- `id`：`method_<英文关键词>`，如 `method_bcs_mft`、`method_greens`
- `description`：一句话概括核心思路
- `assumptions`：核心假设（bullet list）
- `applicability`：适用场景（bullet list）
- `relationship`：与其他 method 的关系（点名 method id）
- `paper_ref_ids`：关联文献

### 5. Data

- 只收录"可引用名称的公开数据集"，不收录物理量（N(E)、Cs(T) 这些留在 step description 中）
- 如果原文不涉及公开数据集，`data` 数组留空

### 6. Paper Refs

- `id`：`ref1` 到 `refN`，与原文参考文献编号对齐
- `citation`：作者 + 标题 + venue/年份（单行字符串）
- `description`：该文献在综述中扮演的角色和具体贡献
- `is_cluster_member`：是否属于该工作流聚类的成员文献

## 内容质量标准

1. **不压缩**：抽取的内容应与原文保持同等丰富度，公式、列表、参数细节不丢失
2. **LaTeX 公式**：使用 `$...$`（行内）和 `$$...$$`（独立行），保留原文公式编号
3. **Markdown 格式**：bullet points 用 `- `，粗体用 `**...**`
4. **引用闭环**：所有 `*_ids` 字段引用的 id 必须在对应池中存在
5. **自然语言优先**：`description` 用自然语言组织，不上演过度形式化

## 常见陷阱

- **把通用工具类别放入 `tools` 数组** → 只放入被明确命名的具体软件
- **把物理量（N(E)、ΔC）放入 `data` 数组** → 物理量留在 step description 中
- **把 step 级字段（工具/参数/坑）塞进 substep** → 这些在原文中是章节级块，不属于任何单个子节
- **结论保留大段套话** → 只保留与 workflow 直接相关的量化结果和判据
- **id 命名不一致** → 全文件统一前缀约定（`step_`、`method_`、`ref`、`c`）