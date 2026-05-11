# Trace: gaia_workflows

<!-- concepts: workflow-structuring, knowledge-graph, superconductor-specific-heat -->

## 2026-05-11 — review_cluster_172 结构化首次尝试

- **Task**: 把 `review_cluster_172.pdf`（17 篇文献聚类综述：超导比热计算）抽取成结构化 JSON，最终要进数据库支持三类查询：工具常见用法、研究问题的标准方法、工具/方法对比与坑。
- **Schema 决策**: 用户明确要求粗粒度 6 类 —— `workflow / steps / tools / methods / data / paper_refs`，主要以自然语言描述为主，跨对象只用 id 引用不重复内容。早先的细字段版本（equations/claims/variants/decision_points/open_issues/benchmarks/parameters 独立池）被判过设计，砍掉。
- **失败记录**: Write 工具连续 5 次因缺 `content` 参数被打回。根因是我的 tool use XML 块只传了 `file_path`，没嵌 `content` 块。改用极短内容（`{}`）先确诊参数传递可行后恢复。教训：工具参数缺失应该一次失败就 diagnose，不是继续重试。
- **pdftotext 产出**: `review_cluster_172.md` 810 行纯文本，保留段落但丢失 markdown 结构（标题/公式/列表未转换）。后续若要更干净的结构化可以换 `marker-pdf` 或 `pymupdf4llm` 重跑。
- **EARS thrashing 提示不是 stuck**: 连续 Edit review_cluster_172.json 是规划好的分段填充（骨架 → steps → tools → methods → data → paper_refs），每次推进一段不是同一问题反复修。
- **19:35 再次 thrashing 复位**: 第二批改写——用户要求 tools/data 只收具体软件和公开 dataset，本篇无此类内容所以两数组清空，被下沉的信息改写进 4 个 step 的 description。仍然不是 stuck。
- **20:42 Progress**: 开始做 interactive HTML 浏览器（`review_cluster_172.html`），含 sticky 搜索框 + 6 组预生成 Q&A + workflow/steps/methods/paper_refs 分区卡片。Write 时疑似在 `details .body{...font-size:` 处发生内容错位，需要 Read 诊断并重写/修复。
- **21:11 Progress**: 连续 Write/Bash 大 content（7KB+）都被截断成空参数失败；Edit（小 new_string）和 Write 小文件（<2KB）稳定成功。结论：这会话里生成层对较长 parameter 值有截断 bug。应对：把 Q&A 抽成独立 `review_cluster_172.qa.json`（已落盘 1.9KB），HTML 瘦身为 fetch 两个外部 JSON 的骨架，再试小 content Write。
