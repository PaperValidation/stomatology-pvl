# PVL Phase 9 — 第二独立队列（GSE85195）横断面佐证

## 目的
选项 B：用101基因有向 z-score 在第二个**独立公开表达队列**上做交叉验证，作为对 *Oral Oncology* 评审意见"缺乏独立恶变验证队列"的回应。本质是**横断面判别**（癌前 vs 癌），**不是前瞻预测**（谁会转化）。

## 队列
- **GSE85195**（Agilent Whole Human Genome 4×44K, GPL6480）：来自龈颊复合体（印度人群）。
  - 15 OPL（oral leukoplakia，癌前）/ 34 早期 OSCC / **1 normal**（n=1，无法做 normal 组比较）。
  - 来源文献：Bhosale PG et al. *Transl Oncol.* 2017;10(3):396-409. doi:10.1016/j.tranon.2017.03.008（已加入参考文献 [23]）。

## 方法（与 PVL_04 完全一致的方向加权 z-score）
- 权重 `w` = GSE26549 进展 DEG 的 logFC（**固定**，仅用发现队列权重，不在本队列重估、不跨平台搬系数）。
- 探针→基因映射：GPL6480.annot（经 REST SOFT 解析为 `GPL6480_probe2gene.tsv`，30,936 带符号）。
- 评分：`score = safe_scale(expr) %*% w`，队列内逐基因 z-score 去除平台量纲。
- 分组来自 `!Sample_title`：OSCC / OPL / Normal。

## 结果
| 指标 | 值 |
|---|---|
| 签名基因映射到 GSE85195 | 83 / 101（82.2%） |
| OPL(n=15) vs OSCC(n=34) 判别 AUC | **0.949**（95% CI 0.891–1.0） |
| 严重度序数（OPL<OSCC）Spearman ρ | **0.717**（P = 6.8×10⁻⁹） |
| 平均评分 OPL / OSCC | −40.8 / +19.5 |

## 诚实定位（关键）
1. **横断面 state 判别，非前瞻预测**：GSE85195 的 OPL 是现患白斑、无转化随访，故只能说明签名区分"癌前 vs 癌"两种状态、并沿严重度轴上升；不能证明"预测个体会转化"。
2. **强独立性证据**：不同平台（Agilent vs 发现队列 Affymetrix）、不同人群/解剖部位（龈颊复合体 vs 原队列）→ 说明签名不是 RNA-seq 或某一人群特有的假象。
3. **局限**：n 较小（OPL 15 / OSCC 34），95% CI 上界触 1.0（精度有限）；normal n=1，未做 normal 组比较；AUC 0.949 高但属小样本横断面，需谨慎解读。
4. 与 GSE227919 互补：GSE227919 给的是"严重度梯度相关"（ρ=0.47），GSE85195 给的是"跨平台跨人群区分癌前/癌状态的判别力"，二者共同把签名从"单 RNA-seq 队列"提升到"双独立队列佐证"，但**前瞻转化预测队列仍是空白**。

## 修复的 bug
初版脚本误用 DEG csv 中**全部基因**（16,147 条）作权重而非限定到 101 基因签名，导致映射率 15987%、结果无效。已改为 `w` 仅取 `sig`（adj.P<0.05 & |logFC|>0.5，101 基因），重跑得上述正确结果。

## 产出文件
- `PVL_09_second_cohort.R`（WSL 脚本，仅读本地文件，无需联网）
- `raw/GSE85195/GSE85195_series_matrix.txt(.gz)`、`raw/GPL6480.soft.txt`、`raw/GPL6480_probe2gene.tsv`
- `results/manuscript_figs/PVL9_ROC_OPL_vs_OSCC.png`（Fig 4A）
- `results/manuscript_figs/PVL9_boxplot_score.png`（Fig 4B）
- `results/manuscript_figs/PVL9_summary.csv`、`PVL9_GSE85195_scores.csv`

## 手稿集成（manuscript_v1.md）
- 摘要：新增 GSE85195 AUC 0.95 / ρ 0.72 一句。
- §2.1 数据来源：新增 GSE85195。
- §2.3：说明对 GSE85195 套用相同盲法、属横断面判别。
- §3.7 新增"Independent cross-sectional corroboration"小节 + Fig 4A/B。
- 表 3：第二队列佐证汇总。
- §4.1 讨论、§4.5 局限：将"独立证据仅 GSE227919"更正为"GSE227919 + GSE85195，但皆横断面"。
- 数据可用性：加入 GSE85195。
- 图注 Fig 4、参考文献 [23]。
