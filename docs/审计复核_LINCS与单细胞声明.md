# 审计复核记录：LINCS 层与单细胞层声明

日期：2026-10-18
复核对象：`manuscript/submission_OralOncology_v1.md` 与 `manuscript/manuscript_v1.md`
方法：从原始/中间数据独立重算，不采信既有结论

---

## 一、727 这个数字的溯源

### 结论：可复现，但手稿措辞两处不准

手稿原文（§3.4）：
> A single-direction LINCS query returned 727 nominally significant compounds, whereas the two-direction reversal intersection narrowed to eight high-confidence candidates.

重跑 `PVL_05_lincs.R` 的两个单方向 Enrichr 查询（adj.P < 0.05），结果落盘于：

| 查询 | 库 | 显著 Term | 唯一 Drug（按 Term 末段解析） |
|---|---|---|---|
| progression-**up**（26 基因） | `LINCS_L1000_Chem_Pert_down` | 8 | 8 |
| progression-**down**（60 基因） | `LINCS_L1000_Chem_Pert_up` | 720 | 244 |

- 简单相加：720 + 8 = **728**
- 两查询共有的 Term：**1**（`LJP008 HEPG2 24H-GSK-1059615-1.11`）
- 并集去重：720 + 8 − 1 = **727** ← 与手稿数字精确一致

### 措辞问题（数值本身无误，描述需修正）

1. 写为 “**single-direction** query returned 727” —— 727 是**两个方向合并去重**的结果，不是单方向。单方向分别只有 8 和 720。
2. 写为 “727 ... **compounds**” —— 727 是 **Term 条目数**（化合物 × 细胞系 × 剂量 × 时长的组合）。按化合物去重后仅 **244** 个。

建议表述：
> Across the two directional LINCS queries (adj. P < 0.05), 727 compound–cell-line–dose entries corresponding to 244 unique compounds were nominally significant in the protective direction; requiring reversal in both directions narrowed this to eight compounds.

数据文件：`F:\Data Temp\PVL\LINCS_single_direction_up_geneset.csv`（3889 条）、
`F:\Data Temp\PVL\LINCS_single_direction_down_geneset.csv`（12923 条）

---

## 二、8 个双向候选可复现性

用 `PVL_05_lincs.R` 的 `build_both()` 逻辑（`Drug` 解析：Term 末段按 `-` 切分，取 `p[2:(n-1)]` 拼接）重算两方向 adj.P < 0.05 的 Drug 交集：

```
计算得到：AS-601245, AZD-7762, GSK-1059615, HG-6-64-01, PF-477736, ZSTK-474, sorafenib, torin-1
归档结果：AS-601245, AZD-7762, GSK-1059615, HG-6-64-01, PF-477736, ZSTK-474, sorafenib, torin-1
差异：在交集但不在归档 = 0；在归档但不在交集 = 0
```

**8 个候选零差异复现。**

---

## 三、"eight high-confidence candidates" 的措辞问题

`PVL5_lincs_reversal_candidates.csv` 的 `consensus_confirmed` 列 8 行全为 `FALSE`；
`phase5_run.log` 明写 `consensus-confirmed: 0`。

### 根因已查明：两重原因叠加

**原因 1（解析器 bug）**：`qry()` 的 `Drug` 提取器假定 Term 形如
`LJP006 HA1E 24H-ZSTK-474-0.12`（末段含 `-` 分隔的化合物名）。
而 `LINCS_L1000_Chem_Pert_Consensus_Sigs` 库的 Term 形如 **`Afatinib Up`**，
末段无 `-`，解析器退化为返回末段本身 → 全库 1413 个 Term 的 `Drug` 值全部变成 `'Up'`。
因此 `grepl(nn(drug), ...)` 永远匹配不上，`consensus_confirmed` 必然全 FALSE。

**原因 2（库本身覆盖不足）**：该 consensus 库共 1,413 个 Term，**全部为 `Up`，`Down` 条目为 0**。
保护性复原方向需要的是 "化合物下调 progression-up 基因"（即 Down 方向），
在这个库里**原理上不可能被确认**。8 个候选中仅 `sorafenib` 以 `Sorafenib Up` 形式存在，
其余 7 个（ZSTK-474 / PF-477736 / AZD-7762 / HG-6-64-01 / GSK-1059615 / AS-601245 / torin-1）
在该库中不存在。

### 建议
- 「high-confidence」不能保留，因为 consensus 层实际未提供任何确认；
- 改为中立表述：「eight compounds showed significant reversal in both directions」，
  并在 Limitations 说明 consensus 层因库覆盖与方向单一性未能提供独立确认。

---

## 四、手稿 §3.3 单细胞声明被数据否定

手稿原文（`submission_OralOncology_v1.md` 第 130 行 / `manuscript_v1.md` 第 108 行）：
> ...the most advanced block (bulky hyperkeratosis) showed enrichment of Treg, plasma and mast-cell states...

### 病灶块对应（GEO GSE196296 元数据 `histopathological diagnosis`）

| 别名 | 病理注释 | 细胞数 |
|---|---|---|
| LesionA | keratosis of unknown significance (KUS) | 8,815 |
| LesionB | corrugated orthohyperkeratosis | 8,814 |
| **LesionC** | **bulky hyperkeratotic epithelial proliferation**（最进展块） | 10,914 |

### 三类细胞状态的逐块占比（marker 定义，`SC_marker_state_by_lesion_stats.csv`）

| 状态 | A (KUS) | B (corrugated) | **C (bulky)** | Cramér's V | 判定 |
|---|---|---|---|---|---|
| Treg（CD3D⁺CD3E⁺PTPRC⁺ 且 FOXP3⁺ 且 IL2RA/CTLA4/TNFRSF18 任一） | 0.726% | 0.862% | **1.164%** | 0.0195 | 方向一致，但效应量可忽略 |
| plasma（MZB1⁺JCHAIN⁺） | **5.865%** | 1.917% | 2.804% | 0.0895 | **峰值在 A，方向相反** |
| mast（TPSAB1⁺TPSB2⁺CPA3⁺） | 0.363% | **1.021%** | 0.687% | 0.0312 | **峰值在 B，方向相反** |

补充：`celltype` 24 类注释中**不含** Treg / plasma / mast 任何细分
（仅有 `T_cells`、`B_cell`、`NK_cell` 等大类），故该声明原本就无法从 celltype 直接得出。

### 同段其余两句成立

| 声明 | 数据 |
|---|---|
| keratinocyte-dominant | Keratinocytes 63.07% / 55.16% / 73.15%；加上皮合计 68.93% ✅ |
| discernible T-cell infiltration | T_cells 12.65% / 21.09% / 13.66%（三块均 >12%）✅ |

### 建议
删除或重写该分句。若保留，只能写成经得起核查的形式，例如：
> Treg-associated transcripts were modestly enriched in the bulky block (1.16% vs 0.73–0.86%; Cramér's V = 0.02), whereas plasma-cell and mast-cell signals peaked in the KUS and corrugated blocks respectively.

数据文件：`results/manuscript_figs/SC_celltype_by_lesion_{counts,props,stats}.csv`、
`results/manuscript_figs/SC_marker_detection_by_lesion.csv`、
`results/manuscript_figs/SC_marker_state_by_lesion_stats.csv`

---

## 五、附：发现的历史产物问题

1. **`PVL5_lincs_promoting.csv`（192 行）是旧版单方向产物**，非现行脚本输出。
   其列结构为 `Term, Drug, reversal_score, combined_score, n_dir, overlap_genes`（`n_dir` 全为 1），
   现行脚本写的是 `Drug, scoreA, termsA, ..., consensus_confirmed`。
   192 恰等于 `down × Chem_Pert_down` 单方向 adj.P < 0.05 的 Term 数 → 应标注为 legacy 单方向结果。

2. **Enrichr API v3 已不返回 `Odds.Ratio` 字段**（现为 9 元组），
   重跑结果的该列为空，已在表中如实留空而非填充伪值。

---

## 六、本轮修复的数据产物缺陷（提交 ec2ac69）

| 缺陷 | 修复 | 验证 |
|---|---|---|
| `signature_101genes.csv` 的 `signif` 列沿用 \|logFC\|>1.0 阈值（仅 13 行 yes） | 按文档阈值 `adj.P<0.05 & \|logFC\|>0.5` 重算 | 101 行全为 yes |
| `PVL6_summary.csv` 的 `KW_p` 数值下溢写为 0 | 写为 `1.127e-11`，并补 `KW_chisq=50.418`、`KW_df=2` | 读取正确 |
| `Fig4_second_cohort.png` 为 700px/120dpi 插值放大 | 从 `PVL9_GSE85195_scores.csv` 以 matplotlib 重绘 | 2274×1092 @300dpi；AUC 复算 0.949020 与归档一致 |
| 缺 83/101 映射基因名单 | 按 `PVL_09` 同逻辑重建 `PVL9_mapped_genes.csv` | 83 行，与归档 `sig_genes_mapped=83` 一致 |
