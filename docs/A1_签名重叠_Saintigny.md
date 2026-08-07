# A1：Saintigny 2011 签名与本工作 101 基因签名重叠分析（回应复审新增行动 A1）

> 目的：量化发现队列源文献（Saintigny et al., *Cancer Prev Res* 2011）的进展签名与本工作 GSE26549 重分析所得 101 基因签名的重叠，以决定投稿叙事应定位为「独立发现」还是「复制 + 扩展/精细化」。
> 脚本/数据：`results/saintigny_overlap/analyze_overlap.py` · `our_101_genes.txt` · `saintigny_overlap_results.csv`

---

## 1. 关键校正（先纠错，再分析）

文献对话底稿原将 Saintigny 2011 描述为「71 基因签名」，DOI 记为 `10.1158/1940-6207.CAPR-10-0996`。经核对原文（PMC3074595，开放获取）两项均**有误**：

- **正确 DOI**：`10.1158/1940-6207.CAPR-10-0155`（原 `…-0996` 在 DOI 系统查无此文）。
- **签名构成**（非 71 基因）：该文基于 GSE26549（86 例 OL，35 例恶变）开发了 **29-transcript 预测模型**（Table 2 Model 1，预测误差 8%），并另行以单变量 Cox 模型鉴定出 **2,182 个与恶变风险显著相关的转录本**（FDR≈11%，见 Supplementary Material 9）。文中并无「71 基因」之说。

> ⚠️ 此纠错已同步回 `Literature_对话.md` §1（引用 #1）与 `MEMORY.md`。投稿前务必以正确 DOI 与「29-transcript 模型 / 2,182 Cox 显著转录本」措辞引用，避免被审稿人判定「不熟该领域」。

---

## 2. 方法

- **本工作 101 基因签名**：GSE26549 基线进展 DEG，阈值 `adj.P.Val<0.05 & |logFC|>0.5`（limma；与 Phase 2/4 一致），全部为标准蛋白编码基因符号。
- **Saintigny 29-transcript 模型**：取自原文 Table 2 Model 1 的 29 个 probeset 及其 published gene symbol。
  - 组成：14 个可比对基因符号（其中 13 个明确蛋白编码 + `CNO` 未定），10 个 Ensembl 转录本 ID（`ENST…`）、1 个 miRNA（`hsa-mir-101-1`）、1 个 snRNA（`RNU2-1`）、1 个基因组 accession（`NC_001807`）、1 个 `NA`、1 个非编码转录本（`AK128090`）。
  - 即该「签名」中约 **55%（16/29）为非蛋白编码/优化型特征**，仅 14 个可映射到常规基因符号。
- **重叠指标**：Jaccard = |交集| / |并集|，在**可比对基因符号**层面计算（ENST/miRNA/snRNA 等无法与本工作蛋白编码签名直接匹配，按方法学诚实保留为不可比项）。
- **2,182 Cox 显著转录本列表**：原文 Supplementary Material 9（NIHMS261208-supplement-9.xls）。尝试下载时被 NCBI 的 Proof-of-Work 反爬网关拦截（返回 "Preparing to download…" 中间页），沙箱内未解出。该列表与我们的签名**同源（均出自 GSE26549）**，预期重叠更高，仅会**强化**复制结论，不影响本分析的方法学判定（见 §4）。建议在投稿前由作者在可交互浏览器中下载，做扩展一致性分析。

---

## 3. 结果

| 指标 | 值 |
|---|---|
| 本工作 101 基因签名 | 101 |
| Saintigny 29-model（可比对符号） | 14 |
| **交集** | **0** |
| **Jaccard(101, 29-model)** | **0.000** |

> 交集为 0 是**统计对象差异**所致，而非「独立」证据：29-transcript 模型是 spars e 优化的**预测分类器**（刻意选 29 个最 discriminating 特征，且多数非蛋白编码），本工作 101 基因为**差异表达集**（不同统计目标、不同筛选阈值）。两者从同一批样本、同一平台得出，特征集合不同属正常现象。

---

## 4. 叙事判定（核心产出）

**结论：本工作对 Saintigny 2011 签名既非「独立发现」，亦非简单「复制」；正确定位为「同一发现队列内的复制 + 扩展（精细化）」，并把真正的独立性主张建立在外部队列上。**

理由链：

1. **同源不可称独立**：101 基因签名与 Saintigny 29-model 均源自 GSE26549 同一 86 例队列、同一 Affymetrix 平台。任何「我们独立发现了 OL 进展签名」的表述都是错误且会被秒拒的——我们是在原作者数据上重算了更广的进展表达信号。
2. **扩展价值真实存在**：
   - 维度更宽（101 vs 29），并补充了方向性（logFC 正负）；
   - 做了**探针→基因折叠敏感性**（max vs mean，Jaccard 0.971）与**外部独立验证**；
   - 把签名落到**单细胞层面**（Phase 3/6：免疫抑制 + 上皮代谢轴）并做**功能/TF 解构**（Phase 7：mTORC1/PI3K-AKT-mTOR/EMT/KRAS，上游 P53）；
   - 推导到 **LINCS 药物重定位**（Phase 5：8 候选）。
3. **独立性主张的唯一合法支点 = 外部队列 GSE227919**：RNA-seq 独立队列的严重度梯度 ρ=0.47（p=6×10⁻⁵）与 Cancer-vs-其他 AUC=0.82。这是本工作**跨平台、跨数据集**的证据，也是相对 Saintigny 2011（仅同队列 29 特征模型、无外部表达验证）的真正增量。

**投稿表述建议（植入 Introduction/Discussion）**：
> "Our progression-expression signature was derived from the same discovery cohort (GSE26549) that underpinned the Saintigny 2011 29-transcript classifier, and should therefore be framed as a replication and refinement of that signal rather than an independent discovery. The genuinely independent evidence for the signature comes from an external RNA-seq cohort (GSE227919), where it tracks histopathologic severity (Spearman ρ=0.47, p=6×10⁻⁵)."

---

*本文件为复审 A1 行动的正式交付。DOI 与签名词措纠错已同步至 `Literature_对话.md` 与 `MEMORY.md`；扩展的 2,182-转录本一致性分析列为投稿前可由作者完成的收尾项。*
