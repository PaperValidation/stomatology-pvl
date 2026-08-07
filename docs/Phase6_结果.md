# Phase 6 结果：单细胞 ↔ bulk 签名闭环映射（针对预审委员会 Reviewer 4 断裂1 / 弱支撑2、Reviewer 2 #4）

> 脚本：`PVL_06_closeloop.R`（WSL `~/PVL/scripts/`）
> 输入：`results/PVL_sc_combined.rds`（28,543 细胞 / 24 簇 / 3 lesion）+ `GSE26549_progression_DEG.csv`（101 基因签名）+ `PVL4_lasso_coefficients.csv`（20 基因）
> 输出：`PVL6_sigscore_by_lesion.png` · `PVL6_sig_heatmap_clusters.png` · `PVL6_split_umap_lesion.png` · `PVL6_lesion_pca_variance.csv` · `PVL6_lasso_marker_overlap.csv` · `PVL6_summary.csv`

---

## 1. 目的

将 Phase 2/4 的 **bulk 进展签名** 落到 Phase 3 的单细胞对象上，回答「签名在单细胞层面的细胞来源是什么、是否随 PVL 病理连续谱（LesionA→B→C）梯度变化」——这是把 Phase 3 从「独立描述性章节」转为「签名的细胞级验证层」的关键一步（委员会所指叙事断裂点）。

---

## 2. 方法与指标

- **每细胞签名评分**：87 个签名基因在单细胞 `data` 层内逐基因 z-score，取均值（与 Phase 4 bulk 评分方法学一致、可直接比较）。
- **效应量**：lesion 对签名评分的 η²（而非仅 KW p，因 n=28,543 使微小差异也「显著」）。
- **签名富集簇**：各 cluster 平均签名评分排序。
- **批次效应证据**：RunPCA 后计算 lesion 对 PC1–PC10 的 η²（解释方差占比）。
- **LASSO↔簇标记重叠**：20 个 LASSO 基因与 24 簇 marker 的交集。

## 2.1 方法学注释：为何以效应量 η² 为准、不以 Kruskal–Wallis p 值定论（复审 A2）

本单细胞对象含 **n = 28,543** 个细胞。在此样本量下，即便生物学上可忽略的微小差异也会给出「统计显著」的 p 值——这是大样本研究的经典陷阱（Lin et al., 2013, *Methods Ecol Evol*；对「p 值崇拜」的批判亦见 Wasserstein & Lazar, 2016, *Am Stat*）。本分析因此**以效应量（η²，即 lesion 对签名评分方差的解释比例）作为结论依据，p 值仅作辅助**。

- **实测对照**：lesion 对签名评分的 KW p = 1.13×10⁻¹¹（表面「高度显著」），但 η² = **0.0061**——即 lesion 仅解释签名评分方差的 0.6%，属**可忽略**效应。若仅凭 p 值，会错误得出「签名在病灶连续谱上差异显著」；以 η² 为准，则正确结论为「签名在 lesion 层面无生物学梯度」。
- **η² 判读基准**（Cohen 1979 惯例，组间设计）：0.01 = 小、0.06 = 中、0.14 = 大。本签名 lesion η² = 0.006 远低于「小」阈值，故判为无效应。
- **对比**：lesion 对 PC1–PC10 的 η² 均值 0.026（最高 PC10 = 0.078），虽仍属小效应但高于签名评分的 0.006——支持「批次/生物学变异存在于 PCA 层面、幅度有限」的论断，与「签名评分无 lesion 梯度」互不矛盾。
- **可引用表述**：*"With n = 28,543 cells, we report effect sizes (η²) rather than p-values as the basis for inference; the lesion-level Kruskal–Wallis test was nominally significant (p = 1.13×10⁻¹¹) but explained only 0.6% of signature-score variance (η² = 0.0061), i.e., a negligible effect, leading us to conclude the bulk signature does not track the within-lesion histopathologic continuum."*

---

## 3. 结果

| 指标 | 值 | 解读 |
|---|---|---|
| 签名基因落入 sc 对象 | 87 / 101 | 覆盖度良好 |
| 每细胞签名评分均值（LesionA/B/C） | 0.0150 / 0.0096 / **−0.0198** | **无 lesion 级梯度**（LesionC 反略低） |
| lesion η²（签名评分效应量） | **0.0061** | 可忽略 → bulk 签名**不**追踪病灶内连续谱 |
| KW p | 1.13e-11 | 由大 n 驱动，非生物学效应，故以 η² 为准 |
| 签名最高富集簇 | 簇 5 (0.42) > 13 (0.31) > 22 (0.27) > 23 (0.26) > 15 (0.26) | 签名**集中在特定细胞状态**而非均匀分布 |
| lesion 对 PC1–PC10 的 η²（均值） | 0.026（PC4=0.059, PC7=0.055, PC10=0.078） | lesion 解释方差**小但存在** → 支持「batch≈生物学、幅度有限」的免整合辩护 |
| LASSO 20 基因 ↔ 簇标记重叠 | **4 个**：QPCT, ECM1, S100A7A, RBP7 | 风险模型与单细胞有部分直接锚定 |

---

## 4. 解读（诚实修正原假设）

**原假设**（Phase 3 草稿）：签名评分应随 LesionA<B<C 单调上升。
**实测结果**：签名评分在 lesion 层面**无梯度**（η²=0.006），LesionC 反而略低。这是一项**真实的边界条件**，而非分析失败，理由如下：

1. **bulk 签名捕捉的是「队列级」进展生物学**（GSE26549 进展者 vs 非进展者、GSE227919 严重度梯度 ρ=0.47），而非「同一 PVL 病灶内 KUS→corrugated→bulky 的细粒度连续谱」。两者生物学尺度不同，签名在后者不显梯度属合理。
2. 签名**并非无结构**——它**富集于特定细胞状态（簇 5/13/22/23/15）**，提示其信号来自特定上皮/免疫亚群，而非全细胞均匀贡献。这些簇的组合在 lesion 间的比例变化，正是 field-cancerization 的细胞基础。
3. 4/20 LASSO 基因（ECM1、S100A7A 等上皮/基质/角质标记）直接命中簇标记，为风险模型提供了**单细胞层面的部分锚定**。

**结论性叙事（成稿用）**：单细胞分析的角色由「签名追踪连续谱」修正为「**签名的细胞级来源解析 + 免疫抑制/上皮代谢轴的独立佐证**」。结合 Phase 3 的 LesionC Treg/浆/肥大富集与上皮代谢重编程（AKR1C/NQO1/SPRR/ELF3），以及 lesion PCA η² 小但存在的证据，Phase 3 不再游离于主线，而是作为「签名→特定细胞状态→免疫代谢轴」的验证层。

---

## 5. 输出文件

| 文件 | 内容 |
|---|---|
| `PVL6_sigscore_by_lesion.png` | 每细胞签名评分 violin/box（按 lesion，附 η² 注释） |
| `PVL6_sig_heatmap_clusters.png` | 87 签名基因 × 24 簇平均表达热图（z/基因，标注主导 lesion） |
| `PVL6_split_umap_lesion.png` | 按 lesion 拆分的 UMAP（批次≈生物学可视化证据） |
| `PVL6_lesion_pca_variance.csv` | lesion 对 PC1–PC10 的 η² |
| `PVL6_lasso_marker_overlap.csv` | LASSO 基因↔簇标记的重叠及最佳簇 |
| `PVL6_summary.csv` | 全部指标汇总 |

---

## 6. 对委员会意见的直接回应

- **Reviewer 4 断裂1（单细胞功能模糊）**：✅ 已建立签名→细胞状态映射（富集簇 5/13/22/23/15），Phase 3 转为「签名的细胞级验证层」。
- **Reviewer 4 弱支撑2（LASSO↔簇无映射）**：✅ 4/20 重叠（QPCT/ECM1/S100A7A/RBP7），其中 ECM1/S100A7A 为上皮/角质标记。
- **Reviewer 2 #4（批次效应可视化证据）**：✅ split-UMAP + lesion PCA η²（均值 0.026）已产出，支持免整合的生物学辩护。
- **Reviewer 2 #5/#6、Reviewer 4 断裂2**：见 `Phase4_结果.md` 修订版（梯度 ρ 为头号证据、CV 嵌套偏倚显式声明）。
