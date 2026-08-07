# Phase 3 结果 — 单细胞 RNA-seq（GSE196296, 10x Genomics）

> 高危口腔白斑 / PVL 恶变-复发一体化风险转录组学项目
> 数据：GSE196296（3 个 10x 样本 = 同一 PVL 病灶内病理连续谱）
> - LesionA = KUS（角化不全鳞状增生）
> - LesionB = corrugated（皱褶状）
> - LesionC = bulky hyperkeratosis（疣状肥厚角化）
> 分析环境：WSL2 Ubuntu-24.04 / R 4.3.3 / Seurat v5 / SingleR + celldex（Human Primary Cell Atlas）

---

## 1. 方法

| 步骤 | 参数 / 说明 |
|---|---|
| 读取 | `readMM` 解析 matrix.mtx.gz + barcodes + features（基因 SYMBOL，重复名 `make.unique`）|
| 合并 | 3 样本 `merge`（按 lesion 加 cell-barcode 前缀）|
| 质控 QC | `nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 25` |
| 标准化 | `NormalizeData` → `FindVariableFeatures(3000)` → `ScaleData` |
| 降维 | `RunPCA(50)` → `RunUMAP(dims=1:30)` |
| 聚类 | `FindNeighbors(dims=1:30)` → `FindClusters(resolution=0.6)` |
| 注释 | SingleR（celldex `HumanPrimaryCellAtlasData`，`label.main`，经 SOCKS5 代理从 ExperimentHub 缓存取得）|
| 标记基因 | `FindAllMarkers`（仅阳性、`min.pct=0.1`、`logfc.threshold=0.25`，**限定 variable features + `max.cells.per.ident=200` 加速**）|

**关键工程修复（Seurat v5 / WSL 环境）：**
- 合并后对象默认生成**按样本拆分的多层**（`data.LesionA/B/C`），`GetAssayData(slot=)` 在 v5 已废弃且对多层报错 → 用 `JoinLayers()` 合并为单一 `data` 层，并以 `LayerData(combined, layer="data")` 取稀疏矩阵传给 SingleR。
- 后台任务须由 Bash 工具 `run_in_background` 托管（WSL 内 `nohup &` 会被包装进程退出时回收）。
- 单细胞进程与另一项目（Dermatology-Phyllodes-tumor 的 Bioconductor 安装）同时跑会抢内存；标记发现步骤用 variable-feature 限定 + 下采样以避免长时间 Wilcoxon 循环。

---

## 2. 结果

### 2.1 规模
- 原始：LesionA 9,200 / LesionB 9,381 / LesionC 11,738 cells × 33,538 features
- QC 后：**28,543 cells**（LesionA 8,815 / LesionB 8,814 / LesionC 10,914）
- 聚类：**24 个 clusters**

### 2.2 整体细胞类型组成（SingleR, Primary Cell Atlas）

| 细胞类型 | 数量 | 占比 | 注释 |
|---|---:|---:|---|
| Keratinocytes 角质形成细胞 | 18,387 | 64.4% | 病灶主体（鳞状上皮）|
| T_cells T 细胞 | 4,460 | 15.6% | 主要免疫浸润 |
| Tissue_stem_cells 组织干细胞 | 1,538 | 5.4% | 多为基底/增殖型角质形成细胞（参考库无口腔特异性，属 SingleR 误注）|
| Epithelial_cells 上皮细胞 | 1,267 | 4.4% | |
| Chondrocytes 软骨样 | 552 | 1.9% | 参考库误注（实为部分鳞状/间质）|
| Fibroblasts 成纤维细胞 | 466 | 1.6% | |
| B_cell B 细胞 | 369 | 1.3% | |
| Endothelial_cells 内皮 | 362 | 1.3% | |
| Monocyte 单核细胞 | 347 | 1.2% | |
| Smooth_muscle_cells 平滑肌 | 242 | 0.8% | 肌成纤维/血管周 |
| NK_cell NK 细胞 | 230 | 0.8% | |
| DC 树突状 | 106 | 0.4% | |
| Macrophage 巨噬 | 128 | 0.4% | |
| Neutrophils 中性粒 | 18 | 0.06% | |

**解读**：典型口腔鳞状上皮病灶 + 明显慢性炎细胞浸润（T 细胞为主，伴 B/浆、NK、髓系），间质成分较少。SingleR 的 `Tissue_stem_cells`/`Chondrocytes` 标签在口腔角质形成细胞背景下可靠性低，**以 cluster 标记基因驱动的注释为准**（见 2.4）。

### 2.3 聚类 × 病变组成（病理连续谱信号）

> 因 LesionC 细胞数更多（10,914 vs A/B ≈8,815），下表用**簇内各 lesion 占比**判读，避免抽样量偏差。阈值：C≥60% 视为 LesionC 偏向；C≤10% 视为 A/B 偏向；中间为均衡。

| Cluster | 主导 lesion（簇内占比）| 代表标记 | 推定身份 |
|---|---|---|---|
| 0 | A/B（C 0.3%）| SPINK6, CKB | 角质形成细胞（分化）|
| 1 | **LesionC（99%）** | AKR1C2, AKR1C3, NQO1 | 上皮代谢/应激态 |
| 2 | A/B/C 均衡 | KLRD1, CTSW, ZNF683 | NK 细胞 |
| 3 | A/B（C 1%）| IL26, TMIGD2, CCR10 | 组织驻留 T（Th17 样）|
| 4 | A/B（C 0.7%）| COL17A1, IL20, LAMB4 | 基底角质形成细胞 |
| 5 | A/B（C 12%）| STMN2, ANGPTL7, ADH1B | 混合（上皮/神经样）|
| 6 | **LesionC（90%）** | FOXP3, LAIR2, RTKN2 | **Treg** |
| 7 | 均衡（C 29%）| IL22, IL17A, ICOS | **Th17** |
| 8 | **LesionC（99%）** | SOST, CAMK2B | SOST⁺ 细胞 |
| 9 | LesionC 富集（C 67%）| SELE, APLNR | 内皮 |
| 10 | A/B（C 3%）| IGHM, IGLC3, CLDN3 | B 细胞 |
| 11 | A/B（C 0.5%）| DSC1, HEPHL1 | 角质形成细胞（桥粒）|
| 12 | **LesionC（99%）** | ELF3, SCD | 上皮（转录因子驱动）|
| 13 | LesionC 富集（C 52%）| POSTN, TNN, F2RL2 | 成纤维/间质 |
| 14 | A/B（C 0.1%）| KRT75, SPINK6, SPINK7 | 分化角质形成细胞 |
| 15 | A/B（C 1%）| MMP10, MMP3, MEDAG | 上皮/间质重塑 |
| 16 | 均衡（C 43%）| CSF1R, LILRB5, FPR3 | 巨噬/髓系 |
| 17 | A/B（C 24%）| MYH11, GPR20 | 肌成纤维/平滑肌 |
| 18 | **LesionC（91%）** | IGLL5, IGHG2 | **浆细胞** |
| 19 | B/C（C 32%）| TPSAB1, TPSB2 | **肥大细胞** |
| 20 | **LesionC（99%）** | SPRR2B, SPRR3, FAM3D | 分化角质形成（cornified envelope）|
| 21 | LesionC 富集（C 68%）| MS4A1, BLK | B 细胞 |
| 22 | 均衡 | CCL21 | 淋巴管内皮/间质 |
| 23 | LesionC 富集（C 53%）| NRXN1, PMP2 | 神经/胶质样 |

**连续谱结论**：KUS（A）与 corrugated（B）共享一套上皮–免疫程序（clusters 0,3,4,10,11,14,15 几乎不含 LesionC）；而 bulky hyperkeratosis（C）拥有**一组独特的细胞状态**——上皮代谢/分化重编程（1,12,20）、**Treg 富集（6）**、**浆细胞浸润（18）**、SOST⁺（8）、内皮（9）、成纤维（13）、B（21）。与“早期病变更均一、晚期疣状肥厚伴随免疫重塑（Treg/浆/肥大）”的 PVL 恶变前场癌变（field cancerization）假说一致。

### 2.4 UMAP 可视化
- `sc_umap.png` / `sc_umap.pdf`：按 lesion / cluster / SingleR celltype 三视图。
- `PVL_sc_combined.rds`：完整注释对象（含 `lesion`、`seurat_clusters`、`celltype` 元数据），供 Phase 4 pseudobulk 复用。

---

## 3. 生物学解读

1. **病灶以角质形成细胞为绝对主体**（≈70% 含上皮祖/上皮），符合口腔白斑/PVL 的鳞状本质。
2. **慢性炎微环境显著**：T 细胞（含 Th17 cluster 3/7、Treg cluster 6）、B/浆、NK、髓系共同存在——炎症被公认是 PVL 恶变的重要驱动。
3. **病理连续谱伴随免疫重塑**：最晚期 LesionC 相对富集 Treg + 浆细胞 + 肥大细胞 + 特定上皮分化程序，提示晚期进入“免疫抑制/特化炎症”状态，与恶变风险升高方向一致。
4. **上皮代谢重编程**：AKR1C2/3、NQO1、SCD、ELF3 在 LesionC 上皮簇高表达，指向脂质/氧化还原代谢与增殖重编程，可作为 Phase 4 风险签名候选。

---

## 4. 方法学局限与注意

1. **未做批次整合**：3 个独立 10x 文库可能存在文库特异性批次效应，部分“LesionC 特异”簇或含技术成分。因三样本本就是**同一病灶的不同病理阶段**，batch≈biology，直接整合可能抹去要研究的梯度；后续应采用**阶段感知（stage-aware）整合**（SCTransform + 限定 stage 的 anchors）以分离生物信号与批次。
2. **抽样不均衡**：LesionC 细胞数更多，绝对计数会高估 C；本报告已用簇内占比校正，但 pseudobulk 差异分析需按 lesion 均衡抽样或加权。
3. **SingleR 为粗注释**：Primary Cell Atlas 无口腔特异性参考，`Tissue_stem_cells`/`Chondrocytes` 等标签不可靠；**cluster 标记驱动注释为金标准**。
4. **标记发现限定于 variable features + 下采样**：足以定义 cluster 身份，但若需全转录本 marker 需放开（更慢、更占内存）。
5. **尚缺空间验证**：连续谱结论最好用空间转录组/多重免疫荧光在原位确认 C 区域确为晚期病变区。

---

## 5. 与 Phase 2 衔接 / 下一步

- **Phase 4（pseudobulk + LASSO 风险模型）**：以 `PVL_sc_combined.rds` 按 lesion/cluster 做 pseudobulk，与 Phase 2 的 GSE26549 进展签名（建议宽松阈值 `adj.P<0.05 & |logFC|>0.5` ≈ 101 基因）交叉，构建 LASSO 风险评分；注意独立外部验证 + 校准 + DCA。
- **Phase 5（LINCS 药物重定位）**：对 Phase 2/3 候选基因做转录组匹配化合物预测，写作 *transcriptome-matched candidate compounds*（非治疗推荐）。
- **Phase 3b（可选）**：生成 lesion × celltype 堆叠条形图（需从 RDS 取元数据，避免与竞争进程抢内存时再做）；阶段感知整合以强化梯度结论。

---

## 6. 输出文件清单（`results/`）

| 文件 | 说明 |
|---|---|
| `PVL_sc_combined.rds` | 完整 Seurat 对象（注释后）|
| `sc_umap.png` / `sc_umap.pdf` | UMAP 三视图 |
| `sc_celltype_counts.csv` | SingleR 细胞类型计数 |
| `sc_cluster_lesion_composition.csv` | 聚类 × 病变 列联表 |
| `sc_cluster_markers.csv` | 24 簇标记基因（6,776 行）|
| `phase3_run.log` | 运行日志 |
