# PVL 项目 Phase 2 结果小结（bulk DEG + 通路富集）

> 运行环境：WSL Ubuntu-24.04 / R 4.3.3 / 经 SOCKS5 代理联网。
> 脚本：`PVL_02_bulk.R`（WSL 内 `~/PVL/scripts/`）。结果 CSV/PNG 已同步至 `results/`。

## 一、管线验证（关键结论）
GSE227919 **Cancer vs Control** 的 MSigDB Hallmark 富集呈教科书级癌症信号，证明整条流程生物学合理：

| Hallmark（上调） | 校正 p | Combined Score |
|---|---|---|
| E2F Targets | 2.6e-45 | 905 |
| G2-M Checkpoint | 7.2e-43 | 800 |
| Interferon Gamma Response | 1.2e-24 | 288 |
| Interferon Alpha Response | 4.4e-23 | 466 |
| Myc Targets V1 | 1.8e-21 | 225 |
| Mitotic Spindle / EMT / mTORC1 / TNFα-NFκB | — | — |

GO BP 上调顶端为 DNA 复制、有丝分裂 G2/M 转换、MCM 家族——典型增殖程序。与"正常→低危→高危→癌变"梯度方向一致。

## 二、各队列 DEG 结果
| 队列 / 对比 | 样本 | 上调 DEG | 下调 DEG | 备注 |
|---|---|---|---|---|
| **GSE26549** 进展 vs 未进展（基线活检） | 62（21 prog / 41 non） | 6（\|logFC\|>1） | 7 | 严格阈值仅 13；adj.P<0.05 共 454，\|logFC\|>0.5 共 101 |
| **GSE227919** Cancer vs Control | 66（9 Cancer / 19 Control） | 2986 | 368 | 强信号 |
| **GSE227919** 进展 vs 稳定 | 27（10 prog / 17 stable） | 0 | 0 | 严格阈值下无显著 DEG（该梯度队列进展对比动力不足，记为阴性/弱信号） |
| **GSE227919** 梯度趋势（Control<HkNR<Dysplasia<Cancer） | 66 | — | — | 5007 个基因呈显著单调趋势 |

GSE26549 进展 top 基因（按 adj.P）：FMO2（↓）、TIMP2（↓）、ZNF667-AS1（↓）、RNASE10（↑）等。

## 三、富集产出（每对比 up/down 各一份）
- `enrich_GO_Biological_Process_2021_*.csv/.png`
- `enrich_KEGG_2021_Human_*.csv/.png`
- `enrich_MSigDB_Hallmark_2020_*.csv/.png`
- DEG 主表：`GSE26549_progression_DEG.csv`、`GSE227919_Cancer_vs_Control_DEG.csv`、`GSE227919_progression_DEG.csv`、`GSE227919_gradient_trend_DEG.csv`
- 火山图：`*_volcano.png`

## 四、方法学重要说明（务必记录）
1. **富集改用 `enrichR`（Enrichr API）**：本机 R 4.3.3 + GCC 13 下 `fgsea` 因 Boost 头 `boost::math::pow` 冲突无法源码编译；且 Posit PPM 仅提供 CRAN 二进制、**不提供 Bioconductor 二进制**，故 clusterProfiler/DOSE/fgsea 族暂时装不上。enrichR 为 CRAN 二进制（已装），结果可引用（Kuleshov et al., Nat Methods 2016）。后续若升级 R 到 4.4（Bioc 3.19，CellChat 也可用）可切回 clusterProfiler。
2. **`getGPL` 已在新版 GEOquery 移除**：GSE26549 探针→基因符号映射改用本地 `GPL6244.annot.gz`（已下载至 raw/），在 R 内解析。
3. **GEO 矩阵表头样本 ID 带双引号**：`parse_gse_matrix` 已 `gsub` 去引号，否则与临床表 GSM 无法匹配。
4. GSE227919 表达矩阵为 **原始非log（TPM 量级）值、空格分隔、基因为符号**（非 Ensembl、非 log2）；GSE26549 为芯片归一化值。

## 五、下一步
- **Phase 3 单细胞**：GSE196296 三个样本（LesionA/B/C，同一 PVL 病灶内 KUS→corrugated→bulky 连续谱）的 10x counts 已全部下载（raw/ 下 9 个文件）。Seurat 已装，可直接做 QC→整合→聚类→注释。
- **Phase 4 风险模型**：GSE26549 候选签名建议用宽松阈值（adj.P<0.05 & \|logFC\|>0.5 ≈ 101 基因，或 top N）送入 LASSO，而非仅 13 个严格 DEG。
- **Phase 5 LINCS 药物重定位**：基于 Phase 2/3 的上调基因做转录组匹配。
