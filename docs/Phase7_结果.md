# Phase 7 结果：签名功能解构（GO/KEGG/Hallmark/上游 TF）

> 脚本：`PVL_07_signature_decon.R`（WSL `~/PVL/scripts/`）
> 方法：对 Phase 2 的 **101 基因进展签名**（GSE26549 `adj.P<0.05 & |logFC|>0.5`）做富集；因本机 fgsea/clusterProfiler 装不上，改用 `enrichR`（Enrichr API，Kuleshov et al., Nat Methods 2016），数据库 GO BP / KEGG / MSigDB Hallmark / ChEA_2022 / ENCODE。
> 输出：`results/PVL7_*.csv` · `PVL7_summary.csv`

---

## 1. 真实富集结果（签名层面）

| 数据库 | adj.P<0.05 项数 | 顶端项（adj.P / Combined Score / 命中基因） |
|---|---|---|
| GO BP | **3** | steroid catabolic process (0.029, CS=556; SPP1;HSD17B6;CYP19A1)；androgen metabolic process (0.046, CS=296; SPP1;HSD17B6;CYP19A1)；negative regulation of proteolysis (0.029, CS=202; ECM1;CAMK2N1;LRRK2;TIMP2) |
| KEGG | **0** | 无 |
| MSigDB Hallmark | **1** | UV Response Dn (0.023, CS=52.6; IGFBP5;PRKAR2B;HAS2;MET;SYNE1) |
| ChEA_2022（上游 TF） | **0** | 无命中 |
| ENCODE（上游 TF） | **0** | 无命中 |

`PVL7_summary.csv` 核验：`sig_n=101`、`GO_adjp05=3`、`KEGG_adjp05=0`、`Hallmark_adjp05=1`、`ChEA_TF_adjp05=0`、`ENCODE_TF_adjp05=0`。

---

## 2. 解读

- **101 基因签名层面机制解构基本为空**：严格 FDR<0.05 下仅 GO 3 项 + Hallmark 1 项可达标，上游 TF（ChEA/ENCODE）**零命中**。这与签名规模小（101 基因）、多重检验偏严直接相关——小基因集做富集本就难达显著。
- 幸存的少数术语与签名已知方向一致：SPP1（进展上调）、ECM1/TIMP2/IGFBP5（进展下调）均出现在 GO/Hallmark 命中基因中，说明富集并非随机噪声，但**不足以构成连贯的通路模块**。
- ⚠️ **重要纠偏（防误述）**：初稿/修订计划曾称"签名 Hallmark 命中 mTORC1 / PI3K-AKT-mTOR / EMT / KRAS，且上游 TF 顶端为 P53"——**该说法不成立**。那些通路来自 **Phase 2 的 bulk GSE227919 Cancer-vs-Control（2986 DEG）富集**，而非本 101 基因签名；两处富集对象不同，不可混为一谈。机制性叙事必须锚定在 *bulk 癌-对照层面*，不能写作「签名→机制→候选」闭环。

---

## 3. 局限

1. 签名集仅 101 基因，富集检验功效低（FDR 严），小集很难产出显著通路。
2. `enrichR` 依赖 Enrichr 在线 API，术语命名/版本随库更新而变。
3. 上游 TF 零命中不排除真实生物学调控，仅说明在 101 基因集层面不可检。

---

## 4. 输出文件

| 文件 | 内容 |
|---|---|
| `PVL7_summary.csv` | 各库达标项数汇总 |
| `PVL7_GO_top.csv` / `PVL7_Hallmark_top.csv` / `PVL7_*_fdr.csv` | 各库顶端项与 FDR 表 |
| `PVL7_enrichment_raw.rds` | 原始富集对象 |

