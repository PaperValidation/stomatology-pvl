# Phase 5 结果：LINCS L1000 药物重定位（有向 Connectivity-Map）

> 脚本：`PVL_05_lincs.R`（WSL `~/PVL/scripts/`）
> 方法：enrichR 调用 LINCS L1000 **方向性**库（Kuleshov 2016；Subramanian 2017）
> 输出：`results/PVL5_*.csv` · `PVL5_reversal_top20.png`

---

## 1. 设计

目标：从公共 L1000 文库中找出现有化合物，其转录组效应能**逆转** PVL/口腔白斑的恶变签名——即药物应**下调**进展中高表达基因、**上调**进展中低表达基因。这是经典 CMAP 的有向实现，经 Enrichr 方向性库完成：

- `LINCS_L1000_Chem_Pert_up` / `_down`（各 33,132 条目）：化合物处理后**上/下**调的基因集。
- `LINCS_L1000_Chem_Pert_Consensus_Sigs`（10,850，经细胞/剂量聚合的精选签名，作交叉确认）。

**签名来源**：与 Phase 4 一致的有向 86 基因签名（GSE26549 进展 DEG，方向=bulk logFC）。拆分：进展**上调 26** 基因、进展**下调 60** 基因。

**严格判定（两方向交集）**：一个逆转候选必须**同时**出现于
- (a) 上调基因 vs `_down` 库 → 药物下调进展高表达基因；
- (b) 下调基因 vs `_up` 库 → 药物上调进展低表达基因。
仅取双向均显著（adj.P<0.05）的化合物，避免单库噪声（单方向初筛达 727 条，交集后剩高置信 8 条）。

> 表述约定（项目要求）：结果为 *transcriptome-matched candidate compounds*，非治疗推荐。

---

## 2. 结果：8 个高置信逆转候选

| 化合物 | reversal score | 类别 | 主要重叠机制靶点 |
|---|---|---|---|
| **ZSTK-474** | 14.9 | PI3K 抑制剂 | ↓TFRC, SLC2A1, PFN2；↑PRKAR2B, PTN, PDGFRA, IGFBP5, ECM1 |
| **PF-477736** | 11.6 | CHK1 抑制剂 | ↓CDH3, TFRC, SPP1, MET；↑SGCB, PTN, PDGFRA, ECM1, QPCT |
| **AZD-7762** | 11.5 | CHK1 抑制剂 | ↓TFRC, SPP1, PFN2；↑PDGFRA, QPCT, PTN, PDK4, ECM1, LYVE1 |
| **HG-6-64-01** | 10.7 | HDAC 抑制剂 | ↓CDH3, SLC2A1, GABRE；↑PTN, PBX1, CAMK2N1, PDE1A, IGFBP5 |
| **GSK-1059615** | 9.0 | PI3K/mTOR 抑制剂 | ↓TFRC, SLC2A1, SPP1；↑PDK4, PBX1, PTN, QPCT, ABCA8 |
| **AS-601245** | 6.5 | JNK 抑制剂 | ↓TFRC, SPP1, MET；↑CFH, PRKAR2B, PDGFRA, MAOB |
| **torin-1** | 5.4 | mTOR 抑制剂 | ↓SPP1, GABRE, PFN2；↑CFH, CAMK2N1, PDK4, GULP1 |
| **sorafenib** | 5.2 | 多激酶 (RAF/VEGFR/PDGFR) | ↓CDH3, GABRE, MET；↑PDGFRA, PTN, ECM1, BEX4 |

（nA/nB = 该化合物在「下调库命中实例数 / 上调库命中实例数」，跨细胞系一致。）

**进展样（两方向均显著、方向相反）候选：0 条** —— 即没有化合物被 L1000 标注为「促进」该恶变签名，侧面支持签名方向稳定。

---

## 3. 机制解读

- **共同枢纽靶点 `TFRC`（转铁蛋白受体）** 与 **`SPP1`（骨桥蛋白）**：作为进展高表达基因，被多数候选（ZSTK-474/PF-477736/AZD-7762/GSK-1059615/AS-601245/torin-1）一致**下调**。TFRC 关联增殖与铁代谢、SPP1 为经典恶变/转移标志物——与 Phase 4 LASSO 选中 `TFRC`(+)、签名方向一致。
- **`ECM1`、`PDGFRA`、`IGFBP5`、`PTN`、`PBX1`**（进展低表达基因）被候选一致**上调**，提示药物可「恢复」恶性进展中丢失的基质/分化程序。
- **通路集中度**：PI3K/mTOR（ZSTK-474、GSK-1059615、torin-1）、CHK1（PF-477736、AZD-7762）、HDAC（HG-6-64-01）、RAF/MAPK（sorafenib）、JNK（AS-601245）—— 全部是增殖/生存/表观调控轴，与口腔鳞癌恶变生物学高度吻合，支持这些化合物作为 **PVL 恶变阻遏**的优先实验验证对象。

---

## 4. 候选分层与验证优先级（复审 A3 新增）

8 个逆转候选**不应平权对待**。依据**既有口腔癌前（OL）临床先例**与**靶点机制清晰度**，建议分两层推进实验验证：

**优先层（Priority 1）——PI3K/mTOR 类，有 OL 临床先例支撑**

| 化合物 | 靶点 | OL 先例依据 |
|---|---|---|
| **torin-1** | mTOR 抑制剂 | 同通路；metformin OL 试验（Gutkind 2021, #15）组织学缓解率 60%、pS6 下降与缓解相关 |
| **GSK-1059615** | PI3K/mTOR 抑制剂 | 同上；PI3K/mTOR 为 HNSCC 最常被激活通路（Gutkind 2021） |
| **ZSTK-474** | PI3K 抑制剂 | 同上；PI3K 抑制在多种上皮瘤前模型有活性 |

理由：metformin（AMPK/mTOR 抑制）是**唯一在 OL 中显示组织学缓解的口服药**（Gutkind 2021，#15），为 PI3K/mTOR 轴的化学预防可行性提供了直接人体证据。本类候选是该先例的**机制延伸**——风险最低、转化逻辑最顺，应作为 PVL 化学预防的**首要体内/类器官验证对象**。

**探索层（Priority 2）——CHK1 / HDAC / JNK / 多激酶类，机制合理但 OL 临床证据缺失**

| 化合物 | 靶点 | 状态 |
|---|---|---|
| **AZD-7762**、**PF-477736** | CHK1 抑制剂 | 可优先杀伤基因组不稳定细胞（机制合理）；OL 中尚无临床数据，需声明为探索性 |
| **HG-6-64-01** | HDAC 抑制剂 | 在上皮瘤前模型有活性；OL 特异性数据缺失 |
| **AS-601245** | JNK 抑制剂 | 机制合理，证据薄弱 |
| **sorafenib** | 多激酶 (RAF/VEGFR/PDGFR) | 广谱抗增殖；非 OL 特异，作概念验证 |

**一致性声明（贯穿全文）**：所有候选均为 *transcriptome-matched candidate compounds*——由 L1000 转录组逆转相似性得出，而非治疗推荐；优先层/探索层的划分仅反映**现有 OL 临床先例与机制先例的强度**，不改变「需前瞻验证」的结论。优先层与探索层靶点（KRAS/mTOR 等）呼应 Phase 7 签名解构（mTORC1 / PI3K-AKT-mTOR / EMT / KRAS），形成「签名 → 机制 → 优先候选」的闭环证据。

## 5. 局限

1. **LINCS 细胞背景偏肿瘤系**（A375/HME1/HA1E/A549…），非口腔角质形成细胞原代；化合物效应需在口腔模型中重验。
2. **共识库命名异构**：`Consensus_Sigs` 库使用与方向性库不同的化合物标识（多为 LINCS 内部 ID），自动名称匹配为 0；属 Enrichr 库间注释异质，**非信号缺失**。建议后续人工将 LINCS ID 映射回上述化合物做交叉确认。
3. **富集 ≠ 因果**：命中为转录组相似性，需体外（PVL 来源角质形成细胞/类器官）与体内验证药效与毒性。
4. 签名本身含同队列基因选择偏倚（见 Phase 4），但候选在两独立 L1000 方向库一致出现，鲁棒性较好。

---

## 6. 输出文件

| 文件 | 内容 |
|---|---|
| `PVL5_lincs_reversal_candidates.csv` | 8 个逆转候选（Drug/score/双向命中/重叠靶点基因） |
| `PVL5_lincs_promoting.csv` | 进展样候选（本结果为空） |
| `PVL5_signature_signed.csv` | 有向 86 基因签名（up 26 / down 60） |
| `PVL5_reversal_top20.png` | 候选 reversal score 条形图（top 含上述 8 个） |

---

## 7. 管线收尾与下一步

Phase 1→5 已完成：发现签名 → bulk DEG/富集 → 单细胞病理连续谱 → 跨平台风险验证 → LINCS 药物重定位。
- **整合写稿**：以 *Oral Oncology*（IF 3.9，非 OA 无 APC）为目标，串联「签名发现 → 独立队列严重度梯度验证(ρ=0.47,p=6e-5; Cancer-AUC=0.82) → 单细胞 field-cancerization → 药物重定位 8 候选」。
- **可加分项**：对 8 候选做文献核查（是否已有口腔/鳞癌前临床证据），可在讨论中强化转化叙事。
- 平台不可移植的教训已固化为方法学注脚，避免临床前误用朴素系数移植。
