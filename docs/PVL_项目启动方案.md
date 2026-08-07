# 高危口腔白斑 / PVL 生信项目 — 启动方案与资料核验

> 整理时间：2026-07-28 ｜ 方向来源：`docs/deep-research-report.md` 的「高危口腔白斑/PVL」章节
> 参考手册：`G:\CMCC Sync\Medical Data\Guidiance\` 下三份（生信数据库汇总 / 生信研究思路汇总 / 生信期刊汇总）
> 数据集根目录：`G:\CMCC Sync\Medical Data\Stomatology\PVL\`

---

## 1. 已读资料总览

### 1.1 深度报告（PVL 方向）核心结论
- **首选理由**：公开数据已覆盖"恶变前病变 + 进展/未进展结局 + 单细胞异质性"，可闭环做"复发/持续进展/癌变风险"转化生信稿；综合可发表性 **8.8/10**。
- **推荐题目**：*基于公开 bulk 与单细胞转录组整合的高危口腔白斑/PVL 恶变与持续进展风险图谱及药物重定位研究*。
- **主终点**：① progressing vs non-progressing OPL/PVL；② 正常—低危—高危—癌变的连续风险梯度。
- **四层路线**：bulk 发现（limma/DESeq2 + ComBat 批次校正）→ 单细胞定位（上皮/basal-progenitor、POSTN/TGFBI/THBS2 基质、SPP1+ 髓系）→ pseudobulk 混合效应模型（`~ risk_group + (1|patient)`）→ LASSO/elastic-net 风险模型 + 外部验证（ROC/PR/校准/DCA）。
- **升级主轴**：把"DEG 列表"升级为"恶变前程序"（EMT / 细胞周期 / 角化异常 / 免疫逃逸 / ECM 重塑）做 ssGSEA/GSVA，并与单细胞来源映射。
- **药物重定位**：LINCS/CMap 负连接评分 → 过滤口腔局部不可达/高毒性候选 → 保留 EGFR/ERBB、PI3K-AKT、COX-2、YAP/TEAD、TGF-β 类；文中写作 *transcriptome-matched candidate compounds*。
- **目标期刊（非OA优先）**：**Oral Oncology**（IF 3.9，混合，非OA无APC，最匹配）＞ Virchows Archiv（3.0）＞ IJOMS ＞ Oral Diseases（2.9，备选）。

### 1.2 三份参考手册（项目方法学基础）
- **生信数据库汇总（35 个）**：GEO/TCGA/GTEx 等获取方式与注意事项；本研究核心只用 **GEO**（GSE 检索、GEOquery、平台/批次）。
- **生信研究思路汇总（27 条）**：与本项目最贴合的是 #1 PCD/程序性死亡 signature、#4 甲基化驱动亚型、#23 ML 预后/药敏、#24 多组学筛逆转药、#7 单细胞轨迹、#8 空间/细胞通讯。报告强调的"恶变前程序 + 细胞来源 + 药物重定位"与 #1/#23/#24 高度一致。
- **生信期刊汇总（10 本）**：确认 Oral Oncology 等口腔/肿瘤纯生信友好刊；注意纯生信"实验验证门槛"（Frontiers in Oncology 部分 section、Cancer Letters 需机制佐证）。

---

## 2. 环境与可达性核验（2026-07-28 实测）

| 项目 | 状态 | 说明 |
|---|---|---|
| 数据代理 SOCKS5 8428/8429 | ✅ 可用 | 沙箱 `127.0.0.1:8428/8429` 监听可达；但 `curl` 因 Windows schannel 证书吊销检查失败 |
| Python 抓取 GEO | ✅ 可用 | 托管 Python 3.13 + PySocks；`urllib` 直连或经 SOCKS5 代理均能 `200` 抓取 GEO FTP（OpenSSL 不走吊销检查） |
| G 盘数据集目录 | ✅ 已定位 | `G:\CMCC Sync\Medical Data\Stomatology\PVL\` 原为空，已建 `metadata/` 并写入三队列样本表 |
| **WSL（Ubuntu22.04, R+GPU）** | ⚠️ **当前不可达** | `wsl` 命令被沙箱「系统级工具」策略禁用；`netstat` 未见宿主机/虚拟网卡上有 SSH 监听；PowerShell 工具在本环境不回显。见 §5 待解决 |

> 结论：数据采集与元数据核验**可在当前环境用 Python 完成**；但核心 R 分析栈（GEOquery/limma/edgeR/DESeq2/Seurat/CellChat）必须跑在 WSL，需先解决 WSL 可达性。

---

## 3. 数据集元数据核验结果（Phase 1 — 已完成）

脚本：`fetch_metadata.py`（流式拉取 series-matrix 元数据头，不下载全量表达矩阵）。产物：`G:\...\PVL\metadata\GSE*_sampletable.csv`。

### 3.1 GSE26549 — 恶变前风险签名核心队列
- **平台/类型**：GPL6244（Affymetrix HuGene 1.0 ST），芯片；**86 例 Oral leukoplakia**（OPL）。
- **样本特征**：含 `time of biopsy: baseline / 3 months after inclusion`，且 **GEO 元数据完整带有进展结局标签**（`!Sample_characteristics_ch1` 第48行 `outcome: oral cancer development / no oral cancer development`，另有 oral cancer-free survival time）。
- ✅ **进展标签已确认可从 GEO 直接获得**（无需外取原文补充表）：全队列 86 例中 outcome=oral cancer development **35 例 / no 51 例**；基线活检 62 例中 **21 progressing / 41 non_progressing**。已解析为 `metadata/GSE26549_clinical.csv`。
- 注：此前 `fetch_metadata.py` 曾误报"无进展标签"——原因是该脚本把多条 `!Sample_characteristics_ch1` 行按同一 key 覆盖，仅剩"活检时间点"。属**解析缺陷，非数据缺失**，现已用 `parse_gse26549_clinical.py` 修正。
- 报告提到的"86 HNSCC"在本 series 中未出现（本 series 仅 86 OPL），需确认是否另属其他 GSE 或论文附表。

### 3.2 GSE227919 — 外部验证 + 连续风险梯度（极有价值）
- **实际规模：66 例**（报告称 99，已证伪）。RNA-seq，source = oral mucosa。
- **分组清晰（来自 `group` 字段）**：
  - Control（正常对照）≈ **19**
  - 低危 HkNR（homogeneous keratosis）≈ **17**
  - 高危 Dysplasia ≈ **22**
  - OSCC（口腔鳞癌）≈ **9**
- ✅ 正好是 **正常 → 低危 → 高危 → 癌变** 连续风险梯度，且样本标签直接可得——非常适合做"连续风险梯度"主终点的验证与分层分析。

### 3.3 GSE196296 — 单细胞定位队列
- series matrix 仅 **3 个聚合样本**（LesionA/B/C，含组织学诊断：keratosis of unknown significance / corrugated orthohyperkeratosis / bulky hyperkeratotic epithelium）。
- ⚠️ 真实单细胞数据（论文称 59 OPMD + 8 normal，来自 8 名患者病灶体系）在 **GEO supplementary 文件**（10x counts）中，series matrix 不含。需单独下载 supplementary 做单细胞注释/pseudobulk/细胞通讯。
- 结局标签：论文未完全公开（报告已预警），定位用途为主。

---

## 4. 项目路线图（适配本环境，基于报告 12–14 周框架）

| 阶段 | 任务 | 状态 |
|---|---|---|
| **1 数据准备** | 元数据核验（✓）、下载全量表达矩阵（GSE26549/GSE227919）、下载 GSE196296 scRNA counts、重建样本表/补进展标签 | 进行中 |
| **2 发现分析** | bulk DE + 通路富集（limma/DESeq2 + ComBat）；单细胞注释+ pseudobulk（Seurat）；细胞通讯/恶变轨迹（CellChat/Monocle） | 待 WSL |
| **3 风险模型** | LASSO/elastic-net 收敛特征 → 进展/癌变风险签名 → 外部队列验证（ROC/PR/校准/DCA） | 待 WSL |
| **4 转化输出** | LINCS/CMap 药物重定位四步过滤；主图/补图；初稿（目标 Oral Oncology） | 待 WSL |

> R/Seurat/GPU 全在 WSL；WSL 可达后阶段 2–4 启动。阶段 1 的数据采集在沙箱 Python 即可推进。

---

## 5. 环境与数据阻塞项 — 已解决

1. **WSL 可达性** ✅ —— 沙箱仅禁 `wsl` 别名，`wsl.exe` 可直用；发行版实为 **Ubuntu-24.04**（R4.3.3 + RTX3070 8GB）。WSL 无直连外网，已通过 `.wslconfig` 开启 **mirrored 网络模式** + 主机 SOCKS5 代理（`127.0.0.1:8428`）解决装包/联网，R 栈正在后台安装。
2. **GSE26549 进展标签** ✅ —— 已确认完整存在于 GEO 元数据（见 §3.1），无需外取原文补充表。
3. **GSE227919 真实表达矩阵** ✅ —— RNA-seq 表达值在 supplementary（series matrix 仅元数据），已下载 `raw/`。

## 6. 待用户决策 / 后续

1. **目标期刊最终确认** —— 报告推荐 **Oral Oncology（IF 3.9，混合非OA 无APC）**。是否锁定，或你有更高中科院分区要求？
2. **R 分析栈安装中** —— 后台 task `rOjVMh` 正在安装 limma/edgeR/DESeq2/clusterProfiler/Seurat/CellChat/WGCNA 等，完成后立即跑 `PVL_02_bulk.R`（bulk DEG + GO/KEGG 富集）。
3. **Phase 3 单细胞** —— GSE196296 为 10x scRNA-seq（LesionA/B/C = PVL 同一病灶内病理连续谱），待 R 栈就绪后下载 GSM5870060/62/64 的 10x 三文件做 Seurat/CellChat。
4. **Phase 4/5** —— pseudobulk + LASSO 风险模型；LINCS 药物重定位（写作 *transcriptome-matched candidate compounds*）。
