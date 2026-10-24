# Oral Oncology 投稿包 — 上传清单与说明

**生成日期：** 2026-10-24
**稿件：** A validated transcriptomic signature of malignant progression in oral leukoplakia and proliferative verrucous leukoplakia: independent replication, single-cell mapping, and in silico drug repositioning
**目标期刊：** Oral Oncology (Elsevier) · Original Research Article
**投稿入口：** https://submit.elsevier.com/OO

---

## 一、本目录文件与 Editorial Manager 上传对应关系

| 文件 | EM 中的 Item Type | 说明 |
|---|---|---|
| `CoverLetter.docx` | Cover Letter | 7 段 668 词；含贡献、主动自我限制、单位与选题说明、常规声明（md 源文件在 `../CoverLetter.md`，不在上传目录） |
| `TitlePage.docx` | Title Page | 独立文件，含五作者、三单位、ORCID、通讯、字数声明、COI 摘要 |
| `Manuscript.docx` | Manuscript | 正文（A4 / TNR 12pt / 双倍行距 / 带页码），已删 Submission notes 与题名页区块；21 页 |
| `Highlights.docx` | Highlights | 5 条，84/83/80/84/83 字符，均 ≤85。**正文内已移除**，仅此处提交 |
| `Fig1.tif` | Figure | 3464×1420 px @300 dpi |
| `Fig2.tif` | Figure | 1860×1680 px @300 dpi |
| `Fig3.tif` | Figure | 2160×1560 px @300 dpi |
| `Fig4.tif` | Figure | 2274×1092 px @300 dpi |
| `Supplementary_FigS1.tif` | Supplementary Material | 1800×2100 px @300 dpi |
| `Table1.docx` | Table | 8 行 × 7 列，可编辑 |
| `Table2.docx` | Table | 2 行 × 5 列，可编辑 |
| `Table3.docx` | Table | 1 行 × 7 列，可编辑 |

**注意：** 图与表合计 7 个，已达期刊上限（4 图 + 3 表）。`Manuscript.docx` 只含图注与三张表，不含图片。

**上传顺序建议：** Cover Letter → Title Page → Manuscript → Highlights → Table 1–3 → Fig 1–4 → Supplementary Fig S1。

---

## 二、投稿前状态

**全部文件已就绪，可直接上传。** 无需再做转换或清理。

- `Manuscript.docx` 已由 `submission_OralOncology_v1.md` 转换完成，并已删除文末 Submission notes 整块、题名页区块、文件头投稿元信息
- 图表均以独立文件提交，正文不嵌图（表格以原生 Word 表格内嵌于正文，同时另有独立 Table 文件备用）
- md 源文件 `../submission_OralOncology_v1.md` 保留在仓库中作为版本追溯，**不上传**

---

## 三、已核验的合规指标

| 指标 | 本稿 | 期刊上限 | 状态 |
|---|---|---|---|
| 正文词数（Introduction–Conclusion） | 3,162 | ≤3,500 | ✅ |
| 摘要词数 | 217 | 250（Elsevier 惯例） | ✅ |
| 关键词 | 9 | ≤10 | ✅ |
| 图 + 表 | 4 + 3 = 7 | ≤7 | ✅ |
| 参考文献 | 23 | ≤60 | ✅ |
| Highlights | 5 条，最长 84 字符 | 3–5 条，≤85 字符 | ✅ |
| 图件分辨率 | 300 dpi（TIFF, LZW） | ≥300 dpi | ✅ |
| 表格格式 | Word 可编辑表格 | 须可编辑文本 | ✅ |

---

## 四、期刊要求的声明项自检

| 期刊要求 | 本稿位置 | 状态 |
|---|---|---|
| Conflict of Interest Statement（正文前部） | Abstract 后 `## Conflict of Interest` | ✅ 已补 |
| CRediT（逐人标注） | Declarations | ✅ 五作者已填（Ma Keke / Liu Ruiguang / Li Yunpeng / Liu Hanping / Wang Yuqian） |
| Declaration of generative AI（**manuscript preparation process**） | Declarations | ✅ 标题与模板措辞均按期刊固定要求，工具名已填（ChatGPT and WorkBuddy） |
| Funding sources | Declarations | ✅ |
| Data availability | Declarations | ✅ concept DOI 21823583 + v2 23453360 |
| Ethics statement | Declarations | ✅ 公共数据二次分析，免伦理审批 |
| Acknowledgements | 无（不需要） | — |

---

## 五、Cover letter 说明

已撰写 `CoverLetter.docx`（7 段 668 词，md 源文件在 `../CoverLetter.md`）。结构与要点：

1. 投稿意图与稿件类型
2. 问题与文献空白：dysplasia 分级主观、切除未证降低 OSCC 发生率；Saintigny 29-transcript 与 Foy 亚型同源 GSE26549 且缺 RNA-seq 独立验证
3. 三点贡献（均有具体数字）：跨平台验证 / 单细胞定位 / LINCS 8 候选
4. **主动自我限制**（核心段落）：CV-AUC=1.0 标为乐观上界、η²=0.006 作诚实 null、GSE85195 仅横断面、签名级富集稀疏故 mTORC1 限 bulk 层面、Saintigny 明认同源非独立
5. 单位与选题理由：主职眼科 + 公共数据二次分析 + 不声称口腔临床专长 + 选口腔因连续谱最完整
6. 常规声明：原创、未发表未他投、全体同意、无 COI、无资助、数据代码公开、免伦理、AI 已声明
7. 收束：只提供可检验假设，不作临床推荐

**已确认（2026-10-24）：**
- 主编称呼用中性的 "Dear Editor-in-Chief,"（公开信息对 EiC 是 Psyrri 还是 Ferris 有冲突），作者已确认无误
- 无基金资助：`no specific funding` 与手稿 Funding 段一致，作者已确认属实

---

## 六、生成脚本

- `build_submission_package.py`（仓库根目录）：从 `manuscript/submission_OralOncology_v1.md` 与 `results/manuscript_figs/` 自动生成 `TitlePage.docx`、`Highlights.docx`、`Table1–3.docx` 与 5 个 TIFF。稿件内容若再修改，重跑该脚本即可保证图表与表格同步。
- `manuscript/_build_manuscript_docx.py`：生成 `Manuscript.docx`（正文 Word）。
- `manuscript/_verify_manuscript_docx.py`：校验 `Manuscript.docx` 的文本一致性、字体、行距、表格与残留项。

两个脚本均只读取 md 源文件，不修改任何原始数据。
