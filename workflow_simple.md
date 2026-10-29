# 最简版 workflow：TE 外显化候选过滤与敲除验证前评估

**目的**：用最小数据量跑通 E1–E7 每一步，用 IFNAR2-S 阳性对照检验每步实现是否合理，最后产出 5–20 个可做简单敲除验证的靶点。**这是 pilot，不是完整挖掘**。本文件自包含，不依赖其他规范文件；Codex 可直接按本文件 + `codex_plan_simple.md` + `config_simple.yaml` 实现，并包装成可复用 skill。

## 0. 范围：保留什么、砍掉什么

| 保留 | 砍掉（本轮不做） |
|---|---|
| 全事件保留、注释式验证、六态状态 | TCGA 四个癌种 junction 分析 |
| TRAILS 种子全集、≥1bp overlap、双分母 overlap 计算 | 质谱（PXD040962/PXD054305）、Ribo-seq |
| TransDecoder ORF + 50nt NMD 规则 | pseudobulk DE/DTU 统计推断 |
| NetMHCpan 4.1（可选，无许可则 blocked） | 跨数据集整合、批次校正 |
| PRJNA599962 单细胞重跑（CellTypist 注释）为主，Doyle 为备选 | sgRNA 设计与实验（下一轮） |
| liftOver 人→鼠 + unresolved 分类 | |

**通用状态词**（六态，全流程统一）：`supported`（有证据）/ `not_detected`（查了没有）/ `not_testable`（输入不完整，如单外显子无法判 NMD）/ `not_applicable`（不适用）/ `blocked`（外部原因，如许可）/ `failed`（实现或输入错误）。不以缺失冒充阴性。

**事件类型**（E1 分类用）：`cassette_exon`（内部可变外显子，两侧 junction 均在）、`alt_terminal_exon`（替代末端外显子，只需一个进入 junction）、`alt_donor` / `alt_acceptor`（5'/3' 剪接位点移动）、`te_tss`（TE 提供新第一外显子/启动子）、`intron_retention`、`utr_only`（TE overlap 仅落在 UTR）、`internal_overlap`（TE 与外显子重叠但不改变剪接结构）、`unresolved`。

**E5 说明**：E5 = 单细胞数据中逐细胞类型检测（TRAILS 亚群表达归入 E1）。**主路 = PRJNA599962 原始数据重跑 + CellTypist 注释**（条形码环节零障碍：10x 公共体系 + 作者公开代码）；**备选 = Doyle**（两级，见 E5-alt）。

## 1. 数据

| 数据 | 用途 | 入口 | 大小级别 |
|---|---|---|---|
| TRAILS `TRAILS.gtf.gz` + `isoform_info.txt.gz` + `male_PBMC.gtf.gz` | E1 种子全集 + 亚群注释 + 独立验证 | github.com/juninamo/TRAILS | 小（<100 MB） |
| GTEx V9 FLAIR `flair_filter_transcripts.gtf.gz` | E2 正常组织结构检测 | storage.googleapis.com/gtex_analysis_v9/long_read_data/ | 小 |
| **PRJNA599962**（Volden & Vollmers：10x Chromium 3' + R2C2 ONT，2974 PBMC，2 重复；ONT raw + 配对 Illumina fastq） | E5 主路：重跑 + CellTypist | SRA/ENA；作者代码 github.com/rvolden/10xR2C2 | 中（ONT raw ~5600 万 reads + Illumina ~1.8 亿 read pairs） |
| Doyle Zenodo 17341210（最新版本 17372869）：gene/transcript counts + transcript CPM 矩阵 + Ensembl v113 参考 GFF3 | E5 备选 tier 1 | zenodo.org/records/17341210 | tar 1.1 GB（矩阵解压后各 8–42 GB，流式处理） |
| Doyle 原始数据 **PRJNA1432415**（2 run，各 3 flow cell 合并，约 4.9 亿 raw reads） | E5 备选 tier 2 | SRA/ENA；作者管线 github.com/UK-SBCoA-EbbertLab/single_cell_nextflow_pipeline | 大 |
| GENCODE v44 GTF + GRCh38 FASTA + hg38 rmsk | 参考与 TE overlap | GENCODE/UCSC | FASTA ~1 GB |
| hg38→mm39 chain | E6 | UCSC | 小 |

**已核验事实（2026-10-29）**：① Doyle Zenodo **没有**作者 extended annotation GTF，BambuTx 新转录本无坐标；② Doyle 解复用代码公开（PIPseeker + Java Levenshtein 拯救），但 PIPseeker 为 Fluent/Illumina 商业软件，whitelist 随其分发，可获取性未确认；③ PRJNA599962 用 10x 公共体系，作者解复用代码公开且 whitelist 从数据自推断——**无外部依赖，故为主路**。

## 2. 七步流程

### E1. TRAILS 种子库
- 从 TRAILS GTF 提取全部 transcript exon，与 rmsk TE 求 overlap（≥1 bp 全保留）；
- 每行记录：host gene、TE family/class、事件类型（见第 0 节定义）、`isoform_info` 中的亚群特异列（specific_cell_LR 等）与 male_PBMC 验证列；
- **通过标准**：产出 `events_master.tsv`，行数 = 全部 transcript×event×TE overlap；IFNAR2-S 的 E9-Alu 应被分类为 `alt_terminal_exon`。

### E2. FLAIR v9 检测
- 对每个候选事件的 required junction（cassette 两个、terminal 一个），在 FLAIR GTF 中查精确匹配（同 chr/strand、坐标一致）；
- 输出 `e2_flair_match`（exact/partial/absent）。FLAIR GTF 无表达量，**只做结构检测**，不判表达阴性；
- **IFNAR2-S 预期**：E9-Alu 的进入 junction 在 FLAIR GTF 中 exact match。

### E3. ORF / NMD 预测
- `gffread` 提取转录本 FASTA → TransDecoder（`-m 40`，保留 partial）；
- NMD：complete ORF 且 stop 距最后一个 exon junction >50 nt → `predicted_NMD`；single-exon → `not_applicable`；不完整 → `not_testable`；
- **IFNAR2-S 预期**：E9-Alu isoform 产生完整 ORF（诱饵受体，非 NMD）。

### E4. NetMHCpan 4.1（可选）
- 仅当本地已装并有许可：对 TE 编码/跨 junction 肽段枚举 8–11 mer，EL rank <0.5 strong、<2.0 weak；
- 无许可/无 allele 列表 → 整步 `blocked`，**不阻塞 E5–E7**，也不影响靶点初筛；
- **IFNAR2-S 预期**：能跑通则输出预测表；跑不通只标 blocked，不算失败。

### E5. PRJNA599962 重跑 + CellTypist 逐细胞类型检测（主路）

数据：Volden & Vollmers 2022，10x Chromium 3' 全长 cDNA 分流为 Illumina（短读，~5000 分子/细胞）与 R2C2 ONT（长读，~4000 consensus reads/细胞），2974 个 PBMC、2 个重复。**设计要点：用 Illumina 做细胞类型注释（深度高），把标签经条形码对应转移到 ONT 细胞（99.1% 的 R2C2 细胞出现在 Illumina 数据中），用 ONT reads 做 isoform 级检测。**

1. **V0 取数**：从 SRA/ENA 拉 PRJNA599962 的 ONT raw reads（2 重复）与 Illumina fastq（2 重复）；核对 run↔重复对应；
2. **V1 R2C2 consensus**：C3POa v2.2.2 把 ONT raw reads 转成 consensus reads → Splint-UMI 合并 → 10X-UMI 合并（作者工具 ExtractUMIs/MergeUMIs，github.com/rvolden/10xR2C2）；预期产出 ~1480 万合并后 consensus reads；
3. **V2 解复用**：按作者方法——先数出样本中最高频的 ~1500 个细胞条形码作 whitelist，再按 Levenshtein 距离分配（L1<3 且 L1<L2−1），~81% reads 可分配到单个细胞；输出 per-cell R2C2 fastq；
4. **V3 Illumina 基因计数**：Cell Ranger（或 kb-python/STARsolo）处理 Illumina pairs → gene × cell 矩阵；Seurat QC（作者起点参数：nFeature 200–2500、percent.mt<5）；
5. **V4 CellTypist 注释**：标准化（log1p）→ CellTypist `Immune_All_Low` 模型 → 每细胞标签 + 置信度；**标签经条形码匹配转移到 ONT 细胞**；预期主要类型：CD8+ T、CD4+ T、B、单核、NK 等；
   - **分辨率警告**：每个细胞仅 ~4000–5000 分子，CD4/CD8 边界可能模糊——CellTypist 低置信度细胞标 `low_confidence`，CD8 标签需 CD8A/CD8B marker 共表达复核；最终以"细胞类型 × 置信度"双列输出，不假装分辨率高于数据；
6. **V5 逐细胞类型 isoform 定量**：按 CellTypist 标签把 R2C2 reads 合并为每类型 synthetic bulk（同作者做法）→ Bambu v3.8.1 `--annotation` 传 **TRAILS GTF**（不开 discovery）→ TRAILS isoform × 细胞类型计数矩阵；另跑一遍 Ensembl 113 discovery 作 QC（应得到与论文同量级的 B/T/单核转录本数）；
7. **V6 事件检测矩阵**：输出 `e5_volden_celltype.tsv`：event × celltype × 检测状态（该类型 synthetic bulk 中宿主 TRAILS isoform 的 read 数与检出/未检出）；
- **只描述，不检验**（~3000 细胞、2 重复，深度不足以做统计）；2 个重复一致性作 QC 报告，不合并假装重复；
- **IFNAR2-S 预期**：IFNAR2-S 在 T 细胞（及 B/单核）synthetic bulk 中可检出。

### E5-alt. Doyle 备选（两级，按需启动）

**触发条件**（满足其一才启动）：① 短名单前 20 中过半候选在 Volden 中 T 细胞未检出（可能是深度问题而非真阴性）；② CD4/CD8 区分对候选排序有实质影响且 CellTypist 置信度不足以区分。

- **Tier 1（快，无重跑）**：Doyle Zenodo 作者矩阵（transcript counts，流式分块读）+ Ensembl v113 GFF3 坐标，仅匹配 **Ensembl 注释转录本**的事件；BambuTx 与 TRAILS 新结构 isoform 标 `not_testable`；输出 `e5a_doyle_celltype.tsv`（event × celltype × n_cells_detected × fraction）；1 名供者 2 技术重复合并为 1 个生物学单位，只描述不检验；
- **Tier 2（完整，重跑）**：PRJNA1432415 原始数据（2 run，各 3 flow cell 合并，~4.9 亿 reads）按 Doyle 原文复现：dorado super-accurate（若 POD5）→ pychopper → **作者公开 Nextflow 管线 Step 0**（github.com/UK-SBCoA-EbbertLab/single_cell_nextflow_pipeline：PIPseeker v02.01.04 + Java Levenshtein 分层条形码拯救）→ minimap2 v2.28 splice-aware（MAPQ≥10，去 secondary/supplementary）→ Bambu v3.8.1 双跑（① TRAILS GTF 注释定量；② Ensembl 113 discovery QC，应回收 ~126 novel isoform）→ 合并技术重复 → `e5b_doyle_celltype.tsv`。**PIPseeker/whitelist 不可获取则 Tier 2 blocked**，不猜条形码结构。Tier 2 可给出 CD8 naive/CM/EM 分辨率，是 Volden 分辨率不足时的最终升级路径。

### E6. liftOver 人→鼠
- 对每个事件的 TE exon 用 liftOver（hg38→mm39，minMatch 0.95）；
- 映射成功 → 检查直系同源区域是否含同源 TE exon；unmapped/split → `unresolved`（**不等于 human-specific**）；
- **IFNAR2-S 预期**：E9-Alu 在小鼠同源区域缺失，流程应输出 `syntenic_absence_supported` 或 `unresolved`，不得输出错误的"保守"。

### E7. 宿主基因功能
- 对每个候选的 host gene 打两个标签：① 九个免疫 GO 术语命中——GO:0006955、GO:0006954、GO:0034097、GO:0005125、GO:0001816、GO:0001817、GO:0032606、GO:0060337、GO:0034340（用 GO ID 查现名，不沿用原文可能对调的标签）；② 简单功能注释（膜蛋白/信号/转录因子，来自 isoform_info 已有列）；
- 不做富集分析。

## 3. 靶点排序（简单字典序，无加权打分）

按以下顺序取前 5–20 个：
1. E1 结构明确（cassette 或 terminal，非 UTR-only）；
2. E2 FLAIR exact match；
3. E3 完整 ORF 且非 predicted_NMD（敲掉才可能有效果）；
4. E5 主路 T 细胞（含 CD8 亚群，若 CellTypist 置信度足够）中检出；Doyle 备选的 CD8 检出作为升级证据；
5. E7 免疫 GO 命中优先；
6. E6 unresolved 或 human-specific 优先（新颖性，非硬性）。

**同时输出三条对照行**：IFNAR2-S（阳性）、一个 E2 absent 的事件（阴性流程对照）、一个 UTR-only 事件（分类对照）。

## 4. 简单敲除验证建议（下一轮，仅方向不设计 sgRNA）

- 策略：TE exon 上下游各 2–3 条 crRNA，双 guide 删除 TE exon，PCR+凝胶确认缺失，RT-qPCR 测 JET/canonical junction PSI 确认剪接改变；
- 模型：先 Jurkat 或 HEK293T（易编辑）验证剪接改变，再进 PMA/ionomycin 激活的原代 CD8 看增殖/杀伤表型；
- IFNAR2-S 不需要敲除验证（已发表），但可作为 **qPCR/PSI 实验流程的阳性对照**。

## 5. 明确边界

- 本轮不声称发现功能性靶点；产出是"流程跑通 + 每步合理性经 IFNAR2-S 校验 + 候选短名单"；
- E2 结构检测 ≠ 表达验证；E5 单供者/低深度检测 ≠ 群体结论，CD4/CD8 分辨率受深度限制并显式标注置信度；E4、E5-alt blocked 不影响其余步骤；
- 所有事件保留在 master，排序只生成派生视图；
- E5 主路与备选的计数单位分开报告，不混用；Volden 的 2 个重复作 QC 一致性报告，不充当统计重复。
