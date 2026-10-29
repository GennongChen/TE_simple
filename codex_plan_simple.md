# Codex 实施计划（最简版，自包含）

对应 `workflow_simple.md`（自包含规范，无外部依赖）。顺序脚本，不用 Snakemake；每步一个入口，失败即停并写状态。所有命令经 `mamba run -n seurat5`。**本文件是实现合同**：脚本、目录、测试由 Codex 实现并通过验收后才算完成。目标形态：整个 `simple/` 可被包装成可复用 skill（输入 = config 中 URL/accession，输出 = results/ 下的表）。

## 目录

```
simple/
  config_simple.yaml
  scripts/
    S0_preflight.py        # 环境/参考/数据入口检查
    S1_download.py         # 小文件下载+校验（TRAILS/GTEx GTF）
    S2_catalog.py          # E1: TRAILS×rmsk overlap + 事件分类
    S3_flair_match.py      # E2: FLAIR v9 junction 匹配
    S4_orf_nmd.py          # E3: gffread + TransDecoder + NMD
    S5_netmhcpan.py        # E4: 可选，无许可则 blocked
    V0_fetch_sra.py        # E5 主路: PRJNA599962 取数（ONT raw + Illumina fastq）
    V1_r2c2_consensus.py   # E5 主路: C3POa + Splint/10X-UMI 合并
    V2_demux.py            # E5 主路: 数据推断 whitelist + Levenshtein 分配
    V3_illumina_counts.py  # E5 主路: Cell Ranger/kb-python 基因矩阵 + Seurat QC
    V4_celltypist.py       # E5 主路: CellTypist 注释 + 标签转移到 ONT 细胞
    V5_bambu_quant.R       # E5 主路: 按类型合成 bulk → Bambu(TRAILS GTF) + QC 跑
    V6_event_matrix.R      # E5 主路: event × celltype 检测矩阵
    A1_doyle_tier1.R       # E5 备选 tier 1: Doyle 作者矩阵 ID 匹配（流式）
    A2_doyle_tier2_*.py/R  # E5 备选 tier 2: PRJNA1432415 重跑（取数/预处理/解复用/比对/Bambu/合并）
    S7_liftover.py         # E6: liftOver + 分类
    S8_function_rank.py    # E7 + 排序 + 输出
  refs/  data/  results/{events_master.tsv, evidence/, candidate_shortlist.tsv, control_rows.tsv}  tests/
```

## 每步合同

| 步骤 | 输入 | 输出 | 通过标准 |
|---|---|---|---|
| S0 | config | preflight.json | YAML 可解析；gffread/TransDecoder/bedtools/liftOver/R 包存在；E5 主路依赖（C3POa、cellranger 或 kb-python、celltypist+模型）缺失只标 E5 blocked；缺 NetMHCpan 只标 E4 blocked |
| S1 | URL 列表 | data/ + sha256 | 每文件校验和记录；>100 MB 需人工确认 |
| S2 | TRAILS GTF + rmsk | events_master.tsv | 全 overlap 保留；IFNAR2-S E9-Alu 分类为 alt_terminal_exon |
| S3 | events_master + FLAIR GTF | evidence/e2.tsv | junction 精确匹配；IFNAR2-S exact |
| S4 | GTF→FASTA | evidence/e3.tsv | ORF+NMD 列齐全；IFNAR2-S 完整 ORF 非 NMD |
| S5 | 肽段 + alleles | evidence/e4.tsv 或 blocked | 无许可时 blocked，不阻塞后续 |
| V0–V2 | PRJNA599962 ONT raw | per-cell R2C2 fastq | consensus ~1480 万；解复用率 ~81% 量级；显著偏低 = 参数问题 |
| V3–V4 | Illumina fastq | gene×cell 矩阵 + celltypist 标签 | QC 过滤后聚类出主要免疫类型；CD8 标签需 CD8A/CD8B 复核；低置信度标 low_confidence |
| V5–V6 | per-cell fastq + 标签 | evidence/e5_volden_celltype.tsv | TRAILS GTF 定量跑通；Ensembl discovery QC 与论文量级一致；IFNAR2-S 在 T 细胞 bulk 检出 |
| A1 | Doyle Zenodo 矩阵 | evidence/e5a_doyle_celltype.tsv 或 skipped | 仅 Ensembl ID 匹配；BambuTx/TRAILS 新结构标 not_testable；流式读矩阵 |
| A2 | PRJNA1432415 | evidence/e5b_doyle_celltype.tsv 或 blocked | 作者 Nextflow 管线 Step 0 解复用；PIPseeker 不可获取 → blocked；discovery QC ~126 novel isoform |
| S7 | TE exon BED + chain | evidence/e6.tsv | unmapped→unresolved；IFNAR2-S 非保守 |
| S8 | 全部 evidence | candidate_shortlist.tsv + control_rows.tsv | 三条对照行正确；短名单 5–20 |

## 关键实现约束

1. **坐标**：GTF 1-based → BED 0-based half-open；overlap 双分母（fraction_exon、fraction_TE）分别算；负链 exon rank 反向。
2. **事件分类**：terminal exon 只要求一个进入 junction；FSM 也保留；UTR-only 不删。
3. **E2**：FLAIR GTF 是结构目录，无表达量；只输出 exact/partial/absent，不写 TPM。
4. **E5 主路（V0–V6）**：
   - 严格按 Volden 原文复现：C3POa v2.2.2 → Splint-UMI 合并 → 10X-UMI 合并（ExtractUMIs/MergeUMIs）→ whitelist 取最高频 ~1500 条形码 → Levenshtein 分配（L1<3 且 L1<L2−1）；
   - **注释用 Illumina、isoform 用 ONT**：CellTypist 跑在 Illumina 基因矩阵上（深度高），标签经条形码匹配转移到 ONT 细胞（99.1% 重叠）；
   - CellTypist 低置信度细胞标 `low_confidence`；CD8 标签必须 CD8A/CD8B marker 共表达复核；输出"细胞类型 × 置信度"双列；
   - V5 按标签合成 per-celltype bulk → Bambu v3.8.1 `--annotation` 传 TRAILS GTF、不开 discovery；另跑 Ensembl 113 discovery 作 QC；两跑分开存放；
   - 2 个重复只作 QC 一致性报告，不合并假装重复、不做统计检验；
   - celltypist（Python 包 + Immune_All_Low 模型）缺失 → E5 标 blocked，不用 marker 硬注释冒充。
5. **E5 备选（A1/A2）**：仅在触发条件满足时运行（workflow E5-alt）；A1 流式分块读 8–42 GB 矩阵，counts 列是 raw integer，CPM 文件只做 presence 描述；A2 按 Doyle 原文（dorado/pychopper/minimap2 v2.28/MAPQ≥10/Bambu v3.8.1），解复用只用作者 Nextflow 管线 Step 0，PIPseeker 不可获取 = blocked；每 flow cell checkpoint。
6. **E6**：liftOver 失败 = `unresolved`，不是 human-specific；不合并不同 chain 输出。
7. **E7**：GO 用 ID 查现名（GO:0006954/0006955 在原文可能对调）；不跑富集。
8. **排序**：字典序（workflow 第 3 节），无加权分数；不足 20 个就输出实际数量；E5 主路 T 细胞检出为排序键，Doyle CD8 检出只作升级证据列。
9. **失败行为**：任何一步失败写 status.json 并停止依赖步骤；E4、E5 备选 blocked 不阻塞其余步骤；不产出全 0 假表。

## 命令

```bash
mamba run -n seurat5 python scripts/S0_preflight.py
mamba run -n seurat5 python scripts/S1_download.py
mamba run -n seurat5 python scripts/S2_catalog.py
mamba run -n seurat5 python scripts/S3_flair_match.py
mamba run -n seurat5 python scripts/S4_orf_nmd.py
mamba run -n seurat5 python scripts/S5_netmhcpan.py   # 可 blocked
# E5 主路
mamba run -n seurat5 python scripts/V0_fetch_sra.py
mamba run -n seurat5 python scripts/V1_r2c2_consensus.py
mamba run -n seurat5 python scripts/V2_demux.py
mamba run -n seurat5 python scripts/V3_illumina_counts.py
mamba run -n seurat5 python scripts/V4_celltypist.py
mamba run -n seurat5 Rscript scripts/V5_bambu_quant.R
mamba run -n seurat5 Rscript scripts/V6_event_matrix.R
# E5 备选（触发条件满足时）
mamba run -n seurat5 Rscript scripts/A1_doyle_tier1.R
mamba run -n seurat5 python scripts/A2_doyle_tier2_run.py
mamba run -n seurat5 python scripts/S7_liftover.py
mamba run -n seurat5 python scripts/S8_function_rank.py
```

## 测试（tests/）

1. **IFNAR2-S 回归**（核心验收）：S2–S8 各步对 IFNAR2-S 的输出符合 workflow 第 2 节预期；任何一步不符 = 该步实现有 bug，修完才能进 S8。
2. 边界 fixture：1bp overlap、负链、terminal vs cassette、NMD 50/51nt、liftOver unmapped、Levenshtein 阈值边界（L1=2 保留 / L1=3 丢弃 / L1=L2−1 歧义）。
3. 空结果合法：某事件 E2 absent 输出 absent 而非删行。
4. E5 QC：解复用率 ~81% 量级；Ensembl discovery QC 的 B/T/单核转录本数与论文 Table 2 同量级（T 细胞 ~14,925 多外显子 isoform 量级）；CellTypist 标签分布与论文三主群一致且 CD8/CD4 可分。

## 验收

- `events_master.tsv` 行数 = TRAILS 全部 transcript×event×TE overlap；
- `control_rows.tsv` 三条对照全部正确；
- `candidate_shortlist.tsv` 5–20 行，每行含 host gene、TE 坐标、event 类型、各步状态；
- E5 备选若启动：相对主路的覆盖增量有报告；
- 不声称已完成敲除实验或发现功能性靶点。

## Skill 包装要求

- 所有外部资源（URL、SRA accession、Zenodo DOI、工具版本、CellTypist 模型名）只从 `config_simple.yaml` 读取，脚本内不硬编码；
- 每个脚本幂等：输入已存在且校验和匹配则跳过；
- `status.json` 记录每步状态（六态），skill 调用方据此决定是否重入。
