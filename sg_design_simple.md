# sgRNA 设计规范：外显化 TE 敲除筛选（SpRY-Cas9 tiling）

依据：Shi et al., Cell 2026（HF1-SpRY-Cas9 全位点 tiling 突变筛选，DOI 10.1016/j.cell.2026.10.002）。本文件自包含，供 Codex 实现设计脚本；输入为 `sg_test.csv`（events master，8,611 行 = 3,039 事件 / 1,411 基因）。

## 1. 靶区定义

- **靶区 = 外显化 TE 本体**：TE ∩ 外显子（`min(te_end, exon_end) − max(te_start, exon_start) + 1`），**不是** TE-CDS 重叠区。CDS overlap 只用于解释编码后果，不限定 guide 位置。
- 长度分布（实测）：中位 132 nt；≥60 nt 占 81%；<20 nt 仅 5.6%。
- **必须用 SpRY-Cas9（HF1-SpRY）**：CDS overlap 内 30% 事件无任何 NGG，SpCas9 无法均匀覆盖；SpRY 近无 PAM 限制，任何 ≥20 nt 区域都可设计。SpRY 单 guide 效率方差大 → 靠每区域多 guide + MAGeCK RRA 聚合弥补。

## 2. Guide 选择规则（按 TE 本体长度分级）

| TE 本体长度 | 占比 | 设计 |
|---|---|---|
| ≥60 nt | 81% | 每链 2 条（2+2），窗口间隔分布，优先剪接 junction 近端 |
| 20–60 nt | 13% | 每链 1 条，置于剪接 junction 近端（破坏剪接 > 破坏外显子中部） |
| <20 nt | 5.6% | ① junction-spanning guide（横跨外显子-内含子边界，SpRY 可精确切 junction，且序列天然唯一）；② 不可行则标 `not_testable`，不硬凑 |

**唯一性过滤（硬门槛）**：事件间 20-mer 有 70% 质量共享（Alu 拷贝同源）。每条 guide 必须全基因组 mappability 过滤：完美或 ≤2 错配的额外位点 ≤2，否则丢弃（多切割 → p53/毒性污染增殖表型，验证时切乱基因组）。**跨 TE/宿主边界或跨剪接 junction 的 guide 天然唯一，最高优先级。**

**每区域 guide 数**：预算允许时 ≥60 nt 区域升到每链 3–4 条（RRA 聚合越强越好）；最终以库大小反推。

## 3. 筛选读出与 coverage（库按 10,000 条 guide 计）

两种读出**分开做、分开算 coverage**，互为补充：增殖读出抓适应度，CD69 读出抓激活/信号（急性窗口，不被增殖缺陷混淆）。

### 3.1 增殖 screen（MYC 原文式 growth screen）

| 参数 | 取值 | 依据 |
|---|---|---|
| 每 guide 覆盖 | **≥500×，推荐 1000×** | MYC 文维持 >1,000×/guide |
| 感染细胞数 | 10,000 × 1,000 = **1×10⁷** | |
| 起始细胞（MOI 0.3，感染率 ~26%） | **~4×10⁷**/重复（500× 时 2×10⁷） | 单 guide tiling 必须 MOI≤0.3 |
| 时间点 | day 0 / 7 / 14 | 原文 day 0/8/20；原代 CD8 用 14 天 |
| 供者 | ≥3 名健康供者 | 生物学重复 |
| 细胞来源 | LR-BC（1–5×10⁸/供者）可满足 | |

### 3.2 CD69 screen（FACS 分选读出）

| 参数 | 取值 | 依据 |
|---|---|---|
| 每 guide 每臂覆盖 | **≥500×** | 二元分选无动态范围累积，单次测量需高表征 |
| 分选臂 | CD69^hi 与 CD69^lo 两臂 + 未分选 input 对照 | MAGeCK 用 input 作基线 |
| 分选细胞总量 | 10,000 × 500 × 2 = **1×10⁷**（回收率 ~50% → 起始 3–5×10⁷） | |
| 刺激 | PMA/离子霉素 4–6 h（CD69 峰值 4–8 h） | 急性窗口 |
| MOI / 供者 | 同 3.1（MOI≤0.3，≥3 供者） | |

### 3.3 预算分配（10,000 guide 内）

- 靶事件：排序前 **300–600 事件** × 4–8 条 ≈ 2.4k–4.8k；
- 宿主基因对照：每基因组成型外显子 2+2 = 4 条（区分 TE 异构体特异 vs 基因水平表型）；
- IFNAR2-S TE exon guide：实验流程阳性对照（已发表位点）；
- NTC ≥200 条 + OTC（intergenic one-target）≥200 条；
- 剩余预算给 junction-spanning guide 与 <20 nt 事件的补救设计。

## 4. 输出与模式标记

- sgRNA 表必须含 `mode` 列：`single_tiling`（单 guide，indel 破坏位点）vs `paired_deletion`（双 guide 删整个 TE exon，需 tRNA-spacer 双 guide 载体）；MAGeCK 分别聚合，不混池分析。
- 每条 guide 输出：坐标（0-based half-open）、链、PAM、CutSite、mappability（全基因组命中数）、RRA 预期聚合组、所属 event_id/gene、设计层级（A=标准 2+2 / B=junction 近端 / C=junction-spanning / D=not_testable）。
- 短/重复导致无法设计的**事件保留在表中**并标状态，不删行。

## 5. 边界

- 本文件是设计规范，不含已设计的 guide 序列；实际设计需全基因组参考（GRCh38 + rmsk）做 mappability，运行时验证 SpRY-Cas9 试剂可得性。
- CD69 读出不能推断增殖表型，反之亦然；两 screen 命中交集优先验证。
