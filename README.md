# PVT1 / miR-152-3p in Colorectal Cancer — Manuscript & Computational Archive

结直肠癌中长链非编码 RNA **PVT1** 上调并通过 ceRNA 机制吸附 **miR-152-3p** 的诊断与机制研究。


> Corresponding author: FU Zexian (付泽娴), School of Medicine, Hebei University of Engineering, Handan, China.
> Email: fuzexian@hebeu.edu.cn

---

## 1. 研究概要 / Summary

- **主线假设（可证伪链）**：PVT1 上调 → ceRNA 吸附 miR-152-3p → "PVT1-high / miR-152-3p-low" 失衡 → 早期高判别力标志物。
- **核心结果**：
  - TCGA 中 PVT1 显著上调（log2FC 2.64；肿瘤-正常 AUC **0.995**）。
  - 三个 GEO/GPL570 队列（统一采用唯一具判别力的探针 **1558290_a_at**）AUC 分别为 **0.987 / 0.994 / 1.000**；四队列随机效应合并 AUC **0.989**（95% CI 0.982–0.993，I²=0%）。
  - 双队列证据轴：GSE39582 PVT1 AUC **0.990** + GSE38389 miR-152-3p AUC **0.740**。
  - PVT1 非泛炎症标志：溃疡性结肠炎（UC）未发炎黏膜 AUC 0.518、发炎黏膜 0.640；进展为肿瘤者合并 AUC 0.697。
- **限制**：所有准确性均为病例-对照（肿瘤组织 vs 非恶性黏膜），筛查效用未确立；阵列估计仅依赖单探针 1558290_a_at，需 qPCR 非探针依赖法确认。

---

## 2. 仓库结构 / Repository structure

```
TG/
├── manuscript/
│   └── 
├── figures/                                         
│   ├──
│   ├──
│   ├──
│   ├──
│   ├── 
│   └── 
├── results/                                          # 73 个处理结果表（CSV/JSON），支撑全部图表与结论
│   ├── AUC_summary.csv                  # 各队列 PVT1 诊断 AUC（单探针 1558290_a_at）
│   ├── N3_meta_result.json              # 诊断 meta 合并 AUC / CI / Q / I² / τ²
│   ├── N3_stage_gradient.csv            # 分期梯度 AUC
│   ├── N3_screening_pp.csv              # 筛查 PPV/特异度投影
│   ├── affy_probe_validity.csv          # 探针效度（Supplementary Table S3 数据源）
│   ├── ibd_*.csv / ibd_meta.json        # IBD 对照队列结果（Supplementary Table S2）
│   ├── *_raw.csv                        # GSE39582/9348/22598 处理后表达矩阵
│   ├── L1000CDS2_*.csv                  # 药物重定位（Supplementary Table S1 数据源）
│   ├── depmap_dependency_CRC.csv        # DepMap 虚拟敲除
│   ├── survival_*.csv, enrichment_*.csv, immune_*.csv, ...
│   └── ...（其余为各分析模块的中间/汇总结果）
├── code/                                             # 45 个分析脚本（Python + R）
│   ├── 
│   ├── Task12_N3_diagnostic_upgrade.py  # 诊断 meta（含连续性校正 + Hanley-McNeil 方差）
│   ├── Task17_probe_fix.py              # GPL570 单探针 1558290_a_at 重写
│   ├── Task14_make_fig8_fig9.py         # 诊断 meta 图
│   ├── Task16_IBD_cohorts.py            # IBD 对照队列
│   ├── analyze_*.py / *.R, Task*_*.py, download_supplement*.R, gdc_*.py ...
│   └── ...
├── docs/                                             # 数据来源与纳入说明
│   ├── IBD队列结果_纳入说明_2026-10-30.md
│   ├── GEO_下载清单_癌前对照.md
│   └── GEO_下载清单_良性癌前与IBD.md
├── requirements.txt
├── .gitignore
└── LICENSE
```

---

## 3. 稿件与补充材料 / Manuscript & Supplementary Information



| 编号 | 内容 | 数据源（位于 `results/`） |
|------|------|---------------------------|
| Fig. S1 | L1000CDS2 药物重定位（逆转/模拟 PVT1-high 特征） | `Fig_L1000CDS2.png`, `L1000CDS2_reverse_summary.csv`, `L1000CDS2_mimic_summary.csv`, `depmap_dependency_CRC.csv`, `binding_sites_effectors.csv` |
| Table S1 | 候选效应子支持证据 | `depmap_dependency_CRC.csv`, `binding_sites_effectors.csv` |
| Table S2 | IBD 对照队列特征与判别 | `ibd_comparisons.csv`, `ibd_meta.json`, `ibd_probe_qc.csv` |
| Table S3 | Affymetrix PVT1 探针效度（5 探针逐一评估） | `affy_probe_validity.csv` |


---

## 4. 数据来源（外部数据库，未纳入本仓库）/ Data sources

原始数据体量大（约 1 GB：TCGA RNA-seq、GEO 原始矩阵、DepMap CRISPR 矩阵），**未纳入本仓库**，请按下列 accession 自行下载：

| 来源 | Accession / 平台 | 用途 |
|------|------------------|------|
| GEO | GSE39582 (GPL570) | PVT1 诊断（主队列，分期） |
| GEO | GSE9348 (GPL570) | PVT1 诊断（早癌） |
| GEO | GSE22598 (GPL570) | PVT1 配对肿瘤-正常；5-aza 去甲基化探索 |
| GEO | GSE38389 | miR-152-3p 配对组织方向（外验证） |
| GEO | GSE11223 (GPL1708) | UC 炎症对照 |
| GEO | GSE236055 (GPL21185) | UC 进展者 vs 非进展者 |
| GEO | GSE151911 (GPL20844) | UC 进展者 vs 非进展者 |
| GEO | GSE38713 (GPL570) | 经预设探针 QC 排除（专一 PVT1 探针处背景水平） |
| TCGA (GDC) | COAD + READ | RNA-seq / miRNA / 临床（PVT1 上调、生存、相关性） |
| DepMap / CCLE | CRISPRGeneDependency / CRISPRGeneEffect | PVT1 虚拟敲除、候选效应子依赖性 |
| L1000CDS2 | — | 药物重定位（逆转 PVT1-high 特征） |

---

## 5. 复现说明 / Reproduction

所有脚本**已改为相对路径**：根目录 `ROOT` 由各脚本自身位置自动推算
（`TG/code/某脚本.py` → 上级目录即 `TG/`），不再硬编码 `D:/AIwork/YYG`。
子目录映射：`data/`、`figures/`、`results/`（原 `ROC_results/`，已统一更名）。
脚本在任意位置运行均能正确定位文件，**无需手动改路径**。

### 5.1 原始数据（需自备才能从头复现）
`data/` 未纳入仓库（约 1 GB）。请按第 4 节 accession 下载原始数据，并放入 `TG/data/` 以下布局：
- `data/TCGA_COAD_RNA.gz`、`data/TCGA_READ_RNA.gz`、`data/TCGA_*_clinical.tsv`、`data/results_tcga.json`（GDC 下载脚本产出）
- `data/GEOcache/<GSE>_series_matrix.txt.gz`（经 `download_supplement*.R` 下载）
- `data/GEOcache/platform_tables/`（经 `fetch_platform_table.py` 下载）
- `data/raw/`、`data/processed/`（经 `Task1_download_datasets.R` 整理）
- `data/DepMap_CRISPRGeneEffect.csv` 等（DepMap 官方下载）

### 5.2 开箱即用的部分
`figures/` 与 `results/` 已包含全部终版配图与处理结果。只读取这两类产物、不依赖原始数据的脚本可直接运行；
需要原始 `data/` 的脚本（如 `gen_draft.py`、`Task12`、`Task16`、`Task17` 及多数 R 下游）只有在 `data/` 就位后才能完整复现。

### 5.3 运行示例
```bash
cd TG
python code/Task12_N3_diagnostic_upgrade.py   # 诊断 meta（读 results，写 results）
python code/Task16_IBD_cohorts.py             # IBD 队列（需 data/GEOcache）
```
Python 依赖见 `requirements.txt`；R 脚本依赖 `limma`、`edgeR`、`Biobase` 等 Bioconductor 包。

---

## 6. 本仓库已排除的文件 / Excluded

- `data/`（约 1 GB 原始矩阵/依赖矩阵）—— 体积超限，按 accession 下载。

---

## 7. 上传 GitHub 快速步骤 / Upload to GitHub

```bash
cd TG
git init
git add .
git commit -m "PVT1/miR-152-3p CRC manuscript M4 + computational archive"
git branch -M main
git remote add origin <your-repo-url>
git push -u origin main
```

---

## 8. 许可 / License

- 稿件、配图与结果数据表：**CC BY 4.0**。
- 分析代码：**MIT License**。
详见 `LICENSE`。
"# -lncRNA-PVT1-CRC" 
[![DOI](https://zenodo.org/badge/1439913245.svg)](https://doi.org/10.5281/zenodo.23812014)
