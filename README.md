# gene-signature-generalization
Code, frozen gene lists (S5–S200), and evaluation pipeline for: Signature Size, Feature Stability, and Cross-Cohort Transportability in Breast Cancer.
# Gene Signature Size and Cross-Cohort Transportability in Breast Cancer

Official repository for:
**Signature Size, Feature Stability, and Cross-Cohort Transportability in Breast Cancer**  
*Mohammed Rozil (2026)*

## Overview
This repository contains the complete execution pipeline, frozen gene lists ($S_5$–$S_{200}$), and evaluation scripts to reproduce the finding that larger nested transcriptomic gene signatures increase the discovery-to-external generalization gap without improving mean external discrimination.

## Public Data Sources
The pipeline processes the following publicly available cohorts:
- **TCGA-BRCA** (Discovery, $n=1,093$, Illumina HiSeq RNA-seq) — [GDC Portal](https://portal.gdc.cancer.gov)
- **GSE20685** (External, $n=327$, Agilent Microarray) — [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20685)
- **GSE74314** (External, $n=220$, NextSeq RNA-seq) — [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74314)
- **GSE136661** (External, $n=214$, NovaSeq RNA-seq) — [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136661)
- **METABRIC** (External, $n=1,980$, Agilent Microarray) — [Synapse syn1688369](https://www.synapse.org/#!Synapse:syn1688369)

## Quickstart & Reproducibility

```bash
# 1. Clone repository
git clone [https://github.com/yourusername/gene-signature-generalization.git](https://github.com/yourusername/gene-signature-generalization.git)
cd gene-signature-generalization

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run nested cross-validation and external evaluation
python scripts/02_nested_cv_elasticnet.py
