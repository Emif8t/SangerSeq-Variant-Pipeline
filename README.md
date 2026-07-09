# SangerSeq-Variant-Pipeline

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)]()

## Overview

The **ASS1 Sanger Variant Pipeline** is an end-to-end Python workflow for the analysis of targeted Sanger sequencing data. The pipeline automates chromatogram processing, sequence alignment, variant detection, HGVS nomenclature generation, functional annotation, quality validation, case–control association analysis, Hardy–Weinberg equilibrium testing, and publication-quality visualization of genetic variants.

The pipeline was developed as part of a PhD study investigating **ASS1 genetic variation in prostate cancer** and is designed to provide a reproducible framework for targeted Sanger sequencing analysis.

---

## Features

✔ Automatic import of ABI chromatogram (.ab1) files

✔ Retrieval of reference sequences from NCBI RefSeq

✔ PCR primer verification

✔ Local sequence alignment

✔ Variant calling

✔ Genotype determination

✔ HGVS cDNA nomenclature generation

✔ Functional annotation using Ensembl Variant Effect Predictor (VEP)

✔ Variant quality validation

✔ Case–control association analysis

✔ Exact Hardy–Weinberg equilibrium testing

✔ Publication-quality mutation plots

✔ Export of publication-ready tables

---

## Pipeline Workflow

```
ABI Chromatograms (.ab1)
            │
            ▼
Manual Chromatogram Inspection
            │
            ▼
Sequence Extraction
            │
            ▼
Reference Retrieval (NCBI RefSeq)
            │
            ▼
Primer Verification
            │
            ▼
Local Pairwise Alignment
            │
            ▼
Variant Calling
            │
            ▼
Genotype Determination
            │
            ▼
HGVS Nomenclature
            │
            ▼
Ensembl VEP Annotation
            │
            ▼
Variant Quality Validation
            │
            ▼
Association Analysis
            │
            ▼
Hardy–Weinberg Equilibrium
            │
            ▼
Publication Figures & Tables
```

---

## Repository Structure

```text
ASS1-Sanger-Variant-Pipeline/
│
├── scripts/
│   ├── Step01_Project_Settings.py
│   ├── Step02_Download_Reference.py
│   ├── ...
│   ├── Step20_Mutation_Plots.py
│
├── data/
│
├── results/
│
├── figures/
│
├── examples/
│
├── requirements.txt
│
├── LICENSE
│
└── README.md
```

---

## Requirements

Python 3.10 or later

Required packages include:

- Biopython
- pandas
- numpy
- scipy
- matplotlib
- requests
- statsmodels
- openpyxl

Install all dependencies using:

```bash
pip install -r requirements.txt
```

---

## Running the Pipeline

Run the scripts sequentially.

Example:

```bash
python Step01_Project_Settings.py

python Step02_Download_Reference.py

...

python Step20_Mutation_Plots.py
```

---

## Outputs

The pipeline generates:

- Quality assessment reports
- Variant tables
- HGVS nomenclature tables
- Ensembl VEP annotations
- Case–control association tables
- Hardy–Weinberg equilibrium tables
- Mutation spectrum plots
- Lollipop plots
- Publication-ready figures

---

## Applications

The workflow can be adapted for:

- Candidate gene studies
- Clinical genetics
- Molecular diagnostics
- Rare disease studies
- Cancer genomics
- Biomarker discovery

---

## Citation

If you use this pipeline in your research, please cite:

> Israel, E. (2026). ASS1 Sanger Variant Pipeline. Version 1.0.

*A Zenodo DOI will be added after the first public release.*

---

## License

This project is distributed under the MIT License.

---

## Acknowledgements

The pipeline was developed during doctoral research in Biochemistry investigating molecular biomarkers for prostate cancer. It integrates publicly available resources from:

- National Center for Biotechnology Information (NCBI)
- Ensembl Variant Effect Predictor (VEP)
- Biopython Project
- SciPy
- Statsmodels
- Matplotlib

---

## Contact

**Emmanuel Israel**

PhD Candidate in Biochemistry

GitHub: *(add your GitHub profile link)*
