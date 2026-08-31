# SangerSeq-Variant-Pipeline

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

**SangerSeq-Variant-Pipeline** is an end-to-end Python workflow for the analysis of targeted Sanger sequencing data.

The pipeline automates key stages of Sanger sequencing analysis, including chromatogram processing, reference sequence handling, sequence alignment, variant detection, genotype determination, HGVS nomenclature generation, functional annotation, quality control, statistical analysis, and final validation.

The pipeline was developed as part of doctoral research investigating **ASS1 genetic variation in prostate cancer** and provides a reproducible framework for targeted Sanger sequencing studies.

---

## Key Features

- Import of ABI Sanger chromatogram (`.ab1`) files
- Reference sequence retrieval and preparation
- Primer and reference sequence verification
- Local sequence alignment
- Variant detection
- Genotype determination
- HGVS cDNA nomenclature generation
- Functional annotation using the Ensembl Variant Effect Predictor (VEP)
- Variant quality assessment
- Carrier identification
- Case-control association analysis
- Exact Hardy-Weinberg equilibrium testing
- Generation of summary tables
- Final variant-table validation
- Reproducible analysis workflow

---

## Pipeline Workflow

The general workflow is:

```text
Sanger chromatograms (.ab1)
            |
            v
     Preprocessing
            |
            v
   Reference preparation
            |
            v
     Sequence alignment
            |
            v
      Variant calling
            |
            v
   Genotype determination
            |
            v
    HGVS nomenclature
            |
            v
     VEP annotation
            |
            v
     Quality control
            |
            +----------------------+
            |                      |
            v                      v
   Association analysis       HWE analysis
            |                      |
            +----------+-----------+
                       |
                       v
              Final variant table
                       |
                       v
                Final validation

Repository Structure
SangerSeq-Variant-Pipeline/
|
├── data/
│   ├── annotation/
│   ├── metadata/
│   ├── raw/
│   └── reference/
|
├── output/
│   ├── association/
│   ├── genotypes/
│   ├── hwe/
│   ├── qc/
│   ├── results/
│   └── variants/
|
├── scripts/
│   ├── alignment.py
│   ├── annotation.py
│   ├── association.py
│   ├── final_table.py
│   ├── genotype.py
│   ├── hgvs.py
│   ├── hwe.py
│   ├── preprocessing.py
│   ├── reference.py
│   ├── validate_final_table.py
│   └── validation.py
|
├── .gitignore
├── config.py
├── LICENSE
├── main.py
├── README.md
├── repository_structure.txt
├── requirements.txt
└── utils.py

Raw sequencing files and generated analysis outputs are excluded from version control where appropriate to protect data integrity, reduce repository size, and support reproducible use with independent datasets.

Requirements

Python 3.10 or later is recommended.

The pipeline uses scientific and bioinformatics Python packages including:

- Biopython
- pandas
- NumPy
- SciPy
- Matplotlib
- Requests
- statsmodels
- openpyxl

Install the required dependencies with:

pip install -r requirements.txt
Running the Pipeline

The main workflow can be initiated using:

python main.py

Individual analysis modules can also be executed directly when required.

For example:

python scripts/preprocessing.py
python scripts/reference.py
python scripts/alignment.py
python scripts/genotype.py
python scripts/hgvs.py
python scripts/annotation.py
python scripts/validation.py
python scripts/association.py
python scripts/hwe.py
python scripts/final_table.py
python scripts/validate_final_table.py

The exact execution order may depend on the configuration and available input files.

Final Variant Table

The pipeline produces a final annotated variant table containing information such as:

- HGVS cDNA nomenclature
- Reference allele
- Alternate allele
- Variant consequence
- Variant impact
- Gene symbol
- Gene information
- Transcript/feature information
- Carrier count
- Variant frequency
- Carrier sample identities

The final table is generated at:

output/results/ASS1_Final_Annotated_Variants.csv
Final Validation

The pipeline includes an independent validation step for the final variant table.

Run:

python scripts/validate_final_table.py

The validation procedure checks:

- Number of final variants
- HGVS variant identities
- Reference and alternate alleles
- Carrier counts
- Carrier sample identities
- Variant frequencies
- VEP annotation completeness
- Variant consequences and impact
- Consistency with the genotype table

A successful validation produces:

OVERALL RESULT: PASS

This provides an additional quality-control step before downstream interpretation or reporting.

Example Validated ASS1 Variants

The current validated analysis identified two ASS1 variants:

| HGVS cDNA            | REF | ALT | Consequence        | Impact | Carrier Count | Variant Frequency |
| -------------------- | --- | --- | ------------------ | ------ | ------------: | ----------------: |
| NM_000050.4:c.783T>C | T   | C   | synonymous_variant | LOW    |             3 |          1.000000 |
| NM_000050.4:c.876T>C | T   | C   | synonymous_variant | LOW    |             2 |          0.666667 |


The carrier identities are determined directly from the confirmed genotype calls rather than from the annotation table.

Analysis Outputs

The pipeline can generate:

- Quality-control summaries
- Genotype tables
- Variant summary tables
- HGVS nomenclature tables
- VEP annotations
- Variant quality reports
- Case-control association results
- Hardy-Weinberg equilibrium results
- Final annotated variant tables
- Reproducibility

The pipeline is designed to support reproducible targeted Sanger sequencing analysis.

Input data, configuration parameters, intermediate files, and generated results can be maintained separately, allowing the workflow to be applied to additional datasets without modifying the core analysis framework.

Raw sequencing data are not included in the public repository.

Applications

The workflow can be adapted for:

- Candidate gene studies
- Cancer genomics
- Molecular biomarker research
- Clinical genetics research
- Molecular diagnostics
- Rare disease studies
- Targeted sequencing studies
- Genetic association studies

Citation

If you use this pipeline in your research, please cite:

Israel, E. (2026). SangerSeq-Variant-Pipeline. Version 1.0.

A DOI will be added following the first public release through Zenodo.

License

This project is distributed under the MIT License.

See the LICENSE file for details.

Acknowledgements

The pipeline integrates publicly available resources and open-source software from:

National Center for Biotechnology Information (NCBI)
Ensembl Variant Effect Predictor (VEP)
Biopython
NumPy
pandas
SciPy
statsmodels
Matplotlib

Author
Emmanuel Israel
PhD fellow in CApIC-ACE

GitHub: https://github.com/Emif8t