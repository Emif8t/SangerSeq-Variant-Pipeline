# SangerSeq-Variant-Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22313559.svg)](https://doi.org/10.5281/zenodo.22313559)

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

**SangerSeq-Variant-Pipeline** is an end-to-end Python workflow for the analysis of targeted Sanger sequencing data.

The pipeline automates key stages of Sanger sequencing analysis, from ABI chromatogram processing and sequence alignment through variant identification, genotype determination, HGVS nomenclature generation, functional annotation, quality control, association analysis, and Hardy-Weinberg equilibrium testing.

The workflow includes:

- ABI chromatogram processing
- Sequencing-read preparation
- Sequencing quality control
- Reference transcript retrieval
- PCR primer verification
- Local sequence alignment
- Alignment processing
- Genotype calling
- High-confidence variant filtering
- Variant summarisation
- HGVS cDNA nomenclature generation
- Functional annotation using Ensembl Variant Effect Predictor (VEP)
- Final annotated variant-table construction
- Carrier identification
- Variant quality control
- Case-control association analysis
- Exact Hardy-Weinberg equilibrium analysis
- Final variant-table validation

The pipeline was developed as part of doctoral research investigating **ASS1 genetic variation in prostate cancer** and provides a reproducible framework for targeted Sanger sequencing analysis.

Although the workflow is designed as a general framework for targeted Sanger sequencing analysis, the current published configuration is **ASS1-focused**.

---

## DOI and Versioning

The SangerSeq Variant Pipeline is archived on Zenodo to provide a
persistent, citable record of each released software version.

The version 1.0.0 release has been archived and assigned the following
DOI:

**DOI:** https://doi.org/10.5281/zenodo.22313559

**Version:** 1.0.0

The version-specific Zenodo record preserves the exact software release
associated with this version. Future releases will be assigned their own
version-specific Zenodo records through the GitHub-Zenodo integration.

# Pipeline Workflow

The complete workflow consists of 16 major stages:

```text
Sanger chromatograms (.ab1)
            |
            v
1. Load ABI files
            |
            v
2. Prepare sequencing reads
            |
            v
3. Quality control
            |
            v
4. Reference retrieval
            |
            v
5. Primer verification
            |
            v
6. Local sequence alignment
            |
            v
7. Alignment processing
            |
            v
8. Genotype calling
            |
            v
9. High-confidence variant filtering
            |
            v
10. Variant summarisation
            |
            v
11. HGVS nomenclature generation
            |
            v
12. VEP functional annotation
            |
            v
13. Final annotated variant table
            |
            v
14. Variant quality control
            |
            v
15. Case-control association analysis
            |
            v
16. Hardy-Weinberg equilibrium analysis
            |
            v
Final results

The complete workflow is orchestrated by main.py.

Key Features
End-to-end analysis

The pipeline integrates multiple stages of targeted Sanger sequencing analysis into a single reproducible workflow.

Automated variant processing

Detected sequence differences are processed through genotype calling, high-confidence filtering, variant summarisation, and HGVS nomenclature generation.

Functional annotation

Variants can be functionally annotated using the Ensembl Variant Effect Predictor (VEP).

The current implementation supports annotation through the Ensembl REST API.

Carrier identification

The pipeline identifies sequencing samples carrying each confirmed variant and calculates carrier counts and variant frequencies.

Independent final-table validation

A dedicated validation script checks the consistency of the final variant table against the genotype and HGVS tables.

Case-control association analysis

The pipeline supports association analysis between identified variants and sample phenotype/group information.

Hardy-Weinberg equilibrium analysis

Hardy-Weinberg equilibrium analysis is performed using the control-group genotype information.

Reproducibility

Important analysis parameters are centralised in config.py, allowing the reference transcript, primer sequences, quality threshold, and annotation method to be explicitly documented.

Requirements
Software requirements

The pipeline requires:

Python 3.10 or later
Internet access for reference-sequence retrieval from NCBI
Internet access when using the Ensembl VEP REST API
A computer capable of running Python and the required scientific-computing packages
Python dependencies

The required Python packages are specified in:

requirements.txt

The pipeline uses packages including:

Biopython
NumPy
pandas
SciPy
Matplotlib
statsmodels
requests
openpyxl
Installation
Step 1: Clone the repository

Clone the repository from GitHub:
git clone https://github.com/Emif8t/SangerSeq-Variant-Pipeline.git 

Move into the project directory:

cd SangerSeq-Variant-Pipeline
Step 2: Create a virtual environment

Creating a virtual environment is recommended to keep the pipeline dependencies isolated.

Windows
python -m venv .venv

Activate the environment:

.venv\Scripts\activate
macOS/Linux
python3 -m venv .venv

Activate the environment:

source .venv/bin/activate
Step 3: Install dependencies

Install the required packages:

pip install -r requirements.txt
Repository Structure

The repository is organised as follows:

SangerSeq-Variant-Pipeline/
│
├── config.py
├── main.py
├── requirements.txt
├── utils.py
├── LICENSE
├── README.md
├── repository_structure.txt
│
└── scripts/
    ├── alignment.py
    ├── annotation.py
    ├── association.py
    ├── final_table.py
    ├── genotype.py
    ├── hgvs.py
    ├── hwe.py
    ├── preprocessing.py
    ├── reference.py
    ├── validate_final_table.py
    └── validation.py

Local analysis directories may include:

data/
├── raw/
├── metadata/
├── annotation/
└── reference/

output/
├── association/
├── genotypes/
├── hwe/
├── qc/
├── results/
├── variants/
└── annotation/

The data/raw/ directory is intended for local sequencing files and should not contain data that is committed to the public repository.

Generated analysis outputs are also intended to remain local unless explicitly required for a reproducible release.

Input Data

Before running the pipeline, users need to prepare the required input files.

Sanger chromatogram files

ABI chromatogram files (.ab1) should be placed in:

data/raw/

For example:

data/raw/
├── Sample01.ab1
├── Sample02.ab1
└── Sample03.ab1

The pipeline automatically searches the configured raw-data directory for ABI chromatograms.

Sample Group Information

For case-control association analysis and Hardy-Weinberg equilibrium analysis, the pipeline requires sample group information.

The expected metadata file is:

data/metadata/Sample_Groups.xlsx

The file should contain sample identifiers and their corresponding phenotype/group classifications.

The sample identifiers must correspond to the sequencing samples used by the pipeline.

For example:

Sample	Group
A5.ab1	Case
A7.ab1	Case
A9.ab1	Control

The exact column structure should remain consistent with the implementation of the association-analysis module.

VEP Annotation

The pipeline uses the Ensembl Variant Effect Predictor (VEP) for functional annotation.

The annotation method is controlled through config.py.

Ensembl REST API

The current implementation supports VEP annotation through the Ensembl REST API.

Set:

ANNOTATION_METHOD = "api"

When this method is selected, the annotation module sends HGVS variants to the Ensembl VEP REST endpoint and retrieves annotation information.

The resulting annotation table is generated locally in the output directory.

For example:

output/annotation/Ensembl_VEP_Annotation.csv

The VEP annotation can include information such as:

Variant class
Allele string
Consequence
Impact
Gene
Gene ID
Gene symbol
HGNC ID
Transcript
Feature
Biotype
Exon
Canonical transcript status
MANE information
RefSeq transcript
HGVSc
HGVSp
Protein ID
Protein position
CDS position
cDNA position
Variant allele
HTTP status

For synonymous variants, HGVSp may contain protein-level notation such as:

ENSP00000253004.6:p.His261=

indicating that the amino acid remains unchanged.

Configuration

Before running the pipeline, open:

config.py

and review the configuration parameters.

## NCBI Email Configuration

The pipeline requires an email address when communicating with NCBI
services.

For security and reproducibility, the email address is not stored directly
in `config.py`. Instead, provide it through the `NCBI_EMAIL` environment
variable.

### Windows PowerShell

```powershell
$env:NCBI_EMAIL="your_real_email@example.com"
```

### macOS/Linux

```bash
export NCBI_EMAIL="your_real_email@example.com"
```

The pipeline will report an error if NCBI_EMAIL has not been configured.

Reference Transcript

The current configuration targets the ASS1 transcript:

REFSEQ_ID = "NM_000050.4"
TRANSCRIPT = "NM_000050.4"

The coding-sequence start is configured as:

CDS_START = 357

These parameters are important for HGVS cDNA nomenclature generation and downstream variant interpretation.

If adapting the pipeline to another gene or transcript, the reference transcript, coding-sequence coordinates, primers, and other gene-specific parameters must be reviewed and appropriately changed.

PCR Primers

The current configuration contains the following PCR primer sequences:

FORWARD_PRIMER = "CAACACCCCTGACATTCTCG"
REVERSE_PRIMER = "ACTTTCCCTTCCACTCGCTC"

The pipeline uses these sequences to verify the expected target amplicon against the retrieved reference sequence.

If analysing another target, these primers must be replaced with the appropriate primers for that target.

Primer sequences should be provided in the expected 5' to 3' orientation.

Sequencing Quality

The minimum Phred-quality threshold is configured as:

MIN_PHRED = 20

This threshold is used during quality assessment and downstream genotype and variant processing.

Users should review this value according to the quality requirements of their sequencing experiment.

Running the Pipeline

Once:

Python is installed
Dependencies are installed
ABI files are placed in data/raw/
Sample-group metadata is prepared
config.py has been reviewed
Internet access is available for NCBI and/or Ensembl services as required

the complete pipeline can be run from the project root directory.

Run:

python main.py

The user does not need to manually execute the individual pipeline modules.

main.py orchestrates the complete workflow from chromatogram loading through final Hardy-Weinberg equilibrium analysis.

Pipeline Steps

When the pipeline is executed, it performs the following stages.

Step 1 — Load ABI chromatograms

ABI sequencing files are loaded from the configured raw-data directory.

data/raw/
Step 2 — Prepare sequencing reads

Sequencing reads are extracted and prepared for downstream analysis.

Step 3 — Quality control

Sequencing quality metrics are calculated and saved to the output directory.

Step 4 — Reference retrieval

The configured RefSeq transcript is retrieved from NCBI.

Step 5 — Primer verification

The forward and reverse primers are checked against the reference sequence to verify the expected target region.

Step 6 — Local sequence alignment

Sequencing reads are aligned to the expected reference sequence.

Step 7 — Alignment processing

The alignments are examined at the nucleotide level to support downstream genotype and variant calling.

Step 8 — Genotype determination

Genotype calls are generated from the sequencing data.

Step 9 — High-confidence variant filtering

Variant calls are filtered according to the configured quality criteria.

Step 10 — Variant summarisation

Detected variants are summarised by position and sequence change.

Step 11 — HGVS nomenclature

Detected variants are converted into HGVS cDNA nomenclature.

Step 12 — Functional annotation

Variants are functionally annotated using Ensembl VEP.

Step 13 — Final annotated variant table

Variant, genotype, carrier, HGVS, and VEP information are combined into the final annotated variant table.

Step 14 — Variant quality control

Quality-control reports are generated for the final variant table.

Step 15 — Association analysis

Case-control association analysis is performed using the supplied sample-group information.

Step 16 — Hardy-Weinberg equilibrium

Hardy-Weinberg equilibrium analysis is performed using control-group genotype information.

Output Files

After successful execution, results are written to the output/ directory.

The principal output structure includes:
output/
│
├── association/
│   └── ASS1_Association_Analysis.csv
│
├── annotation/
│   └── Ensembl_VEP_Annotation.csv
│
├── genotypes/
│   └── Genotype_Table.csv
│
├── hwe/
│   └── ASS1_Exact_HWE_Wigginton2.csv
│
├── qc/
│   ├── Consequence_Summary.csv
│   ├── Final_Table_QC_Copy.csv
│   ├── Impact_Summary.csv
│   ├── Variant_QC_Report.csv
│   └── Variant_Type_Summary.csv
│
├── results/
│   └── Final_Annotated_Variants.csv
│
└── variants/
    ├── HGVS_Table.csv
    ├── HighConfidence_Variants.csv
    └── Variant_Summary.csv

The exact set of generated files may depend on the pipeline configuration and analysis results.

Important Output Files
Genotype table
output/genotypes/Genotype_Table.csv

Contains genotype calls generated from the sequencing reads.

This table is used for downstream carrier identification, variant-frequency calculation, association analysis, and HWE analysis.

HGVS table
output/variants/HGVS_Table.csv

Contains the generated HGVS nomenclature for detected variants.

VEP annotation table
output/annotation/Ensembl_VEP_Annotation.csv

Contains functional annotation retrieved from Ensembl VEP.

Final annotated variant table
output/results/Final_Annotated_Variants.csv

This is the principal final variant table.

It combines variant information with HGVS nomenclature, VEP annotation, genotype-derived carrier information, and variant frequencies.

The table can contain information including:

HGVS cDNA notation
Reference allele
Alternate allele
Variant type
Consequence
Impact
Most severe consequence
Protein HGVS
Protein position
Amino acid change
Codon change
Gene
Gene ID
Transcript
Feature
Assembly
Chromosome
Genomic position
HGVS genomic notation
dbSNP identifiers
Transcript position
cDNA position
CDS position
Carrier count
Variant frequency
Carrier sample identities
VEP HGVSc
VEP HGVSp
Protein ID
MANE information
RefSeq transcript information
Final Variant Validation

The repository includes a dedicated validation script:

scripts/validate_final_table.py

After running the main pipeline, the final table can be independently validated using:

python scripts\validate_final_table.py

The validation checks include:

Number of final variants
Required final-table columns
HGVS variant identities
HGVS completeness
Reference alleles
Alternate alleles
Carrier counts
Carrier sample identities
Variant frequencies
Required genotype-table columns
Genotype carrier consistency
VEP annotation completeness
VEP consequence
VEP impact
VEP gene information
VEP feature information
HGVSc completeness
HGVSp completeness
Numeric field validation

The validation also cross-checks carrier identities between the final variant table and the genotype table.

A successful validation should end with:

FINAL VALIDATION RESULT

OVERALL RESULT: PASS

The final variant table passed all generic validation checks.

VALIDATION COMPLETE

This provides an additional quality-control step between pipeline execution and downstream interpretation.

Reproducibility

For reproducible analyses, users should retain:

The version of the pipeline used
The input ABI files
Sample-group metadata
The configuration used
Reference transcript/accession
Primer sequences
Sequencing quality thresholds
VEP annotation source
VEP annotation date/version where applicable
Python version
Installed package versions
Generated output files

The pipeline configuration is central to reproducibility because parameters such as the reference transcript, primer sequences, quality threshold, coding-sequence coordinates, and annotation method can influence downstream results.

Adapting the Pipeline to Another Gene

The current implementation is configured for ASS1.

To adapt the workflow to another targeted gene, users should review at minimum:

config.py

and the following parameters:

REFSEQ_ID
TRANSCRIPT
CDS_START
FORWARD_PRIMER
REVERSE_PRIMER

The expected input data and any gene-specific assumptions within the analysis modules should also be reviewed.

The reference transcript, coding-sequence coordinates, primer sequences, and annotation settings must correspond to the target gene.

Therefore, although the workflow provides a general framework for targeted Sanger sequencing analysis, the current published configuration should be regarded as an ASS1-focused implementation.

Troubleshooting
No ABI files found

If the pipeline reports that no ABI chromatograms were loaded, check that .ab1 files are located in:

data/raw/

Also ensure that the pipeline is being executed from the repository root:

SangerSeq-Variant-Pipeline/
Reference sequence cannot be retrieved

Check:

Internet connectivity
REFSEQ_ID
NCBI_EMAIL

in config.py.

Primer verification fails

Check that:

FORWARD_PRIMER
REVERSE_PRIMER

correspond to the primers used for the target amplicon and are written in the expected 5' to 3' orientation.

VEP annotation fails

If VEP annotation fails when using the REST API, check:

Internet connectivity
The configured annotation method
The Ensembl REST service availability
The HGVS nomenclature supplied to VEP
The HTTP status reported in the annotation output

The pipeline records the HTTP status associated with VEP requests.

Final validation fails

If:

OVERALL RESULT: FAIL

do not immediately modify the final CSV manually.

Instead, inspect:

output/genotypes/Genotype_Table.csv
output/variants/HGVS_Table.csv
output/annotation/Ensembl_VEP_Annotation.csv
output/results/Final_Annotated_Variants.csv

and review the validation message to determine which component is inconsistent.

Scientific Use

This pipeline is intended as a research tool for targeted Sanger sequencing analysis.

It can support applications including:

Candidate-gene studies
Cancer genomics
Molecular biomarker research
Genetic variation studies
Case-control sequencing studies
Targeted variant discovery
Genotype-phenotype analysis

The pipeline should not be considered a clinical diagnostic system without appropriate clinical validation, regulatory assessment, laboratory quality assurance, and validation according to applicable requirements.

## Citation

If you use the SangerSeq Variant Pipeline in research, please cite the
specific software version used.

### Version 1.0.0

> Israel, E. (2026). *SangerSeq Variant Pipeline* (Version 1.0.0).
> Zenodo. https://doi.org/10.5281/zenodo.22313559

**DOI:** https://doi.org/10.5281/zenodo.22313559

The repository also includes a `CITATION.cff` file containing
machine-readable citation metadata for use with GitHub and other
research software citation systems.

License

This project is distributed under the MIT License.

See:

LICENSE

for the full licence text.

Acknowledgements

The pipeline integrates publicly available resources and open-source scientific software, including:

National Center for Biotechnology Information (NCBI)
Ensembl Variant Effect Predictor (VEP)
Biopython
NumPy
pandas
SciPy
statsmodels
Matplotlib
Contact

Emmanuel Israel

PhD Fellow, CApIC-ACE

GitHub:

https://github.com/Emif8t