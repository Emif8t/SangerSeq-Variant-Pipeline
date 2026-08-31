# SangerSeq-Variant-Pipeline

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

**SangerSeq-Variant-Pipeline** is an end-to-end Python workflow for the analysis of targeted Sanger sequencing data.

The pipeline automates key stages of Sanger sequencing analysis, including:

- ABI chromatogram processing
- Sequencing-read preparation
- Quality control
- Reference transcript retrieval
- PCR primer verification
- Local sequence alignment
- Variant calling
- Genotype determination
- High-confidence variant filtering
- HGVS cDNA nomenclature generation
- Functional annotation using Ensembl Variant Effect Predictor (VEP)
- Variant quality control
- Carrier identification
- Case-control association analysis
- Exact Hardy-Weinberg equilibrium analysis
- Generation of summary tables
- Final annotated variant-table construction
- Final validation of variant and carrier information

The pipeline was developed as part of doctoral research investigating **ASS1 genetic variation in prostate cancer** and provides a reproducible framework for targeted Sanger sequencing analysis.

---

# Pipeline Workflow

The complete workflow consists of 16 major steps:

```text
Sanger chromatograms (.ab1)
            |
            v
      1. Load ABI files
            |
            v
      2. Prepare reads
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
     11. HGVS nomenclature
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
            +------------------------+
            |                        |
            v                        v
     15. Association analysis   16. HWE analysis
            |                        |
            +------------+-----------+
                         |
                         v
                 Final results

The complete workflow can be executed through main.py.

1. Requirements
Software requirements

The pipeline requires:

Python 3.10 or later
Internet access for reference-sequence retrieval from NCBI
Internet access if Ensembl VEP is accessed through the REST API
A computer capable of running Python and the required scientific-computing packages
Python dependencies

The required Python packages are specified in requirements.txt:

Biopython
NumPy
pandas
SciPy
Matplotlib
statsmodels
requests
openpyxl
2. Installation
Step 1: Clone the repository

Clone the repository from GitHub:

git clone https://github.com/Emif8t/SangerSeq-Variant-Pipeline.git

Move into the project directory:

cd SangerSeq-Variant-Pipeline
Step 2: Create a virtual environment

Creating a virtual environment is recommended to keep the pipeline dependencies isolated.

On Windows:

python -m venv .venv

Activate it:

.venv\Scripts\activate

On macOS/Linux:

python3 -m venv .venv

Activate it:

source .venv/bin/activate
Step 3: Install dependencies

Run:

pip install -r requirements.txt
3. Repository Structure

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
├── data/
│   ├── annotation/
│   ├── metadata/
│   ├── raw/
│   └── reference/
│
├── output/
│   ├── association/
│   ├── genotypes/
│   ├── hwe/
│   ├── qc/
│   ├── results/
│   └── variants/
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

The data/raw/ directory is intended for local sequencing files and is excluded from version control.

Generated analysis outputs are also intended to remain local rather than being committed to the public repository.

4. Input Data

Before running the pipeline, users need to prepare the required input files.

4.1 Sanger chromatogram files

ABI chromatogram files (.ab1) should be placed in:

data/raw/

For example:

data/raw/
├── Sample01.ab1
├── Sample02.ab1
├── Sample03.ab1
└── ...

The pipeline automatically searches this directory for ABI chromatograms.

5. Sample Group Information

For case-control analysis and Hardy-Weinberg equilibrium analysis, the pipeline requires sample group information.

The expected file is:

data/metadata/Sample_Groups.xlsx

This file should contain the sample identifiers and their corresponding phenotype/group classifications.

The sample identifiers must correspond to the names used by the sequencing files.

For example:

Sample     Group
A5.ab1     Case
A7.ab1     Case
A9.ab1     Control

The exact column structure should be maintained consistently with the implementation of the association-analysis module.

6. VEP Annotation

The pipeline supports two annotation approaches.

The annotation method is controlled through config.py.

ANNOTATION_METHOD = "web"
Option 1: Existing VEP output

When:

ANNOTATION_METHOD = "web"

the pipeline expects an existing VEP output file at:

data/annotation/VEP_HGVS_OUTPUT.xlsx

This approach is useful when variants have already been submitted to the Ensembl VEP web interface.

Place the resulting VEP annotation file in:

data/annotation/

with the expected filename:

VEP_HGVS_OUTPUT.xlsx
Option 2: Ensembl REST API

The pipeline can also use the Ensembl REST API.

Set:

ANNOTATION_METHOD = "api"

The annotation module will then query the configured Ensembl REST endpoint.

Internet access is required when using this option.

7. Configuration

Before running the pipeline, open:

config.py

and review the configuration parameters.

NCBI email

Provide a valid email address:

NCBI_EMAIL = "your_email@example.com"

Replace the placeholder with your own email address.

For example:

NCBI_EMAIL = "researcher@example.com"

This is used when communicating with NCBI services.

8. Reference Transcript

The pipeline currently uses:

REFSEQ_ID = "NM_000050.4"

and:

TRANSCRIPT = "NM_000050.4"

The current configuration therefore targets the ASS1 transcript NM_000050.4.

The coding-sequence start is configured as:

CDS_START = 357

These parameters are important for HGVS cDNA nomenclature generation.

If adapting the pipeline to another gene or transcript, the reference transcript, coding-sequence coordinates, primers, and other gene-specific parameters must be reviewed and appropriately changed.

9. PCR Primers

The current configuration contains the PCR primer sequences:

FORWARD_PRIMER = "CAACACCCCTGACATTCTCG"

REVERSE_PRIMER = "ACTTTCCCTTCCACTCGCTC"

The pipeline uses these sequences to verify the expected amplicon against the retrieved reference sequence.

If analysing another target, these primers must be replaced with the appropriate primers for that target.

10. Sequencing Quality

The minimum Phred-quality threshold is configured as:

MIN_PHRED = 20

This threshold is used during quality assessment and downstream genotype/variant processing.

Users should review this value according to the quality requirements of their sequencing experiment.

11. Running the Pipeline

Once:

Python is installed
Dependencies are installed
ABI files are placed in data/raw/
Sample-group information is prepared
VEP annotation is prepared when required
config.py has been reviewed

the complete pipeline can be run from the project root directory.

Run:

python main.py

You do not need to manually execute all 16 pipeline modules.

main.py orchestrates the complete workflow from chromatogram loading through final HWE analysis.

12. Pipeline Steps

When the pipeline is executed, it performs the following steps.

Step 1 — Load ABI chromatograms

ABI sequencing files are loaded from:

data/raw/
Step 2 — Prepare sequencing reads

Sequencing reads are extracted and prepared for downstream analysis.

Step 3 — Quality control

Quality metrics are calculated and saved to the output directory.

Step 4 — Reference retrieval

The configured RefSeq transcript is retrieved from NCBI.

Step 5 — Primer verification

The forward and reverse primers are checked against the reference sequence.

Step 6 — Local sequence alignment

Sequencing reads are aligned to the expected reference amplicon.

Step 7 — Alignment processing

The alignments are examined at the nucleotide level.

Step 8 — Genotype determination

Genotype calls are generated from the sequencing data.

Step 9 — High-confidence variant filtering

Variant calls are filtered according to the configured quality threshold.

Step 10 — Variant summarisation

Detected variants are summarised by position and sequence change.

Step 11 — HGVS nomenclature

Variants are converted into HGVS cDNA nomenclature.

Step 12 — Functional annotation

Variants are annotated using Ensembl VEP.

Step 13 — Final variant table

The final annotated variant table is constructed by combining the variant, genotype, and annotation information.

Step 14 — Variant quality control

Quality-control reports are generated for the final variant table.

Step 15 — Association analysis

Case-control association analysis is performed using the supplied sample-group information.

Step 16 — Hardy-Weinberg equilibrium

Hardy-Weinberg equilibrium analysis is performed using the control genotypes.

13. Output Files

After successful execution, results are written to the output/ directory.

The main output structure is:

output/
│
├── QC_Summary.csv
│
├── association/
│   └── ASS1_Association_Analysis.csv
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
│   └── ASS1_Final_Annotated_Variants.csv
│
└── variants/
    ├── HGVS_Table.csv
    ├── HighConfidence_Variants.csv
    └── Variant_Summary.csv
14. Important Output Files
Genotype table
output/genotypes/Genotype_Table.csv

Contains the genotype calls generated from the sequencing reads.

HGVS table
output/variants/HGVS_Table.csv

Contains the generated HGVS nomenclature for detected variants.

Final annotated variant table
output/results/ASS1_Final_Annotated_Variants.csv

This is the principal final variant table.

It combines variant information with VEP annotation and confirmed carrier information.

The table can contain information such as:

HGVS cDNA notation
Reference allele
Alternate allele
Consequence
IMPACT
Gene
Transcript/feature information
HGVS annotation
Carrier count
Variant frequency
Carrier sample identities
15. Final Variant Validation

The repository includes a dedicated validation script:

scripts/validate_final_table.py

After running the main pipeline, the final table can be independently validated using:

python scripts\validate_final_table.py

The validation checks include:

Number of final variants
HGVS variant identities
Reference alleles
Alternate alleles
Carrier counts
Carrier sample identities
Variant frequencies
VEP annotation completeness
VEP consequence
VEP IMPACT
Genotype-table consistency

A successful validation should end with:

OVERALL RESULT: PASS

This provides an additional quality-control step between pipeline execution and downstream interpretation.

16. Example Validation

For example, a successful validation may report:

FINAL VALIDATION RESULT

OVERALL RESULT: PASS

The final variant table contains the correct variants.
The carrier counts and carrier identities are correct.
The final table is consistent with the expected genotype results.

VALIDATION COMPLETE
17. Reproducibility

For reproducible analyses, users should retain:

The version of the pipeline used
The input ABI files
The sample-group metadata
The configuration used
The reference transcript/accession
The VEP annotation source and version/date
The Python version
The installed package versions
The generated output files

The pipeline's configuration is central to reproducibility because parameters such as the reference transcript, primer sequences, quality threshold, CDS start, and annotation method influence downstream analysis.

18. Adapting the Pipeline to Another Gene

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

Therefore, although the workflow provides a general framework for targeted Sanger sequencing analysis, the current published configuration should be regarded as an ASS1-focused implementation.

19. Troubleshooting
No ABI files found

If the pipeline reports that no ABI chromatograms were loaded, check that .ab1 files are located in:

data/raw/

and that the pipeline is being executed from the repository root:

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

VEP annotation file not found

If:

ANNOTATION_METHOD = "web"

make sure the VEP output file exists at:

data/annotation/VEP_HGVS_OUTPUT.xlsx

Alternatively, configure the pipeline to use the Ensembl REST API where appropriate:

ANNOTATION_METHOD = "api"
Validation fails

If:

OVERALL RESULT: FAIL

do not immediately modify the final CSV manually.

Instead, inspect:

output/genotypes/Genotype_Table.csv
output/variants/HGVS_Table.csv
output/results/ASS1_Final_Annotated_Variants.csv

and review the validation message to determine which component is inconsistent.

20. Scientific Use

This pipeline is intended as a research tool for targeted Sanger sequencing analysis.

It can support applications including:

Candidate-gene studies
Cancer genomics
Molecular biomarker research
Genetic variation studies
Case-control sequencing studies
Targeted variant discovery
Genotype-phenotype analysis

The pipeline should not be considered a clinical diagnostic system without appropriate clinical validation, regulatory assessment, and laboratory quality assurance.

21. Citation

If you use this pipeline in research, please cite:

Israel, E. (2026). SangerSeq-Variant-Pipeline. Version 1.0.

A DOI will be added following the first public release through Zenodo.

22. License

This project is distributed under the MIT License.

See:

LICENSE

for the full licence text.

23. Acknowledgements

The pipeline integrates publicly available resources and open-source scientific software, including:

National Center for Biotechnology Information (NCBI)
Ensembl Variant Effect Predictor (VEP)
Biopython
NumPy
pandas
SciPy
statsmodels
Matplotlib

24. Contact
Emmanuel Israel
PhD fellow at CApIC-ACE

GitHub:
https://github.com/Emif8t