"""
Configuration file for the SangerSeq Variant Pipeline.

Edit the parameters below before running the pipeline.
"""

import os


# ============================================================
# INPUT / OUTPUT DIRECTORIES
# ============================================================

DATA_FOLDER = "data"

AB1_FOLDER = os.path.join(
    DATA_FOLDER,
    "raw"
)

OUTPUT_FOLDER = "output"


# ============================================================
# NCBI SETTINGS
# ============================================================

# Email address used when communicating with NCBI services.
# Set the NCBI_EMAIL environment variable before running the pipeline.
#
# The email address is intentionally not stored in the repository.

NCBI_EMAIL = os.getenv("NCBI_EMAIL")


# ============================================================
# REFERENCE TRANSCRIPT
# ============================================================

REFSEQ_ID = "NM_000050.4"


# ============================================================
# PCR AMPLICON
# ============================================================

# Optional explicit amplicon coordinates.
#
# Leave commented if the pipeline should determine the
# amplicon from the PCR primers.

# AMPLICON_START = 1007
# AMPLICON_END = 1402


# ============================================================
# SEQUENCING QUALITY
# ============================================================

MIN_PHRED = 20


# ============================================================
# PCR PRIMERS
# ============================================================

FORWARD_PRIMER = "CAACACCCCTGACATTCTCG"

REVERSE_PRIMER = "ACTTTCCCTTCCACTCGCTC"


# ============================================================
# HGVS SETTINGS
# ============================================================

TRANSCRIPT = REFSEQ_ID

# Coding sequence starts at transcript position 357.
CDS_START = 357


# ============================================================
# ENSEMBL VEP
# ============================================================

ENSEMBL_SERVER = "https://rest.ensembl.org"

ENSEMBL_HEADERS = {
    "Accept": "application/json"
}

REQUEST_TIMEOUT = 30

MAX_RETRIES = 3

REQUEST_DELAY = 0.1


# ============================================================
# ANNOTATION METHOD
# ============================================================

# Default annotation method.
#
# "api" = automatically annotate variants using the
#         Ensembl VEP REST API.
#
# "web" = use an existing VEP web-exported Excel file.
#
# API is recommended as the default because it makes the
# pipeline more reproducible and avoids requiring the user
# to manually upload variants to the VEP website.

ANNOTATION_METHOD = "api"


# ============================================================
# VEP WEB OUTPUT
# ============================================================

# This file is only required when:
#
#     ANNOTATION_METHOD = "web"
#
# It is not required for normal API-based operation.

VEP_OUTPUT_FILE = os.path.join(
    DATA_FOLDER,
    "annotation",
    "VEP_HGVS_OUTPUT.xlsx"
)


# ============================================================
# PHENOTYPE / SAMPLE GROUP INFORMATION
# ============================================================

PHENOTYPE_FILE = os.path.join(
    DATA_FOLDER,
    "metadata",
    "Sample_Groups.xlsx"
)