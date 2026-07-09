"""
Configuration file for the SangerSeq Variant Pipeline.

Edit the parameters below before running the pipeline.
"""

# ======================================================
# INPUT / OUTPUT DIRECTORIES
# ======================================================

AB1_FOLDER = "data/raw"
OUTPUT_FOLDER = "output"

# ======================================================
# NCBI SETTINGS
# ======================================================

NCBI_EMAIL = "your_email@example.com"

# ======================================================
# REFERENCE TRANSCRIPT
# ======================================================

REFSEQ_ID = "NM_000050.4"

# PCR amplicon coordinates (RefSeq cDNA)

#AMPLICON_START = 1007
#AMPLICON_END = 1402

# ======================================================
# QUALITY SETTINGS
# ======================================================

MIN_PHRED = 20

# ======================================================
# PCR PRIMERS
# ======================================================

FORWARD_PRIMER = "CAACACCCCTGACATTCTCG"

REVERSE_PRIMER = "ACTTTCCCTTCCACTCGCTC"


# ======================================================
# HGVS SETTINGS
# ======================================================

TRANSCRIPT = "NM_000050.4"

# Coding sequence starts at transcript position

CDS_START = 357

# ======================================================
# ENSEMBL VEP
# ======================================================

ENSEMBL_SERVER = "https://rest.ensembl.org"

ENSEMBL_HEADERS = {
    "Accept": "application/json"
}

REQUEST_TIMEOUT = 30

MAX_RETRIES = 2

REQUEST_DELAY = 0.05
