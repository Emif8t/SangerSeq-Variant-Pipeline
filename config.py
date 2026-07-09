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

AMPLICON_START = 1007
AMPLICON_END = 1402

# ======================================================
# QUALITY SETTINGS
# ======================================================

MIN_PHRED = 20

# ======================================================
# PCR PRIMERS
# ======================================================

FORWARD_PRIMER = "CAACACCCCTGACATTCTCG"

REVERSE_PRIMER = "ACTTTCCCTTCCACTCGCTC"
