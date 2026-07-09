# =========================================================
# STEP 1 — PROJECT SETTINGS
# =========================================================

# Folder containing ABI files
ab1_folder = r"C:\Users\USER\Desktop\data_analysis\25022026\ASS1\abfiles"

# NCBI requires an email address
Entrez.email = "your_email@example.com"

# Reference transcript
REFSEQ_ID = "NM_000050.4"

# PCR amplicon coordinates (RefSeq cDNA)
AMPLICON_START = 1007
AMPLICON_END   = 1402

print("ASS1 Sanger Pipeline")
print("Reference:", REFSEQ_ID)
