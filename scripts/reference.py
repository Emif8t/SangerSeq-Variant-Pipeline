# =========================================================
# STEP 2 — DOWNLOAD REFSEQ TRANSCRIPT
# =========================================================

print("\nDownloading RefSeq transcript...")

with Entrez.efetch(
    db="nucleotide",
    id=REFSEQ_ID,
    rettype="fasta",
    retmode="text"
) as handle:

    record = SeqIO.read(handle, "fasta")

ref_seq = str(record.seq).upper()

print("Downloaded:", record.id)
print("Transcript length:", len(ref_seq))

# Extract the PCR amplicon

ref_amplicon = ref_seq[AMPLICON_START-1:AMPLICON_END]

print("Amplicon length:", len(ref_amplicon))
