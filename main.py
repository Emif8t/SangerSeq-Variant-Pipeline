"""
main.py

Main entry point for the SangerSeq Variant Pipeline.
"""

from scripts.reference import (

    download_reference,

    extract_amplicon,

    verify_primers

)

from config import *

from scripts.preprocessing import (

    load_ab1_files,

    prepare_reads,

    calculate_qc_metrics,

    save_qc_summary

)

from scripts.alignment import (

    perform_local_alignment,

    walk_alignment

)

from scripts.hgvs import (

    generate_hgvs_table,

    save_hgvs_table

)

from scripts.annotation import (

    annotate_variants,

    save_annotation_table

)

def main():

    # --------------------------------------------
    # Load ABI chromatograms
    # --------------------------------------------

    abi_records = load_ab1_files(

        AB1_FOLDER

    )

    # --------------------------------------------
    # Prepare sequencing reads
    # --------------------------------------------

    processed_reads = prepare_reads(

        abi_records

    )

    # --------------------------------------------
    # Quality-control summary
    # --------------------------------------------

    qc_summary = calculate_qc_metrics(

        processed_reads

    )

    save_qc_summary(

        qc_summary,

        OUTPUT_FOLDER

    )

    print(

        f"Pipeline completed successfully.\n"

        f"Loaded {len(abi_records)} ABI chromatograms."

    )


if __name__ == "__main__":

    main()

# --------------------------------------------
# Download reference transcript
# --------------------------------------------

reference_sequence = download_reference(

    REFSEQ_ID,

    NCBI_EMAIL

)

# --------------------------------------------
# Extract reference amplicon
# --------------------------------------------

reference_amplicon = extract_amplicon(

    reference_sequence,

    AMPLICON_START,

    AMPLICON_END

)

# --------------------------------------------
# Verify PCR primers
# --------------------------------------------

primer_information = verify_primers(

    reference_sequence,

    FORWARD_PRIMER,

    REVERSE_PRIMER

)

# --------------------------------------------
# Local sequence alignment
# --------------------------------------------

alignments = perform_local_alignment(

    processed_reads,

    primer_information["amplicon_sequence"]

)

# --------------------------------------------
# Walk alignment
# --------------------------------------------

alignment_df = walk_alignment(

    alignments,

    reference_sequence,

    primer_information["amplicon_start"],

    MIN_PHRED

)

# --------------------------------------------
# HGVS nomenclature
# --------------------------------------------

hgvs_df = generate_hgvs_table(

    summary_df,

    variant_df,

    TRANSCRIPT,

    CDS_START

)

save_hgvs_table(

    hgvs_df,

    os.path.join(

        OUTPUT_FOLDER,

        "variants"

    )

)

# --------------------------------------------
# Functional annotation
# --------------------------------------------

annotation_df = annotate_variants(

    hgvs_df,

    ENSEMBL_SERVER,

    ENSEMBL_HEADERS,

    REQUEST_TIMEOUT,

    MAX_RETRIES,

    REQUEST_DELAY

)

save_annotation_table(

    annotation_df,

    os.path.join(

        OUTPUT_FOLDER,

        "annotation"

    )

)
