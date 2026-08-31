"""
main.py

Main entry point for the SangerSeq Variant Pipeline.

Pipeline workflow
-----------------
1. Load ABI chromatograms
2. Prepare sequencing reads
3. Perform quality control
4. Download reference transcript
5. Verify PCR primers
6. Perform local sequence alignment
7. Walk alignments
8. Call genotypes
9. Filter high-confidence variants
10. Summarise variants
11. Generate HGVS nomenclature
12. Annotate variants
13. Build final publication-ready table
14. Perform variant quality control
15. Perform association analysis
16. Perform Hardy-Weinberg equilibrium analysis
"""

import os

from config import (
    AB1_FOLDER,
    OUTPUT_FOLDER,
    NCBI_EMAIL,
    REFSEQ_ID,
    MIN_PHRED,
    FORWARD_PRIMER,
    REVERSE_PRIMER,
    TRANSCRIPT,
    CDS_START,
    ENSEMBL_SERVER,
    ENSEMBL_HEADERS,
    REQUEST_TIMEOUT,
    MAX_RETRIES,
    REQUEST_DELAY,
    ANNOTATION_METHOD,
    VEP_OUTPUT_FILE,
    PHENOTYPE_FILE
)

# =========================================================
# PREPROCESSING
# =========================================================

from scripts.preprocessing import (
    load_ab1_files,
    prepare_reads,
    calculate_qc_metrics,
    save_qc_summary
)

# =========================================================
# REFERENCE
# =========================================================

from scripts.reference import (
    download_reference,
    verify_primers
)

# =========================================================
# ALIGNMENT
# =========================================================

from scripts.alignment import (
    perform_local_alignment,
    walk_alignment
)

# =========================================================
# GENOTYPE
# =========================================================

from scripts.genotype import (
    call_genotypes,
    filter_high_confidence_variants,
    summarize_variants,
    save_genotypes,
    save_variants,
    save_variant_summary
)

# =========================================================
# HGVS
# =========================================================

from scripts.hgvs import (
    generate_hgvs_table,
    save_hgvs_table
)

# =========================================================
# ANNOTATION
# =========================================================

from scripts.annotation import (
    annotate_variants,
    save_annotation_table
)

# =========================================================
# FINAL TABLE
# =========================================================

from scripts.final_table import (
    build_final_variant_table,
    save_final_table
)

# =========================================================
# VALIDATION
# =========================================================

from scripts.validation import (
    run_variant_quality_control,
    save_qc_reports
)

# =========================================================
# ASSOCIATION
# =========================================================

from scripts.association import (
    load_sample_groups,
    create_sample_sets,
    run_association_analysis,
    apply_multiple_testing,
    save_association_results
)

# =========================================================
# HWE
# =========================================================

from scripts.hwe import (
    prepare_control_genotypes,
    run_hwe_analysis,
    save_hwe_results
)


# =========================================================
# MAIN PIPELINE
# =========================================================

def main():

    print()
    print("=" * 70)
    print("SangerSeq Variant Pipeline")
    print("=" * 70)
    print()

    # =====================================================
    # 1. LOAD ABI CHROMATOGRAMS
    # =====================================================

    print("[1/16] Loading ABI chromatograms...")

    abi_records = load_ab1_files(
        AB1_FOLDER
    )

    if not abi_records:
        raise ValueError(
            "No ABI chromatograms were loaded. "
            f"Please check the input directory: {AB1_FOLDER}"
        )

    print(
        f"      Loaded {len(abi_records)} ABI chromatograms."
    )

    # =====================================================
    # 2. PREPARE SEQUENCING READS
    # =====================================================

    print()
    print("[2/16] Preparing sequencing reads...")

    processed_reads = prepare_reads(
        abi_records
    )

    if not processed_reads:
        raise ValueError(
            "No sequencing reads were successfully prepared."
        )

    print(
        f"      Prepared {len(processed_reads)} sequencing reads."
    )

    # =====================================================
    # 3. QUALITY CONTROL
    # =====================================================

    print()
    print("[3/16] Calculating sequencing quality metrics...")

    qc_summary = calculate_qc_metrics(
        processed_reads
    )

    save_qc_summary(
        qc_summary,
        OUTPUT_FOLDER
    )

    print(
        "      Quality-control summary saved."
    )

    # =====================================================
    # 4. DOWNLOAD REFERENCE TRANSCRIPT
    # =====================================================

    print()
    print("[4/16] Loading reference transcript...")

    reference_sequence = download_reference(
        REFSEQ_ID,
        NCBI_EMAIL
    )

    if not reference_sequence:
        raise ValueError(
            "Reference sequence could not be retrieved."
        )

    print(
        f"      Reference transcript: {REFSEQ_ID}"
    )

    # =====================================================
    # 5. VERIFY PCR PRIMERS
    # =====================================================

    print()
    print("[5/16] Verifying PCR primers...")

    primer_information = verify_primers(
        reference_sequence,
        FORWARD_PRIMER,
        REVERSE_PRIMER
    )

    if not primer_information:
        raise ValueError(
            "PCR primer verification failed."
        )

    print(
        "      Forward primer found at: "
        f"{primer_information['forward_start']}"
    )

    print(
        "      Reverse primer found at: "
        f"{primer_information['reverse_start']}"
    )

    print(
        "      Amplicon coordinates: "
        f"{primer_information['amplicon_start']}-"
        f"{primer_information['amplicon_end']}"
    )

    # =====================================================
    # 6. LOCAL SEQUENCE ALIGNMENT
    # =====================================================

    print()
    print("[6/16] Performing local sequence alignment...")

    alignments = perform_local_alignment(
        processed_reads,
        primer_information["amplicon_sequence"]
    )

    if not alignments:
        raise ValueError(
            "No sequencing reads could be aligned "
            "to the reference amplicon."
        )

    print(
        f"      Successfully aligned {len(alignments)} reads."
    )

    # =====================================================
    # 7. WALK ALIGNMENT
    # =====================================================

    print()
    print("[7/16] Generating nucleotide-level alignment table...")

    alignment_df = walk_alignment(
        alignments,
        reference_sequence,
        primer_information["amplicon_start"],
        MIN_PHRED
    )

    if alignment_df.empty:
        raise ValueError(
            "The alignment table is empty."
        )

    print(
        f"      Generated {len(alignment_df)} alignment records."
    )

    # =====================================================
    # 8. CALL GENOTYPES
    # =====================================================

    print()
    print("[8/16] Calling genotypes...")

    genotype_df = call_genotypes(
        alignment_df,
        MIN_PHRED
    )

    if genotype_df.empty:
        raise ValueError(
            "No genotype calls were generated."
        )

    save_genotypes(
        genotype_df,
        os.path.join(
            OUTPUT_FOLDER,
            "genotypes"
        )
    )

    print(
        f"      Generated {len(genotype_df)} genotype records."
    )

    # =====================================================
    # 9. FILTER HIGH-CONFIDENCE VARIANTS
    # =====================================================

    print()
    print("[9/16] Filtering high-confidence variants...")

    high_confidence_df, variant_df = (
        filter_high_confidence_variants(
            genotype_df,
            MIN_PHRED
        )
    )

    save_variants(
        variant_df,
        os.path.join(
            OUTPUT_FOLDER,
            "variants"
        )
    )

    print(
        f"      High-confidence genotype calls: "
        f"{len(high_confidence_df)}"
    )

    print(
        f"      High-confidence variant calls: "
        f"{len(variant_df)}"
    )

    # =====================================================
    # 10. SUMMARISE VARIANTS
    # =====================================================

    print()
    print("[10/16] Summarising variants...")

    if variant_df.empty:

        print(
            "      No variants detected."
        )

        summary_df = summarize_variants(
            variant_df
        )

    else:

        summary_df = summarize_variants(
            variant_df
        )

    save_variant_summary(
        summary_df,
        os.path.join(
            OUTPUT_FOLDER,
            "variants"
        )
    )

    print(
        f"      Variant summary contains "
        f"{len(summary_df)} variant positions."
    )

    # =====================================================
    # 11. HGVS NOMENCLATURE
    # =====================================================

    print()
    print("[11/16] Generating HGVS nomenclature...")

    if summary_df.empty:

        print(
            "      No variants available for HGVS generation."
        )

        hgvs_df = generate_hgvs_table(
            summary_df,
            variant_df,
            TRANSCRIPT,
            CDS_START
        )

    else:

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

    print(
        f"      Generated {len(hgvs_df)} HGVS records."
    )

    # =====================================================
    # 12. FUNCTIONAL ANNOTATION
    # =====================================================

    print()
    print("[12/16] Performing functional annotation...")

    if ANNOTATION_METHOD.lower() == "api":

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

    elif ANNOTATION_METHOD.lower() == "web":

        if not os.path.exists(VEP_OUTPUT_FILE):

            raise FileNotFoundError(
                "ANNOTATION_METHOD is set to 'web', but the "
                f"VEP output file was not found:\n"
                f"{VEP_OUTPUT_FILE}"
            )

        print(
            "      Using existing VEP web annotation:"
        )

        print(
            f"      {VEP_OUTPUT_FILE}"
        )

        annotation_df = None

    else:

        raise ValueError(
            "ANNOTATION_METHOD must be either "
            "'api' or 'web'."
        )

    # =====================================================
    # 13. BUILD FINAL PUBLICATION TABLE
    # =====================================================

    print()
    print("[13/16] Building final annotated variant table...")

    final_df = build_final_variant_table(
        hgvs_df,
        genotype_df,
        VEP_OUTPUT_FILE
    )

    final_output_folder = os.path.join(
        OUTPUT_FOLDER,
        "results"
    )

    final_output_file = save_final_table(
        final_df,
        final_output_folder
    )

    print(
        "      Final variant table saved to:"
    )

    print(
        f"      {final_output_file}"
    )

    # =====================================================
    # 14. VARIANT QUALITY CONTROL
    # =====================================================

    print()
    print("[14/16] Performing variant quality control...")

    (
        qc_report,
        consequence_summary,
        impact_summary,
        variant_type_summary
    ) = run_variant_quality_control(
        final_df
    )

    save_qc_reports(
        qc_report,
        consequence_summary,
        impact_summary,
        variant_type_summary,
        final_df,
        os.path.join(
            OUTPUT_FOLDER,
            "qc"
        )
    )

    print(
        "      Variant quality-control reports saved."
    )

    # =====================================================
    # 15. ASSOCIATION ANALYSIS
    # =====================================================

    print()
    print("[15/16] Performing association analysis...")

    phenotype_df = load_sample_groups(
        PHENOTYPE_FILE
    )

    case_samples, control_samples = create_sample_sets(
        phenotype_df
    )

    print(
        f"      Case samples: {len(case_samples)}"
    )

    print(
        f"      Control samples: {len(control_samples)}"
    )

    association_df = run_association_analysis(
        final_df,
        case_samples,
        control_samples
    )

    association_df = apply_multiple_testing(
        association_df
    )

    save_association_results(
        association_df,
        os.path.join(
            OUTPUT_FOLDER,
            "association"
        )
    )

    print(
        "      Association analysis completed."
    )

    # =====================================================
    # 16. HARDY-WEINBERG EQUILIBRIUM
    # =====================================================

    print()
    print("[16/16] Performing Hardy-Weinberg equilibrium analysis...")

    control_df = prepare_control_genotypes(
        genotype_df,
        control_samples
    )

    hwe_df = run_hwe_analysis(
        final_df,
        control_df
    )

    save_hwe_results(
        hwe_df,
        os.path.join(
            OUTPUT_FOLDER,
            "hwe"
        )
    )

    print(
        "      HWE analysis completed."
    )

    # =====================================================
    # COMPLETE
    # =====================================================

    print()
    print("=" * 70)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print()
    print(
        f"ABI chromatograms loaded : {len(abi_records)}"
    )
    print(
        f"Prepared reads            : {len(processed_reads)}"
    )
    print(
        f"Alignment records         : {len(alignment_df)}"
    )
    print(
        f"Genotype records          : {len(genotype_df)}"
    )
    print(
        f"High-confidence variants : {len(variant_df)}"
    )
    print(
        f"Variant positions        : {len(summary_df)}"
    )
    print(
        f"HGVS records              : {len(hgvs_df)}"
    )
    print(
        f"Final table records       : {len(final_df)}"
    )
    print()
    print(
        f"Final output:"
    )
    print(
        f"{final_output_file}"
    )
    print()


# =========================================================
# SCRIPT ENTRY POINT
# =========================================================

if __name__ == "__main__":
    main()