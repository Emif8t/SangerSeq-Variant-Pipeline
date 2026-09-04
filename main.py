"""
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

import pandas as pd

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
    sort_association_results,
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
# Utility functions
# =========================================================

def _require_nonempty(
    dataframe,
    name,
    context
):
    """
    Stop the pipeline when a required DataFrame
    is unexpectedly empty.
    """

    if dataframe is None:
        raise ValueError(
            f"{context}: {name} is None. "
            "The previous pipeline step did not return "
            "a valid DataFrame."
        )

    if not isinstance(dataframe, pd.DataFrame):
        raise TypeError(
            f"{context}: {name} must be a pandas DataFrame, "
            f"got {type(dataframe).__name__}."
        )

    if dataframe.empty:
        raise ValueError(
            f"{context}: {name} is empty. "
            "The pipeline cannot continue because no "
            "records are available."
        )


def _validate_phenotype_data(
    phenotype_df
):
    """
    Validate the phenotype/sample-group table.

    The pipeline is intentionally sample-size agnostic.
    Therefore, this function validates the required
    structure rather than expecting a specific number
    of samples.
    """

    required_columns = {
        "Sample",
        "Group"
    }

    missing_columns = (
        required_columns
        - set(phenotype_df.columns)
    )

    if missing_columns:
        raise ValueError(
            "Phenotype file is missing required columns: "
            f"{sorted(missing_columns)}"
        )

    if phenotype_df.empty:
        raise ValueError(
            "The phenotype/sample-group file is empty."
        )

    phenotype_df = phenotype_df.copy()

    phenotype_df["Sample"] = (
        phenotype_df["Sample"]
        .astype(str)
        .str.strip()
    )

    phenotype_df["Group"] = (
        phenotype_df["Group"]
        .astype(str)
        .str.strip()
    )

    if phenotype_df["Sample"].eq("").any():
        raise ValueError(
            "The phenotype file contains empty sample IDs."
        )

    valid_groups = {
        "Case",
        "Control"
    }

    observed_groups = set(
        phenotype_df["Group"].unique()
    )

    unexpected_groups = (
        observed_groups
        - valid_groups
    )

    if unexpected_groups:
        raise ValueError(
            "Unexpected phenotype groups detected: "
            f"{sorted(unexpected_groups)}. "
            "Expected only 'Case' and 'Control'."
        )

    duplicated_samples = (
        phenotype_df["Sample"]
        [phenotype_df["Sample"].duplicated()]
        .unique()
    )

    if len(duplicated_samples) > 0:
        raise ValueError(
            "Duplicate sample IDs detected in the "
            "phenotype file: "
            f"{sorted(duplicated_samples)}"
        )

    case_count = (
        phenotype_df["Group"]
        .eq("Case")
        .sum()
    )

    control_count = (
        phenotype_df["Group"]
        .eq("Control")
        .sum()
    )

    if case_count == 0:
        raise ValueError(
            "The phenotype file contains no Case samples."
        )

    if control_count == 0:
        raise ValueError(
            "The phenotype file contains no Control samples."
        )

    return phenotype_df


def _validate_hgvs_input(
    summary_df,
    variant_df
):
    """
    Validate the invariant between the variant summary
    and genotype calls.

    Every row in summary_df is generated from variant_df.
    Therefore, each summary variant must have at least
    one exact matching high-confidence genotype call.
    """

    required_summary = {
        "cDNA_Position",
        "REF",
        "ALT"
    }

    required_variants = {
        "cDNA_Position",
        "REF",
        "ALT",
        "Sample",
        "Quality",
        "Alignment_Score"
    }

    missing_summary = (
        required_summary
        - set(summary_df.columns)
    )

    missing_variants = (
        required_variants
        - set(variant_df.columns)
    )

    if missing_summary:
        raise ValueError(
            "HGVS validation failed: summary_df is missing "
            "columns: "
            f"{sorted(missing_summary)}"
        )

    if missing_variants:
        raise ValueError(
            "HGVS validation failed: variant_df is missing "
            "columns: "
            f"{sorted(missing_variants)}"
        )

    def normalise_series(series):
        return (
            series
            .astype(str)
            .str.strip()
            .str.upper()
        )

    summary_positions = pd.to_numeric(
        summary_df["cDNA_Position"],
        errors="coerce"
    )

    variant_positions = pd.to_numeric(
        variant_df["cDNA_Position"],
        errors="coerce"
    )

    summary_ref = normalise_series(
        summary_df["REF"]
    )

    summary_alt = normalise_series(
        summary_df["ALT"]
    )

    variant_ref = normalise_series(
        variant_df["REF"]
    )

    variant_alt = normalise_series(
        variant_df["ALT"]
    )

    unmatched = []

    for index in summary_df.index:

        position = summary_positions.loc[index]
        ref = summary_ref.loc[index]
        alt = summary_alt.loc[index]

        matches = (
            (variant_positions == position)
            & (variant_ref == ref)
            & (variant_alt == alt)
        )

        if not matches.any():

            if pd.notna(position):

                position_text = f"{position:g}"

            else:

                position_text = str(position)

            unmatched.append(
                f"{position_text} {ref}>{alt}"
            )

    if unmatched:
        raise ValueError(
            "HGVS validation failed: the following summary "
            "variants have no exact matching records in "
            "variant_df: "
            + ", ".join(unmatched)
            + ". This indicates an upstream "
            "data-integrity mismatch."
        )


def _print_variant_diagnostics(
    summary_df,
    variant_df
):
    """
    Print concise diagnostics for variant-to-HGVS matching.
    """

    print()
    print(
        "      Variant-to-HGVS diagnostic check"
    )
    print(
        "      --------------------------------"
    )

    if summary_df.empty:

        print(
            "      Variant summary: EMPTY"
        )

        return

    print(
        f"      Variant summary rows : "
        f"{len(summary_df)}"
    )

    print(
        f"      Variant call rows    : "
        f"{len(variant_df)}"
    )

    for _, row in summary_df.iterrows():

        position = row["cDNA_Position"]

        ref = (
            str(row["REF"])
            .strip()
            .upper()
        )

        alt = (
            str(row["ALT"])
            .strip()
            .upper()
        )

        positions = pd.to_numeric(
            variant_df["cDNA_Position"],
            errors="coerce"
        )

        matches = variant_df[
            (
                positions
                == pd.to_numeric(
                    position,
                    errors="coerce"
                )
            )
            & (
                variant_df["REF"]
                .astype(str)
                .str.strip()
                .str.upper()
                == ref
            )
            & (
                variant_df["ALT"]
                .astype(str)
                .str.strip()
                .str.upper()
                == alt
            )
        ]

        print(
            f"      {position} {ref}>{alt}: "
            f"{len(matches)} matching genotype call(s)"
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

    print(
        "[1/16] Loading ABI chromatograms..."
    )

    abi_records = load_ab1_files(
        AB1_FOLDER
    )

    if not abi_records:

        raise ValueError(
            "No ABI chromatograms were loaded. "
            f"Please check the input directory: "
            f"{AB1_FOLDER}"
        )

    print(
        f"      Loaded {len(abi_records)} "
        "ABI chromatograms."
    )

    # =====================================================
    # 2. PREPARE SEQUENCING READS
    # =====================================================

    print()
    print(
        "[2/16] Preparing sequencing reads..."
    )

    processed_reads = prepare_reads(
        abi_records
    )

    if not processed_reads:

        raise ValueError(
            "No sequencing reads were successfully prepared."
        )

    print(
        f"      Prepared {len(processed_reads)} "
        "sequencing reads."
    )

    # =====================================================
    # 3. QUALITY CONTROL
    # =====================================================

    print()
    print(
        "[3/16] Calculating sequencing quality metrics..."
    )

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
    print(
        "[4/16] Loading reference transcript..."
    )

    reference_sequence = download_reference(
        REFSEQ_ID,
        NCBI_EMAIL
    )

    if not reference_sequence:

        raise ValueError(
            "Reference sequence could not be retrieved."
        )

    print(
        f"      Reference transcript: "
        f"{REFSEQ_ID}"
    )

    # =====================================================
    # 5. VERIFY PCR PRIMERS
    # =====================================================

    print()
    print(
        "[5/16] Verifying PCR primers..."
    )

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
    print(
        "[6/16] Performing local sequence alignment..."
    )

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
        f"      Successfully aligned "
        f"{len(alignments)} reads."
    )

    # =====================================================
    # 7. WALK ALIGNMENT
    # =====================================================

    print()
    print(
        "[7/16] Generating nucleotide-level "
        "alignment table..."
    )

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
        f"      Generated {len(alignment_df)} "
        "alignment records."
    )

    # =====================================================
    # 8. CALL GENOTYPES
    # =====================================================

    print()
    print(
        "[8/16] Calling genotypes..."
    )

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
        f"      Generated {len(genotype_df)} "
        "genotype records."
    )

    # =====================================================
    # 9. FILTER HIGH-CONFIDENCE VARIANTS
    # =====================================================

    print()
    print(
        "[9/16] Filtering high-confidence variants..."
    )

    (
        high_confidence_df,
        variant_df
    ) = filter_high_confidence_variants(
        genotype_df,
        MIN_PHRED
    )

    print("\nDEBUG: confirmed variant_df columns:")
    print(variant_df.columns.tolist())

    print("\nDEBUG: confirmed variant_df:")
    print(variant_df.to_string(index=False))

    _require_nonempty(
        high_confidence_df,
        "high_confidence_df",
        "Step 9"
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
    print(
        "[10/16] Summarising variants..."
    )

    if variant_df.empty:

        raise ValueError(
            "Step 10: No high-confidence variants "
            "were detected. The pipeline cannot "
            "generate HGVS nomenclature or a final "
            "variant table. Review sequencing quality, "
            "alignment, genotype calls, and MIN_PHRED."
        )

    summary_df = summarize_variants(
        variant_df
    )

    _require_nonempty(
        summary_df,
        "summary_df",
        "Step 10"
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
    print(
        "[11/16] Generating HGVS nomenclature..."
    )

    _print_variant_diagnostics(
        summary_df,
        variant_df
    )

    _validate_hgvs_input(
        summary_df,
        variant_df
    )

    hgvs_df = generate_hgvs_table(
        summary_df,
        variant_df,
        TRANSCRIPT,
        CDS_START
    )

    _require_nonempty(
        hgvs_df,
        "hgvs_df",
        "Step 11"
    )

    save_hgvs_table(
        hgvs_df,
        os.path.join(
            OUTPUT_FOLDER,
            "variants"
        )
    )

    print(
        f"      Generated {len(hgvs_df)} "
        "HGVS records."
    )

    # =====================================================
    # 12. FUNCTIONAL ANNOTATION
    # =====================================================

    print()
    print(
        "[12/16] Performing functional annotation..."
    )

    if ANNOTATION_METHOD.lower() == "api":

        annotation_df = annotate_variants(
            hgvs_df,
            ENSEMBL_SERVER,
            ENSEMBL_HEADERS,
            REQUEST_TIMEOUT,
            MAX_RETRIES,
            REQUEST_DELAY,
            TRANSCRIPT,
            VEP_OUTPUT_FILE
        )

        save_annotation_table(
            annotation_df,
            os.path.join(
                OUTPUT_FOLDER,
                "annotation"
            )
        )

    elif ANNOTATION_METHOD.lower() == "web":

        if not os.path.exists(
            VEP_OUTPUT_FILE
        ):

            raise FileNotFoundError(
                "ANNOTATION_METHOD is set to 'web', "
                "but the VEP output file was not found:\n"
                f"{VEP_OUTPUT_FILE}"
            )

        print(
            "      Using existing VEP web annotation:"
        )

        print(
            f"      {VEP_OUTPUT_FILE}"
        )

        annotation_df = pd.read_excel(
            VEP_OUTPUT_FILE
        )

    else:

        raise ValueError(
            "ANNOTATION_METHOD must be either "
            "'api' or 'web'."
        )

    _require_nonempty(
        annotation_df,
        "annotation_df",
        "Step 12"
    )

    # =====================================================
    # 13. BUILD FINAL PUBLICATION TABLE
    # =====================================================

    print()
    print(
        "[13/16] Building final annotated variant table..."
    )

    final_df = build_final_variant_table(
        hgvs_df,
        genotype_df,
        annotation_df,
        confirmed_variant_df=variant_df,
        min_phred=MIN_PHRED,
    )

    _require_nonempty(
        final_df,
        "final_df",
        "Step 13"
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
    print(
        "[14/16] Performing variant quality control..."
    )

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
    print(
        "[15/16] Performing association analysis..."
    )

    phenotype_df = load_sample_groups(
        PHENOTYPE_FILE
    )

    phenotype_df = _validate_phenotype_data(
        phenotype_df
    )

    (
        case_samples,
        control_samples
    ) = create_sample_sets(
        phenotype_df
    )

    print(
        f"      Phenotype Case samples: "
        f"{len(case_samples)}"
    )

    print(
        f"      Phenotype Control samples: "
        f"{len(control_samples)}"
    )

    # -----------------------------------------------------
    # Association analysis uses actual genotype data.
    #
    # Samples without a callable genotype at a particular
    # variant are treated as missing rather than as
    # non-carriers. The association module therefore
    # determines the callable denominator separately for
    # each variant.
    # -----------------------------------------------------

    association_df = run_association_analysis(
        final_df,
        genotype_df,
        case_samples,
        control_samples,
        min_phred=MIN_PHRED
    )

    _require_nonempty(
        association_df,
        "association_df",
        "Step 15"
    )

    association_df = apply_multiple_testing(
        association_df
    )

    association_df = sort_association_results(
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
    print(
        "[16/16] Performing Hardy-Weinberg "
        "equilibrium analysis..."
    )

    # -----------------------------------------------------
    # Use only high-confidence genotype calls.
    #
    # HWE is evaluated in controls only. The HWE module
    # further restricts these calls to the supplied control
    # sample set.
    # -----------------------------------------------------

    control_df = prepare_control_genotypes(
        high_confidence_df,
        control_samples
    )

    hwe_df = run_hwe_analysis(
        final_df,
        control_df
    )

    _require_nonempty(
        hwe_df,
        "hwe_df",
        "Step 16"
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
        f"ABI chromatograms loaded : "
        f"{len(abi_records)}"
    )

    print(
        f"Prepared reads            : "
        f"{len(processed_reads)}"
    )

    print(
        f"Alignment records         : "
        f"{len(alignment_df)}"
    )

    print(
        f"Genotype records          : "
        f"{len(genotype_df)}"
    )

    print(
        f"High-confidence genotypes: "
        f"{len(high_confidence_df)}"
    )

    print(
        f"High-confidence variants : "
        f"{len(variant_df)}"
    )

    print(
        f"Variant positions        : "
        f"{len(summary_df)}"
    )

    print(
        f"HGVS records             : "
        f"{len(hgvs_df)}"
    )

    print(
        f"Final table records      : "
        f"{len(final_df)}"
    )

    print()
    print(
        "Final output:"
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