"""
hgvs.py

Functions for generating HGVS nomenclature for
validated Sanger sequencing variants.

This module is gene-agnostic. The transcript accession
and CDS start position are supplied by the pipeline
configuration.
"""

import os

import pandas as pd


# ============================================================
# NORMALISE TRANSCRIPT
# ============================================================

def normalise_transcript(
    transcript: str
) -> str:
    """
    Clean a transcript accession.

    Parameters
    ----------
    transcript : str

    Returns
    -------
    str
    """

    if transcript is None:
        raise ValueError(
            "Transcript accession cannot be None."
        )

    transcript = str(
        transcript
    ).strip()

    if not transcript:
        raise ValueError(
            "Transcript accession cannot be empty."
        )

    return transcript


# ============================================================
# GENERATE HGVS TABLE
# ============================================================

def generate_hgvs_table(
    summary_df: pd.DataFrame,
    variant_df: pd.DataFrame,
    transcript: str,
    cds_start: int
) -> pd.DataFrame:
    """
    Generate HGVS cDNA nomenclature for validated
    sequence variants.

    One row is generated for each unique:

        transcript position
        REF
        ALT

    combination.

    Parameters
    ----------
    summary_df : pandas.DataFrame
        Variant-level summary.

    variant_df : pandas.DataFrame
        High-confidence variant calls.

    transcript : str
        RefSeq transcript accession, for example:

            NM_000050.4

    cds_start : int
        First CDS nucleotide position on the transcript,
        using 1-based transcript coordinates.

    Returns
    -------
    pandas.DataFrame
        HGVS variant table.
    """

    # ========================================================
    # VALIDATE INPUTS
    # ========================================================

    transcript = normalise_transcript(
        transcript
    )

    try:

        cds_start = int(cds_start)

    except (
        TypeError,
        ValueError
    ):

        raise ValueError(
            "cds_start must be an integer."
        )

    if cds_start < 1:

        raise ValueError(
            "cds_start must be >= 1."
        )

    required_summary_columns = [
        "cDNA_Position",
        "REF",
        "ALT",
    ]

    missing_summary = [
        column
        for column in required_summary_columns
        if column not in summary_df.columns
    ]

    if missing_summary:

        raise KeyError(
            "Variant summary is missing required columns:\n"
            + "\n".join(
                f"  - {column}"
                for column in missing_summary
            )
        )

    required_variant_columns = [
        "cDNA_Position",
        "REF",
        "ALT",
    ]

    missing_variant = [
        column
        for column in required_variant_columns
        if column not in variant_df.columns
    ]

    if missing_variant:

        raise KeyError(
            "Variant table is missing required columns:\n"
            + "\n".join(
                f"  - {column}"
                for column in missing_variant
            )
        )

    # ========================================================
    # HANDLE EMPTY INPUT
    # ========================================================

    if summary_df.empty:

        return pd.DataFrame(
            columns=[
                "Transcript",
                "Transcript_Position",
                "HGVS_cDNA_Position",
                "HGVS_cDNA",
                "REF",
                "ALT",
                "Variant_Type",
                "Carrier_Count",
                "Variant_Calls",
                "Mean_Quality",
                "Mean_Alignment_Score",
                "Samples",
            ]
        )

    # ========================================================
    # COPY DATA
    # ========================================================

    summary = summary_df.copy()

    variants = variant_df.copy()

    # ========================================================
    # NORMALISE DATA TYPES
    # ========================================================

    summary["cDNA_Position"] = pd.to_numeric(
        summary["cDNA_Position"],
        errors="coerce"
    )

    variants["cDNA_Position"] = pd.to_numeric(
        variants["cDNA_Position"],
        errors="coerce"
    )

    # Remove invalid positions

    summary = summary[
        summary["cDNA_Position"].notna()
    ].copy()

    variants = variants[
        variants["cDNA_Position"].notna()
    ].copy()

    summary["cDNA_Position"] = (
        summary["cDNA_Position"]
        .astype(int)
    )

    variants["cDNA_Position"] = (
        variants["cDNA_Position"]
        .astype(int)
    )

    # ========================================================
    # NORMALISE REF / ALT
    # ========================================================

    for df in [
        summary,
        variants,
    ]:

        df["REF"] = (
            df["REF"]
            .astype(str)
            .str.strip()
            .str.upper()
        )

        df["ALT"] = (
            df["ALT"]
            .astype(str)
            .str.strip()
            .str.upper()
        )

    # ========================================================
    # BUILD HGVS ROWS
    # ========================================================

    hgvs_rows = []

    for _, row in summary.iterrows():

        transcript_position = int(
            row["cDNA_Position"]
        )

        ref = str(
            row["REF"]
        ).strip().upper()

        # ----------------------------------------------------
        # Obtain ALT alleles directly from the validated
        # genotype calls.
        # ----------------------------------------------------

        matching_variants = variants[
            (
                variants["cDNA_Position"]
                == transcript_position
            )
            &
            (
                variants["REF"]
                == ref
            )
        ].copy()

        alt_values = (
            matching_variants["ALT"]
            .dropna()
            .astype(str)
            .str.strip()
            .str.upper()
            .unique()
            .tolist()
        )

        # Remove invalid / reference ALT values

        alt_values = sorted({
            alt
            for alt in alt_values
            if alt
            and alt != ref
            and alt != "N"
            and alt != "-"
        })

        # ----------------------------------------------------
        # If no validated ALT exists, do not create an HGVS
        # variant record.
        # ----------------------------------------------------

        if not alt_values:
            continue

        # ----------------------------------------------------
        # One HGVS record per ALT
        # ----------------------------------------------------

        for alt in alt_values:

            # ------------------------------------------------
            # Calculate HGVS cDNA coordinate.
            #
            # Transcript positions are 1-based.
            # ------------------------------------------------

            hgvs_position = (
                transcript_position
                - cds_start
                + 1
            )

            if hgvs_position < 1:

                raise ValueError(
                    f"Invalid HGVS coordinate for "
                    f"transcript position "
                    f"{transcript_position}. "
                    f"Calculated cDNA position: "
                    f"{hgvs_position}. "
                    f"Check CDS_START in config.py."
                )

            # ------------------------------------------------
            # Determine variant type.
            #
            # Current Sanger SNV pipeline primarily handles
            # single-nucleotide substitutions.
            #
            # The length check makes the output generic enough
            # to distinguish simple substitutions from
            # insertions/deletions if they appear.
            # ------------------------------------------------

            if len(ref) == 1 and len(alt) == 1:

                variant_type = "SNV"

            elif len(ref) == 1 and len(alt) > 1:

                variant_type = "Insertion"

            elif len(ref) > 1 and len(alt) == 1:

                variant_type = "Deletion"

            else:

                variant_type = "Complex"

            # ------------------------------------------------
            # HGVS cDNA notation
            # ------------------------------------------------

            hgvs = (
                f"{transcript}:c."
                f"{hgvs_position}"
                f"{ref}>{alt}"
            )

            # ------------------------------------------------
            # Obtain genotype-derived summary for this
            # exact variant.
            # ------------------------------------------------

            exact_variant = matching_variants[
                matching_variants["ALT"]
                == alt
            ].copy()

            if exact_variant.empty:
                continue

            # ------------------------------------------------
            # Carrier count
            # ------------------------------------------------

            if "Sample" in exact_variant.columns:

                samples = (
                    exact_variant["Sample"]
                    .dropna()
                    .astype(str)
                    .str.strip()
                    .str.replace(
                        " ",
                        "",
                        regex=False
                    )
                    .loc[
                        lambda x: x.ne("")
                    ]
                    .drop_duplicates()
                    .tolist()
                )

                carrier_count = len(
                    samples
                )

                samples_string = ";".join(
                    sorted(samples)
                )

            elif "Samples" in row:

                samples_string = str(
                    row["Samples"]
                )

                if samples_string in {
                    "",
                    "nan",
                    "None",
                }:

                    samples_string = ""

                samples = [
                    sample.strip()
                    for sample in
                    samples_string.split(";")
                    if sample.strip()
                ]

                carrier_count = len(
                    set(samples)
                )

            else:

                samples_string = ""
                carrier_count = 0

            # ------------------------------------------------
            # Variant call count
            # ------------------------------------------------

            variant_calls = len(
                exact_variant
            )

            # ------------------------------------------------
            # Quality
            # ------------------------------------------------

            if "Quality" in exact_variant.columns:

                mean_quality = pd.to_numeric(
                    exact_variant["Quality"],
                    errors="coerce"
                ).mean()

            elif "Mean_Quality" in row:

                mean_quality = pd.to_numeric(
                    row["Mean_Quality"],
                    errors="coerce"
                )

            else:

                mean_quality = None

            # ------------------------------------------------
            # Alignment score
            # ------------------------------------------------

            if "Alignment_Score" in exact_variant.columns:

                mean_alignment_score = pd.to_numeric(
                    exact_variant["Alignment_Score"],
                    errors="coerce"
                ).mean()

            elif "Mean_Alignment_Score" in row:

                mean_alignment_score = pd.to_numeric(
                    row["Mean_Alignment_Score"],
                    errors="coerce"
                )

            else:

                mean_alignment_score = None

            # ------------------------------------------------
            # Store row
            # ------------------------------------------------

            hgvs_rows.append({

                "Transcript":
                    transcript,

                "Transcript_Position":
                    transcript_position,

                "HGVS_cDNA_Position":
                    int(hgvs_position),

                "HGVS_cDNA":
                    hgvs,

                "REF":
                    ref,

                "ALT":
                    alt,

                "Variant_Type":
                    variant_type,

                "Carrier_Count":
                    carrier_count,

                "Variant_Calls":
                    variant_calls,

                "Mean_Quality":
                    (
                        round(
                            mean_quality,
                            2
                        )
                        if pd.notna(mean_quality)
                        else None
                    ),

                "Mean_Alignment_Score":
                    (
                        round(
                            mean_alignment_score,
                            2
                        )
                        if pd.notna(
                            mean_alignment_score
                        )
                        else None
                    ),

                "Samples":
                    samples_string,

            })

    # ========================================================
    # CREATE DATAFRAME
    # ========================================================

    hgvs_df = pd.DataFrame(
        hgvs_rows
    )

    # ========================================================
    # EMPTY RESULT
    # ========================================================

    if hgvs_df.empty:

        return pd.DataFrame(
            columns=[
                "Transcript",
                "Transcript_Position",
                "HGVS_cDNA_Position",
                "HGVS_cDNA",
                "REF",
                "ALT",
                "Variant_Type",
                "Carrier_Count",
                "Variant_Calls",
                "Mean_Quality",
                "Mean_Alignment_Score",
                "Samples",
            ]
        )

    # ========================================================
    # REMOVE DUPLICATES
    # ========================================================

    hgvs_df = (
        hgvs_df
        .drop_duplicates(
            subset=[
                "HGVS_cDNA"
            ],
            keep="first"
        )
        .reset_index(drop=True)
    )

    # ========================================================
    # SORT
    # ========================================================

    hgvs_df = hgvs_df.sort_values(
        by=[
            "HGVS_cDNA_Position",
            "REF",
            "ALT",
        ]
    ).reset_index(drop=True)

    # ========================================================
    # FINAL VALIDATION
    # ========================================================

    if hgvs_df["HGVS_cDNA"].duplicated().any():

        raise RuntimeError(
            "Duplicate HGVS variants detected after "
            "HGVS generation."
        )

    if (
        hgvs_df["REF"].isna().any()
        or hgvs_df["ALT"].isna().any()
    ):

        raise RuntimeError(
            "HGVS table contains missing REF or ALT values."
        )

    return hgvs_df


# ============================================================
# SAVE HGVS TABLE
# ============================================================

def save_hgvs_table(
    hgvs_df: pd.DataFrame,
    output_folder: str
):
    """
    Save the HGVS variant table.

    Parameters
    ----------
    hgvs_df : pandas.DataFrame

    output_folder : str

    Returns
    -------
    str
        Path to saved file.
    """

    os.makedirs(
        output_folder,
        exist_ok=True
    )

    output_file = os.path.join(
        output_folder,
        "HGVS_Table.csv"
    )

    hgvs_df.to_csv(
        output_file,
        index=False
    )

    return output_file