"""
genotype.py

Functions for genotype determination, quality filtering,
and summarisation of Sanger sequencing variants.

This module is gene-agnostic and does not contain
gene-specific assumptions.
"""

import os
from typing import Tuple

import pandas as pd


# ============================================================
# IUPAC AMBIGUITY CODES
# ============================================================

IUPAC_CODES = {
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def clean_base(value):
    """
    Convert a nucleotide value to an uppercase string.

    Missing values are returned as None.
    """

    if pd.isna(value):
        return None

    value = str(value).strip().upper()

    if value in {"", "NAN", "NONE", "NULL"}:
        return None

    return value


# ============================================================
# GENOTYPE CALLING
# ============================================================

def call_genotypes(
    alignment_df: pd.DataFrame,
    min_phred: int = 20
) -> pd.DataFrame:
    """
    Assign genotype calls from nucleotide-level alignment data.

    Standard bases:
        Reference base -> Homozygous_Reference
        Alternative base -> Homozygous_Variant

    IUPAC ambiguity codes:
        High-quality ambiguity -> Heterozygous
        Low-quality ambiguity -> LowConfidence_Heterozygous

    Unknown base:
        N -> Unknown

    Parameters
    ----------
    alignment_df : pandas.DataFrame
        Nucleotide-level alignment table.

    min_phred : int
        Minimum Phred quality required for a confident
        heterozygous call.

    Returns
    -------
    pandas.DataFrame
        Genotype-level calls.
    """

    required_columns = [
        "Sample",
        "Alignment_Score",
        "Read_Position",
        "Amplicon_Position",
        "cDNA_Position",
        "REF",
        "Observed_Base",
        "Quality",
        "Confidence",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in alignment_df.columns
    ]

    if missing_columns:
        raise KeyError(
            "Alignment table is missing required columns:\n"
            + "\n".join(
                f"  - {column}"
                for column in missing_columns
            )
        )

    genotype_rows = []

    for _, row in alignment_df.iterrows():

        ref = clean_base(row["REF"])
        obs = clean_base(row["Observed_Base"])

        # ----------------------------------------------------
        # Ignore unusable positions
        # ----------------------------------------------------

        if ref is None or obs is None:
            continue

        if obs == "-":
            continue

        quality = pd.to_numeric(
            row["Quality"],
            errors="coerce"
        )

        # ====================================================
        # UNKNOWN BASE
        # ====================================================

        if obs == "N":

            genotype_rows.append({

                "Sample":
                    row["Sample"],

                "Alignment_Score":
                    row["Alignment_Score"],

                "Read_Position":
                    row["Read_Position"],

                "Amplicon_Position":
                    row["Amplicon_Position"],

                "cDNA_Position":
                    row["cDNA_Position"],

                "REF":
                    ref,

                "Observed_Base":
                    obs,

                "ALT":
                    None,

                "Alleles":
                    None,

                "Genotype":
                    None,

                "Zygosity":
                    "Unknown",

                "Is_Variant":
                    False,

                "Quality":
                    quality,

                "Confidence":
                    row["Confidence"],

            })

            continue

        # ====================================================
        # STANDARD NUCLEOTIDES
        # ====================================================

        if obs in {"A", "C", "G", "T"}:

            if obs == ref:

                genotype = 0
                zygosity = "Homozygous_Reference"
                is_variant = False
                alt = None

            else:

                genotype = 2
                zygosity = "Homozygous_Variant"
                is_variant = True
                alt = obs

            genotype_rows.append({

                "Sample":
                    row["Sample"],

                "Alignment_Score":
                    row["Alignment_Score"],

                "Read_Position":
                    row["Read_Position"],

                "Amplicon_Position":
                    row["Amplicon_Position"],

                "cDNA_Position":
                    row["cDNA_Position"],

                "REF":
                    ref,

                "Observed_Base":
                    obs,

                "ALT":
                    alt,

                "Alleles":
                    obs,

                "Genotype":
                    genotype,

                "Zygosity":
                    zygosity,

                "Is_Variant":
                    is_variant,

                "Quality":
                    quality,

                "Confidence":
                    row["Confidence"],

            })

            continue

        # ====================================================
        # IUPAC HETEROZYGOSITY
        # ====================================================

        if obs in IUPAC_CODES:

            alleles = IUPAC_CODES[obs]

            alternative_alleles = [
                allele
                for allele in alleles
                if allele != ref
            ]

            # ------------------------------------------------
            # If the reference is not represented in the
            # IUPAC code, the call cannot safely be interpreted
            # as a conventional REF/ALT heterozygote.
            # ------------------------------------------------

            if ref not in alleles:

                genotype_rows.append({

                    "Sample":
                        row["Sample"],

                    "Alignment_Score":
                        row["Alignment_Score"],

                    "Read_Position":
                        row["Read_Position"],

                    "Amplicon_Position":
                        row["Amplicon_Position"],

                    "cDNA_Position":
                        row["cDNA_Position"],

                    "REF":
                        ref,

                    "Observed_Base":
                        obs,

                    "ALT":
                        None,

                    "Alleles":
                        "/".join(alleles),

                    "Genotype":
                        None,

                    "Zygosity":
                        "Ambiguous",

                    "Is_Variant":
                        False,

                    "Quality":
                        quality,

                    "Confidence":
                        row["Confidence"],

                })

                continue

            alt = (
                alternative_alleles[0]
                if alternative_alleles
                else None
            )

            # ------------------------------------------------
            # Quality determines whether the heterozygous
            # call is considered a confirmed variant.
            # ------------------------------------------------

            if (
                pd.notna(quality)
                and quality >= min_phred
                and alt is not None
            ):

                genotype = 1
                zygosity = "Heterozygous"
                is_variant = True

            else:

                genotype = None
                zygosity = "LowConfidence_Heterozygous"
                is_variant = False

            genotype_rows.append({

                "Sample":
                    row["Sample"],

                "Alignment_Score":
                    row["Alignment_Score"],

                "Read_Position":
                    row["Read_Position"],

                "Amplicon_Position":
                    row["Amplicon_Position"],

                "cDNA_Position":
                    row["cDNA_Position"],

                "REF":
                    ref,

                "Observed_Base":
                    obs,

                "ALT":
                    alt,

                "Alleles":
                    "/".join(alleles),

                "Genotype":
                    genotype,

                "Zygosity":
                    zygosity,

                "Is_Variant":
                    is_variant,

                "Quality":
                    quality,

                "Confidence":
                    row["Confidence"],

            })

            continue

        # ====================================================
        # OTHER / UNRECOGNISED BASE
        # ====================================================

        genotype_rows.append({

            "Sample":
                row["Sample"],

            "Alignment_Score":
                row["Alignment_Score"],

            "Read_Position":
                row["Read_Position"],

            "Amplicon_Position":
                row["Amplicon_Position"],

            "cDNA_Position":
                row["cDNA_Position"],

            "REF":
                ref,

            "Observed_Base":
                obs,

            "ALT":
                None,

            "Alleles":
                None,

            "Genotype":
                None,

            "Zygosity":
                "Unknown",

            "Is_Variant":
                False,

            "Quality":
                quality,

            "Confidence":
                row["Confidence"],

        })

    return pd.DataFrame(genotype_rows)


# ============================================================
# HIGH-CONFIDENCE VARIANT FILTERING
# ============================================================

def filter_high_confidence_variants(
    genotype_df: pd.DataFrame,
    min_phred: int = 20
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter genotype calls to retain high-confidence calls
    and validated variants.

    A call is considered high confidence when:

        1. Quality >= min_phred
        2. Genotype is not missing
        3. Is_Variant is True for variant calls

    Parameters
    ----------
    genotype_df : pandas.DataFrame

    min_phred : int

    Returns
    -------
    high_confidence_df : pandas.DataFrame
        High-quality genotype calls, including reference calls.

    variant_df : pandas.DataFrame
        High-confidence variant calls only.
    """

    required_columns = [
        "Quality",
        "Genotype",
        "Is_Variant",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in genotype_df.columns
    ]

    if missing_columns:
        raise KeyError(
            "Genotype table is missing required columns:\n"
            + "\n".join(
                f"  - {column}"
                for column in missing_columns
            )
        )

    filtered = genotype_df.copy()

    filtered["Quality_Numeric"] = pd.to_numeric(
        filtered["Quality"],
        errors="coerce"
    )

    high_confidence_df = filtered[
        (
            filtered["Quality_Numeric"]
            >= min_phred
        )
        &
        filtered["Genotype"].notna()
    ].copy()

    # --------------------------------------------------------
    # Remove temporary numeric quality field
    # --------------------------------------------------------

    high_confidence_df.drop(
        columns=["Quality_Numeric"],
        inplace=True
    )

    # --------------------------------------------------------
    # Confirmed variants
    # --------------------------------------------------------

    variant_df = high_confidence_df[
        high_confidence_df["Is_Variant"].astype(bool)
    ].copy()

    # --------------------------------------------------------
    # ALT must be present for a valid variant
    # --------------------------------------------------------

    variant_df = variant_df[
        variant_df["ALT"].notna()
    ].copy()

    variant_df = variant_df[
        variant_df["ALT"]
        .astype(str)
        .str.strip()
        .ne("")
    ].copy()

    # --------------------------------------------------------
    # Sort
    # --------------------------------------------------------

    sort_columns = [
        column
        for column in [
            "cDNA_Position",
            "Sample",
            "ALT",
        ]
        if column in variant_df.columns
    ]

    if sort_columns:

        variant_df = variant_df.sort_values(
            by=sort_columns
        ).reset_index(drop=True)

    return high_confidence_df, variant_df

def get_high_confidence_samples(
    high_confidence_df: pd.DataFrame
) -> set:
    """
    Return sample identifiers with at least one
    high-confidence genotype call.

    A high-confidence genotype call is defined by
    filter_high_confidence_variants() as having:

    - Quality >= minimum Phred threshold
    - Non-missing Genotype
    """

    if high_confidence_df.empty:
        return set()

    if "Sample" not in high_confidence_df.columns:
        raise KeyError(
            "High-confidence genotype table must contain "
            "'Sample'."
        )

    samples = (
        high_confidence_df["Sample"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    return set(
        sample
        for sample in samples
        if sample
    )

# ============================================================
# VARIANT SUMMARY
# ============================================================

def summarize_variants(
    variant_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Generate a variant-level summary from confirmed
    high-confidence genotype calls.

    Each unique combination of:

        cDNA_Position
        REF
        ALT

    represents one variant.

    Carrier information is calculated directly from
    genotype calls.

    Parameters
    ----------
    variant_df : pandas.DataFrame

    Returns
    -------
    pandas.DataFrame
    """

    if variant_df.empty:

        return pd.DataFrame(
            columns=[
                "cDNA_Position",
                "REF",
                "ALT",
                "Variant_Calls",
                "Carrier_Count",
                "Mean_Quality",
                "Mean_Alignment_Score",
                "Samples",
            ]
        )

    required_columns = [
        "cDNA_Position",
        "REF",
        "ALT",
        "Sample",
        "Quality",
        "Alignment_Score",
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in variant_df.columns
    ]

    if missing_columns:
        raise KeyError(
            "Variant table is missing required columns:\n"
            + "\n".join(
                f"  - {column}"
                for column in missing_columns
            )
        )

    df = variant_df.copy()

    # --------------------------------------------------------
    # Clean sample identifiers
    # --------------------------------------------------------

    df["Sample_Clean"] = (
        df["Sample"]
        .astype(str)
        .str.strip()
        .str.replace(" ", "", regex=False)
    )

    # --------------------------------------------------------
    # Numeric fields
    # --------------------------------------------------------

    df["Quality_Numeric"] = pd.to_numeric(
        df["Quality"],
        errors="coerce"
    )

    df["Alignment_Score_Numeric"] = pd.to_numeric(
        df["Alignment_Score"],
        errors="coerce"
    )

    # --------------------------------------------------------
    # Group by actual variant
    # --------------------------------------------------------

    summary_df = (
        df
        .groupby(
            [
                "cDNA_Position",
                "REF",
                "ALT",
            ],
            as_index=False
        )
        .agg(
            Variant_Calls=(
                "Sample_Clean",
                "count"
            ),

            Carrier_Count=(
                "Sample_Clean",
                lambda values:
                    values.nunique()
            ),

            Mean_Quality=(
                "Quality_Numeric",
                "mean"
            ),

            Mean_Alignment_Score=(
                "Alignment_Score_Numeric",
                "mean"
            ),

            Samples=(
                "Sample_Clean",
                lambda values:
                    ";".join(
                        sorted(
                            set(
                                value
                                for value in values
                                if value
                            )
                        )
                    )
            ),
        )
    )

    # --------------------------------------------------------
    # Round summary values
    # --------------------------------------------------------

    summary_df["Mean_Quality"] = (
        summary_df["Mean_Quality"]
        .round(2)
    )

    summary_df["Mean_Alignment_Score"] = (
        summary_df["Mean_Alignment_Score"]
        .round(2)
    )

    # --------------------------------------------------------
    # Sort by genomic/transcript position
    # --------------------------------------------------------

    summary_df = summary_df.sort_values(
        by=[
            "cDNA_Position",
            "REF",
            "ALT",
        ]
    ).reset_index(drop=True)

    return summary_df


# ============================================================
# SAVE GENOTYPE TABLE
# ============================================================

def save_genotypes(
    genotype_df: pd.DataFrame,
    output_folder: str
):
    """
    Save complete genotype table.
    """

    os.makedirs(
        output_folder,
        exist_ok=True
    )

    output_file = os.path.join(
        output_folder,
        "Genotype_Table.csv"
    )

    genotype_df.to_csv(
        output_file,
        index=False
    )

    return output_file


# ============================================================
# SAVE HIGH-CONFIDENCE VARIANTS
# ============================================================

def save_variants(
    variant_df: pd.DataFrame,
    output_folder: str
):
    """
    Save high-confidence variant calls.
    """

    os.makedirs(
        output_folder,
        exist_ok=True
    )

    output_file = os.path.join(
        output_folder,
        "HighConfidence_Variants.csv"
    )

    variant_df.to_csv(
        output_file,
        index=False
    )

    return output_file


# ============================================================
# SAVE VARIANT SUMMARY
# ============================================================

def save_variant_summary(
    summary_df: pd.DataFrame,
    output_folder: str
):
    """
    Save variant-level summary.
    """

    os.makedirs(
        output_folder,
        exist_ok=True
    )

    output_file = os.path.join(
        output_folder,
        "Variant_Summary.csv"
    )

    summary_df.to_csv(
        output_file,
        index=False
    )

    return output_file