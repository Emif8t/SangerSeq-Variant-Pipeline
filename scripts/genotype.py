"""
genotype.py

Functions for genotype determination, quality filtering,
and summarization of Sanger sequencing variants.
"""

import os
import pandas as pd


# ======================================================
# IUPAC AMBIGUITY CODES
# ======================================================

IUPAC_CODES = {

    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"]

}


# ======================================================
# GENOTYPE CALLING
# ======================================================

def call_genotypes(
    alignment_df: pd.DataFrame,
    min_phred: int = 20
) -> pd.DataFrame:
    """
    Assign genotypes from aligned nucleotide calls.

    Parameters
    ----------
    alignment_df : pandas.DataFrame

    min_phred : int

    Returns
    -------
    pandas.DataFrame
    """

    genotype_rows = []

    for _, row in alignment_df.iterrows():

        ref = row["REF"]
        obs = row["Observed_Base"]

        if obs == "-" or ref == "-" or pd.isna(obs):
            continue

        obs = obs.upper()

        # ---------------------------------------------
        # Unknown base
        # ---------------------------------------------

        if obs == "N":

            genotype_rows.append({

                "Sample": row["Sample"],
                "Alignment_Score": row["Alignment_Score"],
                "Read_Position": row["Read_Position"],
                "Amplicon_Position": row["Amplicon_Position"],
                "cDNA_Position": row["cDNA_Position"],
                "REF": ref,
                "Observed_Base": obs,
                "ALT": None,
                "Alleles": None,
                "Genotype": None,
                "Zygosity": "Unknown",
                "Is_Variant": False,
                "Quality": row["Quality"],
                "Confidence": row["Confidence"]

            })

            continue

        # ---------------------------------------------
        # Standard nucleotide
        # ---------------------------------------------

        if obs in ["A", "C", "G", "T"]:

            genotype = 0 if obs == ref else 2

            genotype_rows.append({

                "Sample": row["Sample"],
                "Alignment_Score": row["Alignment_Score"],
                "Read_Position": row["Read_Position"],
                "Amplicon_Position": row["Amplicon_Position"],
                "cDNA_Position": row["cDNA_Position"],
                "REF": ref,
                "Observed_Base": obs,
                "ALT": obs,
                "Alleles": obs,
                "Genotype": genotype,
                "Zygosity":

                    "Homozygous_Reference"

                    if genotype == 0

                    else "Homozygous_Variant",

                "Is_Variant": genotype == 2,
                "Quality": row["Quality"],
                "Confidence": row["Confidence"]

            })

            continue

        # ---------------------------------------------
        # IUPAC ambiguity
        # ---------------------------------------------

        if obs in IUPAC_CODES:

            alleles = IUPAC_CODES[obs]

            alt = [a for a in alleles if a != ref]

            alt = alt[0] if alt else ref

            if row["Quality"] >= min_phred:

                genotype = 1

                zygosity = "Heterozygous"

                variant = True

            else:

                genotype = None

                zygosity = "LowConfidence_Heterozygous"

                variant = False

            genotype_rows.append({

                "Sample": row["Sample"],
                "Alignment_Score": row["Alignment_Score"],
                "Read_Position": row["Read_Position"],
                "Amplicon_Position": row["Amplicon_Position"],
                "cDNA_Position": row["cDNA_Position"],
                "REF": ref,
                "Observed_Base": obs,
                "ALT": alt,
                "Alleles": "/".join(alleles),
                "Genotype": genotype,
                "Zygosity": zygosity,
                "Is_Variant": variant,
                "Quality": row["Quality"],
                "Confidence": row["Confidence"]

            })

    return pd.DataFrame(genotype_rows)


# ======================================================
# HIGH-CONFIDENCE FILTERING
# ======================================================

def filter_high_confidence_variants(
    genotype_df: pd.DataFrame,
    min_phred: int = 20
):
    """
    Filter genotype calls based on quality.

    Returns
    -------
    tuple
        high_confidence_df, variant_df
    """

    high_confidence_df = genotype_df[
        genotype_df["Quality"] >= min_phred
    ].copy()

    high_confidence_df = high_confidence_df[
        high_confidence_df["Genotype"].notna()
    ].copy()

    variant_df = high_confidence_df[
        high_confidence_df["Is_Variant"]
    ].copy()

    variant_df = variant_df.sort_values(

        by=["cDNA_Position", "Sample"]

    ).reset_index(drop=True)

    return high_confidence_df, variant_df


# ======================================================
# VARIANT SUMMARY
# ======================================================

def summarize_variants(
    variant_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Generate a summary of validated variants.

    Parameters
    ----------
    variant_df : pandas.DataFrame

    Returns
    -------
    pandas.DataFrame
    """

    summary_df = (

        variant_df

        .groupby(

            ["cDNA_Position", "REF"],

            as_index=False

        )

        .agg(

            Variant_Calls=("Sample", "count"),

            Mean_Quality=("Quality", "mean"),

            Mean_Alignment_Score=("Alignment_Score", "mean"),

            Samples=(

                "Sample",

                lambda x:

                ";".join(

                    sorted(

                        x.unique()

                    )

                )

            )

        )

    )

    summary_df["Carrier_Count"] = (

        summary_df["Samples"]

        .str.split(";")

        .apply(

            lambda x:

            len(set(x))

        )

    )

    summary_df = summary_df.sort_values(

        by=[

            "Carrier_Count",

            "Mean_Quality"

        ],

        ascending=False

    ).reset_index(drop=True)

    return summary_df


# ======================================================
# SAVE FUNCTIONS
# ======================================================

def save_genotypes(
    genotype_df: pd.DataFrame,
    output_folder: str
):
    os.makedirs(output_folder, exist_ok=True)
    genotype_df.to_csv(
        os.path.join(output_folder, "Genotype_Table.csv"),
        index=False
    )


def save_variants(
    variant_df: pd.DataFrame,
    output_folder: str
):
    os.makedirs(output_folder, exist_ok=True)
    variant_df.to_csv(
        os.path.join(output_folder, "HighConfidence_Variants.csv"),
        index=False
    )


def save_variant_summary(
    summary_df: pd.DataFrame,
    output_folder: str
):
    os.makedirs(output_folder, exist_ok=True)
    summary_df.to_csv(
        os.path.join(output_folder, "Variant_Summary.csv"),
        index=False
    )
