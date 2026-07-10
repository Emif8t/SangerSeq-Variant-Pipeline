"""
validation.py

Functions for quality assessment and validation
of the final annotated variant dataset.
"""

import os
import pandas as pd


# ======================================================
# GENERATE QC REPORT
# ======================================================

def generate_qc_report(
    final_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Generate a quality-control report for
    the annotated variant dataset.

    Parameters
    ----------
    final_df : pandas.DataFrame

    Returns
    -------
    pandas.DataFrame
    """

    duplicate_hgvs = final_df.duplicated(
        subset="HGVS_cDNA"
    ).sum()

    duplicate_genomic = final_df.duplicated(
        subset=[
            "Chromosome",
            "Genomic_Position",
            "ALT"
        ]
    ).sum()

    report = pd.DataFrame({

        "Metric":[

            "Total variants",

            "Duplicate HGVS",

            "Duplicate genomic coordinates",

            "Variants with dbSNP rsID",

            "Variants without dbSNP rsID",

            "Missing Gene",

            "Missing HGVS",

            "Missing Consequence",

            "Missing Impact",

            "Missing Genomic Position"

        ],

        "Value":[

            len(final_df),

            duplicate_hgvs,

            duplicate_genomic,

            final_df["dbSNP_ID"].notna().sum(),

            final_df["dbSNP_ID"].isna().sum(),

            final_df["Gene"].isna().sum(),

            final_df["HGVS_cDNA"].isna().sum(),

            final_df["Consequence"].isna().sum(),

            final_df["Impact"].isna().sum(),

            final_df["Genomic_Position"].isna().sum()

        ]

    })

    return report


# ======================================================
# SUMMARIZE CONSEQUENCES
# ======================================================

def summarize_variant_consequences(
    final_df: pd.DataFrame
) -> pd.Series:

    return (

        final_df["Consequence"]

        .fillna("Unknown")

        .value_counts()

    )

# ======================================================
# SUMMARIZE IMPACT
# ======================================================

def summarize_variant_impact(
    final_df: pd.DataFrame
) -> pd.Series:

    return (

        final_df["Impact"]

        .fillna("Unknown")

        .value_counts()

    )


# ======================================================
# SUMMARIZE VARIANT TYPES
# ======================================================

def summarize_variant_types(
    final_df: pd.DataFrame
) -> pd.Series:

    return (

        final_df["Variant_Type"]

        .fillna("Unknown")

        .value_counts()

    )


# ======================================================
# RUN QC
# ======================================================

def run_variant_quality_control(
    final_df: pd.DataFrame
):
    """
    Run all quality-control summaries.

    Returns
    -------
    tuple
    """

    qc_report = generate_qc_report(
        final_df
    )

    consequence_summary = (

        summarize_variant_consequences(
            final_df
        )

    )

    impact_summary = (

        summarize_variant_impact(
            final_df
        )

    )

    variant_type_summary = (

        summarize_variant_types(
            final_df
        )

    )

    return (

        qc_report,

        consequence_summary,

        impact_summary,

        variant_type_summary

    )


# ======================================================
# SAVE QC REPORTS
# ======================================================

def save_qc_reports(

    qc_report,

    consequence_summary,

    impact_summary,

    variant_type_summary,

    final_df,

    output_folder

):
    """
    Save all QC reports.
    """

    os.makedirs(

        output_folder,

        exist_ok=True

    )

    qc_report.to_csv(

        os.path.join(

            output_folder,

            "ASS1_QC_Report.csv"

        ),

        index=False

    )

    consequence_summary.to_csv(

        os.path.join(

            output_folder,

            "ASS1_Consequence_Summary.csv"

        )

    )

    impact_summary.to_csv(

        os.path.join(

            output_folder,

            "ASS1_Impact_Summary.csv"

        )

    )

    variant_type_summary.to_csv(

        os.path.join(

            output_folder,

            "ASS1_VariantType_Summary.csv"

        )

    )

    final_df.to_csv(

        os.path.join(

            output_folder,

            "ASS1_Final_Annotated_Variants_QC.csv"

        ),

        index=False

    )
