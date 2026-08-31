"""
validation.py

Functions for quality control and validation of the
final annotated Sanger sequencing variant table.
"""

import os
import pandas as pd


# =========================================================
# 1. GENERATE QC REPORT
# =========================================================

def generate_qc_report(
    final_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Generate a general quality-control report for the
    final annotated variant table.

    The function adapts to the columns available in the
    final table. It does not require genomic coordinates
    if the VEP output does not provide them.
    """

    final_df = final_df.copy()

    total_variants = len(final_df)

    # -----------------------------------------------------
    # Determine the most appropriate variant identifier
    # -----------------------------------------------------

    if "HGVS_cDNA" in final_df.columns:
        variant_key = ["HGVS_cDNA"]

    elif "Transcript_Position" in final_df.columns:
        variant_key = ["Transcript_Position", "REF", "ALT"]

    elif "HGVS_cDNA_Position" in final_df.columns:
        variant_key = ["HGVS_cDNA_Position", "REF", "ALT"]

    elif "Location" in final_df.columns:
        variant_key = ["Location", "REF", "ALT"]

    else:
        variant_key = None

    # -----------------------------------------------------
    # Duplicate variant detection
    # -----------------------------------------------------

    if variant_key:

        duplicate_variant_count = int(
            final_df.duplicated(
                subset=[
                    column
                    for column in variant_key
                    if column in final_df.columns
                ],
                keep=False
            ).sum()
        )

    else:

        duplicate_variant_count = 0

    # -----------------------------------------------------
    # Missing annotation
    # -----------------------------------------------------

    if "Gene" in final_df.columns:

        missing_gene = int(
            final_df["Gene"]
            .replace("", pd.NA)
            .isna()
            .sum()
        )

    else:

        missing_gene = total_variants

    if "Consequence" in final_df.columns:

        missing_consequence = int(
            final_df["Consequence"]
            .replace("", pd.NA)
            .isna()
            .sum()
        )

    else:

        missing_consequence = total_variants

    if "Impact" in final_df.columns:

        missing_impact = int(
            final_df["Impact"]
            .replace("", pd.NA)
            .isna()
            .sum()
        )

    elif "IMPACT" in final_df.columns:

        missing_impact = int(
            final_df["IMPACT"]
            .replace("", pd.NA)
            .isna()
            .sum()
        )

    else:

        missing_impact = total_variants

    # -----------------------------------------------------
    # HGVS completeness
    # -----------------------------------------------------

    if "HGVS_cDNA" in final_df.columns:

        missing_hgvs = int(
            final_df["HGVS_cDNA"]
            .replace("", pd.NA)
            .isna()
            .sum()
        )

    else:

        missing_hgvs = total_variants

    # -----------------------------------------------------
    # Database identifier completeness
    # -----------------------------------------------------

    if "dbSNP_ID" in final_df.columns:

        variants_with_dbsnp = int(
            (
                final_df["dbSNP_ID"]
                .astype(str)
                .str.strip()
                .replace("", pd.NA)
                .notna()
            ).sum()
        )

    elif "Existing_variation" in final_df.columns:

        variants_with_dbsnp = int(
            final_df["Existing_variation"]
            .astype(str)
            .str.contains(
                r"\brs\d+\b",
                regex=True,
                na=False
            )
            .sum()
        )

    else:

        variants_with_dbsnp = 0

    # -----------------------------------------------------
    # Genotype information
    # -----------------------------------------------------

    if "Zygosity" in final_df.columns:

        zygosity_counts = (
            final_df["Zygosity"]
            .fillna("Unknown")
            .value_counts()
            .to_dict()
        )

    else:

        zygosity_counts = {}

    # -----------------------------------------------------
    # Construct QC report
    # -----------------------------------------------------

    qc_rows = [

        {
            "Metric": "Total variants",
            "Value": total_variants
        },

        {
            "Metric": "Duplicate variant records",
            "Value": duplicate_variant_count
        },

        {
            "Metric": "Variants with missing Gene",
            "Value": missing_gene
        },

        {
            "Metric": "Variants with missing Consequence",
            "Value": missing_consequence
        },

        {
            "Metric": "Variants with missing Impact",
            "Value": missing_impact
        },

        {
            "Metric": "Variants with missing HGVS cDNA",
            "Value": missing_hgvs
        },

        {
            "Metric": "Variants with dbSNP identifier",
            "Value": variants_with_dbsnp
        }

    ]

    # Add zygosity information
    for zygosity, count in zygosity_counts.items():

        qc_rows.append(
            {
                "Metric": f"Zygosity: {zygosity}",
                "Value": count
            }
        )

    return pd.DataFrame(qc_rows)


# =========================================================
# 2. CONSEQUENCE SUMMARY
# =========================================================

def generate_consequence_summary(
    final_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Summarise variant consequences.
    """

    if "Consequence" not in final_df.columns:

        return pd.DataFrame(
            columns=[
                "Consequence",
                "Count"
            ]
        )

    consequence_df = (

        final_df["Consequence"]
        .fillna("Unknown")
        .replace("", "Unknown")
        .value_counts()
        .rename_axis("Consequence")
        .reset_index(name="Count")

    )

    return consequence_df


# =========================================================
# 3. IMPACT SUMMARY
# =========================================================

def generate_impact_summary(
    final_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Summarise predicted variant impacts.
    """

    if "Impact" in final_df.columns:

        impact_column = "Impact"

    elif "IMPACT" in final_df.columns:

        impact_column = "IMPACT"

    else:

        return pd.DataFrame(
            columns=[
                "Impact",
                "Count"
            ]
        )

    impact_df = (

        final_df[impact_column]
        .fillna("Unknown")
        .replace("", "Unknown")
        .value_counts()
        .rename_axis("Impact")
        .reset_index(name="Count")

    )

    return impact_df


# =========================================================
# 4. VARIANT TYPE SUMMARY
# =========================================================

def generate_variant_type_summary(
    final_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Summarise variant types.
    """

    # -----------------------------------------------------
    # Use existing Variant_Type if available
    # -----------------------------------------------------

    if "Variant_Type" in final_df.columns:

        variant_type_df = (

            final_df["Variant_Type"]
            .fillna("Unknown")
            .replace("", "Unknown")
            .value_counts()
            .rename_axis("Variant_Type")
            .reset_index(name="Count")

        )

        return variant_type_df

    # -----------------------------------------------------
    # Otherwise infer from REF and ALT
    # -----------------------------------------------------

    if (
        "REF" not in final_df.columns
        or "ALT" not in final_df.columns
    ):

        return pd.DataFrame(
            columns=[
                "Variant_Type",
                "Count"
            ]
        )

    def classify_variant(row):

        ref = str(row["REF"]).strip()
        alt = str(row["ALT"]).strip()

        if not ref or not alt:
            return "Unknown"

        # Multiple ALT alleles
        if "," in alt:
            return "Multiple_ALTs"

        # SNV
        if len(ref) == 1 and len(alt) == 1:
            return "SNV"

        # Insertion
        if len(alt) > len(ref):
            return "Insertion"

        # Deletion
        if len(alt) < len(ref):
            return "Deletion"

        return "Complex"

    variant_types = (
        final_df
        .apply(classify_variant, axis=1)
        .value_counts()
        .rename_axis("Variant_Type")
        .reset_index(name="Count")
    )

    return variant_types


# =========================================================
# 5. MASTER QC FUNCTION
# =========================================================

def run_variant_quality_control(
    final_df: pd.DataFrame
):
    """
    Run all variant quality-control analyses.

    Returns
    -------
    tuple
        qc_report,
        consequence_summary,
        impact_summary,
        variant_type_summary
    """

    final_df = final_df.copy()

    qc_report = generate_qc_report(
        final_df
    )

    consequence_summary = generate_consequence_summary(
        final_df
    )

    impact_summary = generate_impact_summary(
        final_df
    )

    variant_type_summary = generate_variant_type_summary(
        final_df
    )

    return (
        qc_report,
        consequence_summary,
        impact_summary,
        variant_type_summary
    )


# =========================================================
# 6. SAVE QC REPORTS
# =========================================================

def save_qc_reports(
    qc_report: pd.DataFrame,
    consequence_summary: pd.DataFrame,
    impact_summary: pd.DataFrame,
    variant_type_summary: pd.DataFrame,
    final_df: pd.DataFrame,
    output_folder: str
):
    """
    Save all variant quality-control reports.
    """

    os.makedirs(
        output_folder,
        exist_ok=True
    )

    # -----------------------------------------------------
    # General QC report
    # -----------------------------------------------------

    qc_report.to_csv(
        os.path.join(
            output_folder,
            "Variant_QC_Report.csv"
        ),
        index=False
    )

    # -----------------------------------------------------
    # Consequence summary
    # -----------------------------------------------------

    consequence_summary.to_csv(
        os.path.join(
            output_folder,
            "Consequence_Summary.csv"
        ),
        index=False
    )

    # -----------------------------------------------------
    # Impact summary
    # -----------------------------------------------------

    impact_summary.to_csv(
        os.path.join(
            output_folder,
            "Impact_Summary.csv"
        ),
        index=False
    )

    # -----------------------------------------------------
    # Variant type summary
    # -----------------------------------------------------

    variant_type_summary.to_csv(
        os.path.join(
            output_folder,
            "Variant_Type_Summary.csv"
        ),
        index=False
    )

    # -----------------------------------------------------
    # Final table used for QC
    # -----------------------------------------------------

    final_df.to_csv(
        os.path.join(
            output_folder,
            "Final_Table_QC_Copy.csv"
        ),
        index=False
    )