import os

import numpy as np
import pandas as pd

from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.multitest import multipletests


# =========================================================
# Load sample groups
# =========================================================

def load_sample_groups(
    phenotype_file: str
) -> pd.DataFrame:
    """
    Load sample phenotype information.
    """

    return pd.read_excel(
        phenotype_file
    )


# =========================================================
# Create sample sets
# =========================================================

def create_sample_sets(
    phenotype_df: pd.DataFrame
):
    """
    Create case and control sample sets.
    """

    case_samples = set(
        phenotype_df.loc[
            phenotype_df["Group"] == "Case",
            "Sample"
        ].astype(str)
    )

    control_samples = set(
        phenotype_df.loc[
            phenotype_df["Group"] == "Control",
            "Sample"
        ].astype(str)
    )

    return case_samples, control_samples


# =========================================================
# Analyse individual variant
# =========================================================

def analyse_variant(
    row,
    case_samples,
    control_samples
):
    """
    Analyse a single variant using a case-control
    carrier association framework.

    Calculates:

    - Case carriers
    - Control carriers
    - Case non-carriers
    - Control non-carriers
    - Fisher's exact test P-value
    - Chi-square P-value
    - Odds ratio
    - 95% confidence interval
    """

    # -----------------------------------------------------
    # Samples carrying this variant
    # -----------------------------------------------------

    carriers = set()

    if pd.notna(row["Samples"]):

        carriers = set(
            s.strip()
            for s in str(row["Samples"]).split(";")
        )

    # -----------------------------------------------------
    # Carrier counts
    # -----------------------------------------------------

    case_carriers = len(
        carriers & case_samples
    )

    control_carriers = len(
        carriers & control_samples
    )

    case_noncarriers = (
        len(case_samples) - case_carriers
    )

    control_noncarriers = (
        len(control_samples) - control_carriers
    )

    # -----------------------------------------------------
    # 2 × 2 contingency table
    # -----------------------------------------------------

    table = np.array([
        [case_carriers, case_noncarriers],
        [control_carriers, control_noncarriers]
    ])

    # -----------------------------------------------------
    # Fisher's exact test
    # -----------------------------------------------------

    try:

        _, fisher_p = fisher_exact(table)

    except Exception:

        fisher_p = np.nan

    # -----------------------------------------------------
    # Chi-square test
    # -----------------------------------------------------

    try:

        _, chi_p, _, _ = chi2_contingency(
            table,
            correction=False
        )

    except Exception:

        chi_p = np.nan

    # -----------------------------------------------------
    # Odds ratio + 95% CI
    #
    # Table2x2 applies the Haldane correction internally
    # when zero cells are encountered.
    # -----------------------------------------------------

    try:

        t = Table2x2(table)

        odds_ratio = t.oddsratio

        ci_low, ci_high = (
            t.oddsratio_confint()
        )

    except Exception:

        odds_ratio = np.nan

        ci_low = np.nan

        ci_high = np.nan

    # -----------------------------------------------------
    # Return results
    # -----------------------------------------------------

    return {
        **row.to_dict(),

        "Case_Carriers": case_carriers,

        "Control_Carriers": control_carriers,

        "Case_NonCarriers": case_noncarriers,

        "Control_NonCarriers": control_noncarriers,

        "Odds_Ratio": odds_ratio,

        "CI_Lower": ci_low,

        "CI_Upper": ci_high,

        "Fisher_P": fisher_p,

        "ChiSquare_P": chi_p
    }


# =========================================================
# Run association analysis
# =========================================================

def run_association_analysis(
    final_df,
    case_samples,
    control_samples
):
    """
    Perform case-control association analysis
    for every variant.
    """

    results = []

    for _, row in final_df.iterrows():

        results.append(
            analyse_variant(
                row,
                case_samples,
                control_samples
            )
        )

    return pd.DataFrame(
        results
    )


# =========================================================
# Multiple-testing correction
# =========================================================

def apply_multiple_testing(
    association_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Apply Benjamini-Hochberg FDR correction
    to Fisher's exact test P-values.
    """

    association_df = association_df.copy()

    association_df["FDR_P"] = (
        multipletests(
            association_df["Fisher_P"],
            method="fdr_bh"
        )[1]
    )

    association_df["Significant_FDR"] = (
        association_df["FDR_P"] < 0.05
    )

    return association_df


# =========================================================
# Sort association results
# =========================================================

def sort_association_results(
    association_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Sort variants by Fisher's exact test P-value
    and odds ratio.
    """

    return association_df.sort_values(
        by=[
            "Fisher_P",
            "Odds_Ratio"
        ],
        ascending=[
            True,
            False
        ]
    ).reset_index(
        drop=True
    )


# =========================================================
# Save association results
# =========================================================

def save_association_results(
    association_df: pd.DataFrame,
    output_folder: str
) -> str:
    """
    Save association analysis results.

    Parameters
    ----------
    association_df : pandas.DataFrame
        Association analysis results.

    output_folder : str
        Directory in which to save the results.

    Returns
    -------
    str
        Path to the saved association results file.
    """

    os.makedirs(
        output_folder,
        exist_ok=True
    )

    output_file = os.path.join(
        output_folder,
        "ASS1_Association_Analysis.csv"
    )

    association_df.to_csv(
        output_file,
        index=False
    )

    return output_file