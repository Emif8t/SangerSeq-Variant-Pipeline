import os

import numpy as np
import pandas as pd


# =========================================================
# Exact Hardy-Weinberg Equilibrium
# Wigginton et al. (2005) algorithm
# =========================================================

def exact_hwe_wigginton(
    obs_hom1: int,
    obs_het: int,
    obs_hom2: int
) -> float:
    """
    Calculate the exact Hardy-Weinberg equilibrium
    P-value using the Wigginton et al. (2005) algorithm.
    """

    obs_homc = max(
        obs_hom1,
        obs_hom2
    )

    obs_homr = min(
        obs_hom1,
        obs_hom2
    )

    rare = (
        2 * obs_homr
        + obs_het
    )

    genotypes = (
        obs_het
        + obs_homc
        + obs_homr
    )

    if genotypes == 0:
        return np.nan

    probs = np.zeros(
        rare + 1
    )

    mid = int(
        rare
        * (2 * genotypes - rare)
        / (2 * genotypes)
    )

    if (rare & 1) != (mid & 1):
        mid += 1

    probs[mid] = 1.0

    total = 1.0

    curr_hets = mid

    curr_homr = (
        rare - mid
    ) // 2

    curr_homc = (
        genotypes
        - curr_hets
        - curr_homr
    )

    # -----------------------------------------------------
    # Move downward through possible heterozygote counts
    # -----------------------------------------------------

    while curr_hets >= 2:

        prob = (
            probs[curr_hets]
            * curr_hets
            * (curr_hets - 1)
        ) / (
            4
            * (curr_homr + 1)
            * (curr_homc + 1)
        )

        probs[
            curr_hets - 2
        ] = prob

        total += prob

        curr_hets -= 2

        curr_homr += 1

        curr_homc += 1

    # -----------------------------------------------------
    # Move upward through possible heterozygote counts
    # -----------------------------------------------------

    curr_hets = mid

    curr_homr = (
        rare - mid
    ) // 2

    curr_homc = (
        genotypes
        - curr_hets
        - curr_homr
    )

    while curr_hets <= rare - 2:

        prob = (
            probs[curr_hets]
            * 4
            * curr_homr
            * curr_homc
        ) / (
            (curr_hets + 2)
            * (curr_hets + 1)
        )

        probs[
            curr_hets + 2
        ] = prob

        total += prob

        curr_hets += 2

        curr_homr -= 1

        curr_homc -= 1

    # -----------------------------------------------------
    # Normalize probabilities
    # -----------------------------------------------------

    probs /= total

    p = probs[
        probs <= probs[obs_het]
    ].sum()

    return min(
        1.0,
        p
    )


# =========================================================
# Prepare control genotypes
# =========================================================

def prepare_control_genotypes(
    genotype_df: pd.DataFrame,
    control_samples: set
) -> pd.DataFrame:
    """
    Prepare control genotypes for
    Hardy-Weinberg equilibrium analysis.

    Only control samples are retained and
    low-confidence heterozygotes are removed.
    """

    control_df = genotype_df[
        genotype_df["Sample"].isin(
            control_samples
        )
    ].copy()

    control_df = control_df[
        control_df["Zygosity"]
        != "LowConfidence_Heterozygous"
    ].copy()

    return control_df


# =========================================================
# Calculate HWE for one variant
# =========================================================

def calculate_variant_hwe(
    variant,
    control_df
):
    """
    Calculate exact Hardy-Weinberg equilibrium
    for a single variant.
    """

    pos = variant[
        "Transcript_Position"
    ]

    subset = control_df[
        control_df["cDNA_Position"] == pos
    ]

    # -----------------------------------------------------
    # Genotype counts
    # -----------------------------------------------------

    AA = (
        subset["Genotype"] == 0
    ).sum()

    AB = (
        subset["Genotype"] == 1
    ).sum()

    BB = (
        subset["Genotype"] == 2
    ).sum()

    # -----------------------------------------------------
    # Exact HWE test
    # -----------------------------------------------------

    p = exact_hwe_wigginton(
        AA,
        AB,
        BB
    )

    # -----------------------------------------------------
    # Return result
    # -----------------------------------------------------

    return {

        "Transcript_Position": pos,

        "Controls": (
            AA + AB + BB
        ),

        "AA": AA,

        "AB": AB,

        "BB": BB,

        "Exact_HWE_P": p,

        "HWE_Status":
            (
                "In Equilibrium"
                if (
                    pd.notna(p)
                    and p >= 0.05
                )
                else "Deviation"
            )
    }


# =========================================================
# Run HWE analysis
# =========================================================

def run_hwe_analysis(
    final_df,
    control_df
):
    """
    Perform exact Hardy-Weinberg equilibrium
    analysis for every variant.
    """

    rows = []

    for _, variant in final_df.iterrows():

        rows.append(
            calculate_variant_hwe(
                variant,
                control_df
            )
        )

    return pd.DataFrame(
        rows
    )


# =========================================================
# Merge variant annotation
# =========================================================

def merge_variant_annotation(
    hwe_df,
    final_df
):
    """
    Merge HWE results with variant annotation.
    """

    hwe_df = hwe_df.merge(
        final_df,
        on="Transcript_Position",
        how="left"
    )

    return hwe_df


# =========================================================
# Reorder HWE columns
# =========================================================

def reorder_hwe_columns(
    hwe_df
):
    """
    Reorder HWE results into the
    final reporting structure.
    """

    cols = [

        "Gene",

        "Transcript",

        "Transcript_Position",

        "HGVS_cDNA_Position",

        "HGVS_cDNA",

        "Chromosome",

        "Genomic_Position",

        "REF",

        "ALT",

        "Consequence",

        "Impact",

        "Controls",

        "AA",

        "AB",

        "BB",

        "Exact_HWE_P",

        "HWE_Status"
    ]

    return hwe_df[
        [
            c
            for c in cols
            if c in hwe_df.columns
        ]
    ]


# =========================================================
# Sort HWE results
# =========================================================

def sort_hwe_results(
    hwe_df
):
    """
    Sort HWE results by transcript position.
    """

    return hwe_df.sort_values(
        "Transcript_Position"
    ).reset_index(
        drop=True
    )


# =========================================================
# Complete HWE analysis pipeline
# =========================================================

def generate_hwe_results(
    final_df,
    genotype_df,
    control_samples
):
    """
    Generate the complete HWE results table.
    """

    control_df = prepare_control_genotypes(
        genotype_df,
        control_samples
    )

    hwe_df = run_hwe_analysis(
        final_df,
        control_df
    )

    hwe_df = merge_variant_annotation(
        hwe_df,
        final_df
    )

    hwe_df = reorder_hwe_columns(
        hwe_df
    )

    hwe_df = sort_hwe_results(
        hwe_df
    )

    return hwe_df


# =========================================================
# Save HWE results
# =========================================================

def save_hwe_results(
    hwe_df,
    output_folder
):
    """
    Save exact HWE results as CSV.
    """

    os.makedirs(
        output_folder,
        exist_ok=True
    )

    output_file = os.path.join(
        output_folder,
        "ASS1_Exact_HWE_Wigginton2.csv"
    )

    hwe_df.to_csv(
        output_file,
        index=False
    )

    return output_file
