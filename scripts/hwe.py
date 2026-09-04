import os

import numpy as np
import pandas as pd


# =========================================================
# Exact Hardy-Weinberg Equilibrium
# Wigginton et al. (2005)
# =========================================================


def exact_hwe_wigginton(
    obs_hom1: int,
    obs_het: int,
    obs_hom2: int
) -> float:
    """
    Calculate the exact Hardy-Weinberg equilibrium P-value
    using the Wigginton et al. (2005) algorithm.

    Parameters
    ----------
    obs_hom1 : int
        Number of homozygous reference individuals.

    obs_het : int
        Number of heterozygous individuals.

    obs_hom2 : int
        Number of homozygous alternate individuals.

    Returns
    -------
    float
        Exact HWE P-value.

        Returns NaN if no genotype observations are available.
    """

    if min(
        obs_hom1,
        obs_het,
        obs_hom2
    ) < 0:
        raise ValueError(
            "Genotype counts cannot be negative."
        )

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
        rare + 1,
        dtype=float
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

    p_value = probs[
        probs <= probs[obs_het]
    ].sum()

    return min(
        1.0,
        p_value
    )


# =========================================================
# Prepare control genotypes
# =========================================================


def prepare_control_genotypes(
    genotype_df: pd.DataFrame,
    control_samples: set,
    min_phred: int = 20
) -> pd.DataFrame:
    """
    Prepare high-confidence control genotypes for HWE analysis.

    Only control samples with valid genotype calls and
    sufficient sequencing quality are retained.

    Parameters
    ----------
    genotype_df : pandas.DataFrame
        Genotype calls.

    control_samples : set
        Sample identifiers classified as controls.

    min_phred : int, default=20
        Minimum Phred quality threshold.

    Returns
    -------
    pandas.DataFrame
        Filtered control genotype calls.
    """

    required_columns = {
        "Sample",
        "cDNA_Position",
        "Genotype",
        "Quality"
    }

    missing_columns = (
        required_columns
        - set(genotype_df.columns)
    )

    if missing_columns:
        raise ValueError(
            "Missing required genotype columns: "
            f"{sorted(missing_columns)}"
        )

    control_df = genotype_df[
        genotype_df["Sample"].isin(
            control_samples
        )
    ].copy()

    control_df["Quality"] = pd.to_numeric(
        control_df["Quality"],
        errors="coerce"
    )

    control_df["Genotype"] = pd.to_numeric(
        control_df["Genotype"],
        errors="coerce"
    )

    # Keep only genotype calls meeting the
    # minimum sequencing-quality threshold.
    control_df = control_df[
        control_df["Quality"] >= min_phred
    ].copy()

    # Keep only valid diploid genotype classifications:
    # 0 = homozygous reference
    # 1 = heterozygous
    # 2 = homozygous alternate
    control_df = control_df[
        control_df["Genotype"].isin(
            [0, 1, 2]
        )
    ].copy()

    return control_df


# =========================================================
# Calculate HWE for one variant
# =========================================================


def calculate_variant_hwe(
    variant: pd.Series,
    control_df: pd.DataFrame
) -> dict:
    """
    Calculate exact HWE for a single variant.

    Genotypes are matched using transcript position,
    reference allele and alternate allele where available.

    HWE is calculated using callable control genotypes only.
    """

    position = variant[
        "Transcript_Position"
    ]

    subset = control_df[
        control_df["cDNA_Position"] == position
    ].copy()

    # -----------------------------------------------------
    # Match REF/ALT where these columns are available
    # -----------------------------------------------------

    if (
        "REF" in variant.index
        and "ALT" in variant.index
        and "REF" in control_df.columns
        and "ALT" in control_df.columns
    ):

        ref = variant["REF"]
        alt = variant["ALT"]

        if pd.notna(ref) and pd.notna(alt):

            subset = subset[
                (
                    subset["REF"].astype(str)
                    == str(ref)
                )
                & (
                    subset["ALT"].astype(str)
                    == str(alt)
                )
            ]

    # -----------------------------------------------------
    # Remove duplicate calls per sample
    #
    # Keep the highest-quality genotype call if multiple
    # records exist for the same sample and variant.
    # -----------------------------------------------------

    if not subset.empty:

        subset = subset.sort_values(
            "Quality",
            ascending=False
        )

        subset = subset.drop_duplicates(
            subset=["Sample"],
            keep="first"
        )

    # -----------------------------------------------------
    # Genotype counts
    # -----------------------------------------------------

    AA = int(
        (subset["Genotype"] == 0).sum()
    )

    AB = int(
        (subset["Genotype"] == 1).sum()
    )

    BB = int(
        (subset["Genotype"] == 2).sum()
    )

    controls = (
        AA
        + AB
        + BB
    )

    # -----------------------------------------------------
    # HWE note
    # -----------------------------------------------------

    if controls == 0:
        variant_label = variant.get(
            "HGVS_cDNA",
            "unknown variant"
        )

        print(
            "WARNING: No callable control genotypes were "
            f"available for {variant_label}. "
            "HWE is not estimable."
        )

        hwe_note = (
            "No callable control genotypes available "
            "after quality filtering."
        )
    else:
        hwe_note = ""

    # -----------------------------------------------------
    # Exact HWE test
    # -----------------------------------------------------

    p_value = exact_hwe_wigginton(
        AA,
        AB,
        BB
    )

    # -----------------------------------------------------
    # Interpretation
    # -----------------------------------------------------

    if pd.isna(p_value):

        hwe_status = "Not estimable"

    elif p_value < 0.05:

        hwe_status = "Deviation"

    else:

        hwe_status = "In Equilibrium"

    # -----------------------------------------------------
    # Return result
    # -----------------------------------------------------

    return {
        "Transcript_Position": position,
        "Controls": controls,
        "AA": AA,
        "AB": AB,
        "BB": BB,
        "Exact_HWE_P": p_value,
        "HWE_Status": hwe_status,
        "HWE_Note": hwe_note
    }


# =========================================================
# Run HWE analysis
# =========================================================


def run_hwe_analysis(
    final_df: pd.DataFrame,
    control_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Perform exact HWE analysis for every variant.
    """

    if final_df.empty:

        return pd.DataFrame()

    rows = []

    for _, variant in final_df.iterrows():

        rows.append(
            calculate_variant_hwe(
                variant,
                control_df
            )
        )

    return pd.DataFrame(rows)


# =========================================================
# Create stable variant key
# =========================================================


def add_variant_key(
    df: pd.DataFrame
) -> pd.DataFrame:
    """
    Create a stable variant key from transcript position,
    reference allele and alternate allele.
    """

    df = df.copy()

    required = {
        "Transcript_Position",
        "REF",
        "ALT"
    }

    if required.issubset(df.columns):

        df["Variant_Key"] = (
            df["Transcript_Position"]
            .astype(str)
            + ":"
            + df["REF"].astype(str)
            + ">"
            + df["ALT"].astype(str)
        )

    else:

        df["Variant_Key"] = (
            df["Transcript_Position"]
            .astype(str)
        )

    return df


# =========================================================
# Merge variant annotation
# =========================================================


def merge_variant_annotation(
    hwe_df: pd.DataFrame,
    final_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge HWE results with variant annotation using
    a stable variant key.
    """

    hwe_df = add_variant_key(
        hwe_df
    )

    annotation_df = add_variant_key(
        final_df
    )

    annotation_columns = [
        column
        for column in annotation_df.columns
        if column != "Variant_Key"
    ]

    annotation_df = annotation_df[
        ["Variant_Key"]
        + annotation_columns
    ]

    hwe_df = hwe_df.merge(
        annotation_df,
        on="Variant_Key",
        how="left",
        suffixes=(
            "",
            "_annotation"
        )
    )

    return hwe_df


# =========================================================
# Reorder HWE columns
# =========================================================


def reorder_hwe_columns(
    hwe_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Reorder HWE results into the final reporting structure.
    """

    columns = [

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

        "HWE_Status",

        "HWE_Note"
    ]

    return hwe_df[
        [
            column
            for column in columns
            if column in hwe_df.columns
        ]
    ]


# =========================================================
# Sort HWE results
# =========================================================


def sort_hwe_results(
    hwe_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Sort HWE results by transcript position.
    """

    if hwe_df.empty:

        return hwe_df

    return hwe_df.sort_values(
        "Transcript_Position"
    ).reset_index(
        drop=True
    )


# =========================================================
# Complete HWE analysis pipeline
# =========================================================


def generate_hwe_results(
    final_df: pd.DataFrame,
    genotype_df: pd.DataFrame,
    control_samples: set,
    min_phred: int = 20
) -> pd.DataFrame:
    """
    Generate the complete HWE results table.

    HWE is evaluated in controls only.

    Only genotype calls meeting the minimum Phred
    quality threshold and having valid genotype
    classifications are included.
    """

    control_df = prepare_control_genotypes(
        genotype_df,
        control_samples,
        min_phred=min_phred
    )

    hwe_df = run_hwe_analysis(
        final_df,
        control_df
    )

    if hwe_df.empty:

        return hwe_df

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
    hwe_df: pd.DataFrame,
    output_folder: str
) -> str:
    """
    Save Hardy-Weinberg equilibrium results.
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