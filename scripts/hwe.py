import os
import numpy as np
import pandas as pd

--------------------------------------------------------------
#exact_hwe_wigginton
-------------------------------------------------------------

def exact_hwe_wigginton(
    obs_hom1: int,
    obs_het: int,
    obs_hom2: int
) -> float:
  
    obs_homc = max(obs_hom1, obs_hom2)
    obs_homr = min(obs_hom1, obs_hom2)

    rare = 2 * obs_homr + obs_het
    genotypes = obs_het + obs_homc + obs_homr

    if genotypes == 0:
        return np.nan

    probs = np.zeros(rare + 1)

    mid = int(
        rare * (2 * genotypes - rare) /
        (2 * genotypes)
    )

    if (rare & 1) != (mid & 1):
        mid += 1

    probs[mid] = 1.0
    total = 1.0

    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = genotypes - curr_hets - curr_homr

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

        probs[curr_hets - 2] = prob

        total += prob

        curr_hets -= 2
        curr_homr += 1
        curr_homc += 1

    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = genotypes - curr_hets - curr_homr

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

        probs[curr_hets + 2] = prob

        total += prob

        curr_hets += 2
        curr_homr -= 1
        curr_homc -= 1

    probs /= total

    p = probs[
        probs <= probs[obs_het]
    ].sum()

    return min(1.0, p)


--------------------------------------------------------------
#prepare_control_genotypes
------------------------------------------------------------
def prepare_control_genotypes(
    genotype_df: pd.DataFrame,
    control_samples: set
) -> pd.DataFrame:
    """
    Prepare control genotypes for
    Hardy–Weinberg equilibrium analysis.
    """

    control_df = genotype_df[

        genotype_df["Sample"].isin(
            control_samples
        )

    ].copy()

    control_df = control_df[

        control_df["Zygosity"] !=
        "LowConfidence_Heterozygous"

    ].copy()

    return control_df


-------------------------------------------------------------------------------
#calculate_variant_hwe
------------------------------------------------------------------------------
calculate_variant_hwe(
    variant,
    control_df
)

---------------------------------------------------------------------------
#run_hwe_analysis
--------------------------------------------------------------------

def run_hwe_analysis(

    final_df,

    control_df

):

    rows = []

    for _, variant in final_df.iterrows():

        rows.append(

            calculate_variant_hwe(

                variant,

                control_df

            )

        )

    hwe_df = pd.DataFrame(rows)

    ...

def merge_variant_annotation(

    hwe_df,

    final_df

):

    ...

def reorder_hwe_columns(

    hwe_df

):


--------------------------------------------------------
#save_hwe_results
-------------------------------------------------------

def save_hwe_results(

    hwe_df,

    output_folder

):
