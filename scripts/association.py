import os 
import numpy as np 
import pandas as pd

--------------------------------------------
#load_sample_groups
-------------------------------------------

def load_sample_groups(
    phenotype_file: str
) -> pd.DataFrame:
    """
    Load sample phenotype information.
    """

    return pd.read_excel(
        phenotype_file
    )

----------------------------------------
#create_sample_sets
---------------------------------------

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

-----------------------------------------------
#analyse_variant
----------------------------------------------

analyse_variant(
    row,
    case_samples,
    control_samples
)

------------------------------------------
#run_association_analysis
-----------------------------------------

def run_association_analysis(

    final_df,

    case_samples,

    control_samples

):

    results = []

    for _, row in final_df.iterrows():

        results.append(

            analyse_variant(

                row,

                case_samples,

                control_samples

            )

        )

    return pd.DataFrame(results)


----------------------------------------
#apply_multiple_testing
---------------------------------------

def apply_multiple_testing(
    association_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Apply Benjamini–Hochberg FDR correction.
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

----------------------------------------------
#save_association_results
--------------------------------------------

def save_association_results(

    association_df,

    output_folder

):

    os.makedirs(

        output_folder,

        exist_ok=True

    )

    association_df.to_csv(

        os.path.join(

            output_folder,

            "ASS1_Association_Analysis.csv"

        ),

        index=False

    )
