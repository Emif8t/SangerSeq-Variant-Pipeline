"""
hgvs.py

Functions for generating HGVS nomenclature for
validated sequence variants.
"""

import os
import pandas as pd


# ======================================================
# GENERATE HGVS TABLE
# ======================================================

def generate_hgvs_table(
    summary_df: pd.DataFrame,
    variant_df: pd.DataFrame,
    transcript: str,
    cds_start: int
) -> pd.DataFrame:
    """
    Generate HGVS cDNA nomenclature for all
    validated variants.

    Parameters
    ----------
    summary_df : pandas.DataFrame

    variant_df : pandas.DataFrame

    transcript : str

    cds_start : int

    Returns
    -------
    pandas.DataFrame
    """

    summary_df = summary_df.copy()

    summary_df["HGVS_cDNA_Position"] = (

        summary_df["cDNA_Position"]

        - cds_start

        + 1

    ).astype(int)

    hgvs_rows = []

    for _, row in summary_df.iterrows():

        transcript_position = int(

            row["cDNA_Position"]

        )

        hgvs_position = int(

            row["HGVS_cDNA_Position"]

        )

        ref = row["REF"]

        alt_list = (

            variant_df.loc[

                variant_df["cDNA_Position"]

                == transcript_position,

                "ALT"

            ]

            .dropna()

            .astype(str)

            .unique()

        )

        alt_list = sorted({

            alt

            for alt in alt_list

            if alt != ref

        })

        hgvs_list = [

            f"{transcript}:c.{hgvs_position}{ref}>{alt}"

            for alt in alt_list

        ]

        if len(alt_list) == 1:

            variant_type = "SNV"

        elif len(alt_list) > 1:

            variant_type = "Multiple_ALTs"

        else:

            variant_type = "Reference"

        hgvs_rows.append({

            "Transcript":

                transcript,

            "Transcript_Position":

                transcript_position,

            "HGVS_cDNA_Position":

                hgvs_position,

            "HGVS_cDNA":

                ";".join(hgvs_list),

            "REF":

                ref,

            "ALT":

                ",".join(alt_list),

            "Variant_Type":

                variant_type,

            "Carrier_Count":

                row["Carrier_Count"],

            "Variant_Calls":

                row["Variant_Calls"],

            "Mean_Quality":

                round(

                    row["Mean_Quality"],

                    2

                ),

            "Mean_Alignment_Score":

                round(

                    row["Mean_Alignment_Score"],

                    2

                ),

            "Samples":

                row["Samples"]

        })

    hgvs_df = pd.DataFrame(hgvs_rows)

    hgvs_df = hgvs_df.sort_values(

        "HGVS_cDNA_Position"

    ).reset_index(drop=True)

    return hgvs_df


# ======================================================
# SAVE HGVS TABLE
# ======================================================

def save_hgvs_table(
    hgvs_df: pd.DataFrame,
    output_folder: str
):
    """
    Save HGVS table.

    Parameters
    ----------
    hgvs_df : pandas.DataFrame

    output_folder : str
    """

    os.makedirs(

        output_folder,

        exist_ok=True

    )

    hgvs_df.to_csv(

        os.path.join(

            output_folder,

            "HGVS_Table.csv"

        ),

        index=False

    )
