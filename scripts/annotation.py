"""
annotation.py

Functions for functional annotation of validated
variants using the Ensembl Variant Effect Predictor
(VEP) REST API.
"""

import os
import time
import requests
import pandas as pd

from urllib.parse import quote


# ======================================================
# PARSE VEP RESPONSE
# ======================================================

def parse_vep_response(
    hgvs: str,
    result: dict
) -> dict:
    """
    Extract useful information from a VEP response.

    Parameters
    ----------
    hgvs : str

    result : dict

    Returns
    -------
    dict
    """

    transcript = {}

    if result.get("transcript_consequences"):

        transcript = result["transcript_consequences"][0]

    chromosome = result.get("seq_region_name")

    genomic_position = result.get("start")

    genomic_coordinate = None

    if chromosome and genomic_position:

        genomic_coordinate = (

            f"chr{chromosome}:{genomic_position}"

        )

    return {

        "HGVS_cDNA":

            hgvs,

        "Status":

            "Annotated",

        "Assembly":

            result.get("assembly_name"),

        "Gene":

            transcript.get("gene_symbol"),

        "Gene_ID":

            transcript.get("gene_id"),

        "Transcript":

            transcript.get("transcript_id"),

        "Chromosome":

            chromosome,

        "Genomic_Position":

            genomic_position,

        "End_Position":

            result.get("end"),

        "Strand":

            result.get("strand"),

        "Genomic_Coordinate":

            genomic_coordinate,

        "HGVS_Genomic":

            result.get("hgvsg"),

        "dbSNP_rsID":

            result.get("id"),

        "Variant_Class":

            result.get("variant_class"),

        "Most_Severe_Consequence":

            result.get(

                "most_severe_consequence"

            ),

        "Consequence":

            ",".join(

                transcript.get(

                    "consequence_terms",

                    []

                )

            ),

        "Impact":

            transcript.get("impact"),

        "Protein_HGVS":

            transcript.get("hgvsp"),

        "Protein_Position":

            transcript.get(

                "protein_start"

            ),

        "Amino_Acid_Change":

            transcript.get(

                "amino_acids"

            ),

        "Codon_Change":

            transcript.get(

                "codons"

            ),

        "CDS_Position":

            transcript.get(

                "cds_start"

            ),

        "cDNA_Position":

            transcript.get(

                "cdna_start"

            )

    }


# ======================================================
# ENSEMBL VEP ANNOTATION
# ======================================================

def annotate_variants(
    hgvs_df: pd.DataFrame,
    server: str,
    headers: dict,
    timeout: int = 30,
    retries: int = 2,
    delay: float = 0.05
) -> pd.DataFrame:
    """
    Annotate HGVS variants using the Ensembl
    REST API.

    Parameters
    ----------
    hgvs_df : pandas.DataFrame

    server : str

    headers : dict

    Returns
    -------
    pandas.DataFrame
    """

    annotation_rows = []

    for _, row in hgvs_df.iterrows():

        hgvs_entries = str(

            row["HGVS_cDNA"]

        ).split(";")

        for hgvs in hgvs_entries:

            hgvs = hgvs.strip()

            if not hgvs:

                continue

            url = (

                f"{server}/vep/human/hgvs/"

                f"{quote(hgvs, safe='')}"

            )

            response = None

            for _ in range(retries):

                try:

                    response = requests.get(

                        url,

                        headers=headers,

                        timeout=timeout

                    )

                    if response.status_code == 200:

                        break

                except requests.RequestException:

                    pass

                time.sleep(1)

            if response is None or response.status_code != 200:

                annotation_rows.append({

                    "HGVS_cDNA": hgvs,

                    "Status": "Failed"

                })

                continue

            results = response.json()

            if not results:

                annotation_rows.append({

                    "HGVS_cDNA": hgvs,

                    "Status": "Not_Found"

                })

                continue

            annotation_rows.append(

                parse_vep_response(

                    hgvs,

                    results[0]

                )

            )

            time.sleep(delay)

    annotation_df = pd.DataFrame(

        annotation_rows

    )

    annotation_df = (

        annotation_df

        .drop_duplicates(

            subset="HGVS_cDNA"

        )

        .reset_index(drop=True)

    )

    return annotation_df


# ======================================================
# SAVE TABLE
# ======================================================

def save_annotation_table(
    annotation_df: pd.DataFrame,
    output_folder: str
):
    """
    Save Ensembl annotation table.
    """

    os.makedirs(

        output_folder,

        exist_ok=True

    )

    annotation_df.to_csv(

        os.path.join(

            output_folder,

            "Ensembl_VEP_Annotation.csv"

        ),

        index=False

    )
