"""
annotation.py

Functions for functional annotation of validated
variants using the Ensembl Variant Effect Predictor
(VEP) REST API.

The Ensembl VEP REST API is the default annotation
method. An optional web-exported VEP file can be
supported by the main pipeline when configured.

The parser preserves the principal VEP transcript-level
fields required for downstream final-table generation
and validation, including HGVSc and HGVSp.
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
        Input HGVS cDNA variant.

    result : dict
        Single variant result returned by the
        Ensembl VEP REST API.

    Returns
    -------
    dict
        Parsed VEP annotation.
    """

    # --------------------------------------------------
    # Select transcript consequence
    # --------------------------------------------------

    transcript_consequences = result.get(
        "transcript_consequences",
        []
    )

    transcript = {}

    if transcript_consequences:

        # Prefer canonical transcript where available
        canonical_transcripts = [

            item

            for item in transcript_consequences

            if item.get("canonical") == 1

        ]

        if canonical_transcripts:

            transcript = canonical_transcripts[0]

        else:

            transcript = transcript_consequences[0]

    # --------------------------------------------------
    # Basic genomic information
    # --------------------------------------------------

    chromosome = result.get(
        "seq_region_name"
    )

    genomic_position = result.get(
        "start"
    )

    genomic_end = result.get(
        "end"
    )

    genomic_coordinate = None

    if chromosome and genomic_position:

        if genomic_end and genomic_end != genomic_position:

            genomic_coordinate = (

                f"chr{chromosome}:"
                f"{genomic_position}-"
                f"{genomic_end}"

            )

        else:

            genomic_coordinate = (

                f"chr{chromosome}:"
                f"{genomic_position}"

            )

    # --------------------------------------------------
    # VEP transcript-level HGVS
    # --------------------------------------------------

    hgvsc = transcript.get(
        "hgvsc"
    )

    hgvsp = transcript.get(
        "hgvsp"
    )

    # --------------------------------------------------
    # Consequence
    # --------------------------------------------------

    consequence_terms = transcript.get(
        "consequence_terms",
        []
    )

    if isinstance(
        consequence_terms,
        list
    ):

        consequence = ",".join(
            consequence_terms
        )

    else:

        consequence = str(
            consequence_terms
        )

    # --------------------------------------------------
    # Return annotation
    # --------------------------------------------------

    return {

        # ==============================================
        # Input / status
        # ==============================================

        "HGVS_cDNA":
            hgvs,

        "Status":
            "Annotated",

        # ==============================================
        # Assembly / genomic annotation
        # ==============================================

        "Assembly":
            result.get(
                "assembly_name"
            ),

        "Chromosome":
            chromosome,

        "Genomic_Position":
            genomic_position,

        "End_Position":
            genomic_end,

        "Strand":
            result.get(
                "strand"
            ),

        "Genomic_Coordinate":
            genomic_coordinate,

        "HGVS_Genomic":
            result.get(
                "hgvsg"
            ),

        # ==============================================
        # Variant identification
        # ==============================================

        "dbSNP_rsID":
            result.get(
                "id"
            ),

        "Variant_Class":
            result.get(
                "variant_class"
            ),

        "Allele_String":
            result.get(
                "allele_string"
            ),

        "Most_Severe_Consequence":
            result.get(
                "most_severe_consequence"
            ),

        # ==============================================
        # Gene annotation
        # ==============================================

        "Gene":
            transcript.get(
                "gene_symbol"
            ),

        "Gene_ID":
            transcript.get(
                "gene_id"
            ),

        "Gene_Symbol":
            transcript.get(
                "gene_symbol"
            ),

        "HGNC_ID":
            transcript.get(
                "hgnc_id"
            ),

        "Gene_Symbol_Source":
            transcript.get(
                "gene_symbol_source"
            ),

        # ==============================================
        # Transcript annotation
        # ==============================================

        "Transcript":
            transcript.get(
                "transcript_id"
            ),

        "Feature":
            transcript.get(
                "transcript_id"
            ),

        "Biotype":
            transcript.get(
                "biotype"
            ),

        "Exon":
            transcript.get(
                "exon"
            ),

        "Canonical":
            transcript.get(
                "canonical"
            ),

        # ==============================================
        # MANE annotation
        # ==============================================

        "MANE_Select":
            transcript.get(
                "mane_select"
            ),

        "MANE":
            ",".join(
                transcript.get(
                    "mane",
                    []
                )
            )
            if isinstance(
                transcript.get(
                    "mane",
                    []
                ),
                list
            )
            else transcript.get(
                "mane"
            ),

        "RefSeq_Transcript":
            ",".join(
                transcript.get(
                    "refseq_transcript_ids",
                    []
                )
            )
            if isinstance(
                transcript.get(
                    "refseq_transcript_ids",
                    []
                ),
                list
            )
            else transcript.get(
                "refseq_transcript_ids"
            ),

        # ==============================================
        # HGVS annotation
        # ==============================================

        # IMPORTANT:
        # These are the exact fields expected by
        # downstream validation.

        "HGVSc":
            hgvsc,

        "HGVSp":
            hgvsp,

        # Retain the older descriptive field as well
        # for compatibility with existing pipeline
        # components.

        "Protein_HGVS":
            hgvsp,

        # ==============================================
        # Consequence / impact
        # ==============================================

        "Consequence":
            consequence,

        "Impact":
            transcript.get(
                "impact"
            ),

        # ==============================================
        # Protein annotation
        # ==============================================

        "Protein_ID":
            transcript.get(
                "protein_id"
            ),

        "Protein_Position":
            transcript.get(
                "protein_start"
            ),

        "Protein_End":
            transcript.get(
                "protein_end"
            ),

        "Amino_Acid_Change":
            transcript.get(
                "amino_acids"
            ),

        "Codon_Change":
            transcript.get(
                "codons"
            ),

        # ==============================================
        # Coding / cDNA coordinates
        # ==============================================

        "CDS_Position":
            transcript.get(
                "cds_start"
            ),

        "CDS_End":
            transcript.get(
                "cds_end"
            ),

        "cDNA_Position":
            transcript.get(
                "cdna_start"
            ),

        "cDNA_End":
            transcript.get(
                "cdna_end"
            ),

        # ==============================================
        # Variant allele
        # ==============================================

        "Variant_Allele":
            transcript.get(
                "variant_allele"
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
    retries: int = 3,
    delay: float = 0.1,
    transcript: str = None
) -> pd.DataFrame:
    """
    Annotate validated variants using the Ensembl VEP REST API.

    Parameters
    ----------
    hgvs_df : pandas.DataFrame
        DataFrame containing validated variants and HGVS cDNA
        identifiers.

    server : str
        Ensembl REST API server URL.

    headers : dict
        HTTP headers used for VEP API requests.

    timeout : int
        Maximum number of seconds allowed for each API request.

    retries : int
        Maximum number of retry attempts for failed API requests.

    delay : float
        Delay in seconds between consecutive API requests.

    transcript : str, optional
        RefSeq transcript accession used for the analysis.

    Returns
    -------
    pandas.DataFrame
        Annotated variant table containing VEP consequences,
        transcript information, HGVSc, HGVSp, gene information,
        variant class, genomic coordinates, and other available
        VEP annotations.
    """

    annotation_rows = []

    # --------------------------------------------------
    # Validate input
    # --------------------------------------------------

    if "HGVS_cDNA" not in hgvs_df.columns:

        raise ValueError(
            "Input HGVS table must contain "
            "'HGVS_cDNA' column."
        )

    # --------------------------------------------------
    # Process variants
    # --------------------------------------------------

    for _, row in hgvs_df.iterrows():

        hgvs_entries = str(
            row["HGVS_cDNA"]
        ).split(";")

        for hgvs in hgvs_entries:

            hgvs = hgvs.strip()

            if not hgvs:
                continue

            # ------------------------------------------
            # Build VEP REST URL
            # ------------------------------------------

            encoded_hgvs = quote(
                hgvs,
                safe=""
            )

            url = (
                f"{server.rstrip('/')}"
                f"/vep/human/hgvs/"
                f"{encoded_hgvs}"
            )

            response = None

            # ------------------------------------------
            # API request with retry
            # ------------------------------------------

            for attempt in range(retries):

                try:

                    response = requests.get(

                        url,

                        headers=headers,

                        params={
                            "hgvs": 1,
                            "mane": 1,
                            "canonical": 1,
                            "pick": 1,
                            "variant_class": 1,
                            "protein": 1,
                            "numbers": 1,
                            "xref_refseq": 1
                        },

                        timeout=timeout
                    )

                    if response.status_code == 200:
                        break

                except requests.RequestException:

                    response = None

                # --------------------------------------
                # Retry delay
                # --------------------------------------

                if attempt < retries - 1:
                    time.sleep(1)

            # ------------------------------------------
            # Failed request
            # ------------------------------------------

            if (
                response is None
                or response.status_code != 200
            ):

                status_code = (
                    response.status_code
                    if response is not None
                    else None
                )

                annotation_rows.append({

                    "HGVS_cDNA": hgvs,

                    "Status": "Failed",

                    "HTTP_Status": status_code

                })

                continue

            # ------------------------------------------
            # Parse JSON
            # ------------------------------------------

            try:

                results = response.json()

            except ValueError:

                annotation_rows.append({

                    "HGVS_cDNA": hgvs,

                    "Status": "Invalid_JSON"

                })

                continue

            # ------------------------------------------
            # Empty result
            # ------------------------------------------

            if not results:

                annotation_rows.append({

                    "HGVS_cDNA": hgvs,

                    "Status": "Not_Found"

                })

                continue

            # ------------------------------------------
            # VEP returns a list of variant objects
            # ------------------------------------------

            result = results[0]

            annotation = parse_vep_response(
                hgvs,
                result
            )

            # Preserve HTTP status for
            # troubleshooting/reproducibility.

            annotation["HTTP_Status"] = (
                response.status_code
            )

            annotation_rows.append(
                annotation
            )

            # ------------------------------------------
            # Respect API request delay
            # ------------------------------------------

            time.sleep(delay)

    # --------------------------------------------------
    # Construct annotation DataFrame
    # --------------------------------------------------

    annotation_df = pd.DataFrame(
        annotation_rows
    )

    # --------------------------------------------------
    # Remove duplicate HGVS entries
    # --------------------------------------------------

    if not annotation_df.empty:

        annotation_df = (
            annotation_df
            .drop_duplicates(
                subset="HGVS_cDNA"
            )
            .reset_index(
                drop=True
            )
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
    Save Ensembl VEP annotation table.

    Parameters
    ----------
    annotation_df : pandas.DataFrame

    output_folder : str
    """

    os.makedirs(

        output_folder,

        exist_ok=True

    )

    output_file = os.path.join(

        output_folder,

        "Ensembl_VEP_Annotation.csv"

    )

    annotation_df.to_csv(

        output_file,

        index=False

    )

    return output_file