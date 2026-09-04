"""
annotation.py

Functions for functional annotation of validated variants using
the Ensembl Variant Effect Predictor (VEP) REST API.

The Ensembl VEP REST API is the primary annotation method.

If an API annotation fails, the pipeline can optionally fall back
to an existing VEP web-exported Excel file. This provides a
reproducible recovery mechanism for temporary API/server failures.
"""

import os
import time
from urllib.parse import quote

import pandas as pd
import requests


# ============================================================
# VEP RESPONSE PARSER
# ============================================================

def parse_vep_response(
    hgvs: str,
    result: dict
) -> dict:
    """
    Extract useful information from a VEP API response.

    Parameters
    ----------
    hgvs : str
        Input HGVS cDNA variant.

    result : dict
        Single variant result returned by the Ensembl VEP API.

    Returns
    -------
    dict
        Parsed VEP annotation.
    """

    transcript_consequences = result.get(
        "transcript_consequences",
        []
    )

    transcript = {}

    if transcript_consequences:

        canonical_transcripts = [
            item
            for item in transcript_consequences
            if item.get("canonical") == 1
        ]

        if canonical_transcripts:
            transcript = canonical_transcripts[0]
        else:
            transcript = transcript_consequences[0]

    # --------------------------------------------------------
    # Genomic information
    # --------------------------------------------------------

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

        if (
            genomic_end
            and genomic_end != genomic_position
        ):
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

    # --------------------------------------------------------
    # Transcript-level HGVS
    # --------------------------------------------------------

    hgvsc = transcript.get(
        "hgvsc"
    )

    hgvsp = transcript.get(
        "hgvsp"
    )

    # --------------------------------------------------------
    # Consequence
    # --------------------------------------------------------

    consequence_terms = transcript.get(
        "consequence_terms",
        []
    )

    if isinstance(consequence_terms, list):

        consequence = ",".join(
            consequence_terms
        )

    else:

        consequence = str(
            consequence_terms
        )

    # --------------------------------------------------------
    # Return annotation
    # --------------------------------------------------------

    return {

        "HGVS_cDNA":
            hgvs,

        "Status":
            "Annotated",

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

        "Gene_Symbol":
            transcript.get(
                "gene_symbol"
            ),

        "Gene_ID":
            transcript.get(
                "gene_id"
            ),

        "Transcript_ID":
            transcript.get(
                "transcript_id"
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

        "HGVSc":
            hgvsc,

        "HGVSp":
            hgvsp,

        "Protein_HGVS":
            hgvsp,

        "Consequence":
            consequence,

        "Impact":
            transcript.get(
                "impact"
            ),

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

        "Variant_Allele":
            transcript.get(
                "variant_allele"
            )
    }


# ============================================================
# LOAD WEB VEP FALLBACK
# ============================================================

def load_vep_web_fallback(
    vep_output_file: str
) -> pd.DataFrame:
    """
    Load an existing VEP web-exported Excel annotation table.

    Parameters
    ----------
    vep_output_file : str
        Path to VEP Excel export.

    Returns
    -------
    pandas.DataFrame
        VEP web annotation table.
    """

    if not os.path.exists(
        vep_output_file
    ):
        raise FileNotFoundError(
            "VEP fallback file was not found:\n"
            f"{vep_output_file}"
        )

    fallback_df = pd.read_excel(
        vep_output_file
    )

    if fallback_df.empty:
        raise ValueError(
            "The VEP fallback file is empty:\n"
            f"{vep_output_file}"
        )

    return fallback_df


# ============================================================
# NORMALISE HGVS VALUES
# ============================================================

def normalise_hgvs(
    value
) -> str:
    """
    Normalise an HGVS value for matching.

    Parameters
    ----------
    value : object
        HGVS value.

    Returns
    -------
    str
        Normalised HGVS string.
    """

    if pd.isna(value):
        return ""

    return (
        str(value)
        .strip()
        .replace(" ", "")
    )


# ============================================================
# MATCH FALLBACK ANNOTATION
# ============================================================

def find_web_annotation(
    hgvs: str,
    fallback_df: pd.DataFrame
) -> dict | None:
    """
    Find a matching variant in the VEP web-exported table.

    The function attempts to match against common HGVS columns.

    Parameters
    ----------
    hgvs : str
        HGVS cDNA identifier.

    fallback_df : pandas.DataFrame
        VEP web annotation table.

    Returns
    -------
    dict or None
        Matching annotation or None.
    """

    target = normalise_hgvs(
        hgvs
    )

    if not target:
        return None

    # --------------------------------------------------------
    # Candidate HGVS columns
    # --------------------------------------------------------

    candidate_columns = [
        "HGVS_cDNA",
        "HGVSc",
        "HGVSc",
        "Uploaded_variation",
        "Input",
        "Variant"
    ]

    available_columns = [
        column
        for column in candidate_columns
        if column in fallback_df.columns
    ]

    # --------------------------------------------------------
    # Search each candidate column
    # --------------------------------------------------------

    for column in available_columns:

        normalised_values = (
            fallback_df[column]
            .map(normalise_hgvs)
        )

        matches = fallback_df[
            normalised_values == target
        ]

        if not matches.empty:

            return (
                matches.iloc[0]
                .to_dict()
            )

    # --------------------------------------------------------
    # Additional fallback:
    # search entire row text
    # --------------------------------------------------------

    for _, row in fallback_df.iterrows():

        row_values = [
            normalise_hgvs(value)
            for value in row.tolist()
        ]

        if target in row_values:

            return row.to_dict()

    return None


# ============================================================
# MERGE API AND WEB ANNOTATION
# ============================================================

def merge_web_fallback(
    api_annotation: dict,
    web_annotation: dict,
    http_status
) -> dict:
    """
    Merge a VEP web annotation into a failed API annotation.

    Existing API values are retained where present.
    Web values are used to fill missing fields.

    Parameters
    ----------
    api_annotation : dict
        Existing API annotation record.

    web_annotation : dict
        Matching VEP web annotation.

    http_status : int or None
        HTTP status from API request.

    Returns
    -------
    dict
        Recovered annotation.
    """

    recovered = dict(
        api_annotation
    )

    # --------------------------------------------------------
    # Direct field mappings
    # --------------------------------------------------------

    field_mapping = {

        "Consequence":
            [
                "Consequence",
                "Consequence",
                "consequence"
            ],

        "Impact":
            [
                "Impact",
                "IMPACT",
                "impact"
            ],

        "HGVSc":
            [
                "HGVSc",
                "HGVS_cDNA",
                "HGVSc"
            ],

        "HGVSp":
            [
                "HGVSp",
                "Protein_HGVS",
                "HGVSp"
            ],

        "Protein_HGVS":
            [
                "Protein_HGVS",
                "HGVSp"
            ],

        "Gene_Symbol":
            [
                "Gene_Symbol",
                "SYMBOL",
                "Gene"
            ],

        "Gene_ID":
            [
                "Gene_ID",
                "Gene"
            ],

        "Transcript_ID":
            [
                "Transcript_ID",
                "Feature"
            ],

        "Assembly":
            [
                "Assembly",
                "ASSEMBLY"
            ],

        "Chromosome":
            [
                "Chromosome",
                "Location"
            ],

        "Genomic_Position":
            [
                "Genomic_Position",
                "Position"
            ],

        "Genomic_Coordinate":
            [
                "Genomic_Coordinate",
                "Location"
            ],

        "Protein_ID":
            [
                "Protein_ID",
                "Protein"
            ],

        "Protein_Position":
            [
                "Protein_Position"
            ],

        "Amino_Acid_Change":
            [
                "Amino_Acid_Change",
                "Amino_acids"
            ],

        "Codon_Change":
            [
                "Codon_Change",
                "Codons"
            ],

        "Variant_Allele":
            [
                "Variant_Allele",
                "Allele"
            ]
    }

    # --------------------------------------------------------
    # Fill missing API values from web annotation
    # --------------------------------------------------------

    for target_field, source_fields in (
        field_mapping.items()
    ):

        current_value = recovered.get(
            target_field
        )

        current_missing = (
            current_value is None
            or pd.isna(current_value)
            or str(current_value).strip()
            in {"", "nan", "None"}
        )

        if not current_missing:
            continue

        for source_field in source_fields:

            if source_field not in web_annotation:
                continue

            value = web_annotation.get(
                source_field
            )

            if (
                value is not None
                and not pd.isna(value)
                and str(value).strip()
                not in {"", "nan", "None"}
            ):
                recovered[target_field] = value
                break

    # --------------------------------------------------------
    # Final metadata
    # --------------------------------------------------------

    recovered["HGVS_cDNA"] = api_annotation.get(
        "HGVS_cDNA"
    )

    recovered["Status"] = (
        "Annotated"
    )

    recovered["Annotation_Source"] = (
        "VEP_web_fallback"
    )

    recovered["HTTP_Status"] = (
        http_status
    )

    return recovered


# ============================================================
# ENSEMBL VEP ANNOTATION
# ============================================================

def annotate_variants(
    hgvs_df: pd.DataFrame,
    server: str,
    headers: dict,
    timeout: int = 30,
    retries: int = 3,
    delay: float = 0.1,
    transcript: str = None,
    fallback_file: str = None
) -> pd.DataFrame:
    """
    Annotate validated variants using the Ensembl VEP REST API.

    Failed API annotations are optionally recovered from an
    existing VEP web-exported Excel file.

    Parameters
    ----------
    hgvs_df : pandas.DataFrame
        DataFrame containing validated HGVS variants.

    server : str
        Ensembl REST API server URL.

    headers : dict
        HTTP request headers.

    timeout : int
        Request timeout in seconds.

    retries : int
        Maximum number of API attempts.

    delay : float
        Delay between requests.

    transcript : str, optional
        Preferred transcript.

    fallback_file : str, optional
        Existing VEP web-exported Excel file.

    Returns
    -------
    pandas.DataFrame
        Annotation table.
    """

    if "HGVS_cDNA" not in hgvs_df.columns:

        raise ValueError(
            "Input HGVS table must contain "
            "'HGVS_cDNA'."
        )

    # --------------------------------------------------------
    # Load fallback annotation if available
    # --------------------------------------------------------

    fallback_df = None

    if fallback_file:

        if os.path.exists(
            fallback_file
        ):

            print(
                "      VEP web fallback available:"
            )

            print(
                f"      {fallback_file}"
            )

            fallback_df = (
                load_vep_web_fallback(
                    fallback_file
                )
            )

        else:

            print(
                "      VEP web fallback file not found."
            )

            print(
                f"      {fallback_file}"
            )

    annotation_rows = []

    # ========================================================
    # Process variants
    # ========================================================

    for _, row in hgvs_df.iterrows():

        hgvs_entries = str(
            row["HGVS_cDNA"]
        ).split(";")

        for hgvs in hgvs_entries:

            hgvs = hgvs.strip()

            if not hgvs:
                continue

            print(
                f"      Annotating: {hgvs}"
            )

            # ------------------------------------------------
            # Build API URL
            # ------------------------------------------------

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

            # ------------------------------------------------
            # API request with retry
            # ------------------------------------------------

            for attempt in range(
                1,
                retries + 1
            ):

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

                    if response.status_code >= 500:

                        print(
                            f"      API server error "
                            f"{response.status_code}; "
                            f"retry {attempt}/{retries}"
                        )

                    else:

                        print(
                            f"      API returned "
                            f"HTTP {response.status_code}"
                        )

                except requests.RequestException as exc:

                    print(
                        f"      API request error: "
                        f"{exc}"
                    )

                if attempt < retries:

                    time.sleep(
                        min(
                            2 ** attempt,
                            10
                        )
                    )

            # ------------------------------------------------
            # API succeeded
            # ------------------------------------------------

            if (
                response is not None
                and response.status_code == 200
            ):

                try:

                    results = response.json()

                except ValueError:

                    results = []

                    print(
                        "      API returned invalid JSON."
                    )

                if results:

                    annotation = (
                        parse_vep_response(
                            hgvs,
                            results[0]
                        )
                    )

                    annotation["HTTP_Status"] = (
                        response.status_code
                    )

                    annotation["Annotation_Source"] = (
                        "VEP_API"
                    )

                    annotation_rows.append(
                        annotation
                    )

                    time.sleep(
                        delay
                    )

                    continue

                print(
                    "      API returned an empty result."
                )

            # ------------------------------------------------
            # API failed
            # ------------------------------------------------

            status_code = (
                response.status_code
                if response is not None
                else None
            )

            print(
                f"      API annotation failed "
                f"for {hgvs}."
            )

            # ------------------------------------------------
            # Try VEP web fallback
            # ------------------------------------------------

            if fallback_df is not None:

                web_annotation = (
                    find_web_annotation(
                        hgvs,
                        fallback_df
                    )
                )

                if web_annotation is not None:

                    print(
                        "      Web VEP fallback matched."
                    )

                    recovered = merge_web_fallback(
                        {
                            "HGVS_cDNA": hgvs,
                            "Status": "Failed",
                            "HTTP_Status": status_code
                        },
                        web_annotation,
                        status_code
                    )

                    annotation_rows.append(
                        recovered
                    )

                    time.sleep(
                        delay
                    )

                    continue

                print(
                    "      No matching web VEP "
                    "annotation found."
                )

            # ------------------------------------------------
            # No annotation available
            # ------------------------------------------------

            annotation_rows.append(
                {
                    "HGVS_cDNA": hgvs,
                    "Status": "Failed",
                    "HTTP_Status": status_code,
                    "Annotation_Source": (
                        "VEP_API_failed"
                    )
                }
            )

            time.sleep(
                delay
            )

    # ========================================================
    # Construct DataFrame
    # ========================================================

    annotation_df = pd.DataFrame(
        annotation_rows
    )

    # --------------------------------------------------------
    # Remove duplicate HGVS records
    # --------------------------------------------------------

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


# ============================================================
# SAVE ANNOTATION TABLE
# ============================================================

def save_annotation_table(
    annotation_df: pd.DataFrame,
    output_folder: str
):
    """
    Save VEP annotation table.

    Parameters
    ----------
    annotation_df : pandas.DataFrame
        Annotation results.

    output_folder : str
        Output directory.

    Returns
    -------
    str
        Path to saved file.
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