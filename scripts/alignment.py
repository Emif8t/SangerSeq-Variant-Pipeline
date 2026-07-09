"""
alignment.py

Functions for performing local sequence alignment
and generating nucleotide-level alignment tables.
"""

import pandas as pd
from Bio import pairwise2


# ======================================================
# IUPAC AMBIGUITY CODES
# ======================================================

IUPAC_CODES = {

    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"]

}


# ======================================================
# LOCAL ALIGNMENT
# ======================================================

def perform_local_alignment(
    processed_reads: list,
    reference_amplicon: str
) -> list:
    """
    Perform Smith-Waterman local alignment.

    Parameters
    ----------
    processed_reads : list

    reference_amplicon : str

    Returns
    -------
    list
        Successfully aligned reads.
    """

    alignments = []

    for read in processed_reads:

        results = pairwise2.align.localms(

            reference_amplicon,

            read["sequence"],

            2,

            -2,

            -10,

            -1,

            one_alignment_only=True

        )

        if not results:

            continue

        aln = results[0]

        alignments.append({

            "filename": read["filename"],

            "aligned_ref": aln.seqA,

            "aligned_read": aln.seqB,

            "quality": read["quality"],

            "score": aln.score

        })

    return alignments


# ======================================================
# WALK ALIGNMENT
# ======================================================

def walk_alignment(
    alignments: list,
    reference_sequence: str,
    amplicon_start: int,
    min_phred: int = 20
) -> pd.DataFrame:
    """
    Traverse aligned sequences and generate a
    nucleotide-level alignment table.

    Parameters
    ----------
    alignments : list

    reference_sequence : str

    amplicon_start : int

    min_phred : int

    Returns
    -------
    pandas.DataFrame
    """

    alignment_rows = []

    for sample in alignments:

        filename = sample["filename"]

        aligned_ref = sample["aligned_ref"]

        aligned_read = sample["aligned_read"]

        quality = sample["quality"]

        score = sample["score"]

        first_aligned = None

        for i in range(len(aligned_ref)):

            if aligned_ref[i] != "-" and aligned_read[i] != "-":

                first_aligned = i

                break

        if first_aligned is None:

            continue

        ref_index = (

            sum(

                b != "-"

                for b in aligned_ref[:first_aligned]

            ) - 1

        )

        read_index = (

            sum(

                b != "-"

                for b in aligned_read[:first_aligned]

            ) - 1

        )

        for column in range(

            first_aligned,

            len(aligned_ref)

        ):

            ref_base = aligned_ref[column]

            read_base = aligned_read[column]

            if ref_base != "-":

                ref_index += 1

                amplicon_position = ref_index + 1

                cdna_position = (

                    amplicon_start +

                    ref_index

                )

                reference_base = reference_sequence[
                    cdna_position - 1
                ]

            else:

                amplicon_position = None

                cdna_position = None

                reference_base = "-"

            if read_base != "-":

                read_index += 1

                read_position = read_index + 1

                base_quality = (

                    quality[read_index]

                    if read_index < len(quality)

                    else None

                )

            else:

                read_position = None

                base_quality = None

            if ref_base == "-" and read_base != "-":

                event = "Insertion"

            elif ref_base != "-" and read_base == "-":

                event = "Deletion"

            elif ref_base == read_base:

                event = "Match"

            else:

                event = "Mismatch"

            if read_base == "-":

                confidence = "Deletion"

            elif base_quality is None:

                confidence = "Missing"

            elif read_base in IUPAC_CODES:

                confidence = "Ambiguous_Het"

            elif base_quality >= 30:

                confidence = "High"

            elif base_quality >= min_phred:

                confidence = "Medium"

            else:

                confidence = "Low"

            alignment_rows.append({

                "Sample": filename,

                "Alignment_Score": score,

                "Alignment_Column": column + 1,

                "Reference_Index":

                    ref_index

                    if ref_base != "-"

                    else None,

                "Read_Index":

                    read_index

                    if read_base != "-"

                    else None,

                "Amplicon_Position":

                    amplicon_position,

                "Read_Position":

                    read_position,

                "cDNA_Position":

                    cdna_position,

                "REF":

                    reference_base,

                "Observed_Base":

                    read_base,

                "Quality":

                    base_quality,

                "Confidence":

                    confidence,

                "Variant_Base":

                    read_base

                    if (

                        ref_base != read_base

                        and

                        ref_base != "-"

                        and

                        read_base != "-"

                    )

                    else None,

                "Alignment_Event":

                    event

            })

    return pd.DataFrame(alignment_rows)
