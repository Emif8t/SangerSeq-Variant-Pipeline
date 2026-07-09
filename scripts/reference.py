"""
reference.py

Functions for downloading reference transcripts,
extracting PCR amplicons and verifying primer
binding sites.
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq


# ======================================================
# DOWNLOAD REFSEQ TRANSCRIPT
# ======================================================

def download_reference(
    refseq_id: str,
    email: str
) -> str:
    """
    Download a RefSeq transcript from NCBI.

    Parameters
    ----------
    refseq_id : str
        RefSeq accession number.

    email : str
        Email address required by NCBI Entrez.

    Returns
    -------
    str
        Reference nucleotide sequence.
    """

    Entrez.email = email

    with Entrez.efetch(
        db="nucleotide",
        id=refseq_id,
        rettype="fasta",
        retmode="text"
    ) as handle:

        record = SeqIO.read(handle, "fasta")

    return str(record.seq).upper()


# ======================================================
# EXTRACT PCR AMPLICON
# ======================================================

def extract_amplicon(
    reference_sequence: str,
    start: int,
    end: int
) -> str:
    """
    Extract the PCR amplicon from the
    reference transcript.

    Parameters
    ----------
    reference_sequence : str

    start : int

    end : int

    Returns
    -------
    str
        PCR amplicon sequence.
    """

    return reference_sequence[start - 1:end]


# ======================================================
# VERIFY PCR PRIMERS
# ======================================================

def verify_primers(
    reference_sequence: str,
    forward_primer: str,
    reverse_primer: str
) -> dict:
    """
    Verify primer binding positions on the
    reference transcript.

    Parameters
    ----------
    reference_sequence : str

    forward_primer : str

    reverse_primer : str

    Returns
    -------
    dict
        Primer coordinates and amplicon sequence.
    """

    forward_primer = forward_primer.upper()

    reverse_primer = reverse_primer.upper()

    reverse_rc = str(
        Seq(reverse_primer).reverse_complement()
    )

    forward_start = reference_sequence.find(
        forward_primer
    )

    reverse_start = reference_sequence.find(
        reverse_rc
    )

    if forward_start == -1:

        raise ValueError(
            "Forward primer not found."
        )

    if reverse_start == -1:

        raise ValueError(
            "Reverse primer not found."
        )

    forward_end = (
        forward_start +
        len(forward_primer)
    )

    reverse_end = (
        reverse_start +
        len(reverse_rc)
    )

    amplicon_start = forward_start + 1

    amplicon_end = reverse_end

    amplicon_sequence = reference_sequence[
        amplicon_start - 1:
        amplicon_end
    ]

    return {

        "forward_start": forward_start + 1,

        "forward_end": forward_end,

        "reverse_start": reverse_start + 1,

        "reverse_end": reverse_end,

        "amplicon_start": amplicon_start,

        "amplicon_end": amplicon_end,

        "amplicon_sequence": amplicon_sequence

    }
