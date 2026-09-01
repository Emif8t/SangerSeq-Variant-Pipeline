"""
reference.py

Reference sequence handling functions for the SangerSeq
Variant Pipeline.

This module:

1. Downloads a reference transcript from NCBI RefSeq.
2. Extracts reference amplicons.
3. Verifies forward and reverse PCR primer binding sites.
4. Determines the reference amplicon sequence and coordinates.

The module is gene-agnostic. Gene-specific information such as
the RefSeq accession and PCR primers must be supplied through
the pipeline configuration.
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq


# ============================================================
# 1. DOWNLOAD REFSEQ TRANSCRIPT
# ============================================================

def download_reference(
    refseq_id: str,
    email: str
) -> str:
    """
    Download a nucleotide reference sequence from NCBI RefSeq.

    Parameters
    ----------
    refseq_id : str
        RefSeq accession number, for example:
        NM_000050.4

    email : str
        Email address used by NCBI Entrez.

    Returns
    -------
    str
        Reference nucleotide sequence in uppercase.

    Raises
    ------
    ValueError
        If the accession or email is missing.

    RuntimeError
        If the reference sequence cannot be retrieved.
    """

    if not refseq_id:
        raise ValueError(
            "REFSEQ_ID must be provided."
        )

    if not email:
        raise ValueError(
            "NCBI_EMAIL must be provided."
        )

    # --------------------------------------------------------
    # Configure NCBI Entrez
    # --------------------------------------------------------

    Entrez.email = email

    try:

        with Entrez.efetch(
            db="nucleotide",
            id=refseq_id,
            rettype="fasta",
            retmode="text"
        ) as handle:

            record = SeqIO.read(
                handle,
                "fasta"
            )

    except Exception as exc:

        raise RuntimeError(
            f"Failed to retrieve reference sequence "
            f"for {refseq_id} from NCBI."
        ) from exc

    sequence = str(
        record.seq
    ).upper().strip()

    if not sequence:

        raise RuntimeError(
            f"NCBI returned an empty sequence for "
            f"{refseq_id}."
        )

    # --------------------------------------------------------
    # Basic sequence validation
    # --------------------------------------------------------

    valid_bases = set(
        "ACGTN"
    )

    invalid_bases = set(sequence) - valid_bases

    if invalid_bases:

        raise ValueError(
            f"Reference sequence {refseq_id} contains "
            f"unexpected nucleotide characters: "
            f"{sorted(invalid_bases)}"
        )

    return sequence


# ============================================================
# 2. EXTRACT PCR AMPLICON
# ============================================================

def extract_amplicon(
    reference_sequence: str,
    start: int,
    end: int
) -> str:
    """
    Extract a reference amplicon using 1-based inclusive
    coordinates.

    Parameters
    ----------
    reference_sequence : str
        Full reference nucleotide sequence.

    start : int
        Amplicon start coordinate, 1-based and inclusive.

    end : int
        Amplicon end coordinate, 1-based and inclusive.

    Returns
    -------
    str
        Amplicon sequence.

    Raises
    ------
    ValueError
        If coordinates are invalid.
    """

    if not reference_sequence:

        raise ValueError(
            "Reference sequence is empty."
        )

    if start < 1:

        raise ValueError(
            "Amplicon start must be >= 1."
        )

    if end < start:

        raise ValueError(
            "Amplicon end must be greater than "
            "or equal to amplicon start."
        )

    if end > len(reference_sequence):

        raise ValueError(
            f"Amplicon end ({end}) exceeds reference "
            f"sequence length ({len(reference_sequence)})."
        )

    # Python slicing is zero-based and end-exclusive.
    # Convert from 1-based inclusive coordinates.
    amplicon = reference_sequence[
        start - 1:end
    ]

    if not amplicon:

        raise ValueError(
            "Extracted amplicon is empty."
        )

    return amplicon.upper()


# ============================================================
# 3. NORMALISE PRIMER
# ============================================================

def _normalise_primer(
    primer: str,
    primer_name: str
) -> str:
    """
    Normalise a PCR primer sequence.

    Whitespace is removed and the sequence is converted
    to uppercase.

    Parameters
    ----------
    primer : str
        Primer nucleotide sequence.

    primer_name : str
        Human-readable primer name for error messages.

    Returns
    -------
    str
        Normalised primer sequence.
    """

    if primer is None:

        raise ValueError(
            f"{primer_name} was not provided."
        )

    primer = "".join(
        str(primer).split()
    ).upper()

    if not primer:

        raise ValueError(
            f"{primer_name} is empty."
        )

    valid_bases = set(
        "ACGTN"
    )

    invalid_bases = (
        set(primer)
        - valid_bases
    )

    if invalid_bases:

        raise ValueError(
            f"{primer_name} contains invalid nucleotide "
            f"characters: {sorted(invalid_bases)}"
        )

    return primer


# ============================================================
# 4. FIND PRIMER
# ============================================================

def _find_primer(
    reference_sequence: str,
    primer: str,
    primer_name: str
) -> int:
    """
    Find a primer sequence in a reference sequence.

    Coordinates returned internally are zero-based.

    Parameters
    ----------
    reference_sequence : str

    primer : str

    primer_name : str

    Returns
    -------
    int
        Zero-based start coordinate.

    Raises
    ------
    ValueError
        If the primer cannot be found.
    """

    position = reference_sequence.find(
        primer
    )

    if position == -1:

        raise ValueError(
            f"{primer_name} was not found in the "
            f"reference sequence."
        )

    return position


# ============================================================
# 5. VERIFY PCR PRIMERS
# ============================================================

def verify_primers(
    reference_sequence: str,
    forward_primer: str,
    reverse_primer: str
) -> dict:
    """
    Verify forward and reverse PCR primer binding sites.

    The forward primer is searched directly in the reference
    sequence.

    The reverse primer is reverse-complemented before searching,
    because the primer sequence is normally supplied in the
    5' -> 3' orientation of the oligonucleotide.

    Coordinates returned by this function are 1-based and
    inclusive.

    Parameters
    ----------
    reference_sequence : str
        Full reference transcript sequence.

    forward_primer : str
        Forward PCR primer sequence.

    reverse_primer : str
        Reverse PCR primer sequence in its normal 5' -> 3'
        oligonucleotide orientation.

    Returns
    -------
    dict
        Dictionary containing:

        forward_primer
        reverse_primer
        reverse_primer_reverse_complement
        forward_start
        forward_end
        reverse_start
        reverse_end
        amplicon_start
        amplicon_end
        amplicon_length
        amplicon_sequence

    Raises
    ------
    ValueError
        If either primer cannot be found or if the primer
        orientation is inconsistent.
    """

    if not reference_sequence:

        raise ValueError(
            "Reference sequence is empty."
        )

    reference_sequence = (
        "".join(
            str(reference_sequence).split()
        )
        .upper()
    )

    forward_primer = _normalise_primer(
        forward_primer,
        "Forward primer"
    )

    reverse_primer = _normalise_primer(
        reverse_primer,
        "Reverse primer"
    )

    # --------------------------------------------------------
    # Reverse-complement the reverse primer
    # --------------------------------------------------------

    reverse_primer_rc = str(
        Seq(
            reverse_primer
        ).reverse_complement()
    )

    # --------------------------------------------------------
    # Locate forward primer
    # --------------------------------------------------------

    forward_start_zero = _find_primer(
        reference_sequence,
        forward_primer,
        "Forward primer"
    )

    forward_end_zero = (
        forward_start_zero
        + len(forward_primer)
        - 1
    )

    # --------------------------------------------------------
    # Locate reverse primer
    # --------------------------------------------------------

    reverse_start_zero = _find_primer(
        reference_sequence,
        reverse_primer_rc,
        "Reverse primer reverse-complement"
    )

    reverse_end_zero = (
        reverse_start_zero
        + len(reverse_primer_rc)
        - 1
    )

    # --------------------------------------------------------
    # Ensure correct primer orientation
    # --------------------------------------------------------

    if reverse_start_zero <= forward_start_zero:

        raise ValueError(
            "Invalid primer orientation: the reverse primer "
            "binding site occurs before the forward primer "
            "binding site."
        )

    # --------------------------------------------------------
    # Convert to 1-based inclusive coordinates
    # --------------------------------------------------------

    forward_start = (
        forward_start_zero + 1
    )

    forward_end = (
        forward_end_zero + 1
    )

    reverse_start = (
        reverse_start_zero + 1
    )

    reverse_end = (
        reverse_end_zero + 1
    )

    amplicon_start = forward_start

    amplicon_end = reverse_end

    # --------------------------------------------------------
    # Extract amplicon
    # --------------------------------------------------------

    amplicon_sequence = extract_amplicon(
        reference_sequence,
        amplicon_start,
        amplicon_end
    )

    amplicon_length = len(
        amplicon_sequence
    )

    # --------------------------------------------------------
    # Final validation
    # --------------------------------------------------------

    if amplicon_length <= 0:

        raise ValueError(
            "PCR amplicon has zero length."
        )

    return {

        # ----------------------------------------------------
        # Primer information
        # ----------------------------------------------------

        "forward_primer":
            forward_primer,

        "reverse_primer":
            reverse_primer,

        "reverse_primer_reverse_complement":
            reverse_primer_rc,

        # ----------------------------------------------------
        # Forward primer coordinates
        # ----------------------------------------------------

        "forward_start":
            forward_start,

        "forward_end":
            forward_end,

        # ----------------------------------------------------
        # Reverse primer coordinates
        # ----------------------------------------------------

        "reverse_start":
            reverse_start,

        "reverse_end":
            reverse_end,

        # ----------------------------------------------------
        # Amplicon coordinates
        # ----------------------------------------------------

        "amplicon_start":
            amplicon_start,

        "amplicon_end":
            amplicon_end,

        "amplicon_length":
            amplicon_length,

        # ----------------------------------------------------
        # Reference amplicon
        # ----------------------------------------------------

        "amplicon_sequence":
            amplicon_sequence

    }