"""
alignment.py

Functions for performing local sequence alignment and generating
nucleotide-level alignment tables for the SangerSeq Variant Pipeline.

The module is gene-agnostic. Reference sequence, amplicon sequence,
and amplicon coordinates are supplied by the calling pipeline.

Workflow
--------
1. Perform Smith-Waterman local alignment between each sequencing
   read and the reference PCR amplicon.
2. Walk through the alignment column by column.
3. Map aligned reference bases back to reference transcript/cDNA
   coordinates.
4. Preserve read position and Phred quality information.
5. Classify alignment events and base-call confidence.
"""


import pandas as pd

from Bio import pairwise2


# ============================================================
# 1. IUPAC AMBIGUITY CODES
# ============================================================

IUPAC_CODES = {

    "R": ["A", "G"],  # A/G
    "Y": ["C", "T"],  # C/T
    "S": ["G", "C"],  # G/C
    "W": ["A", "T"],  # A/T
    "K": ["G", "T"],  # G/T
    "M": ["A", "C"],  # A/C

    # Additional standard IUPAC ambiguity codes
    "B": ["C", "G", "T"],  # not A
    "D": ["A", "G", "T"],  # not C
    "H": ["A", "C", "T"],  # not G
    "V": ["A", "C", "G"],  # not T
    "N": ["A", "C", "G", "T"],  # any base

}


# ============================================================
# 2. BASIC SEQUENCE VALIDATION
# ============================================================

def _clean_sequence(
    sequence: str,
    sequence_name: str
) -> str:
    """
    Clean and validate a nucleotide sequence.

    Parameters
    ----------
    sequence : str
        Input nucleotide sequence.

    sequence_name : str
        Name used in error messages.

    Returns
    -------
    str
        Uppercase sequence without whitespace.

    Raises
    ------
    ValueError
        If the sequence is empty or contains invalid characters.
    """

    if sequence is None:

        raise ValueError(
            f"{sequence_name} was not provided."
        )

    sequence = "".join(
        str(sequence).split()
    ).upper()

    if not sequence:

        raise ValueError(
            f"{sequence_name} is empty."
        )

    valid_bases = set(
        "ACGTNRYSWKMBDHV"
    )

    invalid_bases = (
        set(sequence)
        - valid_bases
    )

    if invalid_bases:

        raise ValueError(
            f"{sequence_name} contains invalid nucleotide "
            f"characters: {sorted(invalid_bases)}"
        )

    return sequence


# ============================================================
# 3. LOCAL SEQUENCE ALIGNMENT
# ============================================================

def perform_local_alignment(
    processed_reads: list,
    reference_amplicon: str
) -> list:
    """
    Perform Smith-Waterman local alignment between sequencing
    reads and the reference PCR amplicon.

    Parameters
    ----------
    processed_reads : list
        List of processed sequencing-read dictionaries.

        Each read is expected to contain at least:

            filename
            sequence
            quality

    reference_amplicon : str
        Reference PCR amplicon sequence.

    Returns
    -------
    list
        Successfully aligned sequencing reads.

    Notes
    -----
    The alignment is performed against the reference amplicon,
    not directly against the complete transcript.

    Reference transcript coordinates are reconstructed later
    using the amplicon_start coordinate supplied to walk_alignment().
    """

    reference_amplicon = _clean_sequence(
        reference_amplicon,
        "Reference amplicon"
    )

    alignments = []

    for read in processed_reads:

        # ----------------------------------------------------
        # Validate read structure
        # ----------------------------------------------------

        if not isinstance(
            read,
            dict
        ):

            continue

        if "sequence" not in read:

            continue

        if "filename" not in read:

            continue

        read_sequence = _clean_sequence(
            read["sequence"],
            f"Sequencing read {read.get('filename', '')}"
        )

        quality = read.get(
            "quality",
            []
        )

        # ----------------------------------------------------
        # Smith-Waterman local alignment
        # ----------------------------------------------------

        results = pairwise2.align.localms(

            reference_amplicon,

            read_sequence,

            2,       # match score
            -2,      # mismatch penalty
            -10,     # gap opening penalty
            -1,      # gap extension penalty

            one_alignment_only=True

        )

        if not results:

            continue

        alignment = results[0]

        aligned_ref = alignment.seqA
        aligned_read = alignment.seqB

        # ----------------------------------------------------
        # Store alignment
        # ----------------------------------------------------

        alignments.append({

            "filename":
                read["filename"],

            "aligned_ref":
                aligned_ref,

            "aligned_read":
                aligned_read,

            "quality":
                quality,

            "score":
                alignment.score,

        })

    return alignments


# ============================================================
# 4. FIND FIRST ALIGNED BASE
# ============================================================

def _find_first_aligned_column(
    aligned_ref: str,
    aligned_read: str
):
    """
    Find the first alignment column containing both a reference
    base and a read base.

    Parameters
    ----------
    aligned_ref : str

    aligned_read : str

    Returns
    -------
    int or None
        Alignment-column index.
    """

    for column in range(
        len(aligned_ref)
    ):

        if (
            aligned_ref[column] != "-"
            and
            aligned_read[column] != "-"
        ):

            return column

    return None


# ============================================================
# 5. WALK ALIGNMENT
# ============================================================

def walk_alignment(
    alignments: list,
    reference_sequence: str,
    amplicon_start: int,
    min_phred: int = 20
) -> pd.DataFrame:
    """
    Traverse local alignments and generate a nucleotide-level
    alignment table.

    Parameters
    ----------
    alignments : list
        Output from perform_local_alignment().

    reference_sequence : str
        Complete reference transcript sequence.

    amplicon_start : int
        1-based inclusive coordinate of the first base of the
        PCR amplicon within the reference transcript.

    min_phred : int, default=20
        Minimum Phred score used to distinguish lower-confidence
        from medium/high-confidence base calls.

    Returns
    -------
    pandas.DataFrame
        Nucleotide-level alignment table.

    Coordinate system
    -----------------
    Reference transcript coordinates are 1-based and inclusive.

    For example, if:

        amplicon_start = 100

    then:

        first amplicon base  -> cDNA_Position 100
        second amplicon base -> cDNA_Position 101
        third amplicon base  -> cDNA_Position 102

    Insertions relative to the reference do not receive a
    reference cDNA coordinate because they occur between
    reference positions.
    """

    # ========================================================
    # Validate reference
    # ========================================================

    reference_sequence = _clean_sequence(
        reference_sequence,
        "Reference sequence"
    )

    if amplicon_start < 1:

        raise ValueError(
            "amplicon_start must be >= 1."
        )

    if amplicon_start > len(reference_sequence):

        raise ValueError(
            f"amplicon_start ({amplicon_start}) exceeds "
            f"reference sequence length ({len(reference_sequence)})."
        )

    if min_phred < 0:

        raise ValueError(
            "min_phred cannot be negative."
        )

    # ========================================================
    # Alignment rows
    # ========================================================

    alignment_rows = []

    # ========================================================
    # Process each sequencing read
    # ========================================================

    for sample in alignments:

        filename = sample.get(
            "filename"
        )

        aligned_ref = sample.get(
            "aligned_ref"
        )

        aligned_read = sample.get(
            "aligned_read"
        )

        quality = sample.get(
            "quality",
            []
        )

        score = sample.get(
            "score"
        )

        if aligned_ref is None or aligned_read is None:

            continue

        if len(aligned_ref) != len(aligned_read):

            raise ValueError(
                f"Alignment length mismatch for sample "
                f"{filename}."
            )

        # ----------------------------------------------------
        # Ensure quality is usable
        # ----------------------------------------------------

        if quality is None:

            quality = []

        try:

            quality = list(quality)

        except TypeError:

            quality = []

        # ----------------------------------------------------
        # Find beginning of local alignment
        # ----------------------------------------------------

        first_aligned = _find_first_aligned_column(
            aligned_ref,
            aligned_read
        )

        if first_aligned is None:

            continue

        # ----------------------------------------------------
        # Determine how many reference bases occur before the
        # first aligned column.
        #
        # This allows the alignment to begin inside the PCR
        # amplicon while still mapping correctly to transcript
        # coordinates.
        # ----------------------------------------------------

        reference_bases_before_alignment = sum(

            base != "-"

            for base in aligned_ref[
                :first_aligned
            ]

        )

        # ----------------------------------------------------
        # Determine how many read bases occur before the local
        # alignment.
        # ----------------------------------------------------

        read_bases_before_alignment = sum(

            base != "-"

            for base in aligned_read[
                :first_aligned
            ]

        )

        # ----------------------------------------------------
        # Convert counts to zero-based indices.
        #
        # The first aligned base is incremented before it is
        # recorded below.
        # ----------------------------------------------------

        reference_index = (
            reference_bases_before_alignment - 1
        )

        read_index = (
            read_bases_before_alignment - 1
        )

        # ====================================================
        # Walk alignment column by column
        # ====================================================

        for column in range(

            first_aligned,

            len(aligned_ref)

        ):

            ref_base = aligned_ref[column]

            read_base = aligned_read[column]

            # ------------------------------------------------
            # Reference coordinate
            # ------------------------------------------------

            if ref_base != "-":

                reference_index += 1

                amplicon_position = (
                    reference_index + 1
                )

                # ------------------------------------------------
                # Map amplicon coordinate to transcript coordinate.
                #
                # amplicon_position is 1-based.
                # amplicon_start is 1-based.
                # ------------------------------------------------

                cdna_position = (

                    amplicon_start
                    +
                    amplicon_position
                    - 1

                )

                # ------------------------------------------------
                # Retrieve actual reference base from complete
                # transcript.
                # ------------------------------------------------

                if (
                    cdna_position >= 1
                    and
                    cdna_position <= len(reference_sequence)
                ):

                    reference_base = (
                        reference_sequence[
                            cdna_position - 1
                        ]
                    )

                else:

                    raise RuntimeError(

                        f"Calculated reference coordinate "
                        f"{cdna_position} is outside the "
                        f"reference sequence for sample "
                        f"{filename}."

                    )

            else:

                amplicon_position = None

                cdna_position = None

                reference_base = "-"

            # ------------------------------------------------
            # Read coordinate and quality
            # ------------------------------------------------

            if read_base != "-":

                read_index += 1

                read_position = (
                    read_index + 1
                )

                if (
                    read_index >= 0
                    and
                    read_index < len(quality)
                ):

                    base_quality = quality[
                        read_index
                    ]

                else:

                    base_quality = None

            else:

                read_position = None

                base_quality = None

            # ------------------------------------------------
            # Normalise observed base
            # ------------------------------------------------

            if read_base != "-":

                observed_base = (
                    str(read_base).upper()
                )

            else:

                observed_base = "-"

            # ------------------------------------------------
            # Alignment event
            # ------------------------------------------------

            if (
                ref_base == "-"
                and
                read_base != "-"
            ):

                event = "Insertion"

            elif (
                ref_base != "-"
                and
                read_base == "-"
            ):

                event = "Deletion"

            elif (
                ref_base == read_base
            ):

                event = "Match"

            else:

                event = "Mismatch"

            # ------------------------------------------------
            # Determine confidence
            # ------------------------------------------------

            if read_base == "-":

                confidence = "Deletion"

            elif base_quality is None:

                confidence = "Missing"

            elif observed_base in IUPAC_CODES:

                if base_quality >= min_phred:

                    confidence = "Ambiguous_Het"

                else:

                    confidence = (
                        "LowConfidence_Heterozygous"
                    )

            else:

                try:

                    numeric_quality = float(
                        base_quality
                    )

                except (
                    TypeError,
                    ValueError
                ):

                    numeric_quality = None

                if numeric_quality is None:

                    confidence = "Missing"

                elif numeric_quality >= 30:

                    confidence = "High"

                elif numeric_quality >= min_phred:

                    confidence = "Medium"

                else:

                    confidence = "Low"

            # ------------------------------------------------
            # Determine variant base
            # ------------------------------------------------

            if (

                ref_base != "-"

                and

                read_base != "-"

                and

                observed_base != str(
                    reference_base
                ).upper()

            ):

                variant_base = observed_base

            else:

                variant_base = None

            # ------------------------------------------------
            # Store alignment record
            # ------------------------------------------------

            alignment_rows.append({

                "Sample":
                    filename,

                "Alignment_Score":
                    score,

                "Alignment_Column":
                    column + 1,

                "Reference_Index":
                    (
                        reference_index
                        if ref_base != "-"
                        else None
                    ),

                "Read_Index":
                    (
                        read_index
                        if read_base != "-"
                        else None
                    ),

                "Amplicon_Position":
                    amplicon_position,

                "Read_Position":
                    read_position,

                "cDNA_Position":
                    cdna_position,

                "REF":
                    reference_base,

                "Observed_Base":
                    observed_base,

                "Quality":
                    base_quality,

                "Confidence":
                    confidence,

                "Variant_Base":
                    variant_base,

                "Alignment_Event":
                    event

            })

    # ========================================================
    # Create dataframe
    # ========================================================

    alignment_df = pd.DataFrame(
        alignment_rows
    )

    # ========================================================
    # Ensure predictable column order
    # ========================================================

    expected_columns = [

        "Sample",
        "Alignment_Score",
        "Alignment_Column",
        "Reference_Index",
        "Read_Index",
        "Amplicon_Position",
        "Read_Position",
        "cDNA_Position",
        "REF",
        "Observed_Base",
        "Quality",
        "Confidence",
        "Variant_Base",
        "Alignment_Event",

    ]

    if alignment_df.empty:

        return pd.DataFrame(
            columns=expected_columns
        )

    alignment_df = alignment_df[
        expected_columns
    ]

    # ========================================================
    # Numeric conversion
    # ========================================================

    for column in [

        "Alignment_Score",
        "Alignment_Column",
        "Reference_Index",
        "Read_Index",
        "Amplicon_Position",
        "Read_Position",
        "cDNA_Position",
        "Quality",

    ]:

        alignment_df[column] = pd.to_numeric(
            alignment_df[column],
            errors="coerce"
        )

    return alignment_df