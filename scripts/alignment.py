"""
Alignment and reference-walking functions for the SangerSeq Variant Pipeline.

This module performs local alignment of Sanger reads to the target amplicon
using Bio.Align.PairwiseAligner and walks the resulting coordinate blocks
directly to generate nucleotide-level alignment records.

Important design principles
----------------------------
1. PairwiseAligner is used instead of the deprecated pairwise2 API.
2. Alignment coordinates are taken directly from PairwiseAligner.
3. Original per-base PHRED scores are preserved throughout the pipeline.
4. PHRED scores are mapped using the original Read_Index.
5. Insertions receive the PHRED score of the inserted read base.
6. Deletions have no read base and therefore no PHRED score.
7. The public API is compatible with main.py:
       perform_local_alignment(processed_reads, reference)
       walk_alignment(alignments, reference, transcript_start, min_phred)
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple

import pandas as pd
from Bio.Align import PairwiseAligner


# ---------------------------------------------------------------------------
# IUPAC ambiguity codes
# ---------------------------------------------------------------------------

IUPAC_CODES = {
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


# ---------------------------------------------------------------------------
# Basic sequence helpers
# ---------------------------------------------------------------------------

def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    if sequence is None:
        return ""

    complement = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan",
    )

    return sequence.translate(complement)[::-1]


def normalise_base(base: Any) -> str:
    """
    Convert a base or sequence to uppercase string representation.
    """
    if base is None:
        return ""

    return str(base).strip().upper()


def is_iupac_ambiguous(base: str) -> bool:
    """
    Return True if base is an IUPAC ambiguity code.
    """
    return normalise_base(base) in IUPAC_CODES


def is_standard_base(base: str) -> bool:
    """
    Return True for A, C, G or T.
    """
    return normalise_base(base) in {
        "A",
        "C",
        "G",
        "T",
    }


# ---------------------------------------------------------------------------
# Alignment configuration
# ---------------------------------------------------------------------------

def create_aligner() -> PairwiseAligner:
    """
    Create the local PairwiseAligner.

    Scoring corresponds to the previous pairwise2 localms configuration:

        match       = +2
        mismatch    = -2
        gap open    = -10
        gap extend  = -1
    """
    aligner = PairwiseAligner()

    aligner.mode = "local"

    aligner.match_score = 2
    aligner.mismatch_score = -2

    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1

    return aligner


# ---------------------------------------------------------------------------
# Read metadata helpers
# ---------------------------------------------------------------------------

def _extract_read_field(
    read_data: Dict[str, Any],
    *names: str,
) -> Any:
    """
    Return the first available value from a read dictionary.

    This permits compatibility with either upper- or lower-case field names.
    """
    for name in names:
        if name in read_data:
            return read_data[name]

    return None


# ---------------------------------------------------------------------------
# PairwiseAligner coordinate helpers
# ---------------------------------------------------------------------------

def _alignment_blocks(
    alignment: Any,
) -> List[Tuple[Tuple[int, int], Tuple[int, int]]]:
    """
    Extract PairwiseAligner aligned coordinate blocks.

    Returns
    -------
    list
        Each item is:

            ((reference_start, reference_end),
             (read_start, read_end))

        Coordinates are zero-based and end-exclusive.
    """
    blocks: List[Tuple[Tuple[int, int], Tuple[int, int]]] = []

    aligned = alignment.aligned

    reference_blocks = aligned[0]
    read_blocks = aligned[1]

    for ref_block, read_block in zip(
        reference_blocks,
        read_blocks,
    ):
        ref_start = int(ref_block[0])
        ref_end = int(ref_block[1])

        read_start = int(read_block[0])
        read_end = int(read_block[1])

        blocks.append(
            (
                (ref_start, ref_end),
                (read_start, read_end),
            )
        )

    return blocks


def _reconstruct_gapped_alignment(
    reference: str,
    read: str,
    alignment: Any,
) -> Tuple[str, str]:
    """
    Reconstruct gapped alignment strings for diagnostics/display.

    Variant calling itself does not use these strings. It uses the original
    PairwiseAligner coordinate blocks directly.
    """
    blocks = _alignment_blocks(alignment)

    if not blocks:
        return "", ""

    gapped_reference: List[str] = []
    gapped_read: List[str] = []

    previous_ref_end: Optional[int] = None
    previous_read_end: Optional[int] = None

    for block_number, (
        (ref_start, ref_end),
        (read_start, read_end),
    ) in enumerate(blocks):

        if block_number == 0:
            previous_ref_end = ref_start
            previous_read_end = read_start

        assert previous_ref_end is not None
        assert previous_read_end is not None

        ref_gap = ref_start - previous_ref_end
        read_gap = read_start - previous_read_end

        if ref_gap < 0 or read_gap < 0:
            raise ValueError(
                "Alignment blocks are not monotonic."
            )

        if ref_gap > 0 and read_gap > 0:
            raise ValueError(
                "Invalid alignment: both reference and read advance "
                "between adjacent blocks."
            )

        # Deletion relative to the read.
        if ref_gap > 0:
            gapped_reference.append(
                reference[
                    previous_ref_end:ref_start
                ]
            )

            gapped_read.append(
                "-" * ref_gap
            )

        # Insertion relative to the reference.
        if read_gap > 0:
            gapped_reference.append(
                "-" * read_gap
            )

            gapped_read.append(
                read[
                    previous_read_end:read_start
                ]
            )

        block_ref_length = ref_end - ref_start
        block_read_length = read_end - read_start

        if block_ref_length != block_read_length:
            raise ValueError(
                "Aligned reference/read blocks have unequal lengths."
            )

        gapped_reference.append(
            reference[ref_start:ref_end]
        )

        gapped_read.append(
            read[read_start:read_end]
        )

        previous_ref_end = ref_end
        previous_read_end = read_end

    return (
        "".join(gapped_reference),
        "".join(gapped_read),
    )


# ---------------------------------------------------------------------------
# Quality helpers
# ---------------------------------------------------------------------------

def _get_quality_value(
    quality: Optional[Sequence[Any]],
    read_index: Optional[int],
) -> Optional[float]:
    """
    Return the PHRED quality corresponding to an original read position.
    """
    if quality is None:
        return None

    if read_index is None:
        return None

    if read_index < 0 or read_index >= len(quality):
        return None

    try:
        value = float(quality[read_index])
    except (TypeError, ValueError):
        return None

    if pd.isna(value):
        return None

    return value


def _confidence_from_quality(
    observed_base: str,
    quality: Optional[float],
    min_phred: int,
) -> str:
    """
    Classify a nucleotide according to its PHRED quality.
    """
    observed_base = normalise_base(observed_base)

    if observed_base == "-":
        return "Deletion"

    if quality is None:
        return "Missing"

    if is_iupac_ambiguous(observed_base):
        if quality >= min_phred:
            return "Ambiguous_Het"

        return "LowConfidence_Heterozygous"

    if quality >= 30:
        return "High"

    if quality >= min_phred:
        return "Medium"

    return "Low"


# ---------------------------------------------------------------------------
# Single-read local alignment
# ---------------------------------------------------------------------------

def _align_single_read(
    reference: str,
    read: str,
    filename: str = "",
    quality: Optional[Sequence[Any]] = None,
) -> Dict[str, Any]:
    """
    Perform local alignment for one sequencing read.

    The original PHRED quality sequence is retained in the returned
    dictionary.
    """
    reference = normalise_base(reference)
    read = normalise_base(read)

    if not reference:
        raise ValueError(
            "Reference sequence is empty."
        )

    if not read:
        raise ValueError(
            f"Read sequence is empty: {filename}"
        )

    if quality is not None:
        try:
            quality_length = len(quality)
        except TypeError as exc:
            raise TypeError(
                f"Quality scores for {filename or 'read'} "
                "must be a sequence."
            ) from exc

        if quality_length != len(read):
            raise ValueError(
                f"Quality length ({quality_length}) does not match "
                f"read length ({len(read)}) for "
                f"{filename or 'read'}."
            )

    aligner = create_aligner()

    alignments = aligner.align(
        reference,
        read,
    )

    if len(alignments) == 0:
        raise ValueError(
            f"No local alignment found for "
            f"{filename or 'read'}."
        )

    alignment = alignments[0]

    blocks = _alignment_blocks(alignment)

    if not blocks:
        raise ValueError(
            f"Alignment contains no aligned blocks for "
            f"{filename or 'read'}."
        )

    aligned_reference, aligned_read = (
        _reconstruct_gapped_alignment(
            reference,
            read,
            alignment,
        )
    )

    return {
        "Filename": filename,
        "Reference": reference,
        "Read": read,

        "Aligned_Reference": aligned_reference,
        "Aligned_Read": aligned_read,

        "aligned_ref": aligned_reference,
        "aligned_read": aligned_read,

        "Alignment_Score": float(
            alignment.score
        ),
        "alignment_score": float(
            alignment.score
        ),

        "Alignment_Blocks": blocks,
        "alignment_blocks": blocks,

        "Alignment": alignment,

        # IMPORTANT:
        # Preserve the original per-base PHRED scores.
        "Quality": quality,
        "quality": quality,
    }


# ---------------------------------------------------------------------------
# Public alignment function
# ---------------------------------------------------------------------------

def perform_local_alignment(
    reference: Any,
    read: Any = None,
    filename: str = "",
    quality: Optional[Sequence[Any]] = None,
) -> Any:
    """
    Perform local alignment for one read or multiple reads.

    Supported forms
    ---------------

    Batch form used by main.py:

        perform_local_alignment(
            processed_reads,
            reference_amplicon,
        )

    Single-read form:

        perform_local_alignment(
            reference,
            read,
            filename,
            quality,
        )
    """
    # ------------------------------------------------------------------
    # Batch API
    #
    # main.py supplies:
    #
    #     perform_local_alignment(
    #         processed_reads,
    #         reference_amplicon,
    #     )
    # ------------------------------------------------------------------
    if (
        isinstance(reference, (list, tuple))
        and isinstance(read, str)
    ):
        return align_reads(
            reference=read,
            reads=reference,
        )

    # ------------------------------------------------------------------
    # Single-read API
    # ------------------------------------------------------------------
    return _align_single_read(
        reference=reference,
        read=read,
        filename=filename,
        quality=quality,
    )


# ---------------------------------------------------------------------------
# Alignment record generation
# ---------------------------------------------------------------------------

def _make_alignment_record(
    *,
    filename: str,
    reference: str,
    read: str,
    reference_index: Optional[int],
    read_index: Optional[int],
    quality: Optional[float],
    alignment_column: int,
    min_phred: int,
    alignment_score: float,
) -> Dict[str, Any]:
    """
    Create one nucleotide-level alignment record.
    """
    # Reference information.
    if reference_index is not None:
        reference_base = normalise_base(
            reference[reference_index]
        )
        amplicon_position = reference_index + 1
    else:
        reference_base = "-"
        amplicon_position = None

    # Read information.
    if read_index is not None:
        observed_base = normalise_base(
            read[read_index]
        )
        read_position = read_index + 1
    else:
        observed_base = "-"
        read_position = None

    # Alignment event.
    if reference_index is None:
        event = "Insertion"

    elif read_index is None:
        event = "Deletion"

    elif observed_base == reference_base:
        event = "Match"

    else:
        event = "Mismatch"

    # Variant status.
    is_variant = (
        reference_index is not None
        and read_index is not None
        and observed_base != reference_base
    )

    variant_base = (
        observed_base
        if is_variant
        else None
    )

    # Confidence.
    confidence = _confidence_from_quality(
        observed_base=observed_base,
        quality=quality,
        min_phred=min_phred,
    )

    return {
        "Sample": filename,
        "Filename": filename,

        "Alignment_Column": alignment_column,

        "Reference_Index": reference_index,
        "Read_Index": read_index,

        "Amplicon_Position": amplicon_position,
        "Read_Position": read_position,

        "REF": reference_base,
        "Observed_Base": observed_base,

        "Quality": quality,

        "Alignment_Event": event,

        "Is_Variant": is_variant,
        "Variant_Base": variant_base,

        "Confidence": confidence,

        "Alignment_Score": alignment_score,
    }


# ---------------------------------------------------------------------------
# Walk aligned blocks
# ---------------------------------------------------------------------------

def _walk_aligned_block(
    *,
    filename: str,
    reference: str,
    read: str,
    ref_start: int,
    ref_end: int,
    read_start: int,
    read_end: int,
    quality: Optional[Sequence[Any]],
    alignment_column: int,
    min_phred: int,
    alignment_score: float,
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Walk a contiguous aligned PairwiseAligner block.
    """
    ref_length = ref_end - ref_start
    read_length = read_end - read_start

    if ref_length != read_length:
        raise ValueError(
            "Aligned block has unequal reference/read lengths."
        )

    records: List[Dict[str, Any]] = []

    for offset in range(ref_length):
        reference_index = ref_start + offset
        read_index = read_start + offset

        base_quality = _get_quality_value(
            quality,
            read_index,
        )

        alignment_column += 1

        record = _make_alignment_record(
            filename=filename,
            reference=reference,
            read=read,
            reference_index=reference_index,
            read_index=read_index,
            quality=base_quality,
            alignment_column=alignment_column,
            min_phred=min_phred,
            alignment_score=alignment_score,
        )

        records.append(record)

    return records, alignment_column


# ---------------------------------------------------------------------------
# Walk deletions
# ---------------------------------------------------------------------------

def _walk_deletion(
    *,
    filename: str,
    reference: str,
    read: str,
    ref_start: int,
    ref_end: int,
    alignment_column: int,
    min_phred: int,
    alignment_score: float,
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Generate records for reference bases absent from the read.

    Deletions have no read coordinate and therefore no PHRED quality.
    """
    records: List[Dict[str, Any]] = []

    for reference_index in range(
        ref_start,
        ref_end,
    ):
        alignment_column += 1

        record = _make_alignment_record(
            filename=filename,
            reference=reference,
            read=read,
            reference_index=reference_index,
            read_index=None,
            quality=None,
            alignment_column=alignment_column,
            min_phred=min_phred,
            alignment_score=alignment_score,
        )

        records.append(record)

    return records, alignment_column


# ---------------------------------------------------------------------------
# Walk insertions
# ---------------------------------------------------------------------------

def _walk_insertion(
    *,
    filename: str,
    reference: str,
    read: str,
    read_start: int,
    read_end: int,
    quality: Optional[Sequence[Any]],
    alignment_column: int,
    min_phred: int,
    alignment_score: float,
) -> Tuple[List[Dict[str, Any]], int]:
    """
    Generate records for read bases inserted relative to the reference.

    IMPORTANT
    ---------
    Inserted bases retain the PHRED score corresponding to their original
    read coordinate.
    """
    records: List[Dict[str, Any]] = []

    for read_index in range(
        read_start,
        read_end,
    ):
        base_quality = _get_quality_value(
            quality,
            read_index,
        )

        alignment_column += 1

        record = _make_alignment_record(
            filename=filename,
            reference=reference,
            read=read,
            reference_index=None,
            read_index=read_index,
            quality=base_quality,
            alignment_column=alignment_column,
            min_phred=min_phred,
            alignment_score=alignment_score,
        )

        records.append(record)

    return records, alignment_column


# ---------------------------------------------------------------------------
# Main alignment walker
# ---------------------------------------------------------------------------

def _walk_single_alignment(
    alignment_result: Dict[str, Any],
    transcript_start: int,
    min_phred: int,
    quality: Optional[Sequence[Any]] = None,
) -> pd.DataFrame:
    """
    Walk one alignment result.
    """
    if not isinstance(
        alignment_result,
        dict,
    ):
        raise TypeError(
            "alignment_result must be a dictionary."
        )

    filename = alignment_result.get(
        "Filename",
        "",
    )

    reference = normalise_base(
        alignment_result.get(
            "Reference",
            "",
        )
    )

    read = normalise_base(
        alignment_result.get(
            "Read",
            "",
        )
    )

    alignment_score = alignment_result.get(
        "Alignment_Score",
        alignment_result.get(
            "alignment_score"
        ),
    )

    blocks = alignment_result.get(
        "Alignment_Blocks"
    )

    if blocks is None:
        blocks = alignment_result.get(
            "alignment_blocks"
        )

    if not reference:
        raise ValueError(
            f"Reference sequence is missing for {filename}."
        )

    if not read:
        raise ValueError(
            f"Read sequence is missing for {filename}."
        )

    if alignment_score is None:
        raise KeyError(
            f"Alignment score is missing for {filename}."
        )

    if not blocks:
        raise ValueError(
            f"No alignment blocks found for {filename}."
        )

    alignment_score = float(
        alignment_score
    )

    # The alignment-level quality is authoritative.
    if quality is None:
        quality = alignment_result.get(
            "Quality",
            alignment_result.get(
                "quality"
            ),
        )

    if quality is not None:
        try:
            quality_length = len(quality)
        except TypeError as exc:
            raise TypeError(
                f"Quality scores for {filename or 'read'} "
                "must be a sequence."
            ) from exc

        if quality_length != len(read):
            raise ValueError(
                f"Quality length ({quality_length}) does not match "
                f"read length ({len(read)}) for "
                f"{filename or 'read'}."
            )

    records: List[Dict[str, Any]] = []

    alignment_column = 0

    previous_ref_end: Optional[int] = None
    previous_read_end: Optional[int] = None

    for block_number, block in enumerate(blocks):

        if len(block) != 2:
            raise ValueError(
                f"Invalid alignment block {block_number} "
                f"for {filename}."
            )

        ref_block, read_block = block

        if (
            len(ref_block) != 2
            or len(read_block) != 2
        ):
            raise ValueError(
                f"Invalid coordinate block {block_number} "
                f"for {filename}."
            )

        ref_start = int(
            ref_block[0]
        )
        ref_end = int(
            ref_block[1]
        )

        read_start = int(
            read_block[0]
        )
        read_end = int(
            read_block[1]
        )

        # Coordinate validation.
        if not (
            0
            <= ref_start
            <= ref_end
            <= len(reference)
        ):
            raise ValueError(
                f"Reference coordinates out of bounds for "
                f"{filename}: {ref_start}-{ref_end}."
            )

        if not (
            0
            <= read_start
            <= read_end
            <= len(read)
        ):
            raise ValueError(
                f"Read coordinates out of bounds for "
                f"{filename}: {read_start}-{read_end}."
            )

        # Gap before this aligned block.
        if block_number > 0:

            assert previous_ref_end is not None
            assert previous_read_end is not None

            ref_gap = (
                ref_start
                - previous_ref_end
            )

            read_gap = (
                read_start
                - previous_read_end
            )

            if ref_gap < 0 or read_gap < 0:
                raise ValueError(
                    f"Non-monotonic alignment blocks "
                    f"for {filename}."
                )

            if (
                ref_gap > 0
                and read_gap > 0
            ):
                raise ValueError(
                    f"Invalid alignment gap for "
                    f"{filename}: "
                    f"reference gap={ref_gap}, "
                    f"read gap={read_gap}."
                )

            # Deletion.
            if ref_gap > 0:
                gap_records, alignment_column = (
                    _walk_deletion(
                        filename=filename,
                        reference=reference,
                        read=read,
                        ref_start=previous_ref_end,
                        ref_end=ref_start,
                        alignment_column=alignment_column,
                        min_phred=min_phred,
                        alignment_score=alignment_score,
                    )
                )

                records.extend(
                    gap_records
                )

            # Insertion.
            elif read_gap > 0:
                gap_records, alignment_column = (
                    _walk_insertion(
                        filename=filename,
                        reference=reference,
                        read=read,
                        read_start=previous_read_end,
                        read_end=read_start,
                        quality=quality,
                        alignment_column=alignment_column,
                        min_phred=min_phred,
                        alignment_score=alignment_score,
                    )
                )

                records.extend(
                    gap_records
                )

        # Walk the aligned block.
        block_records, alignment_column = (
            _walk_aligned_block(
                filename=filename,
                reference=reference,
                read=read,
                ref_start=ref_start,
                ref_end=ref_end,
                read_start=read_start,
                read_end=read_end,
                quality=quality,
                alignment_column=alignment_column,
                min_phred=min_phred,
                alignment_score=alignment_score,
            )
        )

        records.extend(
            block_records
        )

        previous_ref_end = ref_end
        previous_read_end = read_end

    # ------------------------------------------------------------------
    # Transcript coordinates
    # ------------------------------------------------------------------
    for record in records:

        amplicon_position = record[
            "Amplicon_Position"
        ]

        if amplicon_position is None:
            record["cDNA_Position"] = None

        else:
            record["cDNA_Position"] = (
                transcript_start
                + amplicon_position
                - 1
            )

    # ------------------------------------------------------------------
    # Final variant flags
    # ------------------------------------------------------------------
    for record in records:

        if record[
            "Alignment_Event"
        ] == "Mismatch":

            if (
                record["Confidence"]
                == "LowConfidence_Heterozygous"
            ):
                record["Is_Variant"] = True

        elif record[
            "Alignment_Event"
        ] == "Insertion":

            record["Is_Variant"] = True
            record["Variant_Base"] = (
                record["Observed_Base"]
            )

        elif record[
            "Alignment_Event"
        ] == "Deletion":

            record["Is_Variant"] = True
            record["Variant_Base"] = "-"

    alignment_df = pd.DataFrame(
        records
    )

    required_columns = {
        "Sample",
        "Alignment_Column",
        "Reference_Index",
        "Read_Index",
        "Amplicon_Position",
        "Read_Position",
        "REF",
        "Observed_Base",
        "Quality",
        "Alignment_Event",
        "Is_Variant",
        "Confidence",
        "Alignment_Score",
        "cDNA_Position",
    }

    missing_columns = (
        required_columns
        - set(alignment_df.columns)
    )

    if missing_columns:
        raise KeyError(
            "Generated alignment table is missing "
            f"required columns: {sorted(missing_columns)}"
        )

    return alignment_df


def walk_alignment(
    alignment_result: Any,
    reference_sequence: Optional[str] = None,
    transcript_start: int = 1,
    min_phred: int = 20,
    quality: Optional[Sequence[Any]] = None,
) -> pd.DataFrame:
    """
    Walk one or more alignment results.

    Compatible with main.py:

        walk_alignment(
            alignments,
            reference_sequence,
            amplicon_start,
            MIN_PHRED,
        )

    Parameters
    ----------
    alignment_result
        A single alignment dictionary or a list of alignment dictionaries.

    reference_sequence
        The full transcript/reference supplied by main.py. The actual
        nucleotide walking uses the amplicon stored in each alignment result.

    transcript_start
        One-based transcript coordinate corresponding to amplicon position 1.

    min_phred
        Minimum PHRED threshold.

    quality
        Optional quality array for a single-read utility call.
    """
    # ------------------------------------------------------------------
    # Batch mode
    # ------------------------------------------------------------------
    if isinstance(
        alignment_result,
        (list, tuple),
    ):
        dataframes: List[pd.DataFrame] = []

        for result in alignment_result:

            walked = _walk_single_alignment(
                alignment_result=result,
                transcript_start=transcript_start,
                min_phred=min_phred,
            )

            if not walked.empty:
                dataframes.append(
                    walked
                )

        if not dataframes:
            return pd.DataFrame(
                columns=[
                    "Sample",
                    "Filename",
                    "Alignment_Column",
                    "Reference_Index",
                    "Read_Index",
                    "Amplicon_Position",
                    "Read_Position",
                    "REF",
                    "Observed_Base",
                    "Quality",
                    "Alignment_Event",
                    "Is_Variant",
                    "Variant_Base",
                    "Confidence",
                    "Alignment_Score",
                    "cDNA_Position",
                ]
            )

        alignment_df = pd.concat(
            dataframes,
            ignore_index=True,
        )

    # ------------------------------------------------------------------
    # Single-read mode
    # ------------------------------------------------------------------
    else:
        alignment_df = _walk_single_alignment(
            alignment_result=alignment_result,
            transcript_start=transcript_start,
            min_phred=min_phred,
            quality=quality,
        )

    # ------------------------------------------------------------------
    # Optional reference consistency check.
    #
    # PairwiseAligner coordinates are relative to the amplicon, not the
    # full transcript. Therefore the supplied full transcript is not used
    # for walking.
    # ------------------------------------------------------------------
    if reference_sequence is not None:
        supplied_reference = normalise_base(
            reference_sequence
        )

        if not supplied_reference:
            raise ValueError(
                "Supplied reference sequence is empty."
            )

    return alignment_df


# ---------------------------------------------------------------------------
# Batch alignment
# ---------------------------------------------------------------------------

def align_reads(
    reference: str,
    reads: Sequence[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Align multiple sequencing reads to a reference.

    The input dictionaries may contain:

        Sequence / sequence
        Filename / filename
        Quality / quality / qualities / phred / PHRED

    The PHRED quality array must have the same length as the read sequence.
    """
    if not isinstance(
        reads,
        (list, tuple),
    ):
        raise TypeError(
            "reads must be a sequence of dictionaries."
        )

    results: List[Dict[str, Any]] = []

    for read_data in reads:

        if not isinstance(
            read_data,
            dict,
        ):
            raise TypeError(
                "Each read must be represented "
                "by a dictionary."
            )

        sequence = _extract_read_field(
            read_data,
            "Sequence",
            "sequence",
        )

        filename = _extract_read_field(
            read_data,
            "Filename",
            "filename",
        )

        quality = _extract_read_field(
            read_data,
            "Quality",
            "quality",
            "qualities",
            "phred",
            "PHRED",
        )

        if sequence is None:
            raise ValueError(
                f"Sequence is missing for "
                f"{filename or 'read'}."
            )

        if quality is None:
            raise ValueError(
                f"Per-base PHRED quality is missing for "
                f"{filename or 'read'}."
            )

        result = _align_single_read(
            reference=reference,
            read=sequence,
            filename=filename or "",
            quality=quality,
        )

        # Preserve other metadata generated during preprocessing.
        for key, value in read_data.items():

            if key not in result:
                result[key] = value

        results.append(result)

    return results


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

def validate_alignment_result(
    alignment_result: Dict[str, Any],
) -> None:
    """
    Validate one PairwiseAligner result.
    """
    if not isinstance(
        alignment_result,
        dict,
    ):
        raise TypeError(
            "alignment_result must be a dictionary."
        )

    reference = normalise_base(
        alignment_result.get(
            "Reference",
            "",
        )
    )

    read = normalise_base(
        alignment_result.get(
            "Read",
            "",
        )
    )

    blocks = alignment_result.get(
        "Alignment_Blocks",
        alignment_result.get(
            "alignment_blocks"
        ),
    )

    if not reference:
        raise ValueError(
            "Alignment result has no reference sequence."
        )

    if not read:
        raise ValueError(
            "Alignment result has no read sequence."
        )

    if not blocks:
        raise ValueError(
            "Alignment result has no coordinate blocks."
        )

    alignment_score = alignment_result.get(
        "Alignment_Score",
        alignment_result.get(
            "alignment_score"
        ),
    )

    if alignment_score is None:
        raise ValueError(
            "Alignment result has no alignment score."
        )

    quality = alignment_result.get(
        "Quality",
        alignment_result.get(
            "quality"
        ),
    )

    if quality is None:
        raise ValueError(
            "Alignment result has no PHRED quality array."
        )

    if len(quality) != len(read):
        raise ValueError(
            "PHRED quality length does not match "
            "read sequence length."
        )

    previous_ref_end = None
    previous_read_end = None

    for block in blocks:

        if len(block) != 2:
            raise ValueError(
                "Invalid alignment block."
            )

        (ref_start, ref_end), (
            read_start,
            read_end,
        ) = block

        if not (
            0
            <= ref_start
            <= ref_end
            <= len(reference)
        ):
            raise ValueError(
                "Reference alignment coordinate "
                "out of bounds."
            )

        if not (
            0
            <= read_start
            <= read_end
            <= len(read)
        ):
            raise ValueError(
                "Read alignment coordinate "
                "out of bounds."
            )

        if previous_ref_end is not None:
            if ref_start < previous_ref_end:
                raise ValueError(
                    "Reference alignment blocks overlap."
                )

        if previous_read_end is not None:
            if read_start < previous_read_end:
                raise ValueError(
                    "Read alignment blocks overlap."
                )

        ref_length = (
            ref_end - ref_start
        )

        read_length = (
            read_end - read_start
        )

        if ref_length != read_length:
            raise ValueError(
                "Aligned PairwiseAligner block has "
                "unequal reference/read lengths."
            )

        previous_ref_end = ref_end
        previous_read_end = read_end


def validate_walked_records(
    records: Sequence[Dict[str, Any]],
    reference: str,
    read: str,
) -> None:
    """
    Validate nucleotide-level records generated by walk_alignment().
    """
    if not records:
        raise ValueError(
            "No alignment records were generated."
        )

    reference = normalise_base(
        reference
    )

    read = normalise_base(
        read
    )

    previous_column = 0

    for record in records:

        column = record.get(
            "Alignment_Column"
        )

        if column != previous_column + 1:
            raise ValueError(
                "Alignment columns are not consecutive."
            )

        previous_column = column

        reference_index = record.get(
            "Reference_Index"
        )

        read_index = record.get(
            "Read_Index"
        )

        # Validate reference coordinate.
        if reference_index is not None:

            if not (
                0
                <= reference_index
                < len(reference)
            ):
                raise ValueError(
                    "Reference index is out of bounds."
                )

            expected_ref = reference[
                reference_index
            ].upper()

            if record["REF"] != expected_ref:
                raise ValueError(
                    "REF does not match reference sequence."
                )

        # Validate read coordinate.
        if read_index is not None:

            if not (
                0
                <= read_index
                < len(read)
            ):
                raise ValueError(
                    "Read index is out of bounds."
                )

            expected_observed = read[
                read_index
            ].upper()

            if (
                record["Observed_Base"]
                != expected_observed
            ):
                raise ValueError(
                    "Observed_Base does not match "
                    "read sequence."
                )

            # Every actual read base should have a quality
            # value in the corrected pipeline.
            if record.get("Quality") is None:
                raise ValueError(
                    "A read base has no PHRED quality value."
                )

        event = record.get(
            "Alignment_Event"
        )

        if event == "Match":

            if (
                record["REF"]
                != record["Observed_Base"]
            ):
                raise ValueError(
                    "Match record contains different bases."
                )

        elif event == "Insertion":

            if reference_index is not None:
                raise ValueError(
                    "Insertion should not have "
                    "a reference index."
                )

        elif event == "Deletion":

            if read_index is not None:
                raise ValueError(
                    "Deletion should not have "
                    "a read index."
                )


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------

def align_and_walk(
    reference: str,
    read: str,
    filename: str = "",
    quality: Optional[Sequence[Any]] = None,
    min_phred: int = 20,
    transcript_start: int = 1,
) -> pd.DataFrame:
    """
    Perform alignment and immediately walk one read.
    """
    alignment_result = _align_single_read(
        reference=reference,
        read=read,
        filename=filename,
        quality=quality,
    )

    validate_alignment_result(
        alignment_result
    )

    alignment_df = _walk_single_alignment(
        alignment_result=alignment_result,
        transcript_start=transcript_start,
        min_phred=min_phred,
    )

    validate_walked_records(
        records=alignment_df.to_dict(
            "records"
        ),
        reference=alignment_result["Reference"],
        read=alignment_result["Read"],
    )

    return alignment_df