"""
preprocessing.py

Functions for reading ABI chromatogram files,
extracting nucleotide sequences, and calculating
sequencing quality metrics.
"""

import os
import pandas as pd
from Bio import SeqIO


# ======================================================
# LOAD ABI FILES
# ======================================================

def load_ab1_files(ab1_folder: str) -> list:
    """
    Load ABI chromatogram (.ab1) files.

    Parameters
    ----------
    ab1_folder : str
        Directory containing ABI chromatograms.

    Returns
    -------
    list
        List of dictionaries containing sequence,
        quality scores and ABI metadata.
    """

    abi_records = []

    ab1_files = sorted(

        os.path.join(ab1_folder, filename)

        for filename in os.listdir(ab1_folder)

        if filename.lower().endswith(".ab1")

    )

    for ab1_file in ab1_files:

        abi_record = SeqIO.read(ab1_file, "abi")

        raw = abi_record.annotations["abif_raw"]

        abi_records.append({

            "filename": os.path.basename(ab1_file),

            "sequence": str(abi_record.seq).upper(),

            "quality": abi_record.letter_annotations["phred_quality"],

            "CLRG1": raw.get("CLRG1", 1),

            "CLRG2": raw.get("CLRG2", len(abi_record.seq))

        })

    return abi_records


# ======================================================
# PREPARE READS
# ======================================================

def prepare_reads(abi_records: list) -> list:
    """
    Prepare sequencing reads for downstream analysis.

    Parameters
    ----------
    abi_records : list

    Returns
    -------
    list
        Processed sequencing reads.
    """

    processed_reads = []

    for read in abi_records:

        quality = read["quality"]

        processed_reads.append({

            "filename": read["filename"],

            "sequence": read["sequence"],

            "quality": quality,

            "CLRG1": read["CLRG1"],

            "CLRG2": read["CLRG2"],

            "length": len(read["sequence"]),

            "mean_quality": round(sum(quality) / len(quality), 2),

            "max_quality": max(quality),

            "min_quality": min(quality)

        })

    return processed_reads


# ======================================================
# QUALITY CONTROL SUMMARY
# ======================================================

def calculate_qc_metrics(processed_reads: list) -> pd.DataFrame:
    """
    Calculate sequencing quality metrics.

    Parameters
    ----------
    processed_reads : list

    Returns
    -------
    pandas.DataFrame
        Quality-control summary table.
    """

    qc_summary = pd.DataFrame({

        "Sample":

            [r["filename"] for r in processed_reads],

        "Length":

            [r["length"] for r in processed_reads],

        "Mean_PHRED":

            [r["mean_quality"] for r in processed_reads],

        "Maximum_PHRED":

            [r["max_quality"] for r in processed_reads],

        "Minimum_PHRED":

            [r["min_quality"] for r in processed_reads],

        "CLRG1":

            [r["CLRG1"] for r in processed_reads],

        "CLRG2":

            [r["CLRG2"] for r in processed_reads]

    })

    return qc_summary


# ======================================================
# SAVE QC SUMMARY
# ======================================================

def save_qc_summary(
    qc_summary: pd.DataFrame,
    output_folder: str
) -> None:
    """
    Save quality-control summary.

    Parameters
    ----------
    qc_summary : pandas.DataFrame

    output_folder : str
    """

    os.makedirs(output_folder, exist_ok=True)

    qc_summary.to_csv(

        os.path.join(

            output_folder,

            "QC_Summary.csv"

        ),

        index=False

    )
