"""
main.py

Main entry point for the SangerSeq Variant Pipeline.
"""

from config import *

from scripts.preprocessing import (

    load_ab1_files,

    prepare_reads,

    calculate_qc_metrics,

    save_qc_summary

)


def main():

    # --------------------------------------------
    # Load ABI chromatograms
    # --------------------------------------------

    abi_records = load_ab1_files(

        AB1_FOLDER

    )

    # --------------------------------------------
    # Prepare sequencing reads
    # --------------------------------------------

    processed_reads = prepare_reads(

        abi_records

    )

    # --------------------------------------------
    # Quality-control summary
    # --------------------------------------------

    qc_summary = calculate_qc_metrics(

        processed_reads

    )

    save_qc_summary(

        qc_summary,

        OUTPUT_FOLDER

    )

    print(

        f"Pipeline completed successfully.\n"

        f"Loaded {len(abi_records)} ABI chromatograms."

    )


if __name__ == "__main__":

    main()
