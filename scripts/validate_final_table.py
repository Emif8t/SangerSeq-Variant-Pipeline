"""
validate_final_table.py

Generic validation of the final Sanger sequencing variant table.

The validation checks:

1. Required files exist.
2. Final table contains records.
3. Required columns are present.
4. HGVS variants are valid and unique.
5. REF / ALT fields are populated.
6. Carrier counts agree with carrier sample lists.
7. Variant frequencies agree with carrier counts.
8. Genotype-table carrier identities agree with the final table.
9. HGVS table and final table contain the same variants.
10. Annotation columns are populated where available.

This script deliberately contains NO gene-specific,
transcript-specific, sample-specific, or variant-specific
assumptions.

The current gene/transcript used by the pipeline is defined
in config.py and is therefore not hard-coded here.
"""

import pandas as pd
from pathlib import Path


# ============================================================
# PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent


FINAL_FILE = (
    PROJECT_ROOT
    / "output"
    / "results"
    / "Final_Annotated_Variants.csv"
)


GENOTYPE_FILE = (
    PROJECT_ROOT
    / "output"
    / "genotypes"
    / "Genotype_Table.csv"
)


HGVS_FILE = (
    PROJECT_ROOT
    / "output"
    / "variants"
    / "HGVS_Table.csv"
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def clean(value):
    """
    Convert a value to a stripped string.

    Missing values are returned as an empty string.
    """

    if pd.isna(value):
        return ""

    return str(value).strip()


# ============================================================
# SAMPLE PARSING
# ============================================================

def clean_sample(value):
    """
    Normalise a sample identifier.
    """

    return clean(value).replace(" ", "")


def parse_samples(value):
    """
    Convert a semicolon-separated sample list into a set.
    """

    if pd.isna(value):

        return set()

    value = str(value).strip()

    if value in {
        "",
        "-",
        "nan",
        "None",
    }:

        return set()

    return {
        clean_sample(sample)
        for sample in value.split(";")
        if clean_sample(sample)
    }


# ============================================================
# CHECK FUNCTIONS
# ============================================================

validation_failed = False


def pass_check(message):

    print(
        f"PASS: {message}"
    )


def fail_check(message):

    global validation_failed

    validation_failed = True

    print(
        f"FAIL: {message}"
    )


# ============================================================
# START
# ============================================================

print()

print(
    "=" * 70
)

print(
    "FINAL VARIANT TABLE VALIDATION"
)

print(
    "=" * 70
)

print()


# ============================================================
# 1. CHECK REQUIRED FILES
# ============================================================

print(
    "Checking required files ..."
)

required_files = {

    "Final variant table":
        FINAL_FILE,

    "Genotype table":
        GENOTYPE_FILE,

    "HGVS table":
        HGVS_FILE,

}


for description, file_path in required_files.items():

    if file_path.exists():

        pass_check(
            f"{description}: {file_path}"
        )

    else:

        fail_check(
            f"Missing {description}: {file_path}"
        )


print()


if validation_failed:

    raise SystemExit(
        "\nValidation stopped because required files are missing."
    )


# ============================================================
# 2. LOAD FILES
# ============================================================

print(
    "Loading tables ..."
)

try:

    final_df = pd.read_csv(
        FINAL_FILE
    )

    genotype_df = pd.read_csv(
        GENOTYPE_FILE
    )

    hgvs_df = pd.read_csv(
        HGVS_FILE
    )

except Exception as exc:

    raise RuntimeError(
        f"Unable to load validation tables: {exc}"
    )


# ============================================================
# CLEAN COLUMN NAMES
# ============================================================

final_df.columns = (
    final_df.columns
    .astype(str)
    .str.strip()
)


genotype_df.columns = (
    genotype_df.columns
    .astype(str)
    .str.strip()
)


hgvs_df.columns = (
    hgvs_df.columns
    .astype(str)
    .str.strip()
)


print()

print(
    f"Final table records: {len(final_df)}"
)

print(
    f"Genotype records: {len(genotype_df)}"
)

print(
    f"HGVS records: {len(hgvs_df)}"
)

print()


# ============================================================
# 3. CHECK FINAL TABLE IS NOT EMPTY
# ============================================================

print(
    "-" * 70
)

print(
    "1. FINAL TABLE CHECK"
)

print(
    "-" * 70
)


if len(final_df) > 0:

    pass_check(
        f"Final table contains {len(final_df)} variant record(s)"
    )

else:

    fail_check(
        "Final variant table is empty"
    )


print()


# ============================================================
# 4. REQUIRED FINAL TABLE COLUMNS
# ============================================================

print(
    "-" * 70
)

print(
    "2. REQUIRED FINAL COLUMNS"
)

print(
    "-" * 70
)


required_final_columns = [

    "HGVS_cDNA",

    "REF",

    "ALT",

    "Carrier_Count",

    "Variant_Frequency",

    "Samples",

]


missing_final_columns = [

    column

    for column in required_final_columns

    if column not in final_df.columns

]


if not missing_final_columns:

    pass_check(
        "All required final-table columns are present"
    )

else:

    fail_check(
        "Missing final-table columns: "
        + ", ".join(missing_final_columns)
    )


print()


# ============================================================
# STOP IF CORE COLUMNS ARE MISSING
# ============================================================

if missing_final_columns:

    raise SystemExit(
        "\nValidation stopped because required final-table "
        "columns are missing."
    )


# ============================================================
# 5. HGVS VALIDATION
# ============================================================

print(
    "-" * 70
)

print(
    "3. HGVS VALIDATION"
)

print(
    "-" * 70
)


final_hgvs = [

    clean(value)

    for value in final_df["HGVS_cDNA"]

]


final_hgvs = [

    value

    for value in final_hgvs

    if value

]


if len(final_hgvs) == len(final_df):

    pass_check(
        "All final variants have HGVS_cDNA values"
    )

else:

    fail_check(
        "One or more final variants have missing HGVS_cDNA"
    )


# ------------------------------------------------------------
# Duplicate HGVS
# ------------------------------------------------------------

if len(final_hgvs) == len(set(final_hgvs)):

    pass_check(
        "HGVS_cDNA values are unique"
    )

else:

    duplicates = (

        pd.Series(final_hgvs)

        .value_counts()

        .loc[lambda x: x > 1]

        .index

        .tolist()

    )

    fail_check(
        "Duplicate HGVS variants found: "
        + ", ".join(duplicates)
    )


# ============================================================
# 6. HGVS TABLE CROSS-CHECK
# ============================================================

print()

print(
    "-" * 70
)

print(
    "4. HGVS TABLE CROSS-CHECK"
)

print(
    "-" * 70
)


if "HGVS_cDNA" not in hgvs_df.columns:

    fail_check(
        "HGVS table does not contain HGVS_cDNA"
    )

else:

    hgvs_table_values = [

        clean(value)

        for value in hgvs_df["HGVS_cDNA"]

        if clean(value)

    ]

    final_hgvs_set = set(
        final_hgvs
    )

    hgvs_table_set = set(
        hgvs_table_values
    )

    missing_from_final = (
        hgvs_table_set
        - final_hgvs_set
    )

    missing_from_hgvs = (
        final_hgvs_set
        - hgvs_table_set
    )

    if not missing_from_final and not missing_from_hgvs:

        pass_check(
            "Final table and HGVS table contain the same HGVS variants"
        )

    else:

        if missing_from_final:

            fail_check(
                "HGVS variants present in HGVS table "
                "but missing from final table: "
                + ", ".join(
                    sorted(missing_from_final)
                )
            )

        if missing_from_hgvs:

            fail_check(
                "HGVS variants present in final table "
                "but missing from HGVS table: "
                + ", ".join(
                    sorted(missing_from_hgvs)
                )
            )


# ============================================================
# 7. REF / ALT VALIDATION
# ============================================================

print()

print(
    "-" * 70
)

print(
    "5. REF / ALT VALIDATION"
)

print(
    "-" * 70
)


for index, row in final_df.iterrows():

    hgvs = clean(
        row["HGVS_cDNA"]
    )

    ref = clean(
        row["REF"]
    )

    alt = clean(
        row["ALT"]
    )


    if ref:

        pass_check(
            f"{hgvs}: REF = {ref}"
        )

    else:

        fail_check(
            f"{hgvs}: missing REF"
        )


    if alt:

        pass_check(
            f"{hgvs}: ALT = {alt}"
        )

    else:

        fail_check(
            f"{hgvs}: missing ALT"
        )


# ============================================================
# 8. CARRIER COUNT VALIDATION
# ============================================================

print()

print(
    "-" * 70
)

print(
    "6. CARRIER COUNT VALIDATION"
)

print(
    "-" * 70
)


for _, row in final_df.iterrows():

    hgvs = clean(
        row["HGVS_cDNA"]
    )

    try:

        carrier_count = int(
            row["Carrier_Count"]
        )

    except (
        ValueError,
        TypeError
    ):

        fail_check(
            f"{hgvs}: invalid Carrier_Count"
        )

        continue


    samples = parse_samples(
        row["Samples"]
    )

    sample_count = len(
        samples
    )


    print()

    print(
        f"Variant: {hgvs}"
    )

    print(
        f"Carrier_Count: {carrier_count}"
    )

    print(
        f"Samples: {sorted(samples)}"
    )


    if carrier_count == sample_count:

        pass_check(
            f"{hgvs}: Carrier_Count matches number of unique samples"
        )

    else:

        fail_check(
            f"{hgvs}: Carrier_Count ({carrier_count}) "
            f"does not match sample count ({sample_count})"
        )


# ============================================================
# 9. VARIANT FREQUENCY VALIDATION
# ============================================================

print()

print(
    "-" * 70
)

print(
    "7. VARIANT FREQUENCY VALIDATION"
)

print(
    "-" * 70
)


# ------------------------------------------------------------
# Determine total number of samples
# ------------------------------------------------------------

if "Sample" in genotype_df.columns:

    all_samples = {

        clean_sample(sample)

        for sample in genotype_df["Sample"]

        if clean_sample(sample)

    }

    total_samples = len(
        all_samples
    )

else:

    total_samples = 0


print(
    f"Total unique samples in genotype table: "
    f"{total_samples}"
)


for _, row in final_df.iterrows():

    hgvs = clean(
        row["HGVS_cDNA"]
    )

    try:

        carrier_count = int(
            row["Carrier_Count"]
        )

        observed_frequency = float(
            row["Variant_Frequency"]
        )

    except (
        ValueError,
        TypeError
    ):

        fail_check(
            f"{hgvs}: invalid numeric frequency fields"
        )

        continue


    if total_samples == 0:

        fail_check(
            f"{hgvs}: unable to determine total sample count"
        )

        continue


    expected_frequency = (
        carrier_count
        / total_samples
    )


    print()

    print(
        f"Variant: {hgvs}"
    )

    print(
        f"Expected frequency: "
        f"{expected_frequency:.6f}"
    )

    print(
        f"Observed frequency: "
        f"{observed_frequency:.6f}"
    )


    if abs(
        observed_frequency
        - expected_frequency
    ) < 1e-9:

        pass_check(
            f"{hgvs}: Variant_Frequency is correct"
        )

    else:

        fail_check(
            f"{hgvs}: Variant_Frequency is incorrect"
        )


# ============================================================
# 10. GENOTYPE TABLE STRUCTURE
# ============================================================

print()

print(
    "-" * 70
)

print(
    "8. GENOTYPE TABLE VALIDATION"
)

print(
    "-" * 70
)


required_genotype_columns = [

    "Sample",

    "cDNA_Position",

    "REF",

    "Observed_Base",

    "ALT",

    "Zygosity",

    "Is_Variant",

]


missing_genotype_columns = [

    column

    for column in required_genotype_columns

    if column not in genotype_df.columns

]


if not missing_genotype_columns:

    pass_check(
        "All required genotype-table columns are present"
    )

else:

    fail_check(
        "Missing genotype-table columns: "
        + ", ".join(
            missing_genotype_columns
        )
    )


# ============================================================
# 11. GENOTYPE → FINAL TABLE CROSS-CHECK
# ============================================================

if not missing_genotype_columns:

    print()

    print(
        "-" * 70
    )

    print(
        "9. GENOTYPE CARRIER CROSS-CHECK"
    )

    print(
        "-" * 70
    )


    if "Transcript_Position" in final_df.columns:

        position_column = (
            "Transcript_Position"
        )

    elif "cDNA_Position" in final_df.columns:

        position_column = (
            "cDNA_Position"
        )

    else:

        position_column = None


    if position_column is None:

        print(
            "INFO: No transcript-position column available "
            "for positional genotype cross-check."
        )

    else:

        for _, row in final_df.iterrows():

            hgvs = clean(
                row["HGVS_cDNA"]
            )

            try:

                position = int(
                    row[position_column]
                )

            except (
                ValueError,
                TypeError
            ):

                fail_check(
                    f"{hgvs}: invalid position"
                )

                continue


            records = genotype_df[

                genotype_df["cDNA_Position"]
                == position

            ].copy()


            if records.empty:

                fail_check(
                    f"{hgvs}: no genotype records at "
                    f"position {position}"
                )

                continue


            records["Sample_Clean"] = (

                records["Sample"]

                .apply(clean_sample)

            )


            records["Observed_Clean"] = (

                records["Observed_Base"]

                .apply(clean)

            )


            records["Is_Variant_Clean"] = (

                records["Is_Variant"]

                .astype(str)

                .str.lower()

                .eq("true")

            )


            expected_alt = clean(
                row["ALT"]
            )


            carrier_records = records[

                records["Is_Variant_Clean"]

                & (

                    records["Observed_Clean"]

                    == expected_alt

                )

            ]


            genotype_samples = set(

                carrier_records[
                    "Sample_Clean"
                ]

            )


            final_samples = parse_samples(
                row["Samples"]
            )


            print()

            print(
                f"Variant: {hgvs}"
            )

            print(
                f"Genotype carriers: "
                f"{sorted(genotype_samples)}"
            )

            print(
                f"Final-table carriers: "
                f"{sorted(final_samples)}"
            )


            if genotype_samples == final_samples:

                pass_check(
                    f"{hgvs}: genotype carriers match "
                    "final-table carriers"
                )

            else:

                fail_check(
                    f"{hgvs}: genotype carriers do not "
                    "match final-table carriers"
                )


# ============================================================
# 12. ANNOTATION VALIDATION
# ============================================================

print()

print(
    "-" * 70
)

print(
    "10. VEP ANNOTATION VALIDATION"
)

print(
    "-" * 70
)


annotation_columns = [

    "Consequence",

    "IMPACT",

    "SYMBOL",

    "Gene",

    "Feature",

    "HGVSc",

    "HGVSp",

]


available_annotation_columns = [

    column

    for column in annotation_columns

    if column in final_df.columns

]


if not available_annotation_columns:

    print(
        "INFO: No standard VEP annotation columns "
        "were found in the final table."
    )

else:

    for column in available_annotation_columns:

        populated = (

            final_df[column]

            .notna()

            .sum()

        )


        if populated == len(final_df):

            pass_check(
                f"{column}: populated for all variants"
            )

        else:

            fail_check(
                f"{column}: only "
                f"{populated}/{len(final_df)} "
                "variants populated"
            )


# ============================================================
# 13. NUMERIC FIELD VALIDATION
# ============================================================

print()

print(
    "-" * 70
)

print(
    "11. NUMERIC FIELD VALIDATION"
)

print(
    "-" * 70
)


numeric_columns = [

    "Carrier_Count",

    "Variant_Frequency",

]


for column in numeric_columns:

    converted = pd.to_numeric(

        final_df[column],

        errors="coerce"

    )


    invalid = converted.isna().sum()


    if invalid == 0:

        pass_check(
            f"{column}: all values are numeric"
        )

    else:

        fail_check(
            f"{column}: {invalid} invalid numeric value(s)"
        )


# ============================================================
# 14. FINAL TABLE DISPLAY
# ============================================================

print()

print(
    "=" * 70
)

print(
    "FINAL VARIANT TABLE"
)

print(
    "=" * 70
)


display_columns = [

    "HGVS_cDNA",

    "REF",

    "ALT",

    "Consequence",

    "IMPACT",

    "Carrier_Count",

    "Variant_Frequency",

    "Samples",

]


available_columns = [

    column

    for column in display_columns

    if column in final_df.columns

]


print()

print(

    final_df[

        available_columns

    ]

    .to_string(
        index=False
    )

)


# ============================================================
# 15. FINAL RESULT
# ============================================================

print()

print(
    "=" * 70
)

print(
    "FINAL VALIDATION RESULT"
)

print(
    "=" * 70
)


if validation_failed:

    print()

    print(
        "OVERALL RESULT: FAIL"
    )

    print()

    print(
        "One or more validation checks failed."
    )

    print(
        "Review the FAIL messages above."
    )

else:

    print()

    print(
        "OVERALL RESULT: PASS"
    )

    print()

    print(
        "The final variant table passed all "
        "generic validation checks."
    )

    print(
        "Carrier counts, frequencies, HGVS records, "
        "and genotype carrier identities are consistent."
    )


print()

print(
    "=" * 70
)

print(
    "VALIDATION COMPLETE"
)

print(
    "=" * 70
)

print()


# ============================================================
# EXIT STATUS
# ============================================================

if validation_failed:

    raise SystemExit(1)

else:

    raise SystemExit(0)