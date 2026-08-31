# ============================================================
# final_table.py
#
# FINAL ANNOTATED VARIANT TABLE
#
# ASS1 Sanger Sequencing Variant Pipeline
#
# ============================================================
#
# PURPOSE
# -------
# This script creates the final annotated variant table using
# three independent sources:
#
# 1. HGVS_Table.csv
#       Defines WHICH variants are the final variants.
#
# 2. Genotype_Table.csv
#       Determines WHICH samples are confirmed carriers.
#
# 3. VEP_HGVS_OUTPUT.xlsx
#       Provides functional and database annotation.
#
#
# IMPORTANT DESIGN PRINCIPLE
# --------------------------
# VEP annotation MUST NOT determine:
#
#   - Carrier_Count
#   - Variant_Frequency
#   - Samples
#
# Those fields are calculated exclusively from the confirmed
# genotype calls.
#
#
# POSITION MAPPING
# ----------------
#
# HGVS_Table.csv:
#
#   Transcript_Position = 1139
#   HGVS_cDNA_Position = 783
#   HGVS_cDNA = NM_000050.4:c.783T>C
#
#   Transcript_Position = 1232
#   HGVS_cDNA_Position = 876
#   HGVS_cDNA = NM_000050.4:c.876T>C
#
#
# Genotype_Table.csv:
#
#   cDNA_Position = 1139
#   cDNA_Position = 1232
#
# Therefore:
#
#   HGVS_Table.Transcript_Position
#                |
#                v
#   Genotype_Table.cDNA_Position
#
# is used for genotype matching.
#
#
# EXPECTED RESULT FROM THE CURRENT DATA
# --------------------------------------
#
# NM_000050.4:c.783T>C
#   Carriers: A5.ab1, A7.ab1, A9.ab1
#   Carrier_Count: 3
#   Frequency: 1.000
#
# NM_000050.4:c.876T>C
#   Carriers: A5.ab1, A9.ab1
#   Carrier_Count: 2
#   Frequency: 0.667
#
# The A7.ab1 call at c.876 is:
#
#   LowConfidence_Heterozygous
#   Is_Variant = False
#
# and is therefore NOT counted.
#
# ============================================================


from pathlib import Path

import pandas as pd
import numpy as np


# ============================================================
# 1. PROJECT PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

HGVS_FILE = (
    PROJECT_ROOT
    / "output"
    / "variants"
    / "HGVS_Table.csv"
)

GENOTYPE_FILE = (
    PROJECT_ROOT
    / "output"
    / "genotypes"
    / "Genotype_Table.csv"
)

VEP_FILE = (
    PROJECT_ROOT
    / "data"
    / "annotation"
    / "VEP_HGVS_OUTPUT.xlsx"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "output"
    / "results"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

OUTPUT_FILE = (
    OUTPUT_DIR
    / "ASS1_Final_Annotated_Variants.csv"
)


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def clean_string(value):
    """
    Convert a value to a clean string.

    Missing values are represented by '-'.
    """

    if pd.isna(value):
        return "-"

    value = str(value).strip()

    if value == "":
        return "-"

    if value.lower() in {
        "nan",
        "none",
        "null",
        "na",
        "n/a"
    }:
        return "-"

    return value


# ------------------------------------------------------------
# Normalise HGVS
# ------------------------------------------------------------

def normalise_hgvs(value):
    """
    Normalise HGVS strings by removing whitespace.

    Example:

        NM_000050.4:c. 783T>C

    becomes:

        NM_000050.4:c.783T>C
    """

    if pd.isna(value):
        return ""

    value = str(value).strip()

    # Remove ALL whitespace
    value = "".join(value.split())

    return value


# ------------------------------------------------------------
# Normalise sample names
# ------------------------------------------------------------

def normalise_sample(value):
    """
    Normalise sample identifiers.

    Example:

        A5. ab1

    becomes:

        A5.ab1
    """

    if pd.isna(value):
        return ""

    value = str(value).strip()

    # Remove accidental spaces
    value = value.replace(" ", "")

    return value


# ------------------------------------------------------------
# Numeric conversion
# ------------------------------------------------------------

def to_numeric(value):
    """
    Safely convert a value to numeric.
    """

    return pd.to_numeric(
        value,
        errors="coerce"
    )


# ------------------------------------------------------------
# Boolean conversion
# ------------------------------------------------------------

def is_true(value):
    """
    Convert common representations of TRUE/FALSE
    into a Python boolean.
    """

    if isinstance(value, bool):
        return value

    if pd.isna(value):
        return False

    value = str(value).strip().lower()

    return value in {
        "true",
        "1",
        "yes",
        "y"
    }


# ============================================================
# 3. START
# ============================================================

print()
print("=" * 70)
print("BUILDING FINAL ANNOTATED VARIANT TABLE")
print("=" * 70)
print()


# ============================================================
# 4. CHECK INPUT FILES
# ============================================================

for file_path, description in [
    (HGVS_FILE, "HGVS table"),
    (GENOTYPE_FILE, "genotype table"),
    (VEP_FILE, "VEP annotation"),
]:

    if not file_path.exists():

        raise FileNotFoundError(
            f"\n{description} not found:\n"
            f"{file_path}\n"
        )


# ============================================================
# 5. LOAD HGVS TABLE
# ============================================================

hgvs_df = pd.read_csv(
    HGVS_FILE
)

hgvs_df.columns = (
    hgvs_df.columns
    .astype(str)
    .str.strip()
)


# ------------------------------------------------------------
# Validate required columns
# ------------------------------------------------------------

required_hgvs_columns = [
    "Transcript",
    "Transcript_Position",
    "HGVS_cDNA",
    "REF",
    "ALT",
]

missing_hgvs = [
    col
    for col in required_hgvs_columns
    if col not in hgvs_df.columns
]

if missing_hgvs:

    raise KeyError(
        "HGVS_Table.csv is missing required columns:\n"
        + "\n".join(
            f"  - {x}"
            for x in missing_hgvs
        )
    )


# ============================================================
# 6. NORMALISE HGVS VARIANTS
# ============================================================

hgvs_df["HGVS_NORMALISED"] = (
    hgvs_df["HGVS_cDNA"]
    .apply(normalise_hgvs)
)


# Remove empty HGVS records
hgvs_df = hgvs_df[
    hgvs_df["HGVS_NORMALISED"] != ""
].copy()


# Remove duplicate final variants
hgvs_df = (
    hgvs_df
    .drop_duplicates(
        subset=["HGVS_NORMALISED"],
        keep="first"
    )
    .reset_index(drop=True)
)


print(
    f"HGVS variants found: {len(hgvs_df)}"
)

print()

print(
    "Final variants defined by HGVS_Table.csv:"
)


for i, row in hgvs_df.iterrows():

    print(
        f"{i + 1}. "
        f"{row['HGVS_NORMALISED']} "
        f"("
        f"{clean_string(row['REF'])}"
        ">"
        f"{clean_string(row['ALT'])}"
        ")"
    )


print()


# ============================================================
# 7. LOAD GENOTYPE TABLE
# ============================================================

genotype_df = pd.read_csv(
    GENOTYPE_FILE
)

genotype_df.columns = (
    genotype_df.columns
    .astype(str)
    .str.strip()
)


required_genotype_columns = [
    "Sample",
    "cDNA_Position",
    "REF",
    "ALT",
    "Is_Variant",
]

missing_genotype = [
    col
    for col in required_genotype_columns
    if col not in genotype_df.columns
]

if missing_genotype:

    raise KeyError(
        "Genotype_Table.csv is missing required columns:\n"
        + "\n".join(
            f"  - {x}"
            for x in missing_genotype
        )
    )


print(
    f"Genotype records loaded: {len(genotype_df)}"
)


# ============================================================
# 8. NORMALISE GENOTYPE DATA
# ============================================================

genotype_df["Sample_Clean"] = (
    genotype_df["Sample"]
    .apply(normalise_sample)
)

genotype_df["cDNA_Position_Numeric"] = (
    to_numeric(
        genotype_df["cDNA_Position"]
    )
)

genotype_df["REF_Clean"] = (
    genotype_df["REF"]
    .apply(clean_string)
    .str.upper()
)

genotype_df["ALT_Clean"] = (
    genotype_df["ALT"]
    .apply(clean_string)
    .str.upper()
)

genotype_df["Is_Variant_Clean"] = (
    genotype_df["Is_Variant"]
    .apply(is_true)
)


# ============================================================
# 9. EXCLUDE LOW-CONFIDENCE / AMBIGUOUS CALLS
# ============================================================

# Start with Is_Variant == True
confirmed_genotypes = genotype_df[
    genotype_df["Is_Variant_Clean"]
].copy()


# ------------------------------------------------------------
# Remove low-confidence ambiguous heterozygous calls
# ------------------------------------------------------------

if "Confidence" in confirmed_genotypes.columns:

    confidence = (
        confirmed_genotypes["Confidence"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    ambiguous_confidence = confidence.isin({
        "ambiguous_het",
        "lowconfidence_heterozygous",
        "low confidence heterozygous",
    })

    confirmed_genotypes = confirmed_genotypes[
        ~ambiguous_confidence
    ].copy()


# ------------------------------------------------------------
# Remove low-confidence zygosity calls
# ------------------------------------------------------------

if "Zygosity" in confirmed_genotypes.columns:

    zygosity = (
        confirmed_genotypes["Zygosity"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    low_confidence_zygosity = (
        zygosity.str.contains(
            "lowconfidence",
            na=False
        )
        |
        zygosity.str.contains(
            "low confidence",
            na=False
        )
    )

    confirmed_genotypes = confirmed_genotypes[
        ~low_confidence_zygosity
    ].copy()


print(
    f"Confirmed genotype calls: "
    f"{len(confirmed_genotypes)}"
)

print()


# ============================================================
# 10. TOTAL NUMBER OF SAMPLES
# ============================================================

all_samples = (
    genotype_df["Sample_Clean"]
    .replace("", np.nan)
    .dropna()
    .drop_duplicates()
    .tolist()
)

total_samples = len(all_samples)

print(
    f"Total samples detected: {total_samples}"
)

print(
    "Samples:",
    ";".join(all_samples)
)

print()


# ============================================================
# 11. BUILD GENOTYPE-BASED CARRIER SUMMARY
# ============================================================

print(
    "Matching HGVS variants to confirmed genotype carriers ..."
)

print()


carrier_records = []


for _, hgvs_row in hgvs_df.iterrows():

    hgvs = hgvs_row["HGVS_NORMALISED"]

    ref = (
        clean_string(
            hgvs_row["REF"]
        )
        .upper()
    )

    alt = (
        clean_string(
            hgvs_row["ALT"]
        )
        .upper()
    )


    # --------------------------------------------------------
    # CRITICAL POSITION MAPPING
    #
    # HGVS_Table:
    #     Transcript_Position = 1139
    #
    # Genotype_Table:
    #     cDNA_Position = 1139
    #
    # Therefore these two columns are matched.
    #
    # DO NOT use HGVS_cDNA_Position here because:
    #
    #     1139 -> c.783
    #     1232 -> c.876
    #
    # --------------------------------------------------------

    transcript_position = to_numeric(
        hgvs_row["Transcript_Position"]
    )


    if pd.isna(transcript_position):

        raise ValueError(
            f"Missing Transcript_Position for "
            f"{hgvs}"
        )


    print(
        f"  Checking genotype position: "
        f"{int(transcript_position)}"
    )


    # --------------------------------------------------------
    # Match genotype position
    # --------------------------------------------------------

    position_matches = confirmed_genotypes[
        confirmed_genotypes[
            "cDNA_Position_Numeric"
        ]
        == transcript_position
    ].copy()


    # --------------------------------------------------------
    # Match REF
    # --------------------------------------------------------

    position_matches = position_matches[
        position_matches["REF_Clean"]
        == ref
    ].copy()


    # --------------------------------------------------------
    # Match ALT
    # --------------------------------------------------------

    position_matches = position_matches[
        position_matches["ALT_Clean"]
        == alt
    ].copy()


    # --------------------------------------------------------
    # Unique confirmed carriers
    # --------------------------------------------------------

    carriers = (
        position_matches[
            "Sample_Clean"
        ]
        .replace("", np.nan)
        .dropna()
        .drop_duplicates()
        .tolist()
    )


    carrier_count = len(carriers)


    # --------------------------------------------------------
    # Carrier frequency
    # --------------------------------------------------------

    if total_samples > 0:

        variant_frequency = (
            carrier_count
            / total_samples
        )

    else:

        variant_frequency = np.nan


    samples_string = ";".join(
        carriers
    )


    print(
        f"{hgvs}: "
        f"{carrier_count} carrier(s) -> "
        f"{samples_string if samples_string else 'None'}"
    )

    print()


    # --------------------------------------------------------
    # Store genotype-derived information
    # --------------------------------------------------------

    carrier_records.append({

        # Internal merge key
        "HGVS_NORMALISED": hgvs,

        # Final displayed fields
        "HGVS_cDNA": hgvs,

        "REF": ref,

        "ALT": alt,

        "Carrier_Count": carrier_count,

        "Variant_Frequency": variant_frequency,

        "Samples": samples_string,

    })


# ============================================================
# 12. CREATE CARRIER DATAFRAME
# ============================================================

carrier_df = pd.DataFrame(
    carrier_records
)


# ============================================================
# 13. VALIDATE CARRIER DATAFRAME
# ============================================================

required_carrier_columns = [
    "HGVS_NORMALISED",
    "HGVS_cDNA",
    "REF",
    "ALT",
    "Carrier_Count",
    "Variant_Frequency",
    "Samples",
]

missing_carrier = [
    col
    for col in required_carrier_columns
    if col not in carrier_df.columns
]

if missing_carrier:

    raise RuntimeError(
        "Internal carrier summary is missing:\n"
        + "\n".join(
            f"  - {x}"
            for x in missing_carrier
        )
    )


# ============================================================
# 14. LOAD VEP ANNOTATION
# ============================================================

print(
    "Loading VEP annotation ..."
)

vep_df = pd.read_excel(
    VEP_FILE
)

vep_df.columns = (
    vep_df.columns
    .astype(str)
    .str.strip()
)


# ============================================================
# 15. FIND VEP HGVS COLUMN
# ============================================================

vep_hgvs_candidates = [
    "#Uploaded_variation",
    "Uploaded_variation",
    "HGVS_cDNA",
]

vep_hgvs_column = None


for col in vep_hgvs_candidates:

    if col in vep_df.columns:

        vep_hgvs_column = col
        break


if vep_hgvs_column is None:

    raise KeyError(
        "Could not find an HGVS column in "
        "VEP_HGVS_OUTPUT.xlsx.\n\n"
        "Expected one of:\n"
        + "\n".join(
            f"  - {x}"
            for x in vep_hgvs_candidates
        )
    )


# ============================================================
# 16. NORMALISE VEP HGVS
# ============================================================

vep_df["HGVS_NORMALISED"] = (
    vep_df[
        vep_hgvs_column
    ]
    .apply(normalise_hgvs)
)


# ============================================================
# 17. DEDUPLICATE VEP ANNOTATIONS
# ============================================================

vep_unique = (
    vep_df[
        vep_df["HGVS_NORMALISED"] != ""
    ]
    .drop_duplicates(
        subset=[
            "HGVS_NORMALISED"
        ],
        keep="first"
    )
    .copy()
)


# ============================================================
# 18. MATCH VEP TO FINAL HGVS VARIANTS
# ============================================================

matched_vep = carrier_df[
    ["HGVS_NORMALISED"]
].merge(
    vep_unique[
        ["HGVS_NORMALISED"]
    ],
    on="HGVS_NORMALISED",
    how="left",
    indicator=True
)


vep_match_count = (
    matched_vep["_merge"]
    .eq("both")
    .sum()
)


print(
    "VEP HGVS matching: "
    f"{vep_match_count} of "
    f"{len(carrier_df)} variants matched."
)

print()


# ============================================================
# 19. MERGE VEP WITH CARRIER SUMMARY
# ============================================================
#
# IMPORTANT:
#
# carrier_df is on the LEFT.
#
# Therefore the genotype-derived fields remain the primary
# fields.
#
# VEP is annotation only.
#
# ============================================================

final_df = carrier_df.merge(
    vep_unique,
    on="HGVS_NORMALISED",
    how="left",
    suffixes=(
        "",
        "_VEP"
    )
)


# ============================================================
# 20. REMOVE CONFLICTING VEP-DERIVED FIELDS
# ============================================================
#
# These fields must NEVER replace the genotype-derived fields.
#
# ============================================================

conflicting_columns = [

    "Carrier_Count_VEP",
    "Variant_Carrier_Count",

    "Variant_Frequency_VEP",

    "Samples_VEP",

    "Variant_Calls",
    "Variant_Samples",

    "Genotype_Calls",
    "Genotype_Sample_Count",

]


for col in conflicting_columns:

    if col in final_df.columns:

        final_df.drop(
            columns=[col],
            inplace=True
        )


# ============================================================
# 21. REMOVE INTERNAL VEP/UPLOAD IDENTIFIERS
# ============================================================

for col in [
    "#Uploaded_variation",
    "Uploaded_variation",
]:

    if col in final_df.columns:

        final_df.drop(
            columns=[col],
            inplace=True
        )


# ============================================================
# 22. REMOVE INTERNAL MERGE KEY
# ============================================================

if "HGVS_NORMALISED" in final_df.columns:

    final_df.drop(
        columns=[
            "HGVS_NORMALISED"
        ],
        inplace=True
    )


# ============================================================
# 23. RESTORE / PROTECT CORE GENOTYPE FIELDS
# ============================================================
#
# This is an additional safety measure.
#
# Recalculate these fields directly from carrier_df after
# annotation so that even an unexpected VEP column can never
# overwrite them.
#
# ============================================================

core_genotype_fields = carrier_df[
    [
        "HGVS_cDNA",
        "Carrier_Count",
        "Variant_Frequency",
        "Samples",
    ]
].copy()


# Remove any accidental versions
for col in [
    "Carrier_Count",
    "Variant_Frequency",
    "Samples",
]:

    if col in final_df.columns:

        final_df.drop(
            columns=[col],
            inplace=True
        )


# Merge the protected genotype information back
final_df = final_df.merge(
    core_genotype_fields,
    on="HGVS_cDNA",
    how="left"
)


# ============================================================
# 24. RESTORE REF / ALT FROM HGVS TABLE
# ============================================================
#
# REF and ALT are defined by HGVS_Table.csv.
#
# ============================================================

hgvs_core = hgvs_df[
    [
        "HGVS_NORMALISED",
        "REF",
        "ALT",
    ]
].copy()

hgvs_core["HGVS_cDNA"] = (
    hgvs_core["HGVS_NORMALISED"]
)

hgvs_core = hgvs_core[
    [
        "HGVS_cDNA",
        "REF",
        "ALT",
    ]
]


# Remove possible duplicate REF/ALT
for col in [
    "REF",
    "ALT",
]:

    if col in final_df.columns:

        final_df.drop(
            columns=[col],
            inplace=True
        )


final_df = final_df.merge(
    hgvs_core,
    on="HGVS_cDNA",
    how="left"
)


# ============================================================
# 25. ADD VARIANT TYPE
# ============================================================

if "Variant_Type" not in final_df.columns:

    final_df["Variant_Type"] = "SNV"


# ============================================================
# 26. REORDER CORE COLUMNS
# ============================================================

preferred_order = [

    # --------------------------------------------------------
    # Variant identity
    # --------------------------------------------------------

    "HGVS_cDNA",

    "REF",

    "ALT",

    "Variant_Type",

    # --------------------------------------------------------
    # Functional annotation
    # --------------------------------------------------------

    "Consequence",

    "IMPACT",

    "SIFT",

    "PolyPhen",

    "HGVSc",

    "HGVSp",

    # --------------------------------------------------------
    # Database annotation
    # --------------------------------------------------------

    "Existing_variation",

    "dbSNP_ID",

    "COSMIC_ID",

    "ClinVar_ID",

    "HGMD_ID",

    "Other_ID",

    # --------------------------------------------------------
    # Gene annotation
    # --------------------------------------------------------

    "SYMBOL",

    "Gene",

    "Transcript",

    "Feature",

    "BIOTYPE",

    "Location",

    # --------------------------------------------------------
    # Position annotation
    # --------------------------------------------------------

    "Transcript_Position",

    "HGVS_cDNA_Position",

    "cDNA_position",

    "CDS_position",

    "Protein_position",

    "Amino_acids",

    "Codons",

    # --------------------------------------------------------
    # CRITICAL GENOTYPE-DERIVED INFORMATION
    # --------------------------------------------------------

    "Carrier_Count",

    "Variant_Frequency",

    "Samples",

    # --------------------------------------------------------
    # Quality
    # --------------------------------------------------------

    "Mean_Quality",

    "Mean_Alignment_Score",

]


existing_preferred = [
    col
    for col in preferred_order
    if col in final_df.columns
]


remaining_columns = [
    col
    for col in final_df.columns
    if col not in existing_preferred
]


final_df = final_df[
    existing_preferred
    + remaining_columns
]


# ============================================================
# 27. FORMAT NUMERIC FIELDS
# ============================================================

if "Carrier_Count" in final_df.columns:

    final_df["Carrier_Count"] = (
        pd.to_numeric(
            final_df["Carrier_Count"],
            errors="coerce"
        )
        .fillna(0)
        .astype(int)
    )


if "Variant_Frequency" in final_df.columns:

    final_df["Variant_Frequency"] = (
        pd.to_numeric(
            final_df["Variant_Frequency"],
            errors="coerce"
        )
    )


# ============================================================
# 28. VALIDATE FINAL NUMBER OF VARIANTS
# ============================================================

if len(final_df) != len(hgvs_df):

    raise RuntimeError(
        "\nFINAL VARIANT COUNT ERROR\n"
        f"HGVS_Table.csv variants: "
        f"{len(hgvs_df)}\n"
        f"Final table variants: "
        f"{len(final_df)}\n"
    )


# ============================================================
# 29. VALIDATE THAT HGVS VARIANTS ARE PRESERVED
# ============================================================

expected_hgvs = (
    hgvs_df["HGVS_NORMALISED"]
    .tolist()
)

final_hgvs = (
    final_df["HGVS_cDNA"]
    .apply(normalise_hgvs)
    .tolist()
)


if expected_hgvs != final_hgvs:

    raise RuntimeError(
        "\nHGVS ORDER/IDENTITY ERROR\n\n"
        f"Expected:\n{expected_hgvs}\n\n"
        f"Final:\n{final_hgvs}\n"
    )


# ============================================================
# 30. VALIDATE CARRIER COUNTS
# ============================================================

print(
    "=" * 70
)

print(
    "FINAL CARRIER SUMMARY"
)

print(
    "=" * 70
)


for _, row in final_df.iterrows():

    hgvs = row["HGVS_cDNA"]

    carrier_count = int(
        row["Carrier_Count"]
    )

    frequency = row[
        "Variant_Frequency"
    ]

    samples = (
        ""
        if pd.isna(row["Samples"])
        else str(row["Samples"])
    )


    if samples in {
        "",
        "-",
        "nan",
    }:

        sample_list = []

    else:

        sample_list = [
            x.strip()
            for x in samples.split(";")
            if x.strip()
        ]


    sample_count = len(
        sample_list
    )


    print(
        f"{hgvs}"
    )

    print(
        f"  REF: {row['REF']}"
    )

    print(
        f"  ALT: {row['ALT']}"
    )

    print(
        f"  Carrier_Count: {carrier_count}"
    )

    print(
        f"  Variant_Frequency: "
        f"{frequency:.3f}"
    )

    print(
        f"  Samples: "
        f"{samples if samples else 'None'}"
    )

    print()


    # --------------------------------------------------------
    # Carrier count must equal number of unique samples
    # --------------------------------------------------------

    if carrier_count != sample_count:

        raise RuntimeError(
            f"\nCARRIER COUNT ERROR for {hgvs}\n"
            f"Carrier_Count = {carrier_count}\n"
            f"Number of Samples = {sample_count}\n"
            f"Samples = {samples}\n"
        )


    # --------------------------------------------------------
    # Frequency validation
    # --------------------------------------------------------

    if total_samples > 0:

        expected_frequency = (
            carrier_count
            / total_samples
        )

        if not np.isclose(
            frequency,
            expected_frequency,
            atol=1e-9
        ):

            raise RuntimeError(
                f"\nFREQUENCY ERROR for {hgvs}\n"
                f"Carrier_Count = {carrier_count}\n"
                f"Total samples = {total_samples}\n"
                f"Expected frequency = "
                f"{expected_frequency}\n"
                f"Observed frequency = "
                f"{frequency}\n"
            )


# ============================================================
# 31. SAVE FINAL CSV
# ============================================================

final_df.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# 32. FINAL OUTPUT INFORMATION
# ============================================================

print(
    "=" * 70
)

print(
    "Final variant table saved to:"
)

print(
    OUTPUT_FILE
)

print()

print(
    f"Final number of variants: "
    f"{len(final_df)}"
)

print()


# ============================================================
# 33. PRINT FINAL CORE TABLE
# ============================================================

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


print(
    final_df[
        [
            col
            for col in display_columns
            if col in final_df.columns
        ]
    ].to_string(
        index=False
    )
)


# ============================================================
# 34. FINAL SUCCESS MESSAGE
# ============================================================

print()

print(
    "=" * 70
)

print(
    "PIPELINE STEP COMPLETED SUCCESSFULLY"
)

print(
    "=" * 70
)

print()