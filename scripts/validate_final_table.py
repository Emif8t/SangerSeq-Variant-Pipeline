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
    / "ASS1_Final_Annotated_Variants.csv"
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
# EXPECTED FINAL VARIANTS
#
# These are defined by HGVS_Table.csv.
# The carrier information is expected to come from
# Genotype_Table.csv.
# ============================================================

EXPECTED = {
    "NM_000050.4:c.783T>C": {
        "hgvs_cdna_position": 783,
        "transcript_position": 1139,
        "ref": "T",
        "alt": "C",
        "carrier_count": 3,
        "samples": {
            "A5.ab1",
            "A7.ab1",
            "A9.ab1",
        },
    },

    "NM_000050.4:c.876T>C": {
        "hgvs_cdna_position": 876,
        "transcript_position": 1232,
        "ref": "T",
        "alt": "C",
        "carrier_count": 2,
        "samples": {
            "A5.ab1",
            "A9.ab1",
        },
    },
}


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def clean(value):
    """Clean a simple value."""
    if pd.isna(value):
        return ""
    return str(value).strip()


def clean_sample(value):
    """Normalise sample names."""
    return clean(value).replace(" ", "")


def parse_samples(value):
    """Convert semicolon-separated samples into a set."""
    if pd.isna(value):
        return set()

    return {
        clean_sample(x)
        for x in str(value).split(";")
        if clean_sample(x)
    }


def pass_check(message):
    print(f"PASS: {message}")


def fail_check(message):
    print(f"FAIL: {message}")


# ============================================================
# START
# ============================================================

print()
print("=" * 70)
print("FINAL VARIANT TABLE VALIDATION")
print("=" * 70)
print()


# ============================================================
# CHECK REQUIRED FILES
# ============================================================

print("Checking required files ...")

required_files = [
    FINAL_FILE,
    GENOTYPE_FILE,
    HGVS_FILE,
]

for file_path in required_files:

    if file_path.exists():
        pass_check(str(file_path))
    else:
        fail_check(f"Missing file: {file_path}")
        raise SystemExit(1)

print()


# ============================================================
# LOAD FILES
# ============================================================

final_df = pd.read_csv(FINAL_FILE)
genotype_df = pd.read_csv(GENOTYPE_FILE)
hgvs_df = pd.read_csv(HGVS_FILE)


# Strip whitespace from column names
final_df.columns = final_df.columns.astype(str).str.strip()
genotype_df.columns = genotype_df.columns.astype(str).str.strip()
hgvs_df.columns = hgvs_df.columns.astype(str).str.strip()


print(f"Final table records: {len(final_df)}")
print(f"Genotype records: {len(genotype_df)}")
print(f"HGVS records: {len(hgvs_df)}")
print()


# ============================================================
# 1. FINAL NUMBER OF VARIANTS
# ============================================================

print("-" * 70)
print("1. NUMBER OF FINAL VARIANTS")
print("-" * 70)

if len(final_df) == 2:
    pass_check("Final table contains exactly 2 variants")
else:
    fail_check(
        f"Expected 2 final variants, found {len(final_df)}"
    )

print()


# ============================================================
# 2. HGVS VARIANT CHECK
# ============================================================

print("-" * 70)
print("2. HGVS VARIANT VALIDATION")
print("-" * 70)

actual_hgvs = {
    clean(x)
    for x in final_df["HGVS_cDNA"]
}

expected_hgvs = set(EXPECTED.keys())

if actual_hgvs == expected_hgvs:
    pass_check(
        "Final HGVS variants exactly match the expected HGVS variants"
    )
else:
    fail_check(
        f"HGVS mismatch. Found: {sorted(actual_hgvs)}"
    )

print()

for hgvs, expected in EXPECTED.items():

    print(f"Variant: {hgvs}")

    rows = final_df[
        final_df["HGVS_cDNA"].astype(str).str.strip()
        == hgvs
    ]

    if rows.empty:
        fail_check("Variant is missing from final table")
        print()
        continue

    row = rows.iloc[0]

    # --------------------------------------------------------
    # REF
    # --------------------------------------------------------

    actual_ref = clean(row["REF"])

    if actual_ref == expected["ref"]:
        pass_check(f"REF = {actual_ref}")
    else:
        fail_check(
            f"REF mismatch: expected {expected['ref']}, "
            f"found {actual_ref}"
        )

    # --------------------------------------------------------
    # ALT
    # --------------------------------------------------------

    actual_alt = clean(row["ALT"])

    if actual_alt == expected["alt"]:
        pass_check(f"ALT = {actual_alt}")
    else:
        fail_check(
            f"ALT mismatch: expected {expected['alt']}, "
            f"found {actual_alt}"
        )

    print()


# ============================================================
# 3. CARRIER VALIDATION
# ============================================================

print("-" * 70)
print("3. CARRIER VALIDATION")
print("-" * 70)

carrier_validation_pass = True

for hgvs, expected in EXPECTED.items():

    rows = final_df[
        final_df["HGVS_cDNA"].astype(str).str.strip()
        == hgvs
    ]

    if rows.empty:
        carrier_validation_pass = False
        continue

    row = rows.iloc[0]

    actual_count = int(row["Carrier_Count"])

    actual_samples = parse_samples(
        row["Samples"]
    )

    expected_samples = {
        clean_sample(x)
        for x in expected["samples"]
    }

    print()
    print(hgvs)

    # --------------------------------------------------------
    # Carrier count
    # --------------------------------------------------------

    print(
        f"Expected carrier count: {expected['carrier_count']}"
    )

    print(
        f"Actual carrier count:   {actual_count}"
    )

    if actual_count == expected["carrier_count"]:
        pass_check("Carrier count is correct")
    else:
        fail_check("Carrier count is incorrect")
        carrier_validation_pass = False

    # --------------------------------------------------------
    # Carrier samples
    # --------------------------------------------------------

    print(
        f"Expected carriers: {sorted(expected_samples)}"
    )

    print(
        f"Actual carriers:   {sorted(actual_samples)}"
    )

    if actual_samples == expected_samples:
        pass_check("Carrier sample identities are correct")
    else:
        fail_check("Carrier sample identities are incorrect")
        carrier_validation_pass = False

    # --------------------------------------------------------
    # Variant frequency
    # --------------------------------------------------------

    actual_frequency = float(
        row["Variant_Frequency"]
    )

    expected_frequency = (
        expected["carrier_count"] / 3
    )

    print(
        f"Expected frequency: {expected_frequency:.6f}"
    )

    print(
        f"Actual frequency:   {actual_frequency:.6f}"
    )

    if abs(actual_frequency - expected_frequency) < 1e-6:
        pass_check("Variant frequency is correct")
    else:
        fail_check("Variant frequency is incorrect")
        carrier_validation_pass = False


# ============================================================
# 4. VEP ANNOTATION VALIDATION
# ============================================================

print()
print("-" * 70)
print("4. VEP ANNOTATION VALIDATION")
print("-" * 70)

required_columns = [
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature",
    "HGVSc",
    "HGVSp",
]

for column in required_columns:

    if column not in final_df.columns:

        fail_check(
            f"Missing annotation column: {column}"
        )

    else:

        populated = (
            final_df[column]
            .notna()
            .sum()
        )

        if populated == len(final_df):
            pass_check(
                f"{column} populated for all variants"
            )
        else:
            fail_check(
                f"{column}: {populated}/{len(final_df)} populated"
            )


# ============================================================
# 5. EXPECTED VEP CONSEQUENCES
# ============================================================

print()
print("-" * 70)
print("5. VEP CONSEQUENCE VALIDATION")
print("-" * 70)

for hgvs in expected_hgvs:

    row = final_df[
        final_df["HGVS_cDNA"].astype(str).str.strip()
        == hgvs
    ]

    if row.empty:
        continue

    row = row.iloc[0]

    consequence = clean(
        row["Consequence"]
    )

    impact = clean(
        row["IMPACT"]
    )

    if consequence == "synonymous_variant":
        pass_check(
            f"{hgvs}: Consequence = synonymous_variant"
        )
    else:
        fail_check(
            f"{hgvs}: unexpected consequence = {consequence}"
        )

    if impact == "LOW":
        pass_check(
            f"{hgvs}: IMPACT = LOW"
        )
    else:
        fail_check(
            f"{hgvs}: unexpected IMPACT = {impact}"
        )


# ============================================================
# 6. GENOTYPE TABLE CROSS-CHECK
# ============================================================

print()
print("-" * 70)
print("6. GENOTYPE TABLE CROSS-CHECK")
print("-" * 70)

required_genotype_columns = [
    "Sample",
    "cDNA_Position",
    "REF",
    "Observed_Base",
    "ALT",
    "Zygosity",
    "Is_Variant",
]

missing = [
    c
    for c in required_genotype_columns
    if c not in genotype_df.columns
]

if missing:

    fail_check(
        f"Missing genotype columns: {missing}"
    )

else:

    for hgvs, expected in EXPECTED.items():

        position = expected["transcript_position"]

        records = genotype_df[
            genotype_df["cDNA_Position"]
            == position
        ].copy()

        print()
        print(hgvs)
        print(
            f"Checking genotype position: {position}"
        )

        if records.empty:

            fail_check(
                f"No genotype records found at position {position}"
            )

            continue

        records["Sample_Clean"] = (
            records["Sample"]
            .apply(clean_sample)
        )

        records["ALT_Clean"] = (
            records["ALT"]
            .apply(clean)
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

        carrier_records = records[
            records["Is_Variant_Clean"]
            &
            (
                records["Observed_Clean"]
                == expected["alt"]
            )
        ]

        genotype_samples = set(
            carrier_records["Sample_Clean"]
        )

        expected_samples = {
            clean_sample(x)
            for x in expected["samples"]
        }

        print(
            "Genotype-table carriers:",
            sorted(genotype_samples)
        )

        print(
            "Expected carriers:",
            sorted(expected_samples)
        )

        if genotype_samples == expected_samples:

            pass_check(
                "Genotype carriers match final carrier list"
            )

        else:

            fail_check(
                "Genotype carriers DO NOT match final carrier list"
            )


# ============================================================
# 7. PRINT FINAL TABLE
# ============================================================

print()
print("=" * 70)
print("FINAL VARIANT TABLE")
print("=" * 70)

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
    c for c in display_columns
    if c in final_df.columns
]

print()
print(
    final_df[available_columns]
    .to_string(index=False)
)


# ============================================================
# FINAL RESULT
# ============================================================

print()
print("=" * 70)
print("FINAL VALIDATION RESULT")
print("=" * 70)

if (
    len(final_df) == 2
    and actual_hgvs == expected_hgvs
    and carrier_validation_pass
):

    print()
    print("OVERALL RESULT: PASS")
    print()
    print(
        "The final variant table contains the correct 2 HGVS variants."
    )
    print(
        "The carrier counts and carrier identities are correct."
    )
    print(
        "The final table is consistent with the expected genotype results."
    )

else:

    print()
    print("OVERALL RESULT: FAIL")
    print()
    print(
        "One or more validation checks failed."
    )

print()
print("=" * 70)
print("VALIDATION COMPLETE")
print("=" * 70)