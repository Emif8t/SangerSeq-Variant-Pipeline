"""
Final annotated variant table construction for the
SangerSeq Variant Pipeline.

This module:

1. Loads VEP annotation results
2. Standardises VEP annotation columns
3. Merges HGVS, genotype and VEP information
4. Splits database identifiers
5. Cleans and standardises the final table
6. Reorders columns for publication
7. Saves the final annotated variant table
"""

import os
import pandas as pd


# =========================================================
# 1. LOAD VEP ANNOTATION
# =========================================================

def load_vep_annotation(
    vep_file: str
) -> pd.DataFrame:
    """
    Load VEP annotation results from an Excel file.

    Parameters
    ----------
    vep_file : str
        Path to the VEP Excel output file.

    Returns
    -------
    pandas.DataFrame
        VEP annotation table.
    """

    if not os.path.exists(vep_file):
        raise FileNotFoundError(
            f"VEP annotation file not found: {vep_file}"
        )

    vep_df = pd.read_excel(
        vep_file
    )

    if vep_df.empty:
        raise ValueError(
            "The VEP annotation file is empty."
        )

    return vep_df


# =========================================================
# 2. PREPARE VEP ANNOTATION
# =========================================================

def prepare_vep_annotation(
    vep_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Prepare and standardise VEP annotation data.

    The function:
    - removes completely empty rows
    - removes completely empty columns
    - strips whitespace from column names
    - standardises missing values
    - preserves the original VEP information

    Parameters
    ----------
    vep_df : pandas.DataFrame
        Raw VEP annotation table.

    Returns
    -------
    pandas.DataFrame
        Prepared VEP annotation table.
    """

    vep_df = vep_df.copy()

    # -----------------------------------------------------
    # Remove completely empty rows and columns
    # -----------------------------------------------------

    vep_df = vep_df.dropna(
        axis=0,
        how="all"
    )

    vep_df = vep_df.dropna(
        axis=1,
        how="all"
    )

    # -----------------------------------------------------
    # Clean column names
    # -----------------------------------------------------

    vep_df.columns = [
        str(column).strip()
        for column in vep_df.columns
    ]

    # -----------------------------------------------------
    # Replace missing values
    # -----------------------------------------------------

    vep_df = vep_df.fillna("")

    return vep_df


# =========================================================
# 3. MERGE VARIANT TABLES
# =========================================================

def merge_variant_tables(
    hgvs_df: pd.DataFrame,
    genotype_df: pd.DataFrame,
    vep_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge HGVS, genotype and VEP annotation tables.

    The preferred merge key is Transcript_Position.
    If this is not available, the function attempts to use
    HGVS_cDNA_Position.

    Parameters
    ----------
    hgvs_df : pandas.DataFrame
        HGVS-annotated variant table.

    genotype_df : pandas.DataFrame
        Genotype/carrier information.

    vep_df : pandas.DataFrame
        Prepared VEP annotation table.

    Returns
    -------
    pandas.DataFrame
        Merged variant table.
    """

    hgvs_df = hgvs_df.copy()
    genotype_df = genotype_df.copy()
    vep_df = vep_df.copy()

    # -----------------------------------------------------
    # Clean column names
    # -----------------------------------------------------

    for df in [hgvs_df, genotype_df, vep_df]:

        df.columns = [
            str(column).strip()
            for column in df.columns
        ]

    # =====================================================
    # STEP 1 — Merge HGVS and genotype information
    # =====================================================

    # If both tables contain Transcript_Position, use it.
    if (
        "Transcript_Position" in hgvs_df.columns
        and "Transcript_Position" in genotype_df.columns
    ):

        genotype_columns = [
            column
            for column in genotype_df.columns
            if column not in hgvs_df.columns
            or column in [
                "Samples",
                "Carrier_Count",
                "Genotype",
                "Zygosity"
            ]
        ]

        genotype_subset = genotype_df[
            [
                column
                for column in genotype_columns
                if column in genotype_df.columns
            ]
        ].copy()

        final_df = hgvs_df.merge(
            genotype_subset,
            on="Transcript_Position",
            how="left",
            suffixes=("", "_genotype")
        )

    else:

        # -------------------------------------------------
        # Alternative: merge using cDNA position
        # -------------------------------------------------

        if (
            "HGVS_cDNA_Position" in hgvs_df.columns
            and "cDNA_Position" in genotype_df.columns
        ):

            final_df = hgvs_df.merge(
                genotype_df,
                left_on="HGVS_cDNA_Position",
                right_on="cDNA_Position",
                how="left",
                suffixes=("", "_genotype")
            )

        else:

            raise KeyError(
                "Unable to merge HGVS and genotype tables. "
                "Expected either 'Transcript_Position' "
                "in both tables or "
                "'HGVS_cDNA_Position'/'cDNA_Position'."
            )

    # =====================================================
    # STEP 2 — Merge VEP annotation
    # =====================================================

    # -----------------------------------------------------
    # Preferred VEP merge: Transcript_Position
    # -----------------------------------------------------

    if (
        "Transcript_Position" in final_df.columns
        and "Transcript_Position" in vep_df.columns
    ):

        vep_merge_columns = [
            column
            for column in vep_df.columns
            if (
                column == "Transcript_Position"
                or column not in final_df.columns
            )
        ]

        final_df = final_df.merge(
            vep_df[vep_merge_columns],
            on="Transcript_Position",
            how="left"
        )

    # -----------------------------------------------------
    # HGVS cDNA position
    # -----------------------------------------------------

    elif (
        "HGVS_cDNA_Position" in final_df.columns
        and "HGVS_cDNA_Position" in vep_df.columns
    ):

        vep_merge_columns = [
            column
            for column in vep_df.columns
            if (
                column == "HGVS_cDNA_Position"
                or column not in final_df.columns
            )
        ]

        final_df = final_df.merge(
            vep_df[vep_merge_columns],
            on="HGVS_cDNA_Position",
            how="left"
        )

    # -----------------------------------------------------
    # If no common annotation key exists
    # -----------------------------------------------------

    else:

        raise KeyError(
            "Unable to merge VEP annotation. "
            "No compatible variant-position column was found."
        )

    return final_df


# =========================================================
# 4. SPLIT EXISTING VARIATION IDENTIFIERS
# =========================================================

def split_variant_identifiers(
    final_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Split the Existing_variation column into
    database-specific identifiers.

    Parameters
    ----------
    final_df : pandas.DataFrame
        Final merged variant table.

    Returns
    -------
    pandas.DataFrame
        Variant table with database-specific
        identifier columns.
    """

    final_df = final_df.copy()

    # -----------------------------------------------------
    # Determine source column
    # -----------------------------------------------------

    if "Existing_variation" in final_df.columns:

        source_column = "Existing_variation"

    elif "dbSNP_rsID" in final_df.columns:

        source_column = "dbSNP_rsID"

    else:

        # No identifier column available.
        # Create empty identifier columns.

        for column in [
            "dbSNP_ID",
            "COSMIC_ID",
            "ClinVar_ID",
            "HGMD_ID",
            "Other_ID"
        ]:

            if column not in final_df.columns:
                final_df[column] = ""

        return final_df

    # -----------------------------------------------------
    # Create identifier columns
    # -----------------------------------------------------

    identifier_columns = [
        "dbSNP_ID",
        "COSMIC_ID",
        "ClinVar_ID",
        "HGMD_ID",
        "Other_ID"
    ]

    for column in identifier_columns:

        final_df[column] = ""

    # -----------------------------------------------------
    # Parse identifiers
    # -----------------------------------------------------

    for index, value in (
        final_df[source_column]
        .fillna("")
        .items()
    ):

        if str(value).strip() == "":
            continue

        identifiers = [
            x.strip()
            for x in str(value).split(",")
            if x.strip()
        ]

        dbsnp = []
        cosmic = []
        clinvar = []
        hgmd = []
        other = []

        for identifier in identifiers:

            identifier = identifier.strip()

            # -------------------------------------------------
            # dbSNP
            # -------------------------------------------------

            if identifier.startswith("rs"):

                dbsnp.append(identifier)

            # -------------------------------------------------
            # COSMIC
            # -------------------------------------------------

            elif identifier.startswith(
                ("COSV", "COSM")
            ):

                cosmic.append(identifier)

            # -------------------------------------------------
            # ClinVar
            # -------------------------------------------------

            elif identifier.startswith(
                ("VCV", "RCV", "CD")
            ):

                clinvar.append(identifier)

            # -------------------------------------------------
            # HGMD
            # -------------------------------------------------

            elif identifier.startswith(
                ("CM", "CI")
            ):

                hgmd.append(identifier)

            # -------------------------------------------------
            # Other databases
            # -------------------------------------------------

            else:

                other.append(identifier)

        final_df.at[
            index,
            "dbSNP_ID"
        ] = ";".join(dbsnp)

        final_df.at[
            index,
            "COSMIC_ID"
        ] = ";".join(cosmic)

        final_df.at[
            index,
            "ClinVar_ID"
        ] = ";".join(clinvar)

        final_df.at[
            index,
            "HGMD_ID"
        ] = ";".join(hgmd)

        final_df.at[
            index,
            "Other_ID"
        ] = ";".join(other)

    return final_df


# =========================================================
# 5. CLEAN & STANDARDISE FINAL TABLE
# =========================================================

def clean_final_table(
    final_df: pd.DataFrame,
    genotype_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Clean and standardise the final annotated
    variant table.

    Parameters
    ----------
    final_df : pandas.DataFrame
        Merged variant table.

    genotype_df : pandas.DataFrame
        Genotype table used to calculate total samples.

    Returns
    -------
    pandas.DataFrame
        Cleaned final variant table.
    """

    final_df = final_df.copy()

    # -----------------------------------------------------
    # Rename duplicated merge columns
    # -----------------------------------------------------

    rename_columns = {
        "Transcript_x": "Transcript",
        "Samples_x": "Samples"
    }

    for old_name, new_name in rename_columns.items():

        if old_name in final_df.columns:

            if (
                rename_columns[old_name]
                not in final_df.columns
            ):

                final_df.rename(
                    columns={
                        old_name:
                        rename_columns[old_name]
                    },
                    inplace=True
                )

    # -----------------------------------------------------
    # Remove duplicated columns
    # -----------------------------------------------------

    duplicate_columns = [
        column
        for column in [
            "Transcript_y",
            "Samples_y"
        ]
        if column in final_df.columns
    ]

    if duplicate_columns:

        final_df.drop(
            columns=duplicate_columns,
            inplace=True
        )

    # -----------------------------------------------------
    # Variant frequency
    # -----------------------------------------------------

    if "Sample" in genotype_df.columns:

        total_samples = (
            genotype_df["Sample"]
            .astype(str)
            .nunique()
        )

    else:

        total_samples = 0

    if (
        total_samples > 0
        and "Variant_Frequency"
        not in final_df.columns
        and "Carrier_Count"
        in final_df.columns
    ):

        final_df["Variant_Frequency"] = (

            final_df["Carrier_Count"]
            /
            total_samples

        ).round(3)

    # -----------------------------------------------------
    # Ensure identifier columns exist
    # -----------------------------------------------------

    identifier_columns = [
        "dbSNP_ID",
        "COSMIC_ID",
        "ClinVar_ID",
        "HGMD_ID",
        "Other_ID"
    ]

    for column in identifier_columns:

        if column not in final_df.columns:

            final_df[column] = ""

    # -----------------------------------------------------
    # Replace missing annotation values
    # -----------------------------------------------------

    annotation_columns = [
        "Gene",
        "Consequence",
        "Impact",
        "SIFT",
        "PolyPhen",
        "ClinVar"
    ]

    for column in annotation_columns:

        if column in final_df.columns:

            final_df[column] = (
                final_df[column]
                .replace("", pd.NA)
                .fillna("Unknown")
            )

    # -----------------------------------------------------
    # Replace missing identifiers
    # -----------------------------------------------------

    id_columns = [
        "Existing_variation",
        "dbSNP_ID",
        "COSMIC_ID",
        "ClinVar_ID",
        "HGMD_ID",
        "Other_ID"
    ]

    for column in id_columns:

        if column in final_df.columns:

            final_df[column] = (
                final_df[column]
                .replace("", pd.NA)
                .fillna("-")
            )

    # -----------------------------------------------------
    # Replace missing sample names
    # -----------------------------------------------------

    if "Samples" in final_df.columns:

        final_df["Samples"] = (
            final_df["Samples"]
            .fillna("")
        )

    # -----------------------------------------------------
    # Sort variants
    # -----------------------------------------------------

    if "Transcript_Position" in final_df.columns:

        final_df = (
            final_df
            .sort_values(
                "Transcript_Position"
            )
            .reset_index(
                drop=True
            )
        )

    return final_df


# =========================================================
# 6. REORDER COLUMNS
# =========================================================

def reorder_columns(
    final_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Reorder the final variant table into a
    logical publication-oriented structure.

    Columns that are not present are simply skipped.

    Parameters
    ----------
    final_df : pandas.DataFrame
        Cleaned final variant table.

    Returns
    -------
    pandas.DataFrame
        Reordered variant table.
    """

    preferred_columns = [

        # -----------------------------------------------
        # Variant identification
        # -----------------------------------------------

        "Gene",
        "Transcript",
        "Transcript_Position",
        "HGVS_cDNA_Position",
        "HGVS_cDNA",

        # -----------------------------------------------
        # Genomic information
        # -----------------------------------------------

        "Chromosome",
        "Genomic_Position",
        "REF",
        "ALT",

        # -----------------------------------------------
        # Variant annotation
        # -----------------------------------------------

        "Consequence",
        "Impact",
        "SIFT",
        "PolyPhen",

        # -----------------------------------------------
        # Database identifiers
        # -----------------------------------------------

        "Existing_variation",
        "dbSNP_ID",
        "COSMIC_ID",
        "ClinVar_ID",
        "HGMD_ID",
        "Other_ID",

        # -----------------------------------------------
        # Genotype information
        # -----------------------------------------------

        "Genotype",
        "Zygosity",
        "Carrier_Count",
        "Variant_Frequency",
        "Samples",

        # -----------------------------------------------
        # Additional annotation
        # -----------------------------------------------

        "ClinVar"

    ]

    # Keep columns that actually exist
    ordered_columns = [
        column
        for column in preferred_columns
        if column in final_df.columns
    ]

    # Keep any additional columns not explicitly listed
    remaining_columns = [
        column
        for column in final_df.columns
        if column not in ordered_columns
    ]

    final_df = final_df[
        ordered_columns + remaining_columns
    ]

    return final_df


# =========================================================
# 7. BUILD FINAL PUBLICATION TABLE
# =========================================================

def build_final_variant_table(
    hgvs_df: pd.DataFrame,
    genotype_df: pd.DataFrame,
    vep_file: str
) -> pd.DataFrame:
    """
    Build the final publication-ready annotated
    variant table.

    Parameters
    ----------
    hgvs_df : pandas.DataFrame
        HGVS-annotated variant table.

    genotype_df : pandas.DataFrame
        Genotype information.

    vep_file : str
        Path to VEP annotation Excel file.

    Returns
    -------
    pandas.DataFrame
        Final annotated variant table.
    """

    # -----------------------------------------------------
    # Load VEP output
    # -----------------------------------------------------

    vep_df = load_vep_annotation(
        vep_file
    )

    # -----------------------------------------------------
    # Standardise VEP columns
    # -----------------------------------------------------

    vep_df = prepare_vep_annotation(
        vep_df
    )

    # -----------------------------------------------------
    # Merge HGVS, genotype and VEP
    # -----------------------------------------------------

    final_df = merge_variant_tables(
        hgvs_df,
        genotype_df,
        vep_df
    )

    # -----------------------------------------------------
    # Split Existing_variation identifiers
    # -----------------------------------------------------

    final_df = split_variant_identifiers(
        final_df
    )

    # -----------------------------------------------------
    # Clean final table
    # -----------------------------------------------------

    final_df = clean_final_table(
        final_df,
        genotype_df
    )

    # -----------------------------------------------------
    # Reorder columns
    # -----------------------------------------------------

    final_df = reorder_columns(
        final_df
    )

    return final_df


# =========================================================
# 8. SAVE FINAL TABLE
# =========================================================

def save_final_table(
    final_df: pd.DataFrame,
    output_folder: str
):
    """
    Save the publication-ready variant table.

    Parameters
    ----------
    final_df : pandas.DataFrame
        Final annotated variant table.

    output_folder : str
        Directory in which to save the table.
    """

    os.makedirs(
        output_folder,
        exist_ok=True
    )

    output_file = os.path.join(
        output_folder,
        "ASS1_Final_Annotated_Variants.csv"
    )

    final_df.to_csv(
        output_file,
        index=False
    )

    return output_file
