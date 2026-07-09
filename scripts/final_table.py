def load_vep_annotation(

    vep_file

):

vep_df = pd.read_excel(

    vep_file

)

return vep_df

#Function 2
def prepare_vep_annotation(

    vep_df

):

return vep_df


#Function 3
def merge_variant_tables(

    hgvs_df,

    genotype_df,

    vep_df

):

return final_df


#Function 4
def reorder_columns(

    final_df

):

desired_columns

sort_values()

return final_df


#Function 5
def build_final_variant_table(

    hgvs_df,

    genotype_df,

    vep_file

):

vep_df = load_vep_annotation(

    vep_file

)

vep_df = prepare_vep_annotation(

    vep_df

)

final_df = merge_variant_tables(

    hgvs_df,

    genotype_df,

    vep_df

)

final_df = reorder_columns(

    final_df

)

return final_df

#Save function
def save_final_table(

    final_df,

    output_folder

):

os.makedirs(

    output_folder,

    exist_ok=True

)

final_df.to_csv(

    os.path.join(

        output_folder,

        "ASS1_Final_Annotated_Variants.csv"

    ),

    index=False

)
