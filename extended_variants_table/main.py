from instructions_provider import InstructionsProvider
from extended_table import ExtendedTable
from db import Db, VariantsDb, ValidationDb
import os
import sys


def get_dbs_dir():
    """
    Prompts the user to enter the path to the databases directory.
    Returns the path as a string.
    """
    usr_input = input("Please enter the path to the databases directory: ")
    if not os.path.isdir(usr_input):
        print(f"Directory \"{usr_input}\" does not exist. Please check the path and try again.")
        sys.exit(1)
    return usr_input



def main():
    dbs_dir = input("Please enter the path to the databases directory: ")
    while (not os.path.isdir(dbs_dir)):
        print(f"Directory \"{dbs_dir}\" does not exist. Please check the path and try again or press 'X' to exit.")
        usr_input = input()
        if usr_input.lower() == 'x':
            sys.exit(0)
        dbs_dir = usr_input
    sys.path.append(dbs_dir)

    instructions_provider = InstructionsProvider(dbs_dir)




    # itereate over the directories in the dbs_dir
    variants_dbs = instructions_provider.get_dbs_names()
    print(f"Databases found in the directory: {variants_dbs}")
    validation_dbs = instructions_provider.get_dbs_names("validation")

    print(f"Validation databases found in the directory: {validation_dbs}")

    annotations_funcs = instructions_provider.get_annotations_names()
    annotations_cols = annotations_funcs.copy()
    print(f"Available annotation funcs: {annotations_funcs}")
    # vep_cols = ["Gene","Feature","Feature_type","Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","STRAND", "VARIANT_CLASS", "SYMBOL", "SYMBOL_SOURCE", "SIFT", "PHENO", "PolyPhen", "HGVSc", "HGVSp", "PhastCons100",  "SWISSPROT", "UNIPARC", "EXON", "IMPACT"]
    
    # if "vep_annotations" in annotations_funcs:
    #     annotations_cols.remove("vep_annotations")
    #     annotations_cols.extend(vep_cols)
    print(f"Annotation columns: {annotations_cols}")
    key_cols = instructions_provider.get_key_columns()
    print(f"Key columns for the new table: {key_cols}")
    # Create variant database instances using map
    variant_db_instances = list(map(lambda db_name: instructions_provider.create_db_instance(db_name), variants_dbs))
    print(f"Created {len(variant_db_instances)} variant database instances.")

    # Create validation database instances using map
    validation_db_instances = list(map(lambda db_name: instructions_provider.create_db_instance(db_name, "validation"), validation_dbs))
    print(f"Created {len(validation_db_instances)} validation database instances.")
    # Print details for each validation database instance



    # Create an ExtendedTable instance
    extended_table = ExtendedTable(key_cols=key_cols, instructions_provider=instructions_provider, variant_dbs=variant_db_instances,ann_funcs=annotations_funcs)
   
    # extended_table = ExtendedTable(key_cols=key_cols, instructions_provider=instructions_provider, variant_dbs=variant_db_instances,validation_dbs=validation_db_instances,ann_funcs=annotations_funcs)


    extended_table.merge_all_dbs()
    extended_table.validate_table()
    print("All databases merged into the extended table.")
    print(f"Extended table: \n{extended_table.table}")

    usr_input = input("Do you want to save the table to a file? (y/n): ")
    if usr_input.lower() == 'y':
        # Save the file in the current pwd
        file_path = os.path.join(os.getcwd(), "Results/processed_variants_table.tsv")
        extended_table.save_table(file_path, "tsv")
    else:
        print("Table not saved.")

if __name__ == "__main__":
    main()