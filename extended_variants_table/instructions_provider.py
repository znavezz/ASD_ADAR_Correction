import pandas as pd
import importlib
import os
import sys
from db import Db, VariantsDb, ValidationDb


class InstructionsProvider:
    def __init__(self, dbs_path: str):
        # Add the base DBs directory to sys.path
        sys.path.append(dbs_path)
        self.dbs_path, self.variants_dbs_path, self.validation_dbs_path = check_dbs_dir_structure(dbs_path)
        self.variants_instructions_map = self.get_db_instructions("default")
        self.validation_instructions_map = self.get_db_instructions("default", "validation")
        
    def get_variants_instructions_map(self):
        return self.variants_instructions_map
    
    def set_variants_instructions_map(self, instructions):
        self.variants_instructions_map = instructions
    
    def set_validation_instructions_map(self, instructions):
        self.validation_instructions_map = instructions
    
    def get_db_instructions(self, db_name: str, db_type: str = "variants"):
        """
        Loads the instructions from the specified database.
        Returns the instructions as a dictionary.
        """
        try:
            # Construct the path to the instructions module.
            instructions_path = os.path.join(self.dbs_path, db_type, db_name, "instructions.py")
            if not os.path.isfile(instructions_path):
                raise FileNotFoundError(f"Instructions file not found at {instructions_path}.")
            
            # Compute the module name relative to the base dbs_path.
            rel_path = os.path.relpath(instructions_path, self.dbs_path)
            module_name = rel_path.replace(os.sep, ".").replace(".py", "")
            
            instructions_module = importlib.import_module(module_name)
            # Get the instructions from the module.
            instructions = getattr(instructions_module, "instructions", None)
            if instructions is None:
                raise ImportError(f"Instructions not found in module '{module_name}'.")
            # Check if the instructions are a dictionary.
            if not isinstance(instructions, dict):
                raise TypeError(f"Instructions for '{db_name}' must be a dictionary.")
            
            return instructions
        except ModuleNotFoundError as e:
            raise ImportError(f"Failed to import instructions module: {e}")
        
    def get_final_instructions(self, db_name: str, db_type: str = "variants"):
        db_instructions = self.get_db_instructions(db_name, db_type)
        default_instructions = (self.variants_instructions_map.copy() if db_type == "variants"
                               else self.validation_instructions_map.copy())

        def recursive_merge(default, override):
            merged = default.copy()
            for key, value in override.items():
                # If both the default value and the override value are dictionaries,
                # then recursively merge them.
                if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                    merged[key] = recursive_merge(merged[key], value)
                else:
                    # Otherwise, override the default with the db-specific value.
                    merged[key] = value
            return merged

        merged_instructions = recursive_merge(default_instructions, db_instructions)
        return merged_instructions


    def get_annotations_names(self):
        return list(self.variants_instructions_map["annotations"].keys())
    
    def get_dbs_names(self, db_type: str = "variants"):
        """
        Returns a list of database names in the specified directory.
        """
        dbs_path = self.variants_dbs_path if db_type == "variants" else self.validation_dbs_path
        if not os.path.isdir(dbs_path):
            raise NotADirectoryError(f"Path '{dbs_path}' is not a directory.")
        
        dbs_names = [name for name in os.listdir(dbs_path)
                     if os.path.isdir(os.path.join(dbs_path, name)) and name != "default" and not name.startswith("__")]
        return dbs_names
    
    def get_key_columns(self):
        """
        Returns the key columns for the database.
        """
        key_cols = self.variants_instructions_map.get("key_cols")
        if key_cols is None:
            raise ValueError("Key columns not found in the instructions.")
        return key_cols

    def create_db_instance(self, db_name: str, db_type: str = "variants"):
        """
        Creates a Db instance for the specified database.
        """
        try:
            if db_type == "variants":
                dbs_path = self.variants_dbs_path
            elif db_type == "validation":
                dbs_path = self.validation_dbs_path
            else:
                raise ValueError(f"Invalid db_type '{db_type}'. Must be 'variants' or 'validation'.")
            # Construct the path to the database directory
            db_dir = os.path.join(dbs_path, db_name)
            # Look for a file named variants_table with any extension
            potential_files = [f for f in os.listdir(db_dir) if f.startswith("variants_table")]
            if not potential_files:
                raise FileNotFoundError(f"Database file not found in {db_dir}. Expected a file starting with 'variants_table'.")
            
            # Use the first matching file
            db_path = os.path.join(db_dir, potential_files[0])
            if not os.path.isfile(db_path):
                raise FileNotFoundError(f"Database file not found at {db_path}.")
            
            # Load the instructions for the specified db.
            instructions_map = self.get_final_instructions(db_name, db_type)
            
            # Create a Db instance.
            db_instance = (VariantsDb(
                db_path=db_path,
                instructions=instructions_map,
            ) if db_type == "variants"
            else ValidationDb(
                db_path=db_path,
                instructions=instructions_map,
            ))

            return db_instance
        except Exception as e:
            print(f"Failed to create database instance for '{db_name}': {e}")
            return None



if __name__ == "__main__":
    dbs_path = "/Users/znave/Desktop/DBs"
    instructions_provider = InstructionsProvider(dbs_path)
    # Example usage
    print(instructions_provider.get_variants_instructions_map())
    
    # Example DataFrame
    data = pd.DataFrame({
        "ref": ["G", "A", "T", "T"],
        "alt": ["A", "C", "G", "C"]
    })
    instructions_provider.set_variants_db_instructions("default")
    
    # Apply pre-processing
    processed_data = instructions_provider.variants_instructions_map["pre_processor"](data)

    annotations = instructions_provider.variants_instructions_map["annotations"]
    print(f"Instructions map: {instructions_provider.get_variants_instructions_map()}")
    # Apply annotation functions
    annotations["isADARFixable"]["compute_function"](processed_data)
    annotations["isApoBecFixable"]["compute_function"](processed_data)
    
    print(processed_data)


def check_dbs_dir_structure(dbs_path: str):
    """
    Checks if the dbs_path contains the expected directory structure.
    The expected structure is:
    - dbs_path
        - variants
            - default
                - instructions.py
            - <other_db>
                - instructions.py
        - validation (optional, but if exists must follow this format)
            - default
                - instructions.py
            - <other_db>
                - instructions.py
    """
    if not os.path.isdir(dbs_path):
        raise NotADirectoryError(f"Path '{dbs_path}' is not a directory.")
    
    variants_path = os.path.join(dbs_path, "variants")
    if not os.path.isdir(variants_path):
        raise NotADirectoryError(f"Path '{variants_path}' is not a directory.")
    default_instructions_path = os.path.join(variants_path, "default", "instructions.py")
    if not os.path.isfile(default_instructions_path):
        raise FileNotFoundError(f"Default instructions file not found at {default_instructions_path}.")
            
    validation_path = os.path.join(dbs_path, "validation")
    if os.path.isdir(validation_path):
        default_validation_instructions_path = os.path.join(validation_path, "default", "instructions.py")
        if not os.path.isfile(default_validation_instructions_path):
            raise FileNotFoundError(f"Default validation instructions file not found at {default_validation_instructions_path}.")
    else:
        validation_path = None
    return (dbs_path, variants_path, validation_path)
