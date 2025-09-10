import pandas as pd
from db import Db, VariantsDb, ValidationDb
from instructions_provider import InstructionsProvider


class ExtendedTable:
    def __init__(self ,key_cols: list[str], instructions_provider: InstructionsProvider, variant_dbs: list[VariantsDb] | None = None, validation_dbs: list[ValidationDb] | None = None, ann_funcs: list[str] | None = None) -> None:
        self.table: pd.DataFrame = pd.DataFrame()
        self.key_cols: list[str] = key_cols
        self.instructions_provider: InstructionsProvider = instructions_provider
        self.variant_dbs: list[VariantsDb] = variant_dbs if variant_dbs is not None else []
        self.validation_dbs: list[ValidationDb] = validation_dbs if validation_dbs is not None else []
        self.ann_funcs: list[str] = ann_funcs if ann_funcs is not None else []

    def upload_table(self, file_path: str) -> None:
        """
        Uploads a table from a file path.
        """
        self.table = pd.read_csv(file_path)
        print(f"Table uploaded from {file_path}.")
    
    def _order_with_keys_first(self, cols: pd.Index | list[str]) -> list[str]:
        """
        Return `cols` as a list with the key columns leading and everything
        else afterwards (in the order they originally appeared).

        Assumes every key column is present in `cols`.
        """
        cols = list(cols)                        # keep order as given
        return self.key_cols + [c for c in cols if c not in self.key_cols]


    def create_basic_table(self) -> None:
        """
        Initialise an empty extended table that has
        - the key columns,
        - one indicator column per variant DB,
        - one indicator column per validation DB.
        Annotation columns will appear automatically later.
        """
        variants_dbs_names   = [db.name for db in self.variant_dbs]
        validation_dbs_names = [db.name for db in self.validation_dbs]

        base_cols = self.key_cols + variants_dbs_names + validation_dbs_names
        # just in case a key appears twice, drop duplicates while keeping order
        base_cols = list(dict.fromkeys(base_cols))

        self.table = pd.DataFrame(columns=base_cols)
        print(f"Basic table created with columns: {self.table.columns.tolist()}")



    def register_db(self, db: Db) -> None:
        """
        Adds a database to the list of databases.
        """
        if isinstance(db, VariantsDb):
            self.variant_dbs.append(db)
            print(f"Variant database {db.name} added.")
        elif isinstance(db, ValidationDb):
            self.validation_dbs.append(db)
            print(f"Validation database {db.name} added.")
        else:
            raise ValueError("Unsupported database type.")
        
    def merge_db(self, db: Db) -> None:
        """
        Merges a VariantsDb into the extended table by splitting its variants into:
        - Existing variants: those already in the extended table.
        - New variants: those not present in the extended table.
        
        For existing variants, only the indicator column for this db is updated to 1.
        For new variants, annotation columns are computed and new rows are appended.
        """
        if not isinstance(db, VariantsDb):
            raise ValueError("Only VariantsDb can be merged into the extended table.")
        

        if db.df is None:
            db.upload_db()
        
        # Pre-process the db to get a standardized DataFrame.
        db.pre_process()
        
        # Define the indicator column name for this database.
        indicator_col: str = db.name
        
        # Ensure the indicator column exists in the extended table.
        # If the extended table is empty, then all rows from db.df are new.
        if self.table.empty:
            self.create_basic_table()
            
        # Ensure the indicator column exists in the extended table.
        if indicator_col not in self.table.columns:
            self.table[indicator_col] = 0

        # Use a left merge (with indicator) to identify which rows in db.df are already present.
        merge_df = pd.merge(
            db.df[db.key_cols],
            self.table[self.key_cols],
            on=self.key_cols,
            how='left',
            indicator=True
        )
        
        # 'both' indicates rows already in self.table; 'left_only' are new variants.
        existing_df = merge_df[merge_df['_merge'] == 'both']
        new_df = merge_df[merge_df['_merge'] == 'left_only'].drop(columns=['_merge'])
        print(f"Found {len(existing_df)} existing variants and {len(new_df)} new variants in database '{db.name}'.")
        
        # Update the indicator column for existing variants.
        # Build a set of keys (tuples) for the existing variants.
        existing_keys = {tuple(row) for row in existing_df[self.key_cols].to_numpy()}
        # Create a boolean mask for self.table rows whose key is in the existing_keys set.
        mask = self.table[self.key_cols].apply(lambda row: tuple(row) in existing_keys, axis=1)
        self.table.loc[mask, indicator_col] = 1
        
        # For new variants, compute annotation values and set the indicator.
        new_df[indicator_col] = 1

        # for ann_func in self.ann_funcs:
        #     db.instructions["annotations"][ann_func]["compute_function"](new_df)

        for ann_func in self.ann_funcs:
            fn = db.instructions["annotations"][ann_func]["compute_function"]
            result = fn(new_df)
            if result is not None:        # VEP returns a DataFrame
                new_df = result
        # 1) Union of columns from both frames
        all_cols = self.table.columns.union(new_df.columns)

        # 2) Put key columns first
        all_cols = self._order_with_keys_first(all_cols)

        # 3) Reindex both frames to that ordered union
        self.table = self.table.reindex(columns=all_cols)
        new_df     = new_df.reindex(columns=all_cols)

        # now every column is present in both; NaN wherever data is missing
        
        # Append the new variants to the extended table.
        self.table = pd.concat([self.table, new_df], ignore_index=True)
        
        print(f"Database '{db.name}' merged: {len(existing_df)} existing variants updated and {len(new_df)} new variants added.")
        # Clear the DataFrame in the db instance to free up memory.
        db.df = None

    def merge_all_dbs(self) -> None:
        """
        Merges all registered VariantsDbs into the extended table.
        """
        for db in self.variant_dbs:
            self.merge_db(db)
        print("All databases merged into the extended table.")

    

    def save_table(self, file_path: str, file_format:str="csv") -> None:
        """
        Saves the extended table to a file in the specified format.
        """
        if file_format == "csv":
            self.table.to_csv(file_path, index=False)
            print(f"Table saved to {file_path} in CSV format.")
        elif file_format == "xlsx":
            self.table.to_excel(file_path, index=False)
            print(f"Table saved to {file_path} in Excel format.")
        elif file_format == "tsv":
            self.table.to_csv(file_path, sep='\t', index=False)
            print(f"Table saved to {file_path} in TSV format.")
        else:
            raise ValueError("Unsupported file format. Supported formats are: csv, xlsx, tsv.")

    def validate_table(self) -> None:
        """
        Validates the extended table against registered validation databases.
        This method checks for any discrepancies or errors in the data.
        """
        for db in self.validation_dbs:
            db.validate(self.table)
            print(f"Table validated against {db.name}.")

    def get_table(self) -> pd.DataFrame:
        """
        Returns the current state of the extended table.
        """
        return self.table
    def __clear_table(self) -> None:
        """
        Clears the current state of the extended table.
        """
        self.table = pd.DataFrame()
        self.variant_dbs = []
        self.validation_dbs = []
        print("Extended table cleared.")