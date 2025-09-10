import gzip
import io
import pandas as pd
from .pre_process import pre_process

def validate(db, df) -> pd.DataFrame:
    """
    Validates the DataFrame against the default database.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing variant data.
        
    Returns:
        pd.DataFrame: DataFrame with additional columns for validation results.
    """
    # Check if the DataFrame contains the required columns
    required_cols = db.key_cols
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
        
    
    
    
    # Create a validation column initialized to 0
    # Use merge between db.db and df to check for matches
    df[db.name] = 0
    # Merge with the default database
    # Merge the DataFrame with the database DataFrame
    merged_df = df.merge(db.df, on=db.key_cols, how='left', indicator=True)
    # Update the validation column based on the merge result
    df.loc[merged_df['_merge'] == 'both', db.name] = 1
    # Drop the merge indicator column
    merged_df.drop(columns=['_merge'], inplace=True)
    # Return the updated DataFrame with validation results
    return df


def upload_vcf(db_path):
    """
    Upload the database from the specified path.
    
    Parameters:
        db_path (str): Path to the database file.
        
    Returns:
        pd.DataFrame: DataFrame containing the uploaded data.
    """
    # Original VCF header columns (note the leading '#' for the first column)
    col_names = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    # Modified names: rename "#CHROM" to "chr" and lowercase the rest
    modified_names = ["chr"] + [col.lower() for col in col_names[1:]]
    
    # Open the gzipped file and filter out lines that start with "##"
    with gzip.open(db_path, 'rt') as f:
        filtered_lines = [line for line in f if not line.startswith('##')]
    
    # Create a file-like object from the filtered lines
    filtered_data = io.StringIO("".join(filtered_lines))
    
    # Read the CSV without the usecols argument
    df = pd.read_csv(filtered_data,
                     sep='\t',
                     header=0,
                     names=modified_names,
                     )
    
    return df

def run_validation(db, df):
    """
    Run the validation process on the DataFrame.
    
    Parameters:
        db (object): Database object containing validation information.
        df (pd.DataFrame): DataFrame containing variant data.
        
    Returns:
        pd.DataFrame: DataFrame with validation results.
    """
    db.upload_db()
    # Pre-process the DataFrame
    db.pre_process()
    # Validate the DataFrame against the database  
    validated_df = validate(db, df)
    
    return validated_df
instructions = {
    "key_cols": ["chr", "pos", "ref", "alt"],
    "description": "Default database for validation",
    "validator": validate,
    "pre_processor": pre_process,
    "upload_function": upload_vcf,
    "validator": run_validation,
}
