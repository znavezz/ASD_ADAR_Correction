import pandas as pd
import os
import subprocess

def pre_process(data: pd.DataFrame) -> pd.DataFrame:
    """
    Default pre-processing function for the data.
    """
    return data

def lift_over(df: pd.DataFrame) -> pd.DataFrame:
    """
    Function to perform liftOver from hg19 to hg38.
    This is a placeholder function; actual implementation would depend on the liftOver tool.
    """
    df['start'] = df['pos'].astype(int) - 1  # Convert pos to int and create start column
    df['end'] = df['start'] + df['ref'].str.len()  # Calculate end position based on ref length
    df['chr'] = 'chr' + df['chr'] 
    df = df[['chr', 'start', 'end', 'ref', 'alt']].copy()  # Reorder columns to match expected format
    # 5) Convert 'chr' to string type
    df['chr'] = df['chr'].astype(str)
    # 6) Convert 'start' and 'end' to int type
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)

    # Ensure the DataFrame has the expected columns
    expected_columns = ['chr', 'start', 'end', 'ref', 'alt']
    for col in expected_columns:
        if col not in df.columns:
            raise ValueError(f"Missing expected column: {col}")
    
    tmp_file_path = "tmp_lifted_over_variants.tsv"
    df.to_csv(tmp_file_path, sep='\t', index=False, header=False)  # Save to a temporary file for liftOver

    # run ../../liftOverToHg38.sh
    lift_over_script = os.path.join(os.path.dirname(__file__), '../liftOverToHg38.sh')
    chain_file = os.path.join(os.path.dirname(__file__), '../../.liftOver/hg19ToHg38.over.chain.gz')
    if os.path.exists(lift_over_script):
        if os.path.exists(chain_file):
            try:
                result_path = subprocess.run(['bash', lift_over_script, tmp_file_path, chain_file], check=True, capture_output=True, text=True).stdout.strip()
                if not result_path: 
                    raise ValueError("LiftOver script did not return a valid file path.")
                # run the liftOver script with the temporary file and chain file and save the returned file path

                print("LiftOver completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during LiftOver: {e}")
        else:
            print(f"Chain file not found at {chain_file}")
    else:
        print(f"LiftOver script not found at {lift_over_script}")

    # Load the lifted over variants back into a DataFrame
    lifted_df = pd.read_csv(result_path, sep='\t', header=None, names=expected_columns)
    os.remove(tmp_file_path)  # Clean up the temporary file
    os.remove(result_path)  # Clean up the result file
    # Ensure the lifted DataFrame has the expected columns
    for col in expected_columns:
        if col not in lifted_df.columns:
            raise ValueError(f"Lifted DataFrame is missing expected column: {col}")
        
    print(" DataFrame head:")
    print(df[['chr', 'start', 'end', 'ref', 'alt']].head())  # Debugging output
    print(" Lifted DataFrame head:")
    print(lifted_df[['chr', 'start', 'end', 'ref', 'alt']].head())  # Debugging output
    # Return the lifted DataFrame to the original format
    lifted_df['chr'] = lifted_df['chr'].str.replace('chr', '', regex=False)
    lifted_df['start'] = (lifted_df['start'] + 1).astype(str)  # Convert start back to 1-based index
    # rename 'start' to 'pos'
    lifted_df.rename(columns={'start': 'pos'}, inplace=True)

    return lifted_df[['chr', 'pos', 'ref', 'alt']].reset_index(drop=True)