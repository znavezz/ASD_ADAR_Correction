import gzip
import io
import os
import pandas as pd
import subprocess

def upload_vep_results_file(vep_results_path: str) -> pd.DataFrame:
    """
    Function to upload VEP results file and parse it into a DataFrame.
    Parameters:
        vep_results_path (str): Path to the VEP results file (gzipped or plain text).
    Returns:
        pd.DataFrame: DataFrame containing the parsed VEP results.
    """

    # Column names for VEP output
    col_names = ["#Uploaded_variation", "Location", "Allele", "Gene", "Feature",
                 "Feature_type", "Consequence", "cDNA_position", "CDS_position",
                 "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra"]

    # Choose opener based on file extension
    if vep_results_path.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    # Read and filter out header lines
    with opener(vep_results_path, 'rt') as f:
        filtered = [l for l in f if not l.startswith('##')]

    vep_df = pd.read_csv(
        io.StringIO(''.join(filtered)),
        sep='\t',
        header=0,
        names=col_names
    )

    # Split the "#Uploaded_variation" column into chr / pos / ref_alt
    tmp = vep_df["#Uploaded_variation"].str.split(":", expand=True)
    tmp.columns = ["chr", "pos", "ref", "alt"]


    # Assemble parsed fields
    vep_df = pd.concat([
        vep_df.drop(columns=["#Uploaded_variation"]),
        tmp[["chr", "pos", "ref", "alt"]],
    ], axis=1)

    # Prefix "chr" if not already present
    # vep_df["chr"] = vep_df["chr"].astype(str).apply(lambda x: x if x.lower().startswith('chr') else 'chr' + x)

    return vep_df


def get_from_extra(key: str, extras: str) -> str:
    """
    Given a key and a string of extras, retrieve the value
    associated with the key from the extras string.
    """
    extras_dict = dict(item.split("=", 1) for item in extras.split(";") if "=" in item)
    return extras_dict.get(key, pd.NA)


def vep_annotations(vep_results_path: str, df: pd.DataFrame) -> None:
    """
    Annotate `df` in‑place by adding columns from the VEP output, 
    merging on ["chr","pos","ref","alt"].  Only rows present in both
    will be kept.

    Parameters:
        vep_results_path (str): Path to the VEP file.
        df (pd.DataFrame): DataFrame to annotate (modified in place).

    Returns:
        None
    """
    # 1) Parse the VEP file exactly as before
    vep_df = upload_vep_results_file(vep_results_path)

    # 2) Pull out the extras into their own columns
    keys_from_extras = [
        "STRAND", "VARIANT_CLASS", "SYMBOL", "SYMBOL_SOURCE",
        "SIFT", "PHENO", "PolyPhen", "HGVSc", "HGVSp", "PhastCons46", "SWISSPROT", "UNIPARC", "EXON", "IMPACT"
    ]
    for key in keys_from_extras:
        vep_df[key] = vep_df["Extra"].apply(lambda x: get_from_extra(key, x))
    vep_df.drop(columns=["Extra"], inplace=True)

    # 3) Ensure join‐keys are strings in both frames
    # normalize types & whitespace
    for col in ["chr","pos","ref","alt"]:
        df[col]   = df[col].astype(str).str.strip()
        vep_df[col]= vep_df[col].astype(str).str.strip()

    # find which input keys have no match in VEP
    in_keys  = set(df.set_index(["chr","pos","ref","alt"]).index)
    out_keys = set(vep_df.set_index(["chr","pos","ref","alt"]).index)

    missing = in_keys - out_keys
    print(f"{len(missing)} variants in your input never appeared in VEP output.")
    print(list(missing)[:10])

    # Check duplicates in your original df
    dups_in_df = df.duplicated(subset=["chr","pos","ref","alt"]).sum()
    print(f"Duplicates in input df: {dups_in_df}")

    # Check duplicates in your VEP output
    dups_in_vep = vep_df.duplicated(subset=["chr","pos","ref","alt"]).sum()
    print(f"Duplicates in VEP df: {dups_in_vep}")



    # 4) Inner merge to keep only matching rows
    merged = df.merge(
        vep_df,
        on=["chr","pos","ref","alt"],
        how="inner",
        sort=False
    )
    print(f"VEP annotations: {len(merged)} rows after merge")
    print(f"VEP annotations: {len(vep_df)} rows in VEP file")
    print(f"VEP annotations: {len(df)} rows in input DataFrame")
    print(f"VEP annotations: {len(merged) / len(df):.2%} of input rows matched")
    print(f"VEP annotations: {len(merged) / len(vep_df):.2%} of VEP rows matched")

    # 5) In‑place mutation: clear out old df, then refill from merged
    df.drop(df.index, inplace=True)
    for col in merged.columns:
        df[col] = merged[col].values

    # now `df` contains exactly the inner‑merged table



def isADARFixable(df: pd.DataFrame) -> pd.Series:
    """
    Function to determine if a variant is ADAR fixable.
    """
    return (df["ref"].str.upper() == "G") & (df["alt"].str.upper() == "A")


def isApoBecFixable(df: pd.DataFrame) -> pd.Series:
    """
    Function to determine if a variant is ApoBec fixable.
    """
    return (df["ref"].str.upper() == "T") & (df["alt"].str.upper() == "C")

def add_hg19_column(df: pd.DataFrame,
                    genome_file: str = "/home/alu/aluguest/Nave_Oded_Project/resources/hg19.fa"
                   ) -> pd.DataFrame:
    """
    For each row in the DataFrame, extract the hg19 reference base at 'pos'
    (1-based) using bedtools getfasta, and add the result as a new column 'hg19'.

    Assumptions:
      - df has columns: 'chr', 'pos', 'ref', 'alt'
      - genome_file is the path to an indexed hg19 FASTA (*.fa + .fai)
      - bedtools is installed and in your PATH.

    Returns:
      A new DataFrame with an added column 'hg19' containing the extracted base.
    """
    hg19_sequences = []
    for idx, row in df.iterrows():
        chrom = row["chr"]
        pos   = int(row["pos"])
        # bedtools uses 0-based, half-open intervals:
        bed_start = pos - 1
        bed_stop  = pos      # gives you exactly one base

        bed_line = f"{chrom}\t{bed_start}\t{bed_stop}\n"
        try:
            result = subprocess.run(
                ["bedtools", "getfasta",
                 "-fi", genome_file,
                 "-bed", "-", "-fo", "-"],
                input=bed_line.encode("utf-8"),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True
            )
            fasta_output = result.stdout.decode("utf-8")
            # skip the '>' header line, join and uppercase the rest
            seq = "".join(fasta_output.splitlines()[1:]).upper()
        except subprocess.CalledProcessError as e:
            seq = None
            print(f"Error fetching {chrom}:{bed_start}-{bed_stop}:", 
                  e.stderr.decode("utf-8"))

        hg19_sequences.append(seq)

    df["hg19"] = hg19_sequences
    return df


if __name__ == "__main__":
    # Example usage
    current_dir = os.path.dirname(os.path.abspath(__file__))
    vep_results_path = os.path.join(current_dir, "vep_test.txt")

    df = pd.DataFrame({
        "chr": ["1", "1"],
        "pos": ["880107", "880850"],
        "ref": ["C", "G"],
        "alt": ["A", "A"]
    })

    annotated_df = vep_annotations(vep_results_path, df)
    print(annotated_df)
