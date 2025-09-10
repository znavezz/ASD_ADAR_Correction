import gzip
import io
import os, tempfile
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
    with opener(vep_results_path, 'rt', encoding="utf-8",errors="replace") as f:
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


def vep_annotations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Run VEP via the shell wrapper and merge the annotations into `df` in-place.
    """
    # ── 1.  dump the DataFrame rows we want VEP to see ──────────────────────
    cur_dir  = os.path.dirname(os.path.abspath(__file__))
    # unique TEMP **input** file
    with tempfile.NamedTemporaryFile(
            dir=cur_dir, suffix=".tsv", delete=False, mode="w") as fh:
        df.to_csv(fh.name, sep="\t", index=False, header=False)
        tmp_path = fh.name

    # unique **output** file inside ./vep_results
    out_dir = os.path.join(cur_dir, "vep_results")
    os.makedirs(out_dir, exist_ok=True)
    out_fh  = tempfile.NamedTemporaryFile(dir=out_dir,
                                          prefix="vep_",
                                          suffix=".txt",
                                          delete=False)
    out_fh.close()
    out_path = out_fh.name

    vep_script = os.path.join(cur_dir, "vep_ann.sh")
    if not os.path.exists(vep_script):
        raise FileNotFoundError(f"VEP script not found: {vep_script}")

    # ── 2.  run the wrapper and capture ONLY its stdout (the results path) ──
    try:
        completed = subprocess.run(
            ["bash", vep_script, tmp_path, out_path],
            check=True, capture_output=True, text=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"VEP wrapper failed (exit {e.returncode}).\nSTDERR:\n{e.stderr}"
        ) from None
    finally:
        if os.path.isfile(tmp_path):
            # clean up temp input file if it exists
            # (it may not if the script failed before writing it)
            print(f"Removing temporary input file: {tmp_path}")
            # this is a temp file, so we can safely remove it
            os.remove(tmp_path)          # clean up temp input no matter what

    vep_results_path = completed.stdout.strip()   # wrapper echoes path
    if not os.path.isfile(vep_results_path):
        raise RuntimeError(
            f"VEP finished but results file not found: {vep_results_path}"
        )

    print(f"✓ VEP results: {vep_results_path}")

    # ── 3.  parse VEP output ------------------------------------------------
    vep_df = upload_vep_results_file(vep_results_path)

    # ── 4.  NEW: split the `Extra` column into all its individual tags ──────
    extras_wide = (
        vep_df["Extra"]
        .str.split(';').explode()
        .str.split('=', n=1, expand=True)
        .rename(columns={0: "key", 1: "value"})
        .reset_index()                               # <─ keeps original row-index in a column
        .pivot_table(
            index="index",                           # group by that saved index
            columns="key",
            values="value",
            aggfunc="first"
        )
    )

    vep_df = vep_df.drop(columns=["Extra"]).join(extras_wide)

    # ── 5.  merge VEP data back into the caller's DataFrame ─────────────────
    KEYS = ["chr", "pos", "ref", "alt"]

    # ensure the join keys are clean strings on both sides
    for c in KEYS:
        df[c]     = df[c].astype(str).str.strip()
        vep_df[c] = vep_df[c].astype(str).str.strip()

    # merged = df.merge(vep_df, on=KEYS, how="inner", sort=False)

    # # in-place replace df with the merged content
    # df.drop(df.index, inplace=True)
    # df[merged.columns] = merged
    merged = df.merge(vep_df, on=KEYS, how="inner", sort=False)

    # hand back a fresh, consolidated frame ─ no fragmentation
    return merged.copy()

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

