#!/usr/bin/env python3
"""
2nd_run.py – append genome reference sequence + ADAR/APOBEC flags
               to a copy of the input TSV file

The updated table is written as 2nd_run.csv in the Results dir in the same parent directory as this script.
"""

import sys
import os
import pandas as pd
from pyfaidx import Fasta, FetchError


# ────────────────────────── helpers ─────────────────────────────────────────
def check_columns(df: pd.DataFrame, columns: list[str]) -> None:
    """Raise ValueError if any required column is missing."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")


def add_genome_ref_column(df: pd.DataFrame, fasta_path: str, genome_version: str="hg38") -> None:
    """Create column hg38 or hg19 in *df* by looking up reference sequence in FASTA."""
    if genome_version not in ["hg19", "hg38"]:
        raise ValueError(f"Invalid genome_version: {genome_version}.\n Expected 'hg19' or 'hg38'.")
    check_columns(df, ["chr", "pos", "ref"])

    try:
        fasta = Fasta(fasta_path)
    except Exception as e:
        raise RuntimeError(f"Could not open FASTA: {e}") from None

    seqs: list[str] = []
    for chrom, pos, ref in zip(df["chr"], df["pos"].astype(int), df["ref"]):
        chrom = "chr" + str(chrom).lstrip("chr")
        try:
            seq = fasta[chrom][pos - 1 : pos - 1 + len(ref)].seq.upper()
        except (KeyError, FetchError) as e:
            raise RuntimeError(f"FASTA lookup failed for {chrom}:{pos}: {e}") from None
        seqs.append(seq)

    df[genome_version] = seqs



def isADARFixable(df: pd.DataFrame) -> None:
    """Create boolean column is_ADAR_fixable in *df*."""
    check_columns(df, ["STRAND", "ref", "alt"])
    df["is_ADAR_fixable"] = (
        ((df["STRAND"] == "1")  & (df["ref"] == "G") & (df["alt"] == "A")) |
        ((df["STRAND"] == "-1") & (df["ref"] == "C") & (df["alt"] == "T"))
    )


def isAPOBECFixable(df: pd.DataFrame) -> None:
    """Create boolean column is_APOBEC_fixable in *df*."""
    check_columns(df, ["STRAND", "ref", "alt"])
    df["is_APOBEC_fixable"] = (
        ((df["STRAND"] == "1")  & (df["ref"] == "T") & (df["alt"] == "C")) |
        ((df["STRAND"] == "-1") & (df["ref"] == "A") & (df["alt"] == "G"))
    )

def dbs_count(df: pd.DataFrame) -> None:
    """Create a column dbs_count in *df*."""
    dbs = ["db1_TableS1", "db1_TableS3", "db1_TableS4", "db2_TableS2", "db3", "db4_SD1", "db4_SD2", "db4_SD3", "db5_SD1", "Varicarta"]
    check_columns(df, dbs)
    df["dbs_count"] = df[dbs].apply(lambda x: x.notna().sum(), axis=1)

# ────────────────────────── main ────────────────────────────────────────────
def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    if not (len(sys.argv) == 2 or (len(sys.argv) == 3 and sys.argv[2] in ["hg19", "hg38"])):
        sys.exit(f"Usage: {sys.argv[0]} [input.tsv [hg19|hg38]]\n"
                 f"  input.tsv: full path to input file \n"
                 f"  hg19|hg38: reference genome version (default: hg38)")
    in_path  = os.path.abspath(sys.argv[1])
    if len(sys.argv) == 3:
        genome_version = sys.argv[2]
    else:
        genome_version = "hg38"
    print(f"Using genome version: {genome_version}")
    out_path = os.path.join(here, "Results/2nd_run.csv")
    fasta_path = os.path.join(here, "..", "resources", f"{"hg19.fa" if genome_version == "hg19" else "hg38.fa"}")  # adjust if needed

    # ── 1. read original table ────────────────────────────────────────────
    try:
        original = pd.read_csv(in_path, sep="\t", low_memory=False, dtype=str)
    except FileNotFoundError:
        sys.exit(f"❌  Table not found: {in_path}")
    except Exception as e:
        sys.exit(f"❌  Could not read TSV: {e}")

    df = original.copy()  # work on a separate copy

    # ── 2. add new columns ────────────────────────────────────────────────
    try:
        isADARFixable(df)
        isAPOBECFixable(df)
        add_genome_ref_column(df, fasta_path, genome_version)
        # dbs_count(df)

    except Exception as e:
        sys.exit(f"❌  Processing failed: {e}")

    # ── 3. save as CSV ────────────────────────────────────────────────────
    try:
        df.to_csv(out_path, index=False)
    except Exception as e:
        sys.exit(f"❌  Could not write {out_path}: {e}")

    print(f"✓ Wrote updated table to {out_path}")
    print(df.head())


if __name__ == "__main__":
    main()
