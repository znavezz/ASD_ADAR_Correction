#!/usr/bin/env python3
"""
2nd_run.py – append hg38 reference sequence + ADAR/APOBEC flags
               to a copy of hg38_extended_table.tsv

The updated table is written as 2nd_run.csv in the same directory.
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


def add_hg38_column(df: pd.DataFrame, fasta_path: str) -> None:
    """Add a new column ‘hg38’ with reference bases from *fasta_path*."""
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

    df["hg38"] = seqs


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
    in_path  = os.path.join(here, "hg38_extended_table.tsv")
    out_path = os.path.join(here, "2nd_run.csv")
    fasta_path = os.path.join(here, "..", "resources", "hg38.fa")  # adjust if needed

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
        add_hg38_column(df, fasta_path)
        dbs_count(df)

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
