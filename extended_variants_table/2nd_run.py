#!/usr/bin/env python3
"""
2nd_run.py – append genome reference sequence + ADAR/APOBEC flags
               to a copy of the input TSV file

The updated table is written as 2nd_run.csv in the Results dir in the same parent directory as this script.
"""

import sys
import os
import math
import pandas as pd
from pyfaidx import Fasta, FetchError

from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial
from typing import Optional
try:
    from tqdm import tqdm
    _HAS_TQDM = True
except Exception:
    _HAS_TQDM = False


# ────────────────────────── helpers ─────────────────────────────────────────
def _print(msg: str) -> None:
    print(msg, flush=True)

def check_columns(df: pd.DataFrame, columns: list[str]) -> None:
    """Raise ValueError if any required column is missing."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")

# thread-local Fasta store to keep each worker safe
import threading
_thread_local = threading.local()

def _get_fasta_for_thread(fasta_path: str) -> Fasta:
    fa: Optional[Fasta] = getattr(_thread_local, "fasta", None)
    if fa is None or getattr(_thread_local, "fasta_path", None) != fasta_path:
        # Open read-only; pyfaidx handles its own indexing if .fai exists
        _thread_local.fasta = Fasta(fasta_path)
        _thread_local.fasta_path = fasta_path
    return _thread_local.fasta

def _lookup_one(fasta_path: str, genome_version: str, chrom: str, pos: int, ref: str) -> str:
    fa = _get_fasta_for_thread(fasta_path)
    chrom = "chr" + str(chrom).lstrip("chr")
    seq = fa[chrom][pos - 1 : pos - 1 + len(ref)].seq.upper()
    return seq


def add_genome_ref_column(
    df: pd.DataFrame,
    fasta_path: str,
    genome_version: str = "hg38",
    workers: Optional[int] = None,
    chunk_size: int = 10_000
) -> None:
    """
    Create column hg38 or hg19 in *df* by looking up reference sequence in FASTA.
    Parallelized with a thread pool and per-thread Fasta handles.
    """
    if genome_version not in ["hg19", "hg38"]:
        raise ValueError(f"Invalid genome_version: {genome_version}.\n Expected 'hg19' or 'hg38'.")
    check_columns(df, ["chr", "pos", "ref"])

    # Validate FASTA first (open once in main thread)
    try:
        _ = Fasta(fasta_path)
    except Exception as e:
        raise RuntimeError(f"Could not open FASTA: {e}") from None

    n = len(df)
    _print(f"→ FASTA lookups: {n:,} rows | genome={genome_version} | fasta={os.path.basename(fasta_path)}")
    workers = workers or os.cpu_count() or 4
    _print(f"→ Using {workers} worker threads")

    # Prepare inputs
    chroms = df["chr"].astype(str).to_numpy()
    poss   = df["pos"].astype(int).to_numpy()
    refs   = df["ref"].astype(str).to_numpy()

    # Work in chunks to keep memory bounded and keep progress meaningful
    results = [None] * n
    num_chunks = math.ceil(n / chunk_size)
    iterator = range(num_chunks)
    if _HAS_TQDM:
        iterator = tqdm(iterator, desc="FASTA lookups (chunks)", unit="chunk")

    for ci in iterator:
        start = ci * chunk_size
        end = min(start + chunk_size, n)
        batch_len = end - start

        # inner progress for each chunk (optional)
        if _HAS_TQDM:
            pbar = tqdm(total=batch_len, leave=False, unit="row", desc=f"Chunk {ci+1}/{num_chunks}")
        else:
            pbar = None

        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = []
            for i in range(start, end):
                fut = ex.submit(_lookup_one, fasta_path, genome_version, chroms[i], poss[i], refs[i])
                futs.append((i, fut))

            for i, fut in futs:
                try:
                    seq = fut.result()
                except (KeyError, FetchError) as e:
                    raise RuntimeError(f"FASTA lookup failed for {chroms[i]}:{poss[i]}: {e}") from None
                results[i] = seq
                if pbar: pbar.update(1)

        if pbar: pbar.close()

    df[genome_version] = results  # type: ignore[list-item]


def isADARFixable(df: pd.DataFrame) -> None:
    """Create boolean column is_ADAR_fixable in *df*."""
    check_columns(df, ["STRAND", "ref", "alt"])
    _print("→ Computing is_ADAR_fixable …")
    df["is_ADAR_fixable"] = (
        ((df["STRAND"] == "1")  & (df["ref"] == "G") & (df["alt"] == "A")) |
        ((df["STRAND"] == "-1") & (df["ref"] == "C") & (df["alt"] == "T"))
    )

def isAPOBECFixable(df: pd.DataFrame) -> None:
    """Create boolean column is_APOBEC_fixable in *df*."""
    check_columns(df, ["STRAND", "ref", "alt"])
    _print("→ Computing is_APOBEC_fixable …")
    df["is_APOBEC_fixable"] = (
        ((df["STRAND"] == "1")  & (df["ref"] == "T") & (df["alt"] == "C")) |
        ((df["STRAND"] == "-1") & (df["ref"] == "A") & (df["alt"] == "G"))
    )

def dbs_count(df: pd.DataFrame) -> None:
    """Create a column dbs_count in *df*."""
    dbs = ["db1_TableS1", "db1_TableS3", "db1_TableS4", "db2_TableS2", "db3", "db4_SD1", "db4_SD2", "db4_SD3", "db5_SD1", "Varicarta"]
    check_columns(df, dbs)
    _print("→ Counting database hits (dbs_count) …")
    # tqdm for apply progress (coarse)
    if _HAS_TQDM:
        _print("   (this step is vectorized; progress is fast)")
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
    _print(f"Using genome version: {genome_version}")
    # i want to trim the input file name to create the output file name
    base_name = os.path.basename(in_path)
    name, _ = os.path.splitext(base_name)
    out_path = os.path.join(here, "Results", f"{name}_2nd_run.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fasta_path = os.path.join(here, "..", "resources", ("hg19.fa" if genome_version == "hg19" else "hg38.fa"))  # adjust if needed

    # ── 1. read original table ────────────────────────────────────────────
    _print(f"→ Reading table: {in_path}")
    try:
        original = pd.read_csv(in_path, sep="\t", low_memory=False, dtype=str)
    except FileNotFoundError:
        sys.exit(f"❌  Table not found: {in_path}")
    except Exception as e:
        sys.exit(f"❌  Could not read TSV: {e}")
    _print(f"✓ Read {len(original):,} rows and {len(original.columns)} columns")

    df = original.copy()  # work on a separate copy

    # ── 2. add new columns ────────────────────────────────────────────────
    _print("→ Starting annotation steps …")
    try:
        isADARFixable(df)
        isAPOBECFixable(df)
        add_genome_ref_column(df, fasta_path, genome_version, workers=100, chunk_size=20_000)
        # dbs_count(df)
    except Exception as e:
        sys.exit(f"❌  Processing failed: {e}")

    # ── 3. save as CSV ────────────────────────────────────────────────────
    _print(f"→ Writing CSV to: {out_path}")
    try:
        df.to_csv(out_path, index=False)
    except Exception as e:
        sys.exit(f"❌  Could not write {out_path}: {e}")

    _print("✓ Done.")
    _print(f"✓ Wrote updated table to {out_path}")
    _print("Preview:")
    _print(df.head().to_string())


if __name__ == "__main__":
    main()
