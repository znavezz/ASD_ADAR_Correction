#!/usr/bin/env bash
set -euo pipefail

# -------- arguments ----------------------------------------------------------
in="$1"           # input BED/TSV/VCF-like file (hg19)
chain="$2"        # hg19ToHg38.over.chain.gz
out="${3:-${in%.*}_hg38.vcf}"   # optional 3rd arg, else auto-name
unmapped="${out%.*}_unmapped.txt"
# ---------------------------------------------------------------------------

# -------- sanity checks (all to stderr) -------------------------------------
if [[ ! -f "$in"    ]]; then echo >&2 "❌ File not found: $in";    exit 1; fi
if [[ ! -f "$chain" ]]; then echo >&2 "❌ Chain file not found: $chain"; exit 1; fi
if ! command -v liftOver &>/dev/null; then
  echo >&2 "❌ liftOver command not in PATH"; exit 1
fi
# ---------------------------------------------------------------------------

echo >&2 "» Running liftOver ..."
liftOver "$in" "$chain" "$out" "$unmapped"

echo >&2 "✓ liftOver finished             -> $out"
echo >&2 "✓ unmapped variants (optional)  -> $unmapped"

# ------------- final line: **path only**, sent to stdout --------------------
echo "$out"
