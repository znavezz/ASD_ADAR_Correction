#!/usr/bin/env bash
# vep_ann.sh – annotate a 4-column TSV with VEP + full plugin set
# Usage:
#   bash vep_ann.sh <input.tsv>               # writes default canonical_annotated_results.txt
#   bash vep_ann.sh <input.tsv> <output.txt>  # writes to the file you supply
# The script prints the absolute path of the TXT file on stdout.

set -euo pipefail

# ───────────────────────── arguments ────────────────────────────────────────
in_tsv="$1"                       # required
out_txt_arg="${2:-}"              # optional
[[ -z "${in_tsv:-}" ]] && { echo >&2 "Usage: $0 <input.tsv> [output.txt]"; exit 1; }

in_tsv="$(readlink -f "$in_tsv")"
[[ -f "$in_tsv" ]] || { echo >&2 "❌  File not found: $in_tsv"; exit 1; }

# ───────────────────────── locations ────────────────────────────────────────
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out_dir="$script_dir/vep_results"
mkdir -p "$out_dir"

# unique tmp dir for conversion work
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

# ─── decide output/ stats filenames ─────────────────────────────────────────
if [[ -n "$out_txt_arg" ]]; then             # caller supplied its own file
  out_txt="$(readlink -f "$out_txt_arg")"
  mkdir -p "$(dirname "$out_txt")"
else                                         # fall back to legacy name
  out_txt="$out_dir/canonical_annotated_results.txt"
fi
stats_html="${out_txt%.txt}.stats.html"

# genome FASTA on host
fasta_host="/home/alu/aluguest/Nave_Oded_Project/resources/hg38.fa"
fasta_gz="/home/alu/aluguest/Nave_Oded_Project/resources/hg38.fa.gz"

# VEP cache (includes Plugins/) on host
cache_host="$HOME/.vep"

# ───────────────────────── FASTA (decompress if needed) ─────────────────────
if [[ ! -f "$fasta_host" ]]; then
  echo >&2 "Decompressing $(basename "$fasta_gz") …"
  gunzip -c "$fasta_gz" > "$fasta_host"
else
  echo >&2 "Uncompressed FASTA already present."
fi

# ───────────────────────── ensure FATHMM DB container (optional) ────────────
docker inspect --format='{{.State.Running}}' fathmm-db 2>/dev/null |
  grep -q true ||
  { echo >&2 "Starting fathmm-db container …"; docker start fathmm-db >/dev/null; }

# ───────────────────────── convert TSV → minimal VCF ────────────────────────
echo >&2 "Converting $(basename "$in_tsv") to VCF …"

awk 'BEGIN { FS=OFS="\t" }
NR==1 {
  print "##fileformat=VCFv4.2";
  print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT";
  next
}
{
  id=$1":"$2":"$3":"$4;
  print $1,$2,id,$3,$4,".",".",".",".";
}' "$in_tsv" > "$tmp_dir/raw.vcf"

echo >&2 "Sorting & de-duplicating …"
{
  head -n1 "$tmp_dir/raw.vcf"
  tail -n +2 "$tmp_dir/raw.vcf" |
    sort -k1,1 -k2,2n |
    awk 'BEGIN{FS=OFS="\t"} !seen[$1,$2,$4,$5]++'
} > "$tmp_dir/sorted.vcf"

# ───────────────────────── run VEP (Docker) ────────────────────────────────
echo >&2 "Running VEP …"

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$(readlink -f "$cache_host")":/data/vep \
  -v "$tmp_dir":/data/input \
  -v "$out_dir":/data/output \
  -v /home/alu/aluguest/Nave_Oded_Project/resources:/data/fasta \
  ensemblorg/ensembl-vep vep \
      --offline \
      --cache --dir_cache /data/vep \
      --species homo_sapiens --assembly GRCh38 \
      --fasta /data/fasta/hg38.fa \
      --format vcf \
      -i /data/input/sorted.vcf \
      -o /data/output/"$(basename "$out_txt")" \
      --stats_file /data/output/"$(basename "$stats_html")" \
      --everything \
      --sift b --polyphen b --protein \
      --pick --force_overwrite \
      --plugin AlphaMissense,file=/data/vep/Plugins/AlphaMissense_hg38.tsv.gz,cols=all \
      --plugin Blosum62 \
      --plugin Condel,/data/vep/Plugins/Condel/config,b \
      --plugin ClinPred,file=/data/vep/Plugins/ClinPred_hg38_sorted_tabbed.tsv.gz \
      --plugin Enformer,file=/data/vep/Plugins/enformer_grch38.vcf.gz \
      --plugin FlagLRG,/data/vep/Plugins/list_LRGs_transcripts_xrefs.txt \
      --plugin LoFtool,/data/vep/Plugins/LoFtool_scores.txt \
      --plugin LOVD \
      --plugin MaveDB,file=/data/vep/Plugins/MaveDB_variants.tsv.gz \
      --plugin NMD \
      --plugin PhenotypeOrthologous,file=/data/vep/Plugins/PhenotypesOrthologous_homo_sapiens_112_GRCh38.gff3.gz \
      --plugin pLI,/data/vep/Plugins/pLI_values.txt \
      --plugin REVEL,file=/data/vep/Plugins/new_tabbed_revel_grch38.tsv.gz \
      --custom file=/data/vep/homo_sapiens/113_GRCh38/custom/hg38.phastCons100way.bw,short_name=PhastCons100,format=bigwig

echo >&2 "✓  VEP finished: $(basename "$out_txt")"
echo >&2 "✓  Stats file  : $(basename "$stats_html")"

# ───────────────────────── final line (stdout) ─────────────────────────────
echo "$out_txt"
