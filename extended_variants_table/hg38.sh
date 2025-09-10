#!/bin/bash

table_file="hg38_extended_table.tsv"
genome_file="../resources/hg38.fa"
output_file="hg38_extended_table_tmp.tsv"

# sanity checks
[[ -f "$table_file"   ]] || { echo "ERROR: $table_file not found";   exit 1; }
[[ -f "$genome_file"  ]] || { echo "ERROR: $genome_file not found";  exit 1; }

echo "Building $output_file …"

# write header + new column name
head -n1 "$table_file" | awk '{print $0"\thg38"}'            > "$output_file"

# skip header, process each line
tail -n +2 "$table_file" | \
while IFS=$'\t' read -r chr pos ref alt rest; do
  # compute 0-based bed interval
  # add "chr" prefix

    chr=$(echo "$chr" | sed 's/^chr//g')  # First remove any existing chr prefix
    chr="chr$chr"  # Then add chr prefix
    start=$((pos - 1))
    end=$((start + ${#ref}))

    # fetch hg38    
    hg38_seq=$(printf "%s\t%s\t%s\n" "$chr" "$start" "$end" \
                | bedtools getfasta -fi "$genome_file" -bed - -fo - \
                | tail -n +2 | tr -d '\n' | tr 'a-z' 'A-Z')

    # append original + hg38 to output
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$chr" "$pos" "$ref" "$alt" "$rest" "$hg38_seq" >> "$output_file"
done

rm $table_file
mv $output_file $table_file
rm $output_file
echo "Done: saved with hg38 column to $table_file"

