#!/bin/bash

# Author: RICCARDO PIANEZZA
# Purpose: From bed files containing TE annotations, create a phylogeny tree of all the insertions of a TE family.
# Arguments:
# 1) the name of a TE of interest
# 2) the number of minimum bp for an insertion of the TE to be included in the tree
# 3) the number of maximum bp for an insertion of the TE to be included in the tree
# 4) an output folder
# 5) a folder with the genomes (.fa)
# 6) a folder with the bed files (from the RepeatMasker output file you can run the script "oriout2bed.py" to get the bed files)

te="$1"
length="$2"
max_length="$3"
output_folder="$4"
fasta_folder="$5"
bed_folder="$6"

# Define the paths to python scripts and reference library. Remember to change the path based on
# where you have your files saved. You can find all the files on our GitHub repository.
rename_insertions="/Volumes/Temp1/Dmel-stealthTEs/tree/rename-insertions.py"
calculate_stats="/Volumes/Temp1/Dmel-stealthTEs/tree/divergence.R"
ref_library="/Volumes/Temp1/Dmel-stealthTEs/dmel_TE_invasions/ref/teseqs-3scg-dmel.fasta"
BLAST_PATH="/usr/local/Caskroom/miniconda/base/bin/blastn"

mkdir -p "$output_folder/$te/all"
mkdir -p "$output_folder/$te/fl-$length"
mkdir -p "$output_folder/$te/fl-$length/fasta"
grep -A 1 "${te}" "$ref_library" > "$output_folder/$te/$te.fasta"

# Select BED annotation for the TE of interest (argument 1)
for file in "$bed_folder"/*.bed; do
  if [ -f "$file" ]; then
    grep "${te}_" "$file" > "$output_folder/$te/all/$(basename "$file")"
  fi
done

# Remove fragmented insertions based on TE length (arguments 2-3)
for file in "$output_folder/$te"/all/*.bed; do
  output_file="$output_folder/$te/fl-"$length"/$(basename "$file" .bed.bed)"
  awk -F'\t' -v len="$length" -v max_len="$max_length" '($5 > len && $5 < max_len)' "$file" > "$output_file"
done

# Getfasta from bed file annotation
for file in "$output_folder/$te/fl-$length"/*.bed; do
    basename=$(basename "$file" | sed 's/\.bed$//') 
    output="$output_folder/$te/fl-$length/fasta/$basename.fasta"
    bedtools getfasta -fi "$fasta_folder/$basename.fa" -bed "$file" -fo "$output" -s -name
done

# Concatenate FASTA files in a single file
rm "$output_folder/$te/fl-$length/fasta/$te.fasta"
cat "$output_folder/$te/fl-$length/fasta/"*fasta > $output_folder/$te/fl-$length/fasta/$te.fasta

# BLAST
QUERY_FILE="$output_folder/$te/fl-$length/fasta/$te.fasta"
TARGET_FILE="$output_folder/$te/$te.fasta"
OUTPUT_DIR="$output_folder/$te/fl-$length/blast"
FINAL_OUTPUT_FILE="$OUTPUT_DIR/blast_output.txt"
mkdir -p "$OUTPUT_DIR"
rm "$FINAL_OUTPUT_FILE"

while IFS= read -r line; do
  if [[ $line == ">"* ]]; then
    sequence_id="${line:1}"
    continue
  fi

  echo "$line" | $BLAST_PATH -query <(echo -e ">${sequence_id}\n${line}") -subject "$TARGET_FILE" -outfmt "6 qseqid sseqid pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" >> "$FINAL_OUTPUT_FILE"
done < "$QUERY_FILE"

# Perform Multiple Sequence Alignment (MSA)
mkdir -p "$output_folder/$te/fl-$length/tree-files"
muscle -in "$output_folder/$te/fl-$length/fasta/$te.fasta" -out "$output_folder/$te/fl-$length/tree-files/$te.MSA.fasta"

# Convert MSA fasta into phylip format
muscle -in "$output_folder/$te/fl-$length/tree-files/$te.MSA.fasta" -out "$output_folder/$te/fl-$length/tree-files/$te.MSA.phylip" -refine -phyi
python "$rename_insertions" $output_folder/$te/fl-$length/tree-files/$te.MSA.phylip $output_folder/$te/fl-$length/tree-files/$te.MSA.fasta $output_folder/$te/fl-$length/tree-files/$te.MSA.renamed.phylip $te

# Calculate divergence % of each sequence from the consensus sequence
length_te=$(awk '/^>/ { if (seq) print length(seq); seq=""; } /^([^>])/ { seq = seq $0; } END { if (seq) print length(seq); }' "$output_folder/$te/$te.fasta")
Rscript "$calculate_stats" "$FINAL_OUTPUT_FILE" "$output_folder/$te/fl-$length/tree-files/$te.MSA.renamed.phylip_insertion-names.txt" "$length_te" "$OUTPUT_DIR/divergence.tsv"

# Create the tree files
phyml -i $output_folder/$te/fl-$length/tree-files/$te.MSA.renamed.phylip -d nt -m GTR
mv $output_folder/$te/fl-$length/tree-files/$te.MSA.renamed.phylip_phyml_tree.txt $output_folder/$te-$length"_phyml_tree.txt"