#!/bin/bash

# Set the output directory
OUT_DIR="/home/project18/BVSim/empirical/hg19"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Download the TR regions data if not already present
if [ ! -f "$OUT_DIR/hg19.simpleRepeat.bed.gz" ]; then
    wget -P "$OUT_DIR" https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/resources/hg19.simpleRepeat.bed.gz
fi

if [ ! -f "$OUT_DIR/hg19.simpleRepeat.bed.gz.tbi" ]; then
    wget -P "$OUT_DIR" https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/resources/hg19.simpleRepeat.bed.gz.tbi
fi

# Process each chromosome from 1 to 21
for chr in {1..21}; do
    echo "Processing chromosome chr${chr}..."
    
    # Extract the relevant columns and create the BED file
    zcat "$OUT_DIR/hg19.simpleRepeat.bed.gz" | awk -v chr="chr${chr}" 'BEGIN{OFS="\t"} $1 == chr {print $1, $2, $3}' > "$OUT_DIR/windows_TR_chr${chr}.bed"
    
    # Merge overlapping intervals and remove duplicates
    bedtools sort -i "$OUT_DIR/windows_TR_chr${chr}.bed" | bedtools merge -i stdin | awk 'BEGIN{OFS="\t"} {$4="TR"; print}' | uniq > "$OUT_DIR/windows_TR_unique_chr${chr}.bed"
    
    # Create a final BED file with start and end positions
    awk '{print $2 "\t" $3}' "$OUT_DIR/windows_TR_unique_chr${chr}.bed" > "$OUT_DIR/chr${chr}_TR_unique.bed"
    
    # Clean up intermediate files
    rm "$OUT_DIR/windows_TR_chr${chr}.bed" "$OUT_DIR/windows_TR_unique_chr${chr}.bed"
    
    echo "Finished processing chromosome chr${chr}"
done

echo "All chromosomes processed successfully!"