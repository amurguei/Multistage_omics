#!/bin/bash

# Required paths
BED=/lustre1/home/mass/amalia/montipora_mapping/ref/Mcap_transcripts.bed
BAM_DIR=/lustre1/home/mass/amalia/montipora_mapping/data/bam_files
TIN_OUT=/lustre1/home/mass/amalia/montipora_mapping/outputs/tin

# Create output directory if needed
mkdir -p $TIN_OUT
cd $BAM_DIR

# Run TIN on each BAM file
for bam in *.bam; do
    echo "Running TIN on $bam"
    tin.py -i "$bam" -r "$BED" > $TIN_OUT/${bam%.bam}_tin.txt
done

echo "✅ All TIN calculations done."
