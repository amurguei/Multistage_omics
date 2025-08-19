#Log into hive
ssh amalia@hive02.haifa.ac.il

pwd
/lustre1/home/mass/amalia

mkdir montipora_mapping
cd montipora_mapping

#get an interactive session
salloc --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G --time=04:00:00 --partition=hiveunlim
#enter the node manually
srun --pty /bin/bash



#List of samples
| **Sample Name**                        | **Library Name** | **Run Accession (SRR)** | **Spots**  | **Bases** | **Size** | **Date**   |
| -------------------------------------- | ---------------- | ----------------------- | ---------- | --------- | -------- | ---------- |
| Mcap\_Swimming\_Larvae\_183hpf\_1      | AH1              | SRR14864072             | 36,124,077 | 10.8G     | 3.2Gb    | 2021-06-18 |
| Mcap\_Swimming\_Larvae\_183hpf\_2      | AH2              | SRR14864071             | 39,106,593 | 11.7G     | 3.5Gb    | 2021-06-18 |
| Mcap\_Swimming\_Larvae\_183hpf\_3      | AH3              | SRR14864070             | 36,440,630 | 10.9G     | 3.3Gb    | 2021-06-18 |
| Mcap\_Metamorphosed\_Larvae\_231hpf\_1 | AH4              | SRR14864069             | 36,286,952 | 10.9G     | 3.3Gb    | 2021-06-18 |
| Mcap\_Metamorphosed\_Larvae\_231hpf\_2 | AH5              | SRR14864068             | 34,909,249 | 10.5G     | 3.1Gb    | 2021-06-18 |
| Mcap\_Metamorphosed\_Larvae\_231hpf\_3 | AH6              | SRR14864067             | 41,148,636 | 12.3G     | 3.7Gb    | 2021-06-18 |
| Mcap\_Recruits\_5dps\_1                | AH7              | SRR14864066             | 44,802,304 | 13.4G     | 4Gb      | 2021-06-18 |
| Mcap\_Recruits\_5dps\_2                | AH8              | SRR14864065             | 41,700,000 | 12.5G     | 3.7Gb    | 2021-06-18 |
| Mcap\_Recruits\_5dps\_3                | AH9              | SRR14864064             | 49,300,000 | 14.8G     | 4.4Gb    | 2021-06-18 |

#Creating a new conda environment with the newest SRA tools
conda create -n Mcap -c conda-forge -c bioconda sra-tools
conda activate Mcap
#confirm toolkit works
which prefetch
prefetch --version

#Making an executable script in the interactive session
nano download_mcap.sh

#!/bin/bash

# List of SRR accessions
SRR_LIST=(
SRR14864072
SRR14864071
SRR14864070
SRR14864069
SRR14864068
SRR14864067
SRR14864066
SRR14864065
SRR14864064
)

# Output directory for FASTQ files
OUTDIR="/lustre1/home/mass/amalia/montipora_mapping/fastqreads"
mkdir -p "$OUTDIR"

# Download + Convert loop
for SRR in "${SRR_LIST[@]}"; do
    echo "📥 Downloading $SRR..."
    prefetch "$SRR"

    echo "🔄 Converting $SRR to FASTQ..."
    fasterq-dump "$SRR" --threads 4 --progress --outdir "$OUTDIR"
done

echo "All downloads and conversions completed."

chmod +x download_mcap.sh
conda activate Mcap #if not already activated
./download_mcap.sh

#When done
cd fastqreads
md5sum *.fastq > raw_checksum.md5

(Mcap) [I have no name!@bee77 fastqreads]$ md5sum -c raw_checksum.md5 #check it's OK for all files
SRR14864064_1.fastq: OK
SRR14864064_2.fastq: OK
SRR14864065_1.fastq: OK
SRR14864065_2.fastq: OK
SRR14864066_1.fastq: OK
SRR14864066_2.fastq: OK
SRR14864067_1.fastq: OK
SRR14864067_2.fastq: OK
SRR14864068_1.fastq: OK
SRR14864068_2.fastq: OK
SRR14864069_1.fastq: OK
SRR14864069_2.fastq: OK
SRR14864070_1.fastq: OK
SRR14864070_2.fastq: OK
SRR14864071_1.fastq: OK
SRR14864071_2.fastq: OK
SRR14864072_1.fastq: OK
SRR14864072_2.fastq: OK


zgrep -c "@SRR" *fastq #Count number of reads per file
SRR14864064_1.fastq:49259211
SRR14864064_2.fastq:49259211
SRR14864065_1.fastq:41742625
SRR14864065_2.fastq:41742625
SRR14864066_1.fastq:44802304
SRR14864066_2.fastq:44802304
SRR14864067_1.fastq:41148636
SRR14864067_2.fastq:41148636
SRR14864068_1.fastq:34909249
SRR14864068_2.fastq:34909249
SRR14864069_1.fastq:36286952
SRR14864069_2.fastq:36286952
SRR14864070_1.fastq:36440630
SRR14864070_2.fastq:36440630
SRR14864071_1.fastq:39106593
SRR14864071_2.fastq:39106593
SRR14864072_1.fastq:36124077
SRR14864072_2.fastq:36124077


# Read counts are within a reasonable range (35M–49M) — no sample looks problematic.

# Slight dip in reads for Metamorphosed Larvae Rep 2 (SRR14864068) — but still decent depth.

# Spat (Recruits) tend to have higher read counts overall, especially SRR14864064 (Rep 3).

cd /lustre1/home/mass/amalia/montipora_mapping
mkdir ref
mkdir raw

#Getting the genome assembly file
wget -O - http://cyanophora.rutgers.edu/montipora/Montipora_capitata_HIv3.assembly.fasta.gz | gunzip > Montipora_capitata_HIv3.assembly.fasta

#Making symbolic links
ln -s /lustre1/home/mass/amalia/montipora_mapping/Montipora_capitata_HIv3.assembly.fasta /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.assembly.fasta
ln -s /lustre1/home/mass/amalia/montipora_mapping/fastqreads/*.fastq /lustre1/home/mass/amalia/montipora_mapping/raw/

#Install fastqc 
conda install -c bioconda fastqc

mkdir -p /lustre1/home/mass/amalia/montipora_mapping/fastqc
fastqc /lustre1/home/mass/amalia/montipora_mapping/raw/*.fastq -o /lustre1/home/mass/amalia/montipora_mapping/fastqc

#after done: download multiqc
conda install -c bioconda multiqc

#New Folder for MultiQC
mkdir -p /lustre1/home/mass/amalia/montipora_mapping/outputs/multiqc
multiqc /lustre1/home/mass/amalia/montipora_mapping/fastqc/ \
  -o /lustre1/home/mass/amalia/montipora_mapping/outputs/multiqc

#New directory for cleaned reads
cd /lustre1/home/mass/amalia/montipora_mapping/
mkdir cleaned reads

conda install -c bioconda fastp fastqc multiqc samtools

nano fastp_trim.sh

#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --output=fastp_trim_output.txt
#SBATCH --error=fastp_trim_error.txt
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Activate Conda (via bashrc fallback)
source ~/.bashrc
conda activate Mcap

# Define input/output directories
RAW_DIR=/lustre1/home/mass/amalia/montipora_mapping/raw
OUT_DIR=/lustre1/home/mass/amalia/montipora_mapping/cleaned_reads

# Make sure output directory exists
mkdir -p $OUT_DIR

# Move into raw directory
cd $RAW_DIR

# Create array of *_1.fastq files (paired-end)
array1=($(ls *_1.fastq))

# Loop through and run fastp
for i in "${array1[@]}"; do
    base=$(basename ${i} _1.fastq)
    echo "Processing $base..."
    fastp --in1 ${RAW_DIR}/${base}_1.fastq \
          --in2 ${RAW_DIR}/${base}_2.fastq \
          --out1 ${OUT_DIR}/${base}_cleaned_1.fastq \
          --out2 ${OUT_DIR}/${base}_cleaned_2.fastq \
          --detect_adapter_for_pe \
          --qualified_quality_phred 20 \
          --unqualified_percent_limit 10 \
          --cut_right_window_size 5 \
          --cut_right_mean_quality 20 \
          -h ${OUT_DIR}/${base}.fastp.html \
          -j ${OUT_DIR}/${base}.fastp.json
done

echo "All samples processed."

#ctrl S ctrl X

sbatch fastp_trim.sh

cd /lustre1/home/mass/amalia/montipora_mapping/cleaned_reads
grep -c "@SRR" *fastq
SRR14864064_cleaned_1.fastq:42602597
SRR14864064_cleaned_2.fastq:42602597
SRR14864065_cleaned_1.fastq:36494777
SRR14864065_cleaned_2.fastq:36494777
SRR14864066_cleaned_1.fastq:39614540
SRR14864066_cleaned_2.fastq:39614540
SRR14864067_cleaned_1.fastq:35817524
SRR14864067_cleaned_2.fastq:35817524
SRR14864068_cleaned_1.fastq:30543092
SRR14864068_cleaned_2.fastq:30543092
SRR14864069_cleaned_1.fastq:31706262
SRR14864069_cleaned_2.fastq:31706262
SRR14864070_cleaned_1.fastq:31857800
SRR14864070_cleaned_2.fastq:31857800
SRR14864071_cleaned_1.fastq:34565571
SRR14864071_cleaned_2.fastq:34565571
SRR14864072_cleaned_1.fastq:31441149
SRR14864072_cleaned_2.fastq:31441149

# 1. Create a folder to store FastQC results from cleaned reads
mkdir -p /lustre1/home/mass/amalia/montipora_mapping/outputs/fastqc_clean

# 2. Run FastQC on cleaned reads and save output there
fastqc /lustre1/home/mass/amalia/montipora_mapping/cleaned_reads/*.fastq \
  -o /lustre1/home/mass/amalia/montipora_mapping/outputs/fastqc_clean

# 3. Move into that directory
cd /lustre1/home/mass/amalia/montipora_mapping/outputs/fastqc_clean

# 4. Run MultiQC to summarize results
multiqc ./ -o ./  # Output summary will be `multiqc_report.html` here

# Rename the report
mv multiqc_report.html Multiqc_report_clean.html

# Move it to your target folder
mv Multiqc_report_clean.html /lustre1/home/mass/amalia/montipora_mapping/outputs/multiqc/

# Activate Mcap conda environment
conda activate Mcap

# Install Hisat2 (if not already installed)
conda install -c bioconda hisat2

# Create directory to store index
mkdir -p /lustre1/home/mass/amalia/montipora_mapping/data/hisat2

# Build Hisat2 index
hisat2-build /lustre1/home/mass/amalia/montipora_mapping/Montipora_capitata_HIv3.assembly.fasta \
  /lustre1/home/mass/amalia/montipora_mapping/data/hisat2/Mcap_ref

nano hisat2_align_with_logs.sh

#Ran with logs

#!/bin/bash
#SBATCH --job-name=hisat2_mcap_log
#SBATCH --partition=hive7d
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=hisat2_mcap_log_output.txt
#SBATCH --error=hisat2_mcap_log_error.txt

echo "Starting HISAT2 alignment with logging on $(date)"

# Load Conda
source ~/.bashrc
conda activate Mcap

# Confirm tools are available
which hisat2
which samtools

# Set directories
FASTQ_DIR=/lustre1/home/mass/amalia/montipora_mapping/cleaned_reads
HISAT2_INDEX=/lustre1/home/mass/amalia/montipora_mapping/data/hisat2/Mcap_ref
BAM_DIR=/lustre1/home/mass/amalia/montipora_mapping/data/bam_files

mkdir -p $BAM_DIR

# Align all samples
for R1 in $FASTQ_DIR/*_cleaned_1.fastq; do
    R2=${R1/_cleaned_1.fastq/_cleaned_2.fastq}
    BASENAME=$(basename $R1 _cleaned_1.fastq)

    echo "Aligning $BASENAME on $(date)"

    hisat2 -p $SLURM_CPUS_PER_TASK --dta --rf -x $HISAT2_INDEX \
        -1 $R1 -2 $R2 \
        -S $BAM_DIR/${BASENAME}.sam \
        2> $BAM_DIR/${BASENAME}_hisat2.log

    samtools sort -@ $SLURM_CPUS_PER_TASK -o $BAM_DIR/${BASENAME}.bam $BAM_DIR/${BASENAME}.sam
    rm $BAM_DIR/${BASENAME}.sam

    echo "Finished $BASENAME on $(date)"
done

echo "All alignments completed on $(date)"

sbatch hisat2_align_with_logs.sh

#MultiQC on hisat logs
/lustre1/home/mass/amalia/montipora_mapping/data/bam_files/ \
  -o /lustre1/home/mass/amalia/montipora_mapping/outputs/multiqc_hisat2

cd /lustre1/home/mass/amalia/montipora_mapping
mkdir stringtie

conda install -c bioconda stringtie

nano stringtie_mods.sh

#SBATCH --partition=hive7d
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amalia.murgueitio@gmail.com
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=128GB
#SBATCH -t 100:00:00
#SBATCH -D /lustre1/home/mass/amalia/montipora_mapping/data/bam_files

echo "==== Starting StringTie job ====" $(date)

# Load conda and activate Mcap environment
source ~/.bashrc
conda activate Mcap

# Confirm StringTie is available
which stringtie

# Define output directory
OUTPUT_DIR=/lustre1/home/mass/amalia/montipora_mapping/stringtie_out
mkdir -p $OUTPUT_DIR

# Loop through BAM files
for bam in *.bam; do
    base=$(basename "$bam" .bam)
    echo "Processing $base" $(date)

    stringtie "$bam" \
      -p $SLURM_NTASKS_PER_NODE \
      -e \
      --conservative \
      -G /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.genes_fixed.gff3 \
      -o $OUTPUT_DIR/${base}.gtf \
      -A $OUTPUT_DIR/${base}.abund.tab

    echo "Finished $base" $(date)
done

echo "==== All StringTie assemblies complete ====" $(date)

mkdir -p /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie
ls /lustre1/home/mass/amalia/montipora_mapping/data/bam_files/*.gtf > /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/stringtie_mergelist.txt

stringtie --merge -p 8 \
-G /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.gtf \
-o /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/merged.gtf \
/lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/stringtie_mergelist.txt

gffcompare \
-r /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.gtf \
-G \
-o /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/mcap_compare \
/lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/merged.gtf

conda install -c bioconda gffcompare

#gffcompare -r /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.gtf \
#           -G \
#           -o /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/mcap_compare \
#           /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/merged.gtf

# Summary for dataset: merged.gtf
# Query mRNAs     : 54384 in 54185 loci (36023 multi-exon transcripts)
# Reference mRNAs : 54384 in 54185 loci (36023 multi-exon transcripts)

# This means every reference transcript and locus was represented in the merged.gtf output.

# Super-loci with reference transcripts: 54185
# Total union super-loci: 54185

# ===================== Accuracy Metrics =====================
#                     Sensitivity     Precision
# Base level:           100.0%         100.0%
# Exon level:           100.0%         100.0%
# Intron level:         100.0%         100.0%
# Intron chain level:   100.0%         100.0%
# Transcript level:      66.2%          66.2%
# Locus level:           66.4%          66.4%

# Interpretation:
# - Base/Exon/Intron level: Perfect recovery of individual transcript features.
# - Transcript level: 66.2% of full transcripts (exact exon-intron structure) match reference.
# - Locus level: 66.4% of gene loci match reference structure fully.
#   This is still very good, as UTRs and alternative splicing can slightly affect full matches.

# Novelty/Missing features:
# - Missed exons/introns/loci: 0 — Everything in the reference was detected.
# - Novel exons/introns/loci: 0 — Nothing new was predicted, as expected with `-e -G` flags.

# Output:
# 54384 consensus transcripts were successfully written to:
# /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/mcap_compare.annotated.gtf

# Note:
# This was a reference-guided quantification run (not assembly), so novel features are intentionally not detected.

/lustre1/home/mass/amalia/montipora_mapping/data/bam_files/

mkdir -p /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie


ls /lustre1/home/mass/amalia/montipora_mapping/stringtie_out/*.gtf > /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/stringtie_mergelist.txt

stringtie --merge \
  -p 8 \
  -G /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.gtf \
  -o /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/merged.gtf \
  /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/stringtie_mergelist.txt


ls /lustre1/home/mass/amalia/montipora_mapping/data/bam_files/*.gtf | \
awk -F'/' '{fname=$NF; sub(".gtf", "", fname); print fname "\t" $0}' \
> /lustre1/home/mass/amalia/montipora_mapping/stringtie_sample_list.tsv

gffcompare \
-r /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.gtf \
-G \
-o /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/mcap_compare \
/lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/merged.gtf

python /lustre1/home/mass/amalia/.conda/envs/Mcap/bin/prepDE.py \
  -i /lustre1/home/mass/amalia/montipora_mapping/stringtie_sample_list.tsv \
  -g gene_count_matrix.csv \
  -t transcript_count_matrix.csv \

#Alternative to prepDE
 featureCounts \
  -a /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.gtf \
  -t transcript \
  -g gene_id \
  -p \
  -o /lustre1/home/mass/amalia/montipora_mapping/outputs/stringtie/mcap_gene_counts.txt \
  /lustre1/home/mass/amalia/montipora_mapping/data/bam_files/*.bam

  #Additional: 
# [1] Set up a clean RSeQC environment
conda create -n rseqc_env python=3.9 -y
conda activate rseqc_env

# [2] Install dependencies (if not already installed)
conda install -c bioconda rseqc ucsc-gtfToGenePred ucsc-genePredToBed

# [3] Convert GTF to BED12 format (required for TIN)
cd /lustre1/home/mass/amalia/montipora_mapping/ref

gtfToGenePred Montipora_capitata_HIv3.gtf Mcap.genePred
genePredToBed Mcap.genePred Mcap_transcripts_clean.bed

# Optional: verify BED format (should be 12 fields)
awk '{ print NF }' Mcap_transcripts_clean.bed | sort | uniq -c

# [4] Run TIN calculation using RSeQC
mkdir -p /lustre1/home/mass/amalia/montipora_mapping/outputs/tin
cd /lustre1/home/mass/amalia/montipora_mapping/data/bam_files

for bam in *.bam; do
  echo "Running TIN for $bam"
  tin.py -i "$bam" \
         -r /lustre1/home/mass/amalia/montipora_mapping/ref/Mcap_transcripts_clean.bed \
         > /lustre1/home/mass/amalia/montipora_mapping/outputs/tin/${bam%.bam}.tin.xls
done

# [5] Compute average TIN score per sample
cd /lustre1/home/mass/amalia/montipora_mapping/outputs/tin

echo -e "Sample\tTIN" > multiqc_tin.txt
for f in *.tin.xls; do
  sample=$(basename "$f" .tin.xls)
  tin=$(awk 'NR>1 && $5 > 0 {sum += $5; n++} END {if (n > 0) printf "%.4f", sum/n; else print "0"}' "$f")
  echo -e "${sample}\t${tin}" >> multiqc_tin.txt
done
