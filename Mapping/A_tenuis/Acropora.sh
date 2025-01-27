#Download SRA through conda
conda activate
conda install bioconda/label/cf201901::sra-tools

cd /lustre1/home/mass/amalia/acropora/raw_reads/

#Download files in an SRA format
prefetch DRR318299
prefetch DRR318298
prefetch DRR318297
prefetch DRR318296
prefetch DRR318295
prefetch DRR318294
prefetch DRR318293
prefetch DRR318292
prefetch DRR318291
prefetch DRR318290
prefetch DRR318289
prefetch DRR318288

#Loop to extract fastq files
for dir in */; do
    cd "$dir"  # Enter each subdirectory
    for sra_file in *.sra; do
        fasterq-dump --outdir /lustre1/home/mass/amalia/acropora/raw_reads/fastqreads --threads 4 --progress "$sra_file"
    done
    cd ..  # Go back to the parent directory
done

#Once its done, go to the appropriate directory:

md5sum *.fastq > raw_checksum.md5

md5sum -c raw_checksum.md5 #check it's OK for all files
[amalia@hive02 fastqreads]$ md5sum -c raw_checksum.md5 #check it's OK for all files
#DRR318288_1.fastq: OK
#DRR318288_2.fastq: OK
#DRR318289_1.fastq: OK
#DRR318289_2.fastq: OK
#DRR318290_1.fastq: OK
#DRR318290_2.fastq: OK
#DRR318291_1.fastq: OK
#DRR318291_2.fastq: OK
#DRR318292_1.fastq: OK
#DRR318292_2.fastq: OK
#DRR318293_1.fastq: OK
#DRR318293_2.fastq: OK
#DRR318294_1.fastq: OK
#DRR318294_2.fastq: OK
#DRR318295_1.fastq: OK
#DRR318295_2.fastq: OK
#DRR318296_1.fastq: OK
#DRR318296_2.fastq: OK
#DRR318297_1.fastq: OK
#DRR318297_2.fastq: OK
#DRR318298_1.fastq: OK
#DRR318298_2.fastq: OK
#DRR318299_1.fastq: OK
#DRR318299_2.fastq: OK

zgrep -c "@DRR" *fastq.gz #Count number of reads per file
# DRR318288_1.fastq:20047823
# DRR318288_2.fastq:20047823
# DRR318289_1.fastq:23238105
# DRR318289_2.fastq:23238105
# DRR318290_1.fastq:23389634
# DRR318290_2.fastq:23389634
# DRR318291_1.fastq:23224164
# DRR318291_2.fastq:23224164
# DRR318292_1.fastq:19992103
# DRR318292_2.fastq:19992103
# DRR318293_1.fastq:20972670
# DRR318293_2.fastq:20972670
# DRR318294_1.fastq:21060138
# DRR318294_2.fastq:21060138
# DRR318295_1.fastq:17081185
# DRR318295_2.fastq:17081185
# DRR318296_1.fastq:25849471
# DRR318296_2.fastq:25849471
# DRR318297_1.fastq:25210596
# DRR318297_2.fastq:25210596
# DRR318298_1.fastq:27444050
# DRR318298_2.fastq:27444050
# DRR318299_1.fastq:20781111
# DRR318299_2.fastq:20781111

cd /lustre1/home/mass/amalia/acropora/ref/
cd /lustre1/home/mass/amalia/acropora
mkdir data
cd data
mkdir ref
mkdir raw

wget http://aten.reefgenomics.org/download/aten_0.11.maker_post_001.transcripts.fasta.gz

ln -s /lustre1/home/mass/amalia/acropora/ref/aten_0.11.maker_post_001.transcripts.fasta.gz /lustre1/home/mass/amalia/acropora/data/ref/aten_0.11.maker_post_001.transcripts.fasta.gz
ln -s /lustre1/home/mass/amalia/acropora/raw_reads/fastqreads/*.fastq /lustre1/home/mass/amalia/acropora/data/raw/

#Create conda environment for following steps
conda create -n Atenuis
conda activate Atenuis #cd /lustre1/home/mass/amalia/.conda/envs/Atenuis

#Note, download of fast qc didn't work with conda, trying with module load
module load fastqc/0.12.1
module load java/22.0.1 #Needed for fastqc

mkdir /lustre1/home/mass/amalia/acropora/data/fastqc

cd /lustre1/home/mass/amalia/acropora/data/raw/

mkdir /lustre1/home/mass/amalia/acropora/data/fastqc
mv /lustre1/home/mass/amalia/acropora/data/raw/*fastqc* /lustre1/home/mass/amalia/acropora/data/fastqc/

module load multiqc/1.23
multiqc ./

mkdir /lustre1/home/mass/amalia/acropora/data/cleaned_reads/

#Note, had to run fast p with Conda because there was not the module, had to do fastqc with modules because it wouldn't download with conda...
#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --output=fastp_trim_output.txt
#SBATCH --error=fastp_trim_error.txt
#SBATCH --time=01:00:00
#SBATCH --mem=10G

# Activate the conda environment (ensure to use conda activate, assuming Conda is available)
source /lustre1/home/mass/amalia/miniconda3/etc/profile.d/conda.sh  # Adjust if needed
conda activate Atenuis  # Activate the Atenuis environment

# Change to the directory where the raw reads are located
cd /lustre1/home/mass/amalia/acropora/data/raw/

nano trimming.sh
# Create an array of fastq files
array1=($(ls *_1.fastq))  # Adjust the pattern to match your specific files

# Loop through each file in the array and run fastp for trimming
for i in "${array1[@]}"; do
    fastp --in1 ${i} \
          --in2 $(echo ${i} | sed s/_1/_2/) \
          --out1 /lustre1/home/mass/amalia/acropora/data/cleaned_reads/$(basename ${i} _1.fastq)_cleaned_1.fastq \
          --out2 /lustre1/home/mass/amalia/acropora/data/cleaned_reads/$(basename ${i} _1.fastq)_cleaned_2.fastq \
          --detect_adapter_for_pe \
          --qualified_quality_phred 20 \
          --unqualified_percent_limit 10 \
          --cut_right_window_size 5 --cut_right_mean_quality 20 \
          -h /lustre1/home/mass/amalia/acropora/data/cleaned_reads/$(basename ${i} _1.fastq).fastp.html \
          -j /lustre1/home/mass/amalia/acropora/data/cleaned_reads/$(basename ${i} _1.fastq).fastp.json
done

# Optionally, deactivate the conda environment if needed
conda deactivate

sbatch trimming.sh

#Move to the clean reads directory

zgrep -c "@DRR" *fastq
# DRR318288_cleaned_1.fastq:17004065
# DRR318288_cleaned_2.fastq:17004065
# DRR318289_cleaned_1.fastq:19692512
# DRR318289_cleaned_2.fastq:19692512
# DRR318290_cleaned_1.fastq:19965165
# DRR318290_cleaned_2.fastq:19965165
# DRR318291_cleaned_1.fastq:19849816
# DRR318291_cleaned_2.fastq:19849816
# DRR318292_cleaned_1.fastq:17010711
# DRR318292_cleaned_2.fastq:17010711
# DRR318293_cleaned_1.fastq:17380158
# DRR318293_cleaned_2.fastq:17380158
# DRR318294_cleaned_1.fastq:18119428
# DRR318294_cleaned_2.fastq:18119428
# DRR318295_cleaned_1.fastq:14660658
# DRR318295_cleaned_2.fastq:14660658
# DRR318296_cleaned_1.fastq:22157517
# DRR318296_cleaned_2.fastq:22157517
# DRR318297_cleaned_1.fastq:21927844
# DRR318297_cleaned_2.fastq:21927844
# DRR318298_cleaned_1.fastq:23786159
# DRR318298_cleaned_2.fastq:23786159
# DRR318299_cleaned_1.fastq:18936678
# DRR318299_cleaned_2.fastq:18936678

mkdir /lustre1/home/mass/amalia/acropora/data/fastqc_clean
cd /lustre1/home/mass/amalia/acropora/data/fastqc_clean
mv /lustre1/home/mass/amalia/acropora/data/fastqc/*fastqc* ./
multiqc ./

mv /lustre1/home/mass/amalia/acropora/data/fastqc_clean/multiqc_report.html /lustre1/home/mass/amalia/acropora/data/outputs/

#Hisat 2
conda activate Atenuis
conda install -c bioconda hisat2
mkdir /lustre1/home/mass/amalia/acropora/data/hisat2
hisat2-build /lustre1/home/mass/amalia/acropora/ref/aten_final_0.11.fasta /lustre1/home/mass/amalia/acropora/data/hisat2/Aten_ref

nano hisat2.sh 

#!/bin/bash
#SBATCH --job-name="hisat2Aten"
#SBATCH --partition=hiveunlim
#SBATCH -t 100:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amalia.murgueitio@gmail.com
#SBATCH -D /lustre1/home/mass/amalia/acropora/data/
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=128GB

echo "Loading environment by explicitly setting paths" $(date)

# Path to Conda environment binaries
CONDA_ENV_PATH=/lustre1/home/mass/amalia/.conda/envs/Atenuis/bin

# Add Conda environment binaries to PATH
export PATH=$CONDA_ENV_PATH:$PATH

# Activate Conda environment
source activate Atenuis

# Ensure HISAT2 and SAMtools are correctly set
which hisat2
which samtools

# Define working directories
FASTQ_DIR=/lustre1/home/mass/amalia/acropora/data/cleaned_reads
HISAT2_INDEX=/lustre1/home/mass/amalia/acropora/data/hisat2/Aten_ref
BAM_DIR=/lustre1/home/mass/amalia/acropora/data/bam_files/

# Create BAM output directory if it doesn't exist
mkdir -p $BAM_DIR

# Process paired-end reads
for R1 in $FASTQ_DIR/*_cleaned_1.fastq; do
    R2=${R1/_cleaned_1.fastq/_cleaned_2.fastq}
    BASENAME=$(basename $R1 _cleaned_1.fastq)

    echo "Processing: $BASENAME" $(date)

    # Align reads with HISAT2
    hisat2 -p $SLURM_NTASKS_PER_NODE --rf --dta -x $HISAT2_INDEX \
        -1 $R1 -2 $R2 \
        -S $BAM_DIR/${BASENAME}.sam

    # Convert SAM to BAM and sort
    samtools sort -@ $SLURM_NTASKS_PER_NODE -o $BAM_DIR/${BASENAME}.bam $BAM_DIR/${BASENAME}.sam
    rm $BAM_DIR/${BASENAME}.sam

    echo "Completed: $BASENAME" $(date)
done

echo "All jobs finished!" $(date)


sbatch hisat2.sh

#Stringtie
cd /lustre1/home/mass/amalia/acropora/data/

mkdir stringtie
cd stringtie
mkdir bam

cd /lustre1/home/mass/amalia/acropora/data/bam_files/

for bam_file in *.bam; do
    ln -s /lustre1/home/mass/amalia/acropora/data/bam_files/$bam_file /lustre1/home/mass/amalia/acropora/data/stringtie/bam/$bam_file
done

nano stringtie.sh

#!/bin/bash 
#SBATCH --job-name="stringtie"
#SBATCH --export=ALL
#SBATCH --partition=hiveunlim
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amalia.murgueitio@gmail.com
#SBATCH --exclusive
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=128GB
#SBATCH -D /lustre1/home/mass/amalia/acropora/data/stringtie/bam

echo "Loading programs" $(date)

# Path to Conda environment binaries
CONDA_ENV_PATH=/lustre1/home/mass/amalia/.conda/envs/Atenuis/bin

# Add Conda environment binaries to PATH
export PATH=$CONDA_ENV_PATH:$PATH

# Activate Conda environment
source activate Atenuis

# Ensure StringTie is correctly set
which stringtie

# Continue with the rest of the script
echo "Starting assembly!" $(date)

# Specify working directory
F=/lustre1/home/mass/amalia/acropora/data/stringtie/bam

# StringTie reference-guided assembly
array1=($(ls /lustre1/home/mass/amalia/acropora/data/stringtie/bam/*.bam))

for i in "${array1[@]}"; do
    stringtie -A gene_abund/$(basename ${i}).gene_abund.tab -p $SLURM_NTASKS_PER_NODE --rf -e -G /lustre1/home/mass/amalia/acropora/data/stringtie/aten_0.11.maker_post_001.genes.gff -o ${i}.gtf ${i}
    echo "StringTie-assembly-to-ref ${i}" $(date)
done

stringtie --merge -p 8 -G /lustre1/home/mass/amalia/acropora/data/stringtie/aten_0.11.maker_post_001.genes.gff \
-o /lustre1/home/mass/amalia/acropora/data/stringtie/merged.gtf \

gffcompare -r /lustre1/home/mass/amalia/acropora/data/stringtie/aten_0.11.maker_post_001.genes.gff -G -o /lustre1/home/mass/amalia/acropora/data/stringtie/merged_compare /lustre1/home/mass/amalia/acropora/data/stringtie/merged.gtf

# gffcompare v0.12.6 | Command line used:
# gffcompare -r /lustre1/home/mass/amalia/acropora/data/stringtie/aten_0.11.maker_post_001.genes.gff \
#            -G -o /lustre1/home/mass/amalia/acropora/data/stringtie/merged_compare \
#            /lustre1/home/mass/amalia/acropora/data/stringtie/merged.gtf
#

# Summary for dataset: /lustre1/home/mass/amalia/acropora/data/stringtie/merged.gtf

# Query mRNAs: 30326 transcripts analyzed, found in 30326 unique loci.
#               Of these, 26695 are multi-exon transcripts (contain more than one exon).
#               (There are no loci with multiple transcripts. On average, ~1 transcript per locus.)
# Reference mRNAs: 30326 transcripts exist in the reference GFF, distributed across 30326 loci.
#                  (Of these, 26695 are multi-exon reference transcripts.)

# Super-loci with reference transcripts: 30326
#                 (This indicates that all loci were matched to reference "super-loci".)

# Sensitivity and Precision Overview:
#                     -----------------| Sensitivity | Precision |
# Base level:             100.0%     |   100.0%
# Exon level:             100.0%     |   100.0%
# Intron level:           100.0%     |   100.0%
# Intron chain level:     100.0%     |   100.0%
# Transcript level:        88.0%     |    88.0%
# Locus level:             88.0%     |    88.0%
# 
#                 - Base level: Compares nucleotide positions directly.
#                 - Exon level: Evaluates exon boundaries between query and reference.
#                 - Intron level: Matches intron positions.
#                 - Intron chain level: Considers exact matching of multiple consecutive introns.
#                 - Transcript and locus levels are more stringent measures of accuracy.

# Additional Matching Metrics:
#     Matching intron chains: 26695
#        (Exact matches of intron combinations within transcripts.)
#     Matching transcripts: 26695
#        (Entire transcript matches observed between query and reference.)
#     Matching loci: 26695
#        (Transcript loci matched perfectly between the query and reference.)

# Missed Elements:
#     Missed exons: 0/216455  (0.0%)
#        (No exons in the reference annotation were missed by the query data.)
#     Novel exons: 0/216455  (0.0%)
#        (No additional or unexpected exons were identified in the query data.)
#     Missed introns: 0/186129  (0.0%)
#        (No introns from the reference were missed in the query set.)
#     Novel introns: 0/186129  (0.0%)
#        (No new introns were found.)
#     Missed loci: 0/30326   (0.0%)
#        (All reference loci were successfully matched to query loci.)
#     Novel loci: 0/30326   (0.0%)
#        (No new loci were detected.)

# Total Output:
#     Total union super-loci across all input datasets: 30326
#         (Union of overlapping regions from all query loci and the reference.)
#     30326 out of 30326 consensus transcripts written to:
#         /lustre1/home/mass/amalia/acropora/data/stringtie/merged_compare.annotated.gtf
#     (No transcripts discarded as redundant during merging.)

chmod +x /lustre1/home/mass/amalia/.conda/envs/Atenuis/bin/prepDE.py

# Sample name    Path to the corresponding GTF file from StringTie how the file must be formatted, for the merge step is just paths
#DRR318288    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318288.bam.gtf
#DRR318289    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318289.bam.gtf
#DRR318290    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318290.bam.gtf
#DRR318291    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318291.bam.gtf
#DRR318292    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318292.bam.gtf
#DRR318293    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318293.bam.gtf
#DRR318294    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318294.bam.gtf
#DRR318295    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318295.bam.gtf
#DRR318296    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318296.bam.gtf
#DRR318297    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318297.bam.gtf
#DRR318298    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318298.bam.gtf
#DRR318299    /lustre1/home/mass/amalia/acropora/data/stringtie/bam/DRR318299.bam.gtf


#This didn't work
python /lustre1/home/mass/amalia/.conda/envs/Atenuis/bin/prepDE.py -g /lustre1/home/mass/amalia/acropora/outputs/atenuis_gene_count_matrix.csv -i /lustre1/home/mass/amalia/acropora/data/stringtie/Sample_list_A_tenuis2.txt


#Note: I used feature counts since it didn't work with prepDE: 

featureCounts -a /lustre1/home/mass/amalia/acropora/data/stringtie/aten_0.11.maker_post_001.genes.gff -t gene -o /lustre1/home/mass/amalia/acropora/data/outputs/gene_counts.txt -g ID -p /lustre1/home/mass/amalia/acropora/data/bam_files/*.bam
