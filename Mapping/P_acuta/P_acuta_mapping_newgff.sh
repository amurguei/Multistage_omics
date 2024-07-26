#Preliminary steps

#Connecting to the cluster. In Power Shell
ssh -l amalia.murgueitiocal ssh3.hac.uri.edu
#Insert private key

#Set directory
cd /data/putnamlab/Amalia_Murgueitio/M_capitata/Mapping_P_acuta

#Copy raw reads (specific files, not the folder) from other folder in the cluster to my own folder

cp /data/putnamlab/erin_chille/BSF_3_Stage/Pacu/raw_reads/*.fastq.gz /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/raw_reads

#Upload reference genome

cd /data/putnamlab/Amalia_Murgueitio/M_capitata/Mapping_P_acuta

curl -o- http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv2.assembly.fasta.gz | gunzip > Pocillopora_acuta_HIv2.assembly.fasta

#md5sum command

cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/raw_reads

md5sum *.fastq.gz > raw_checksum.md5
md5sum -c raw_checksum.md5 #check it's OK for all files
#SRR3051863.fastq.gz: OK
#SRR3051864.fastq.gz: OK
#SRR3051865.fastq.gz: OK
#SRR3051866.fastq.gz: OK
#SRR3051867.fastq.gz: OK
#SRR3051868.fastq.gz: OK
#SRR3051869.fastq.gz: OK
#SRR3051870.fastq.gz: OK
#SRR3051871.fastq.gz: OK

zgrep -c "@SRR" *fastq.gz #Count number of reads per file
#SRR3051863.fastq.gz:10693684
#SRR3051864.fastq.gz:13818643
#SRR3051865.fastq.gz:14071891
#SRR3051866.fastq.gz:11895244
#SRR3051867.fastq.gz:11908459
#SRR3051868.fastq.gz:12332327
#SRR3051869.fastq.gz:14306489
#SRR3051870.fastq.gz:12539938
#SRR3051871.fastq.gz:13575827

#Create symbolic links for raw reads and reference genome

/data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/raw_reads/SRR3051866.fastq.gz
ln -s /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/raw_reads/*.fastq.gz /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/raw/
ln -s /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/Pocillopora_acuta_HIv2.assembly.fasta /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/ref/

#Install Miniconda 3. Omit if already installed
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh #accept terms and select correct location

#Conda create
conda create -n Pacunew
conda activate Pacunew

#Download packages to conda environment, press y after each installation to accept the terms
conda install fastqc
conda install multiqc
conda install fastp
conda install hisat2
conda install samtools

#Run fast qc in raw reads
cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/raw/
fastqc ./*.fastq.gz
mkdir /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/fastqc_raw
cd fastqc_raw
mv /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/raw/*fastqc* ./
multiqc ./

#Move multi qc report to new folder cleaned_reads
mkdir /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/cleaned_reads
mv  /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/fastqc_raw/multiqc_report.html  /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/outputs/multiqc_report.html

#Move to directory with raw reads
cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/raw_reads

#Trimming

array1=($(ls *.fastq.gz)) # Make an array of sequences to trim

for i in "${array1[@]}"; do
    fastp --in1 ${i} \
          --in2 $(echo ${i} | sed s/_R1/_R2/) \
          --out1 /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/cleaned_reads/$(basename ${i} _R1_001.fastq.gz)_cleaned_R1.fastq.gz \
          --out2 /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/cleaned_reads/$(basename ${i} _R1_001.fastq.gz)_cleaned_R2.fastq.gz \
          --detect_adapter_for_pe \
          --qualified_quality_phred 20 \
          --unqualified_percent_limit 10 \
          --cut_right_window_size 5 --cut_right_mean_quality 20 \
          -h /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/cleaned_reads/$(basename ${i} _R1_001.fastq.gz).fastp.html \
          -j /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/cleaned_reads/$(basename ${i} _R1_001.fastq.gz).fastp.json
done

#Move to clean reads directory
cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/cleaned_reads

#Count read per file in clean reads
zgrep -c "@SRR" *fastq.gz

#SRR3051863.fastq.gz_cleaned_R1.fastq.gz:10264955
#SRR3051863.fastq.gz_cleaned_R2.fastq.gz:10264955
#SRR3051864.fastq.gz_cleaned_R1.fastq.gz:13283506
#SRR3051864.fastq.gz_cleaned_R2.fastq.gz:13283506
#SRR3051865.fastq.gz_cleaned_R1.fastq.gz:13516600
#SRR3051865.fastq.gz_cleaned_R2.fastq.gz:13516600
#SRR3051866.fastq.gz_cleaned_R1.fastq.gz:10844976
#SRR3051866.fastq.gz_cleaned_R2.fastq.gz:10844976
#SRR3051867.fastq.gz_cleaned_R1.fastq.gz:10876439
#SRR3051867.fastq.gz_cleaned_R2.fastq.gz:10876439
#SRR3051868.fastq.gz_cleaned_R1.fastq.gz:11208816
#SRR3051868.fastq.gz_cleaned_R2.fastq.gz:11208816
#SRR3051869.fastq.gz_cleaned_R1.fastq.gz:13672578
#SRR3051869.fastq.gz_cleaned_R2.fastq.gz:13672578
#SRR3051870.fastq.gz_cleaned_R1.fastq.gz:11947382
#SRR3051870.fastq.gz_cleaned_R2.fastq.gz:11947382
#SRR3051871.fastq.gz_cleaned_R1.fastq.gz:12974476
#SRR3051871.fastq.gz_cleaned_R2.fastq.gz:12974476

mkdir /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/fastqc_clean
cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/fastqc_clean
mv /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/cleaned_reads/*fastqc* ./
multiqc ./

Note: I changed the name of the multiqc report in Cyberduck to multiqc_report_cleaned.html

mv /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/fastqc_clean/multiqc_report_cleaned.html /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/outputs/

mkdir /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2
cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2

ln -s /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/cleaned_reads/*fastq* ./

#Hisat2
cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/ref/
hisat2-build -f /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/ref/Pocillopora_acuta_HIv2.assembly.fasta ./Pacu_ref

mv /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2/*.ht2 /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/ref/

nano PACUHISAT2

#!/bin/bash
#SBATCH --job-name="hisat2Pacu"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amalia.murgueitio@gmail.com
#SBATCH -D /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2/
#SBATCH --exclusive
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=500GB
echo "Loading programs" $(date)
module load HISAT2/2.2.1-foss-2019b #Alignment to reference genome: HISAT2
module load SAMtools/1.16.1-GCC-11.3.0 #Preparation of alignment for assembly: SAMtools

# Specify working directory
F=/data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2/
# Aligning paired-end reads

array1=($(ls $F/*cleaned_R1.fastq.gz*))

# Full path to the HISAT2 index folder
HISAT2_INDEX_FOLDER=/data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/ref/

# This then makes it into a bam file
# And then also sorts the bam file because Stringtie takes a sorted file for input
# And then removes the sam file because I don't need it anymore
for i in "${array1[@]}"; do
    echo "Processing file: $i"
    hisat2 -p $SLURM_NTASKS_PER_NODE --rf --dta -q -x "${HISAT2_INDEX_FOLDER}/Pacu_ref" -1 "${i}" -2 "$(echo ${i} | sed 's/_R1/_R2/')" -S "${i%.fastq.gz}.sam"
    echo "hisat2 finished for $i"
    samtools sort -@ $SLURM_NTASKS_PER_NODE -o "${i%.fastq.gz}.bam" "${i%.fastq.gz}.sam"
    echo "Sorted BAM file: ${i%.fastq.gz}_sorted.bam"
    echo "${i}_bam"
    rm "${i%.fastq.gz}.sam"
    echo "Removed SAM file: ${i%.fastq.gz}.sam"
    echo "HISAT2 PE ${i}" $(date)
done

#Execute the job
chmod u+x PacuHISAT2.sh
sbatch ./PacuHISAT2.sh

mkdir /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie
cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie


mkdir /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2/bam_july

mv /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2/*_cleaned_R1.bam /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2/bam_july/

cd /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie
mkdir bam_july
ln -s /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/hisat2/bam_july/*.bam /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie/bam_july

nano Stringtie_july

#!/bin/bash
#SBATCH --job-name="stringtie"
#SBATCH -t 100:00:00
#SBATCH --export=ALL
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amalia.murgueitio@gmail.com
#SBATCH --exclusive
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=128GB
#SBATCH -q putnamlab
#SBATCH -D /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie

echo "Loading programs" $(date)
module load StringTie/2.1.4-GCC-9.3.0 #Transcript assembly: StringTie
module load GffCompare/0.12.1-GCCcore-8.3.0 #Transcript assembly QC: GFFCompare
module load Python/2.7.18-GCCcore-9.3.0 #Python for prepde.py script
module list

#StringTie reference-guided assembly

#Prior to the Stringtie steps, upload the fixed GFF file available in: 
#https://github.com/amurguei/Multistage_omics/blob/main/Mapping/P_acuta/Pocillopora_acuta_HIv2.genes_fixed.gff3.gz
#and I fixed with this code: https://github.com/amurguei/Multistage_omics/blob/main/Mapping/P_acuta/fix_gff_format_P_acuta.Rmd
echo "Starting assembly!" $(date)

#Specify working directory
F=/data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie

#StringTie reference-guided assembly
#Has the R1 in array1 because of the naming convention in the former script. However, these BAM files contain both forward and reverse reads.
array1=($(ls ${F}/bam_july/*_cleaned_R1.bam))

for i in "${array1[@]}"; do
        stringtie -A gene_abund/$(basename ${i}).gene_abund.tab -p $SLURM_NTASKS_PER_NODE --rf -e -G Pocillopora_acuta_HIv2.genes_fixed.gff3 -o ${i}.gtf ${i}
        #mv /ref-guided-gtfs/${i}.gtf
        echo "StringTie-assembly-to-ref ${i}" $(date)
done

sbatch ./Stringtie_july

stringtie --merge -p 8 -G /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie/Pocillopora_acuta_HIv2.genes_fixed.gff3 -o /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie/stringtiejuly_merged.gtf /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie/sample_list_Pacu_july.txt

gffcompare -r /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie/Pocillopora_acuta_HIv2.genes_fixed.gff3 -G -o /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie/merged_july /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie/stringtiejuly_merged.gtf

interactive
module purge
module load Python/2.7.18-GCCcore-9.3.0
python2 /opt/software/StringTie/2.1.4-GCC-9.3.0/bin/prepDE.py -g /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/outputs/Pacu_gene_count_matrix.csv -i /data/putnamlab/Amalia_Murgueitio/P_acuta/Mapping_P_acuta/data/stringtie/sample_list_Pacu_july_prepDE.txt