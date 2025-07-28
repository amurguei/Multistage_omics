#Hive: 
ssh amalia@hive02.haifa.ac.il
cd /lustre1/home/mass/amalia/orthofinder4sp
mkdir calcification

#Copy protein files from to subdirectory
cp /lustre1/home/mass/amalia/orthofinder4sp/{Spis.genome.annotation.pep.longest.fa,Pocillopora_acuta_HIv2.genes.pep.faa,Montipora_capitata_HIv3.genes.pep.faa,Acropora_tenuis_0.11.maker_post_001.proteins.fasta} /lustre1/home/mass/amalia/orthofinder4sp/calcification/
 cp

 nano orthofinder_calc

#Run orthofinder
#!/bin/bash
#SBATCH --job-name="Orthofinder_calc"
#SBATCH --partition=hiveunlim
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amalia.murgueitio@gmail.com
#SBATCH --output=orthofinder_calc_%j.out
#SBATCH --error=orthofinder_calc_%j.err
#SBATCH -D /lustre1/home/mass/amalia/orthofinder4sp/calcification

# Load Conda and activate environment
source ~/.bashrc
conda activate orthofinder

# Run OrthoFinder on the calcification folder
orthofinder -t 12 -a 12 -f /lustre1/home/mass/amalia/orthofinder4sp/calcification

# Deactivate Conda environment
conda deactivate

#Download to my PC the 4 sp + calc genes files:
scp -r amalia@hive02.haifa.ac.il:/lustre1/home/mass/amalia/orthofinder4sp/calcification/OrthoFinder/Results_Jul22_1/Orthogroups "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/Orthogroups/calcification"

#Now doing Spis only..

cp /lustre1/home/mass/amalia/orthofinder4sp/calcification/'Latest_AmilRS_AdiTakeuchi_SpisPeled_SpisDrake_SpisManju_extras.FIX (5).fasta' \
   /lustre1/home/mass/amalia/orthofinder4sp/Spis.genome.annotation.pep.longest.fa \
   /lustre1/home/mass/amalia/orthofinder4sp/calcification/spis_only/

nano spis_calc.sh

#Run orthofinder

#!/bin/bash
#SBATCH --job-name="Orthofinder_calc"
#SBATCH --partition=hiveunlim
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=amalia.murgueitio@gmail.com
#SBATCH --output=orthofinder_calc_%j.out
#SBATCH --error=orthofinder_calc_%j.err
#SBATCH -D /lustre1/home/mass/amalia/orthofinder4sp/calcification/spis_only

# Load Conda and activate environment
source ~/.bashrc
conda activate orthofinder

# Run OrthoFinder on the calcification folder
orthofinder -t 12 -a 12 -f /lustre1/home/mass/amalia/orthofinder4sp/calcification/spis_only

# Deactivate Conda environment
conda deactivate


