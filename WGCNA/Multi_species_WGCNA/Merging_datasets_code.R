library(dplyr)
library(tidyverse)
library(readr)

setwd("C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA")

# Load the Orthogroups.tsv file
orthogroups <- read.delim("Orthogroups.tsv")

# Load the Orthogroups_SingleCopyOrthologues.txt file. Note: this approach didn't work, keeping it for the records
single_copy_orthogroups <- read.delim("Orthogroups_SingleCopyOrthologues.txt", header = FALSE)

# Check column names and first few rows
head(orthogroups)
#1  OG0000000
#2  OG0000001
#3  OG0000002
#4  OG0000003
#5  OG0000004
#6  OG0000005
head(single_copy_orthogroups)
#1 N0.HOG0000109
#2 N0.HOG0000235
#3 N0.HOG0000321
#4 N0.HOG0000326
#5 N0.HOG0000335
#6 N0.HOG0000338

#There is a naming discrepancy, so the first thing I will do is to remove the "NO.HO" from the single_copy ortologues
# Remove the "N0.HOG" prefix from the single_copy_orthogroups column
single_copy_orthogroups$V1 <- sub("^N0\\.H", "", single_copy_orthogroups$V1)

# Now filter the orthogroups to include only those present in the single_copy_orthogroups
filtered_orthogroups <- orthogroups[orthogroups$Orthogroup %in% single_copy_orthogroups$V1, ]
#There's something strange about this approach, not everything is actually a single_copy_ortologue... Trying it manually

library(dplyr)

# Assuming you have a data frame `orthogroups` with these dataset names:
# Orthogroup, Acropora_tenuis_0.11.maker_post_001.proteins, 
# Montipora_capitata_HIv3.genes.pep, Pocillopora_acuta_HIv2.genes.pep,
# Spis.genome.annotation.pep.longest

library(dplyr)

# Check a few rows in the dataset to confirm the structure and spot any potential issues:
head(orthogroups)

# 1. Check the count of non-NA values in each species column to diagnose any unexpected issues
orthogroups_counts <- orthogroups %>%
  summarise(
    Acropora_count = sum(!is.na(Acropora_tenuis_0.11.maker_post_001.proteins)),
    Montipora_count = sum(!is.na(Montipora_capitata_HIv3.genes.pep)),
    Pocillopora_count = sum(!is.na(Pocillopora_acuta_HIv2.genes.pep)),
    Spis_count = sum(!is.na(Spis.genome.annotation.pep.longest))
  )

# Print counts to check if the dataset has genes in the expected columns
print(orthogroups_counts)

library(dplyr)

# Clean gene list for each species and handle commas
orthogroups_clean <- orthogroups %>%
  mutate(
    # Split genes by commas, remove whitespace around each gene, and ensure non-empty genes
    Acropora_clean = sapply(strsplit(as.character(Acropora_tenuis_0.11.maker_post_001.proteins), ","), function(x) unique(trimws(x[x != ""]))),
    Montipora_clean = sapply(strsplit(as.character(Montipora_capitata_HIv3.genes.pep), ","), function(x) unique(trimws(x[x != ""]))),
    Pocillopora_clean = sapply(strsplit(as.character(Pocillopora_acuta_HIv2.genes.pep), ","), function(x) unique(trimws(x[x != ""]))),
    Spis_clean = sapply(strsplit(as.character(Spis.genome.annotation.pep.longest), ","), function(x) unique(trimws(x[x != ""]))))

# Count the number of unique genes for each species
orthogroups_clean <- orthogroups_clean %>%
  mutate(
    Acropora_count = sapply(Acropora_clean, length),
    Montipora_count = sapply(Montipora_clean, length),
    Pocillopora_count = sapply(Pocillopora_clean, length),
    Spis_count = sapply(Spis_clean, length)
  )

# Filter orthogroups that have exactly one gene in each dataset
filtered_orthogroups <- orthogroups_clean %>%
  filter(Acropora_count == 1 & Montipora_count == 1 & Pocillopora_count == 1 & Spis_count == 1)

# View the filtered orthogroups
head(filtered_orthogroups)

# Finalize the filtered list for orthogroups that meet the condition
final_orthogroups <- orthogroups %>%
  filter(Orthogroup %in% filtered_orthogroups$Orthogroup)

# View the final filtered orthogroups
head(final_orthogroups)

# Read in the gene count matrices for the four species
species1_counts <- read.csv("gene_count_matrix_Acropora.csv", row.names = 1)
species2_counts <- read.csv("Mcap_gene_count_matrix.csv", row.names = 1)
species3_counts <- read.csv("Pacu_gene_count_matrix_newGFF.csv", row.names = 1)
species4_counts <- read.csv("4-Spis-GeneCountMatrix.csv", row.names = 1)

# View the row names of the gene count matrices to identify discrepancies
head(rownames(species1_counts))
#[1] "aten_0.1.m1.1" "aten_0.1.m1.2" "aten_0.1.m1.3" "aten_0.1.m1.4" "aten_0.1.m1.5" "aten_0.1.m1.6"
head(rownames(species2_counts))
#[1] "Montipora_capitata_HIv3___RNAseq.g42319.t1" "Montipora_capitata_HIv3___TS.g49315.t1b"    "Montipora_capitata_HIv3___TS.g49315.t1a"   
#[4] "Montipora_capitata_HIv3___TS.g11988.t1"     "Montipora_capitata_HIv3___RNAseq.g29525.t1" "Montipora_capitata_HIv3___RNAseq.g45785.t1"
head(rownames(species3_counts))
#[1] "Pocillopora_acuta_HIv2___RNAseq.g23616.t1" "Pocillopora_acuta_HIv2___RNAseq.g4197.t1"  "Pocillopora_acuta_HIv2___RNAseq.g12137.t1"
#[4] "Pocillopora_acuta_HIv2___RNAseq.g4145.t1"  "Pocillopora_acuta_HIv2___RNAseq.g4784.t1"  "Pocillopora_acuta_HIv2___TS.g3788.t1"
head(rownames(species4_counts))
#[1] "SpisGene3669" "SpisGene3668" "SpisGene3666" "SpisGene3665" "SpisGene3664" "SpisGene3663"

#Fixing naming inconsistencies. 

# Add ".m1" suffix to gene names in species 1
rownames(species1_counts) <- paste0(rownames(species1_counts), ".m1")

# Proceed with analysis after modifying row names
# Check the updated row names of species 1
head(rownames(species1_counts))

#let's check if fixed: 
# Extract gene names from final_orthogroups (assuming they are in a column like 'Acropora_tenuis_0.11.maker_post_001.proteins')
final_gene_names_acropora <- unique(final_orthogroups$Acropora_tenuis_0.11.maker_post_001.proteins)

# Check which gene names in species1 match final_orthogroups
matching_genes_acropora <- intersect(rownames(species1_counts), final_gene_names_acropora)

# Check if the gene names match, and see how many match
length(matching_genes_acropora)
head(matching_genes_acropora)  # Display a few of the matching gene names

#now let's remove one sample per species in species 1, since all the others have n=3 per life stage.
# Remove unwanted samples from species1_counts
species1_counts <- species1_counts[, !(colnames(species1_counts) %in% c("DRR318292", "DRR318296", "DRR318290"))]

# Check if there are 3 samples left per life stage (adjust depending on column structure)
table(gsub("_[^_]+$", "", colnames(species1_counts)))
#   Chr DRR318288 DRR318289 DRR318291 DRR318293 DRR318294 DRR318295 DRR318297 DRR318298 DRR318299       End    Length     Start    Strand 
#1         1         1         1         1         1         1         1         1         1         1         1         1         1 

#Checking all good for Montipora

# Step 1: Extract gene names for Montipora from final_orthogroups
final_gene_names_montipora <- unique(final_orthogroups$Montipora_capitata_HIv3.genes.pep)

# Step 2: Check which gene names in species2 (Montipora) match final_orthogroups
matching_genes_montipora <- intersect(rownames(species2_counts), final_gene_names_montipora)

# Step 3: Check how many gene names match and display a few of them
length(matching_genes_montipora)
#7129
head(matching_genes_montipora)
#head(matching_genes_montipora)
#[1] "Montipora_capitata_HIv3___RNAseq.g36841.t1" "Montipora_capitata_HIv3___RNAseq.g15942.t1" "Montipora_capitata_HIv3___TS.g7970.t2"     
#[4] "Montipora_capitata_HIv3___RNAseq.g41686.t1" "Montipora_capitata_HIv3___RNAseq.42542_t"   "Montipora_capitata_HIv3___RNAseq.g9640.t1" 

#Looking good, now Pocillopora

# Step 1: Extract gene names for Pocillopora from final_orthogroups
final_gene_names_pocillopora <- unique(final_orthogroups$Pocillopora_acuta_HIv2.genes.pep)

# Step 2: Check which gene names in species3 (Pocillopora) match final_orthogroups
matching_genes_pocillopora <- intersect(rownames(species3_counts), final_gene_names_pocillopora)

# Step 3: Check how many gene names match and display a few of them
length(matching_genes_pocillopora)
#7129
head(matching_genes_pocillopora)
#[1] "Pocillopora_acuta_HIv2___RNAseq.g27841.t1" "Pocillopora_acuta_HIv2___RNAseq.g14011.t1" "Pocillopora_acuta_HIv2___RNAseq.g7479.t3b"
#[4] "Pocillopora_acuta_HIv2___RNAseq.g3340.t1"  "Pocillopora_acuta_HIv2___RNAseq.g25357.t1" "Pocillopora_acuta_HIv2___RNAseq.g28377.t1"
#Looking good!

#Now need to modify Stylophora. Stylophora's problem is the gene count matrix says "SpisGene", rather than "Spis"

# Step 1: Extract the gene number part of the gene names in species4_counts
rownames(species4_counts) <- sub("SpisGene", "Spis", rownames(species4_counts))

# Step 2: Extract gene names from final_orthogroups for Stylophora (without 'Gene')
final_gene_names_stylophora <- unique(final_orthogroups$Spis.genome.annotation.pep.longest)

# Step 3: Check which gene names in species4 (Stylophora) match final_orthogroups
matching_genes_stylophora <- intersect(rownames(species4_counts), final_gene_names_stylophora)

# Step 4: Check how many gene names match and display a few of them
length(matching_genes_stylophora)
#[1] 6349
head(matching_genes_stylophora)
#[1] "Spis3665"  "Spis3663"  "Spis3660"  "Spis16164" "Spis17697" "Spis17349"

#So seems like this fixed most of the issue, but still some mismatch, let's investigate further

# Step 1: Identify unmatched genes
unmatched_genes_stylophora <- setdiff(rownames(species4_counts), final_gene_names_stylophora)

# Step 2: Check how many unmatched genes there are
length(unmatched_genes_stylophora)

# Step 3: Display a few of the unmatched gene names
head (unmatched_genes_stylophora)

# Remove ".t1", ".t2", ".t3" from gene names in final_orthogroups
final_orthogroups_no_suffixes <- gsub("\\.t[123]$", "", final_orthogroups$Spis.genome.annotation.pep.longest)
# Add this column to your final_orthogroups dataframe
final_orthogroups$final_orthogroups_no_suffixes <- final_orthogroups_no_suffixes
# Check the first few genes after removing suffixes
head(final_orthogroups_no_suffixes)

# Extract gene names from the species count matrix (Stylophora in this case)
genes_in_count_matrix <- rownames(species4_counts)


# Step 1: Extract Montipora orthogroups (in case final_orthogroups file stores them in a specific column for this species)
montipora_orthogroups_column <- final_orthogroups$Montipora_capitata_HIv3.genes.pep  # Adjust if necessary
montipora_orthogroups <- unique(montipora_orthogroups_column)

# Step 2: Extract Stylophora orthogroups and remove suffixes from the gene names (remove ".t1", ".t2", etc.)
stylophora_orthogroups <- unique(final_orthogroups$final_orthogroups_no_suffixes) 
stylophora_orthogroups_no_suffixes <- gsub("\\.t\\d+$", "", stylophora_orthogroups)

# Step 3: Identify the orthogroups in Montipora that are NOT present in Stylophora
missing_orthogroups_stylophora <- setdiff(montipora_orthogroups, stylophora_orthogroups_no_suffixes)

# Create a dataframe for species1_counts with gene IDs as a column
species1_counts_df <- data.frame(GeneID = rownames(species1_counts), species1_counts)
