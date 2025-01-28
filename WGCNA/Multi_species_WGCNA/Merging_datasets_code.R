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

# Add this column to your final_orthogroups dataframe
final_orthogroups$final_orthogroups_no_suffixes <- final_orthogroups_no_suffixes
View(final_orthogroups)
View(species4_counts)
head(matching_genes)
#Check how many matched genes there are
length(matching_genes)
#7128

#Now I'm changing the name of final_orthogroups_no_suffixes to something more intuitive
# Rename the column in the final_orthogroups dataframe
colnames(final_orthogroups)[colnames(final_orthogroups) == "final_orthogroups_no_suffixes"] <- "Spis.genome.annotation.pep.longest_no_suffixes"

# Verify the change
colnames(final_orthogroups)

#Let's see now if I can merge the Acropora tenuis data orthogroups with the 

# Step 1: Prepare for the join
# Create a dataframe for species1_counts with gene IDs as a column
species1_counts_df <- data.frame(GeneID = rownames(species1_counts), species1_counts)

#Removing the additional columns from this dataset: 
# Remove columns 1 to 6 from species1_counts_df
species1_counts_df <- species1_counts_df[, -c(2:6)]

# Verify the updated dataframe
head(species1_counts_df)


# Step 3: Left join using GeneID as key
merged_data_acropora <- merge(species1_counts_df, final_orthogroups, by.x = "GeneID", by.y = "Acropora_tenuis_0.11.maker_post_001.proteins", all.x = TRUE)

# Step 4: Check the results
head(merged_data_acropora)

# Overwrite merged_data_acropora to keep only rows where Orthogroup is not NA
merged_data_acropora <- merged_data_acropora[!is.na(merged_data_acropora$Orthogroup), ]

# Check the filtered results
head(merged_data_acropora)

# Reorder the columns to move Orthogroup to the first position
merged_data_acropora <- merged_data_acropora[, c("Orthogroup", setdiff(names(merged_data_acropora), "Orthogroup"))]

# Check the updated dataframe
head(merged_data_acropora)

# Rename GeneID column to GeneID_Atenuis
colnames(merged_data_acropora)[colnames(merged_data_acropora) == "GeneID"] <- "GeneID_Atenuis"

# Check the updated dataframe
head(merged_data_acropora)

#Let's move the other columns to the front

# Identify the columns to be moved
columns_to_move <- c("Montipora_capitata_HIv3.genes.pep", 
                     "Pocillopora_acuta_HIv2.genes.pep", 
                     "Spis.genome.annotation.pep.longest", 
                     "Spis.genome.annotation.pep.longest_no_suffixes")

# Reorder columns
merged_data_acropora <- merged_data_acropora[, c(
  colnames(merged_data_acropora)[1:2],  # Keep the first two columns
  columns_to_move,                      # Add the selected columns here
  setdiff(colnames(merged_data_acropora)[3:ncol(merged_data_acropora)], columns_to_move) # Add the remaining columns
)]

# Check the updated dataframe
head(merged_data_acropora)

#Amazing, now let's try joining Montipora too

# Step 1: Create a dataframe for species2_counts with gene IDs as a column
species2_counts_df <- data.frame(GeneID_Mcap = rownames(species2_counts), species2_counts)

# Step 2: Left join species2_counts with final_orthogroups
merged_data_acropora_montipora <- merge(
  merged_data_acropora,                          # Start with Acropora merged data
  species2_counts_df,                            # Add Montipora data
  by.x = "Montipora_capitata_HIv3.genes.pep",    # Match Montipora column from final_orthogroups
  by.y = "GeneID_Mcap",                          # Use renamed Montipora GeneID column
  all.x = TRUE                                   # Keep only rows with Acropora orthogroups
)

# Step 3: Remove rows where Orthogroup is NA
merged_data_acropora_montipora <- merged_data_acropora_montipora[!is.na(merged_data_acropora_montipora$Orthogroup), ]

# Step 4: Reorder columns to move "Orthogroup" and Montipora gene data
columns_to_move <- c("Orthogroup", "GeneID_Mcap", colnames(species2_counts)) # Bring in species2 columns to the front
merged_data_acropora_montipora <- merged_data_acropora_montipora[, c(
  columns_to_move,                                      # Orthogroup and species2 data
  setdiff(colnames(merged_data_acropora_montipora), columns_to_move)  # The rest of the columns
)]

# Step 5: Check results
head(merged_data_acropora_montipora)

# Step 6 Bring Orthogroup to the first column and rename the GeneID column for Montipora
merged_data_acropora_montipora <- merged_data_acropora_montipora %>%
  select(Orthogroup, everything()) %>%
  rename(GeneID_Mcap = Montipora_capitata_HIv3.genes.pep)

# Step 7 Check the merging results
head(merged_data_acropora_montipora)

# Reorder the columns: Orthogroup, GeneID_Atenuis, GeneID_Mcap, then the rest
merged_data_acropora_montipora <- merged_data_acropora_montipora %>%
  select(Orthogroup, GeneID_Atenuis, GeneID_Mcap, everything())

# Check the final results
head(merged_data_acropora_montipora)

# Pocillopora merger

# Step 1: Prepare Pocillopora data for the merge
# Create a dataframe for Pocillopora counts with GeneIDs as a column
species3_counts_df <- data.frame(GeneID_Pacu = rownames(species3_counts), species3_counts)

# Step 2: Left join Pocillopora data with the merged_data_acropora_montipora that has orthogroups
merged_data_acropora_montipora_pocillopora <- merge(
  merged_data_acropora_montipora,                # Start with the merged Acropora and Montipora data
  species3_counts_df,                            # Add Pocillopora data
  by.x = "Pocillopora_acuta_HIv2.genes.pep",     # Pocillopora column from final_orthogroups
  by.y = "GeneID_Pacu",                          # Use renamed Pocillopora GeneID column
  all.x = TRUE                                   # Keep only rows with the existing orthogroups
)


# Step 4: Rename the Pocillopora gene column if necessary
colnames(merged_data_acropora_montipora_pocillopora)[which(colnames(merged_data_acropora_montipora_pocillopora) == "Pocillopora_acuta_HIv2.genes.pep")] <- "GeneID_Pacu"

# Step 5: Reorder columns to have Orthogroup, GeneID_Atenuis, GeneID_Mcap, GeneID_Pacu
merged_data_acropora_montipora_pocillopora <- merged_data_acropora_montipora_pocillopora %>%
  select(Orthogroup, GeneID_Atenuis, GeneID_Mcap, GeneID_Pacu, everything())

# Check the final results
head(merged_data_acropora_montipora_pocillopora)

#Now let's try Spis and see if the merger allows to identify the missing gene
# Step 1: Prepare Stylophora (Species 4) data
species4_counts_df <- data.frame(GeneID_Spis = rownames(species4_counts), species4_counts)

# Step 2: Left join Stylophora data with the merged_data_acropora_montipora_pocillopora that has orthogroups
merged_data_acropora_montipora_pocillopora_spis <- merge(
  merged_data_acropora_montipora_pocillopora,                      # Start with the merged Acropora, Montipora, and Pocillopora data
  species4_counts_df,                                              # Add Stylophora data
  by.x = "Spis.genome.annotation.pep.longest_no_suffixes",          # Stylophora column from final_orthogroups
  by.y = "GeneID_Spis",                                            # Use renamed Stylophora GeneID column
  all.x = TRUE                                                     # Keep only rows with existing orthogroups
)

# Step 3: Rename the Stylophora gene column if necessary
colnames(merged_data_acropora_montipora_pocillopora_spis)[which(colnames(merged_data_acropora_montipora_pocillopora_spis) == "Spis.genome.annotation.pep.longest_no_suffixes")] <- "GeneID_Spis"

# Step 4: Reorder columns to have Orthogroup, GeneID_Atenuis, GeneID_Mcap, GeneID_Pacu, GeneID_Spis
merged_data_acropora_montipora_pocillopora_spis <- merged_data_acropora_montipora_pocillopora_spis %>%
  select(Orthogroup, GeneID_Atenuis, GeneID_Mcap, GeneID_Pacu, GeneID_Spis, everything())

# Check the final results
head(merged_data_acropora_montipora_pocillopora_spis)

#AMAZING, now let's take a look 

# Check for NAs in the whole dataset
sum(is.na(merged_data_acropora_montipora_pocillopora_spis))

# Check how many NAs are there per column
colSums(is.na(merged_data_acropora_montipora_pocillopora_spis))

# Check the rows where there are NAs (if any)
na_rows <- merged_data_acropora_montipora_pocillopora_spis[!complete.cases(merged_data_acropora_montipora_pocillopora_spis), ]
head(na_rows)

# Check for NAs in the whole dataset
sum(is.na(merged_data_acropora_montipora_pocillopora_spis))

# Check how many NAs are there per column
colSums(is.na(merged_data_acropora_montipora_pocillopora_spis))

# Check the rows where there are NAs (if any)
na_rows <- merged_data_acropora_montipora_pocillopora_spis[!complete.cases(merged_data_acropora_montipora_pocillopora_spis), ]
head(na_rows)

#I was able to identify the NA :):):) I'm going to fix it manually since it's so much easier than through code
#Orthogroup       GeneID_Atenuis                                GeneID_Mcap                           GeneID_Pacu  GeneID_Spis Spis.genome.annotation.pep.longest
#422  OG0007064 aten_0.1.m1.17048.m1 Montipora_capitata_HIv3___RNAseq.g12926.t1 Pocillopora_acuta_HIv2___TS.g21030.t2 Spis11155.t4                       Spis11155.t4
#DRR318288 DRR318289 DRR318291 DRR318293 DRR318294 DRR318295 DRR318297 DRR318298 DRR318299 AH1 AH2 AH3 AH4 AH5 AH6  AH7  AH8  AH9 SRR3051863 SRR3051864 SRR3051865
#422       761       948       564       944      1745       664      5536      5255      4362 728 874 890 989 761 899 1982 2156 2093        248        319        411
#SRR3051866 SRR3051867 SRR3051868 SRR3051869 SRR3051870 SRR3051871 SRR14333319 SRR14333320 SRR14333321 SRR14333322 SRR14333323 SRR14333324 SRR14333325 SRR14333326
#422       1847       2322        927        166        139        172          NA          NA          NA          NA          NA          NA          NA          NA
#SRR14333327
#422          NA

#Now we are going to save everything I need from here
# Save the merged dataset to a CSV file
write.csv(merged_data_acropora_montipora_pocillopora_spis, 
          "merged_data_acropora_montipora_pocillopora_spis.csv", 
          row.names = FALSE)


# Save the final_orthogroups dataframe to a CSV file
write.csv(final_orthogroups, 
          "final_orthogroups.csv", 
          row.names = FALSE)

