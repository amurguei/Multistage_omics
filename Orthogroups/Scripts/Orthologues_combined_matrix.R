# Install if necessary
# install.packages("dplyr")
# install.packages("tidyverse")
# install.packages("readr")

library(dplyr)
library(tidyverse)
library(readr)

setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/Orthogroups")

# Load the Orthogroups.tsv file
orthogroups <- read.delim("Orthogroups.tsv")

# Load the Orthogroups_SingleCopyOrthologues.txt file
single_copy_orthogroups <- read.delim("Orthogroups_SingleCopyOrthologues.txt", header = FALSE)


# Filter the orthogroups to include only single-copy orthologues
filtered_orthogroups <- orthogroups[orthogroups$Orthogroup %in% single_copy_orthogroups$V1, ]

# Read in the gene count matrices for the three species
species1_counts <- read.csv("Mcap_gene_count_matrix.csv", row.names = 1)
species2_counts <- read.csv("Pacu_gene_count_matrix_newGFF.csv", row.names = 1)
species3_counts <- read.csv("4-Spis-GeneCountMatrix.csv", row.names = 1)

# Remove isoform suffixes (e.g., .t1, .t2) from gene names in species3_counts
rownames(species3_counts) <- gsub("\\.t[0-9]+$", "", rownames(species3_counts))


# Rename columns to more convenient names
colnames(filtered_orthogroups) <- c("Orthogroup", "Mcap", "Pacu", "Spis")

# Now, create a dataframe with the orthologue name and the corresponding genes
single_copy_genes <- filtered_orthogroups %>%
  select(Orthogroup, Mcap, Pacu, Spis) %>%
  rename(Species1 = Mcap, Species2 = Pacu, Species3 = Spis)
# Remove the existing 'Spis' prefix if present
single_copy_genes$Species3 <- gsub("^Spis", "", single_copy_genes$Species3)

# Add 'SpisGene' prefix to the cleaned gene names
single_copy_genes$Species3 <- paste0("SpisGene", single_copy_genes$Species3)
# Remove isoform suffixes from gene names in single_copy_genes
single_copy_genes$Species3 <- gsub("\\.t[0-9]+$", "", single_copy_genes$Species3)


# Subset gene count matrices using the filtered single_copy_genes dataframe
species1_subset <- species1_counts[rownames(species1_counts) %in% single_copy_genes$Species1, ]
species2_subset <- species2_counts[rownames(species2_counts) %in% single_copy_genes$Species2, ]
species3_subset <- species3_counts[rownames(species3_counts) %in% single_copy_genes$Species3, ]

# Rename rownames in each subset to Orthogroup
rownames(species1_subset) <- single_copy_genes$Orthogroup[match(rownames(species1_subset), single_copy_genes$Species1)]
rownames(species2_subset) <- single_copy_genes$Orthogroup[match(rownames(species2_subset), single_copy_genes$Species2)]
rownames(species3_subset) <- single_copy_genes$Orthogroup[match(rownames(species3_subset), single_copy_genes$Species3)]

# Combine the matrices by row names (Orthogroup)
combined_matrix <- bind_rows(
  species1_subset %>% mutate(Species = "Species1"),
  species2_subset %>% mutate(Species = "Species2"),
  species3_subset %>% mutate(Species = "Species3")
)

# Reshape to have genes as rows and samples as columns
combined_matrix <- combined_matrix %>%
  pivot_longer(-c(Species), names_to = "Sample", values_to = "Counts") %>%
  pivot_wider(names_from = Sample, values_from = Counts)
