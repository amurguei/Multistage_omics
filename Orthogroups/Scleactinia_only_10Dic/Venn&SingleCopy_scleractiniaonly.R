# 18 November 2024
#Based on the work of Matthew Haas https://github.com/MatthewHaas/genome_assembly/blob/master/venn_diagram_orthogroups_with_rice_relatives.R#L1
# Purpose of this code is to make a venn diagram for overlapping Orthogroups for the NWR genome paper
# This script uses the updated OrthoFinder results that includes Z. latifolia and other Oryza species.
# Done in RStudio

# Set working directory
setwd("C:/Users/amurg/Downloads")

# Load required packages
library(data.table)
library(VennDiagram)
# Since we wanted to add commas as the separator for thousands, I modified the venn.diagram function.
# We need to use source() to load the modified function.
source("modified_venn_diagram.R")

# Read in data
x <- fread("Orthogroups_scleractinia_only.tsv")

# Find all rows for which each species has at least one gene present in the orthogroup
Montipora_capitata_HIv3.genes.pep <- x[Montipora_capitata_HIv3.genes.pep > 0]
Pocillopora_acuta_HIv2.genes.pep <- x[Pocillopora_acuta_HIv2.genes.pep > 0]
Spis.genome.annotation.pep.longest <- x[Spis.genome.annotation.pep.longest > 0]

# Get the orthogroup IDS for each species.
Montipora_capitata_HIv3.genes.pep_orthogroups <- Montipora_capitata_HIv3.genes.pep$Orthogroup
Pocillopora_acuta_HIv2.genes.pep_orthogroups <- Pocillopora_acuta_HIv2.genes.pep$Orthogroup
Spis.genome.annotation.pep.longest_orthogroups <- Spis.genome.annotation.pep.longest$Orthogroup

# Save data
save(Montipora_capitata_HIv3.genes.pep_orthogroups, Pocillopora_acuta_HIv2.genes.pep_orthogroups, Spis.genome.annotation.pep.longest_orthogroups, file="orthogroup_names_for_venn_diagram_update_trial.Rdata")

colors= c("hotpink1", "darkmagenta", "cyan4")

# Load the necessary library
library(VennDiagram)

# Define the orthogroup data for the three species
venn.diagram(
  x = list(
    Montipora_capitata = Montipora_capitata_HIv3.genes.pep_orthogroups,
    Pocillopora_acuta = Pocillopora_acuta_HIv2.genes.pep_orthogroups,
    Spis = Spis.genome.annotation.pep.longest_orthogroups
  ),
  category.names = c(
    expression(atop(italic("Montipora capitata"), plain(13071))),
    expression(atop(italic("Pocillopora acuta"), plain(12957))),
    expression(atop(italic("Stylophora pistillata"), plain(11621)))
  ),
  filename = "venn_diagram_coral_species.png",
  output = TRUE,
  imagetype = "png",
  scaled = FALSE,
  col = "black",
  fill = colors,
  cat.col = colors,
  cat.cex = 0.95,
  cat.dist = 0.35,
  margin = 0.30,
  euler.d = FALSE
)


#Another version with numbers and bigger circles
colors = c("hotpink1", "darkmagenta", "cyan4")

# Ensure the input data for orthogroups are treated as lists
Montipora_capitata_orthogroups <- unique(Montipora_capitata_HIv3.genes.pep_orthogroups)
Pocillopora_acuta_orthogroups <- unique(Pocillopora_acuta_HIv2.genes.pep_orthogroups)
Stylophora_pistillata_orthogroups <- unique(Spis.genome.annotation.pep.longest_orthogroups)

# Check lengths
length_Montipora <- length(Montipora_capitata_orthogroups)
length_Pocillopora <- length(Pocillopora_acuta_orthogroups)
length_Stylophora <- length(Stylophora_pistillata_orthogroups)

# Generate the Venn diagram
venn.diagram(
  x = list(
    Montipora_capitata = Montipora_capitata_orthogroups,
    Pocillopora_acuta = Pocillopora_acuta_orthogroups,
    Stylophora_pistillata = Stylophora_pistillata_orthogroups
  ),
  category.names = c(
    expression(atop(italic("Montipora capitata"), plain(format(length_Montipora, big.mark = ",")))),
    expression(atop(italic("Pocillopora acuta"), plain(format(length_Pocillopora, big.mark = ",")))),
    expression(atop(italic("Stylophora pistillata"), plain(format(length_Stylophora, big.mark = ","))))
  ),
  filename = "venn_diagram_mod.png",
  output = TRUE,
  imagetype = "png",
  scaled = FALSE,
  col = "black",
  fill = colors,
  cat.col = colors,
  cat.cex = 0.95,
  cat.dist = 0.35,
  margin = 0.30,
  euler.d = FALSE
)

str(Montipora_capitata_HIv3.genes.pep_orthogroups)
str(Pocillopora_acuta_HIv2.genes.pep_orthogroups)
str(Spis.genome.annotation.pep.longest_orthogroups)


#Single copy orthologues
# Load necessary libraries
library(tidyverse)

# Load the Orthogroups.tsv file
orthogroups <- read.delim("Orthogroups_scleractinia_only.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
x <- fread("C:/Users/amurgueitio/Documents/Multistage_Omics/Orthogroups/Orthogroups.tsv")

library(dplyr)
library(stringr)

# Step 1: Clean and split genes for each species
cleaned_orthogroups <- orthogroups %>%
  mutate(
    Montipora_clean = str_split(Montipora_capitata_HIv3.genes.pep, ",\\s*"),
    Pocillopora_clean = str_split(Pocillopora_acuta_HIv2.genes.pep, ",\\s*"),
    Spis_clean = str_split(Spis.genome.annotation.pep.longest, ",\\s*")
  )

# Step 2: Filter rows where each species has exactly one gene
single_copy <- cleaned_orthogroups %>%
  filter(
    lengths(Montipora_clean) == 1, # Only one gene for Montipora
    lengths(Pocillopora_clean) == 1, # Only one gene for Pocillopora
    lengths(Spis_clean) == 1        # Only one gene for Spis
  )

# Step 3: Extract single-copy genes into new columns
single_copy <- single_copy %>%
  mutate(
    Montipora_gene = sapply(Montipora_clean, `[`, 1),
    Pocillopora_gene = sapply(Pocillopora_clean, `[`, 1),
    Spis_gene = sapply(Spis_clean, `[`, 1)
  ) %>%
  select(Orthogroup, Montipora_gene, Pocillopora_gene, Spis_gene) # Select only relevant columns

# Step 4: Save the filtered single-copy orthogroups to a file
write.table(single_copy, "Single_Copy_Orthogroups_scleractinia.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Print a summary of the results
cat("Filtered single-copy orthogroups:\n")
print(single_copy)
cat("\nNumber of single-copy orthogroups:", nrow(single_copy), "\n")


# Remove rows where any gene column has an empty value
filtered_single_copy <- single_copy %>%
  filter(
    !is.na(Montipora_gene) & Montipora_gene != "",
    !is.na(Pocillopora_gene) & Pocillopora_gene != "",
    !is.na(Spis_gene) & Spis_gene != ""
  )

# Save the cleaned dataset
write.table(filtered_single_copy, "Cleaned_Single_Copy_Orthogroups.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Save just the orthogroup names to a text file
writeLines(filtered_single_copy$Orthogroup, "Single_Copy_Orthogroups.txt")

# Print a summary
cat("Cleaned single-copy orthogroups dataset saved as 'Cleaned_Single_Copy_Orthogroups.tsv'\n")
cat("Orthogroup names saved as 'Single_Copy_Orthogroups.txt'\n")
cat("\nNumber of cleaned single-copy orthogroups:", nrow(filtered_single_copy), "\n")
