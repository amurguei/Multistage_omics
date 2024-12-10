# Load necessary libraries
library(tidyverse)

# Load the Orthogroups.tsv file
orthogroups <- read.delim("Orthogroups.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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
write.table(single_copy, "Single_Copy_Orthogroups.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

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
