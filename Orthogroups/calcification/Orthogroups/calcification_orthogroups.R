# Load the readr package
library(readr)
library(tidyverse)

setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/Orthogroups/calcification/Orthogroups")

# Read the Biomin toolkit path file 
annotation_data <- read_csv("C:/Users/amurg/Downloads/Biomineralization Toolkit_2_2021-nov-24.xlsx - data.csv")

# Define file path
file_path <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/Orthogroups/calcification/Orthogroups/Orthogroups.tsv"

# Read the file
orthogroups <- read_tsv(file_path)

# Preview data
glimpse(orthogroups)

#Change the name
orthogroups <- orthogroups %>%
  rename(calcification_gene = `Latest_AmilRS_AdiTakeuchi_SpisPeled_SpisDrake_SpisManju_extras.FIX (5)`)


orthogroups <- orthogroups %>%
  filter(!is.na(calcification_gene) & calcification_gene != "")
#93

# Count how many gene IDs are in each row's calcification_gene column
orthogroups_one_calc_gene <- orthogroups %>%
  filter(!is.na(calcification_gene) & calcification_gene != "") %>%
  filter(sapply(strsplit(calcification_gene, ",| "), function(x) length(x[x != ""])) == 1)

# Clean and rename relevant columns
annotation_data <- annotation_data %>%
  rename(
    calcification_gene = `qyery id`,
    product_stylophora_NCBI = `product stylophora NCBI`
  ) %>%
  select(calcification_gene, product_stylophora_NCBI) %>%
  filter(!is.na(calcification_gene) & calcification_gene != "")


# Join to full orthogroups
orthogroups_annotated <- orthogroups %>%
  left_join(annotation_data, by = "calcification_gene")

# Join to filtered one-gene orthogroups
orthogroups_one_calc_gene_annotated <- orthogroups_one_calc_gene %>%
  left_join(annotation_data, by = "calcification_gene")

#Now filtering by genes present in single copy in 2 sp.
# Define species columns explicitly
species_columns <- c(
  "Pocillopora_acuta_HIv2.genes.pep",
  "Spis.genome.annotation.pep.longest",
  "Montipora_capitata_HIv3.genes.pep",
  "Acropora_tenuis_0.11.maker_post_001.proteins"
)

# Filter orthogroups with a calcification gene AND at least 2 species with a single gene copy
orthogroups_single_copy_2plus <- orthogroups_one_calc_gene_annotated %>%
  rowwise() %>%
  mutate(
    single_copy_species_count = sum(sapply(across(all_of(species_columns)), function(x) {
      genes <- unlist(strsplit(x, ",| "))
      genes <- genes[genes != ""]  # remove empty strings
      length(genes) == 1
    }))
  ) %>%
  ungroup() %>%
  filter(single_copy_species_count >= 2) %>%
  select(-single_copy_species_count)

