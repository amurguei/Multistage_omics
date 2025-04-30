# Load required libraries
library(readr)
library(dplyr)
library(stringr)

# ---- 1. Load your gcount_raw data ----
gcount_raw <- read_csv("gcount_raw.csv")

# ---- 2. Load the Spis annotation file ----
spis_annot <- read_csv("Spis_annot.csv")

# ---- 3. Fix gene_id in Spis_annot: remove "Gene" from SpisGene1234 -> Spis1234
spis_annot <- spis_annot %>%
  mutate(gene_id_fixed = str_replace(gene_id, "Gene", ""))

# ---- 4. Merge with annotations ----
merged_data <- gcount_raw %>%
  left_join(spis_annot, by = c("GeneID_Spis" = "gene_id_fixed"))


# ---- 5. Save the new merged file ----
write.csv(merged_data, "gcount_raw_with_annotations.csv", row.names = FALSE)
