library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Load Data
data <- read_csv("GitHub/Multistage_omics/GPCRS/GPCRS_data.csv")


# Convert to long format
data_long <- data %>%
  pivot_longer(cols = starts_with("Larva") | starts_with("Meta") | starts_with("Spat"),
               names_to = "sample", values_to = "expression")

# Extract species information from column names
data_long <- data_long %>%
  mutate(species = case_when(
    str_detect(sample, "A_tenuis") ~ "A_tenuis",
    str_detect(sample, "M_capitata") ~ "M_capitata",
    str_detect(sample, "P_acuta") ~ "P_acuta",
    str_detect(sample, "S_pistillata") ~ "S_pistillata"
  ))

# Log transform expression (optional, uncomment if needed)
# data_long <- data_long %>% mutate(expression = log1p(expression))

# Compute Z-scores per species and gene
data_z <- data_long %>%
  group_by(Orthogroup, species) %>%
  mutate(z_score = scale(expression)) %>%
  pivot_wider(names_from = sample, values_from = z_score)

# Convert to matrix for heatmap
gene_matrix <- as.matrix(data_z[, -c(1,2)])
rownames(gene_matrix) <- data_z$Orthogroup

# Define heatmap colors
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Plot heatmap
Heatmap(gene_matrix, name = "Z-score", col = col_fun,
        cluster_rows = TRUE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = TRUE)

