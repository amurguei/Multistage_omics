# Load required libraries
library(DESeq2)
library(sva)
library(pheatmap)
library(tibble)
library(dplyr)


# ---- 1. Load your raw count matrix ----
# Assuming GPCRS_data is already in the environment as your raw count matrix
# We'll use it directly instead of reloading

GPCRS_data <- read_csv("GitHub/Multistage_omics/GPCRS/GPCRS_data.csv")

# ---- 2. Remove non-count columns if present ----
# Keep only numeric columns (filter out metadata-like columns)
GPCRS_data <- GPCRS_data[, sapply(GPCRS_data, is.numeric)]

# ---- 3. Save original NA positions and replace with 0s ----
na_mask <- is.na(GPCRS_data)  # Save where NAs were originally
GPCRS_data[na_mask] <- 0

# ---- 4. Prepare metadata from column names ----
sample_names <- colnames(GPCRS_data)

# Extract stage and species from sample names
stage <- sub("([A-Za-z]+)[0-9]+_.*", "\\1", sample_names)
species <- sub(".*_([A-Za-z]+_[a-z]+|Spis|Pacu)$", "\\1", sample_names)

# Create metadata dataframe
sample_info <- data.frame(
  sampleID = sample_names,
  stage = stage,
  species = species,
  group = paste(species, stage, sep = "_"),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sampleID

# ---- 5. Create DESeq2 object ----
dds <- DESeqDataSetFromMatrix(countData = round(GPCRS_data),
                              colData = sample_info,
                              design = ~ 1)  # No design needed for VST

# ---- 6. Apply ComBat-seq ----
combat_counts <- ComBat_seq(counts = assay(dds), batch = sample_info$species)

# ---- 7. Create DESeq2 object from corrected counts ----
dds_combat <- DESeqDataSetFromMatrix(countData = combat_counts,
                                     colData = sample_info,
                                     design = ~ 1)

# ---- 8. Apply VST ----
vst_mat <- assay(vst(dds_combat, blind = TRUE))


# ---- 7. Apply VST directly on the ComBat-corrected matrix ----
# Use varianceStabilizingTransformation instead of vst wrapper
dds_combat <- DESeqDataSetFromMatrix(countData = combat_counts,
                                     colData = sample_info,
                                     design = ~ 1)

dds_combat <- estimateSizeFactors(dds_combat)

# ---- 7b. Inspect size factors ----
size_factors <- sizeFactors(dds_combat)
print(size_factors)


vst_obj <- varianceStabilizingTransformation(dds_combat, blind = TRUE)
vst_mat <- assay(vst_obj)


# ---- 8a. Z-score genes (rows) across all samples (global) ----
z_scores_global <- t(scale(t(vst_mat)))
z_scores_global[na_mask] <- NA

# ---- 8b. Z-score genes (rows) within species ----
z_scores_species <- vst_mat
for (sp in unique(sample_info$species)) {
  sp_cols <- sample_info$sampleID[sample_info$species == sp]
  z_scores_species[, sp_cols] <- t(scale(t(vst_mat[, sp_cols])))
}
z_scores_species[na_mask] <- NA

# ---- 9. Reorder samples by stage and species ----
ordered_samples <- sample_info %>% 
  mutate(stage = factor(stage, levels = c("Larva", "Meta", "Spat")),
         species = factor(species, levels = c("A_tenuis", "M_capitata", "P_acuta", "S_pistillata"))) %>%
  arrange(stage, species)

z_scores_ordered <- z_scores_global[, ordered_samples$sampleID]
sample_info <- sample_info[ordered_samples$sampleID, ]

# ---- 10. Heatmap with fixed sample order and grouped annotation ----
pheatmap(z_scores_ordered,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = sample_info[, c("species", "stage", "family")],
         fontsize_col = 7,
         angle_col = 45,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         na_col = "gray",
         main = "Z-score across samples (ordered)")



# Load required libraries
library(DESeq2)
library(sva)
library(pheatmap)
library(tibble)
library(dplyr)

# ---- 1. Load your raw count matrix ----
# Assuming GPCRS_data is already in the environment as your raw count matrix
# We'll use it directly instead of reloading

# ---- 2. Remove non-count columns if present ----
# Keep only numeric columns (filter out metadata-like columns)
GPCRS_data <- GPCRS_data[, sapply(GPCRS_data, is.numeric)]

# ---- 3. Save original NA positions and replace with 0s ----
na_mask <- is.na(GPCRS_data)  # Save where NAs were originally
GPCRS_data[na_mask] <- 0

# ---- 4. Prepare metadata from column names ----
sample_names <- colnames(GPCRS_data)

# Extract stage and species from sample names
stage <- sub("([A-Za-z]+)[0-9]+_.*", "\\1", sample_names)
species <- sub(".*_([A-Za-z]+_[a-z]+|Spis|Pacu)$", "\\1", sample_names)

# Create metadata dataframe
sample_info <- data.frame(
  sampleID = sample_names,
  stage = stage,
  species = species,
  group = paste(species, stage, sep = "_"),
  stringsAsFactors = FALSE
)

rownames(sample_info) <- sample_info$sampleID

# ---- 5. Create DESeq2 object ----
dds <- DESeqDataSetFromMatrix(countData = round(GPCRS_data),
                              colData = sample_info,
                              design = ~ 1)  # No design needed for VST

# ---- 6. Apply ComBat-seq ----
combat_counts <- ComBat_seq(counts = assay(dds), batch = sample_info$species)

# ---- 7. Apply VST directly on the ComBat-corrected matrix ----
dds_combat <- DESeqDataSetFromMatrix(countData = combat_counts,
                                     colData = sample_info,
                                     design = ~ 1)
dds_combat <- estimateSizeFactors(dds_combat)
vst_obj <- varianceStabilizingTransformation(dds_combat, blind = TRUE)
vst_mat <- assay(vst_obj)

# ---- 8. Z-score genes (rows) across all samples (global) ----
z_scores_global <- t(scale(t(vst_mat)))
z_scores_global[na_mask] <- NA

# ---- 9. Reorder samples by stage and species ----
ordered_samples <- sample_info %>% 
  mutate(stage = factor(stage, levels = c("Larva", "Meta", "Spat")),
         species = factor(species, levels = c("A_tenuis", "M_capitata", "P_acuta", "S_pistillata"))) %>%
  arrange(stage, species)

z_scores_ordered <- z_scores_global[, ordered_samples$sampleID]
sample_info <- sample_info[ordered_samples$sampleID, ]

# ---- 10. Retrieve gene family info from rownames if embedded ----
gene_names <- rownames(z_scores_ordered)
gene_family <- sub(".*GPCRS Family ([0-9]+).*", "Family_\\1", gene_names)
annotation_row <- data.frame(Family = gene_family)
rownames(annotation_row) <- gene_names

# ---- 11. Heatmap with fixed sample order and gene family annotation ----
pheatmap(z_scores_ordered,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = sample_info[, c("species", "stage")],
         annotation_row = annotation_row,
         fontsize_col = 7,
         angle_col = 45,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         na_col = "gray",
         main = "Z-score across samples (ordered)")
--------------------------------------------------------------------------------------------

#This is cool but M. capitata is too dominant
#Z scores per sample
  
  
# Load required libraries
library(DESeq2)
library(sva)
library(pheatmap)
library(tibble)
library(dplyr)

# ---- 1. Load your raw count matrix ----
GPCRS_data <- read_csv("GitHub/Multistage_omics/GPCRS/GPCRS_data.csv")

GPCRS_data <- as.data.frame(GPCRS_data)  # Convert to classic data.frame
# Use Orthogroup as rownames (gene ID), and save GPCR family info

rownames(GPCRS_data) <- GPCRS_data$Orthogroup

GPCR_family_info <- GPCRS_data[, "GPCRS Family", drop = FALSE]

# ---- 2. Remove non-count columns if present ----
# Keep only numeric columns (filter out metadata-like columns)
GPCRS_data_numeric <- GPCRS_data[, sapply(GPCRS_data, is.numeric)]

# ---- 3. Save original NA positions and replace with 0s ----
na_mask <- is.na(GPCRS_data_numeric)  # Save where NAs were originally
GPCRS_data_numeric[na_mask] <- 0

# ---- 4. Prepare metadata from column names ----
sample_names <- colnames(GPCRS_data_numeric)

# Extract stage and species from sample names
stage <- sub("([A-Za-z]+)[0-9]+_.*", "\\1", sample_names)
species <- sub(".*_([A-Za-z]+_[a-z]+|Spis|Pacu)$", "\\1", sample_names)

# Create metadata dataframe
sample_info <- data.frame(
  sampleID = sample_names,
  stage = stage,
  species = species,
  group = paste(species, stage, sep = "_"),
  stringsAsFactors = FALSE
)

rownames(sample_info) <- sample_info$sampleID

# ---- 5. Create DESeq2 object ----
dds <- DESeqDataSetFromMatrix(countData = round(GPCRS_data_numeric),
                              colData = sample_info,
                              design = ~ 1)  # No design needed for VST

# ---- 6. Apply ComBat-seq ----
combat_counts <- ComBat_seq(counts = assay(dds), batch = sample_info$species)

# ---- 7. Apply VST directly on the ComBat-corrected matrix ----
dds_combat <- DESeqDataSetFromMatrix(countData = combat_counts,
                                     colData = sample_info,
                                     design = ~ 1)
dds_combat <- estimateSizeFactors(dds_combat)
vst_obj <- varianceStabilizingTransformation(dds_combat, blind = TRUE)
vst_mat <- assay(vst_obj)

# ---- 8. Z-score genes (rows) across all samples (global) ----
z_scores_global <- t(scale(t(vst_mat)))
z_scores_global[na_mask] <- NA

# ---- 9. Reorder samples by stage and species ----
ordered_samples <- sample_info %>% 
  mutate(stage = factor(stage, levels = c("Larva", "Meta", "Spat")),
         species = factor(species, levels = c("A_tenuis", "M_capitata", "P_acuta", "S_pistillata"))) %>%
  arrange(stage, species)

z_scores_ordered <- z_scores_global[, ordered_samples$sampleID]
sample_info <- sample_info[ordered_samples$sampleID, ]

# ---- 10. Reattach saved GPCR family info ----
gene_names <- rownames(z_scores_ordered)

# Assume GPCR_family_info is in the same order as the rows of GPCRS_data
# and GPCRS_data_numeric still has correct rownames matching gene_names
rownames(GPCR_family_info) <- rownames(GPCRS_data_numeric)

# Match to the z-score matrix
gene_family <- GPCR_family_info[gene_names, , drop = FALSE]

# Replace missing family assignments with 'Unknown'
gene_family$`GPCRS Family`[is.na(gene_family$`GPCRS Family`)] <- "Unknown"
annotation_row <- data.frame(Family = gene_family$`GPCRS Family`)
rownames(annotation_row) <- gene_names

# ---- 11. Heatmap with fixed sample order and gene family annotation ----
pheatmap(z_scores_ordered,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = sample_info[, c("species", "stage")],
         annotation_row = annotation_row,
         fontsize_col = 7,
         angle_col = 45,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         na_col = "gray",
         main = "Z-score across samples (ordered)")

-----------------------------------------------------------------------------------------------
#Here I'm trying Z_scores by species  
  
# Load required libraries
library(DESeq2)
library(sva)
library(pheatmap)
library(tibble)
library(dplyr)

# ---- 1. Load your raw count matrix ----
GPCRS_data <- read_csv("GitHub/Multistage_omics/GPCRS/GPCRS_data.csv")
GPCRS_data <- as.data.frame(GPCRS_data)  # Convert to classic data.frame

# ---- 1. Load your raw count matrix ----
# Use Orthogroup as rownames (gene ID), and save GPCR family info
rownames(GPCRS_data) <- GPCRS_data$Orthogroup

GPCR_family_info <- GPCRS_data[, "GPCRS Family", drop = FALSE]

# ---- 2. Remove non-count columns if present ----
# Keep only numeric columns (filter out metadata-like columns)
GPCRS_data_numeric <- GPCRS_data[, sapply(GPCRS_data, is.numeric)]

# ---- 3. Save original NA positions and replace with 0s ----
na_mask <- is.na(GPCRS_data_numeric)  # Save where NAs were originally
GPCRS_data_numeric[na_mask] <- 0

# ---- 4. Prepare metadata from column names ----
sample_names <- colnames(GPCRS_data_numeric)

# Extract stage and species from sample names
stage <- sub("([A-Za-z]+)[0-9]+_.*", "\\1", sample_names)
species <- sub(".*_([A-Za-z]+_[a-z]+|Spis|Pacu)$", "\\1", sample_names)

# Create metadata dataframe
sample_info <- data.frame(
  sampleID = sample_names,
  stage = stage,
  species = species,
  group = paste(species, stage, sep = "_"),
  stringsAsFactors = FALSE
)

rownames(sample_info) <- sample_info$sampleID

# ---- 5. Create DESeq2 object ----
dds <- DESeqDataSetFromMatrix(countData = round(GPCRS_data_numeric),
                              colData = sample_info,
                              design = ~ 1)  # No design needed for VST

# ---- 6. Apply ComBat-seq ----
combat_counts <- ComBat_seq(counts = assay(dds), batch = sample_info$species)

# ---- 7. Apply VST directly on the ComBat-corrected matrix ----
dds_combat <- DESeqDataSetFromMatrix(countData = combat_counts,
                                     colData = sample_info,
                                     design = ~ 1)
dds_combat <- estimateSizeFactors(dds_combat)
vst_obj <- varianceStabilizingTransformation(dds_combat, blind = TRUE)
vst_mat <- assay(vst_obj)

# ---- 8. Z-score genes (rows) within species ----
z_scores_species <- vst_mat
for (sp in unique(sample_info$species)) {
  sp_cols <- sample_info$sampleID[sample_info$species == sp]
  z_scores_species[, sp_cols] <- t(scale(t(vst_mat[, sp_cols])))
}
z_scores_species[na_mask] <- NA

# ---- 9. Reorder samples by stage and species ----
ordered_samples <- sample_info %>% 
  mutate(stage = factor(stage, levels = c("Larva", "Meta", "Spat")),
         species = factor(species, levels = c("A_tenuis", "M_capitata", "P_acuta", "S_pistillata"))) %>%
  arrange(stage, species)

z_scores_ordered <- z_scores_species[, ordered_samples$sampleID]

# ---- 10. Reattach saved GPCR family info ----
gene_names <- rownames(z_scores_ordered)
gene_family <- GPCR_family_info[gene_names, , drop = FALSE]
gene_family$`GPCRS Family`[is.na(gene_family$`GPCRS Family`)] <- "Unknown"
annotation_row <- data.frame(Family = gene_family$`GPCRS Family`)
rownames(annotation_row) <- gene_names

dev.new()  # opens a new graphics device window

# ---- 11. Heatmap with fixed sample order and gene family annotation ----
pheatmap(z_scores_ordered,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = sample_info[, c("species", "stage")],
         annotation_row = annotation_row,
         fontsize_col = 7,
         angle_col = 45,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         na_col = "gray",
         main = "Z-score across samples (ordered)")


#Heatmap with some tweaks: 

# ---- 11. Heatmap with fixed sample order and gene family annotation ----

# Order genes by family (optional: can cluster within families later)
family_order <- c("Family_1", "Family_2", "Family_3", "Family_others")
annotation_row$Family <- factor(annotation_row$Family, levels = family_order)
z_scores_ordered <- z_scores_ordered[order(annotation_row$Family), ]
annotation_row <- annotation_row[rownames(z_scores_ordered), , drop = FALSE]

# Optional: open in separate window
dev.new()

# Plot
pheatmap(z_scores_ordered,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = sample_info[, c("species", "stage")],
         annotation_row = annotation_row,
         fontsize_col = 7,
         angle_col = 45,
         cluster_rows = FALSE,  # turn off row clustering to preserve family grouping
         cluster_cols = FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         na_col = "gray",
         main = "Species-normalized expression across life stages by GPCR family")

#Additional stuff: adding vertical labels for families

# Plot
# Plot all families in one heatmap
heatmap_obj <- pheatmap(z_scores_ordered,
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        annotation_col = sample_info[, c("species", "stage")],
                        fontsize_col = 7,
                        angle_col = 45,
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                        na_col = "gray",
                        main = "Species-normalized expression across life stages by GPCR family",
                        silent = TRUE)

grid.newpage()
grid.draw(heatmap_obj$gtable)

# ---- 12. Add vertical labels for each GPCR family ----
library(grid)
fam_levels <- c("Family 1", "Family 2", "Family 3", "Family others")
fam_positions <- sapply(fam_levels, function(fam) {
  which(annotation_row$Family == fam) |> range() |> mean()
})

for (i in seq_along(fam_levels)) {
  grid.text(label = fam_levels[i],
            x = unit(0.5, "lines"),
            y = unit(1 - (fam_positions[i] / nrow(z_scores_ordered)), "npc"),
            just = c("left", "center"),
            rot = 90,
            gp = gpar(fontsize = 10, fontface = "bold"))
}

