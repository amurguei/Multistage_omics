# Load packages
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library("tidyverse")
library("genefilter")
library("DESeq2")
library("RColorBrewer")
library("WGCNA")
library("flashClust")
library("gridExtra")
library("ComplexHeatmap")
library("goseq")
library("clusterProfiler")
library("pheatmap")
library("magrittr")
library("vegan")
library("factoextra")
library("dplyr")
library("dendsort")
library("NbClust")
library("simplifyEnrichment")
library("factoextra")
library("VennDiagram")
library("patchwork")
library(sva)


#setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff")

setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs")

#Data--------------------------------------------------------------------------------------
#treatmentinfo <- read_csv("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/5-Mcap-SampleInfo.csv")
library(readr)
treatmentinfo <- read_csv("5-Mcap-SampleInfo.csv")
#gcount <- as.data.frame(read.csv("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/Mcap_transcript_count_matrix.csv", row.names="gene_id"), colClasses = double, header=TRUE)
gcount <- as.data.frame(read.csv("Mcap_transcript_count_matrix.csv", row.names="gene_id"), colClasses = double, header=TRUE)


#Quality filter gene counts-----------------------------------------------------------

# Set filter values for PoverA
# Smallest sample size per treatment is 3, so 3/9 (9 samples) is 0.33 (rounded to 2 decimal places)
# This means that 3 out of 9 (0.33) samples need to have counts over 10.
# So P=33 percent of the samples have counts over A=10.
#there's something up with this code :( )
filt <- filterfun(pOverA(0.33, 10))

# Create filter for the counts data
gfilt <- genefilter(gcount, filt)

# Identify genes to keep by count filter
gkeep <- gcount[gfilt,]

# Identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#How many rows do we have before and after filtering?
nrow(gcount) #Before 54384
nrow(gcount_filt) #After 25741

# Normalize our read counts using VST-normalization in DESeq2
# Construct the DESeq2 dataset

#Merge Set group as a factor.
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I","II","III"))

#Create a DESeqDataSet design from gene count matrix and labels. Here we set the design to look at 
#any differences in gene expression across samples attributed to life stage

# Creating an alignment rate table
alignment_table <- tribble(
  ~SampleID, ~SRR, ~AlignRate,
  "AH1", "SRR14864072", 73.0,
  "AH2", "SRR14864071", 80.8,
  "AH3", "SRR14864070", 82.1,
  "AH4", "SRR14864069", 79.9,
  "AH5", "SRR14864068", 71.9,
  "AH6", "SRR14864067", 73.0,
  "AH7", "SRR14864066", 47.3,
  "AH8", "SRR14864065", 65.0,
  "AH9", "SRR14864064", 63.4
)

# Merge alignment rates into your existing sample metadata
treatmentinfo <- treatmentinfo %>%
  left_join(select(alignment_table, SampleID, AlignRate), by = c("sampleID" = "SampleID"))


#Correction based on TINS
# Load TIN data
tin_data <- read.table("multiqc_tin.txt", header = TRUE)

# Rename column for clarity
colnames(tin_data) <- c("SRR", "TIN")

tin_data <- left_join(tin_data, alignment_table, by = "SRR")  # Adds SampleID column

treatmentinfo <- left_join(treatmentinfo, tin_data, by = c("sampleID" = "SampleID"))

treatmentinfo$TIN_scaled <- scale(treatmentinfo$TIN)

tin_summary <- treatmentinfo %>%
  group_by(timepoint) %>%
  summarise(
    mean_TIN = mean(TIN, na.rm = TRUE),
    sd_TIN = sd(TIN, na.rm = TRUE),
    min_TIN = min(TIN, na.rm = TRUE),
    max_TIN = max(TIN, na.rm = TRUE),
    n = n()
  )
print(tin_summary)

ggplot(treatmentinfo, aes(x = timepoint, y = TIN, fill = timepoint)) +
  geom_boxplot() +
  labs(title = "TIN by Life Stage (Timepoint)", x = "Timepoint", y = "TIN") +
  theme_minimal()
ggsave("TIN_by_timepoint_boxplot.pdf", width = 6, height = 4)


#Note, this only includes covariates in the model for VST, for an actual correction there's need to residualize the data

# Function to run DESeq2 + VST + PCA
run_correction_pca <- function(counts, coldata, design_formula, label) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design_formula)
  dds <- DESeq(dds)
  vst_data <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vst_data, intgroup = "timepoint", returnData = TRUE, ntop = 25741)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  spread <- range(pca_data$PC1)
  
  p <- ggplot(pca_data, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
    geom_point(size = 5) +
    xlab(paste0("PC1: ", percent_var[1], "%")) +
    ylab(paste0("PC2: ", percent_var[2], "%")) +
    coord_fixed() +
    theme_classic() +
    ggtitle(label)
  
  list(
    plot = p,
    variance = percent_var[1],
    pc1_range = spread,
    pca_data = pca_data
  )
}

# Prepare metadata and counts
# Reorder metadata to match gcount columns
colnames(gcount) <- colnames(gcount_filt)
treatmentinfo <- treatmentinfo %>%
  filter(sampleID %in% colnames(gcount)) %>%
  arrange(match(sampleID, colnames(gcount)))

# Create spat_status variable for ComBat-Seq
if (!"spat_status" %in% colnames(treatmentinfo)) {
  treatmentinfo <- treatmentinfo %>%
    mutate(spat_status = ifelse(timepoint == "III", "spat", "no_spat"))
}


# Standard design
results <- list()

# 1. No correction
gdds <- DESeqDataSetFromMatrix(gcount_filt, treatmentinfo, design = ~timepoint)
results$raw <- run_correction_pca(gcount_filt, treatmentinfo, ~timepoint, "Raw")

# 2. Alignment rate corrected
treatmentinfo$AlignRate_scaled <- scale(treatmentinfo$AlignRate)
results$align <- run_correction_pca(gcount_filt, treatmentinfo, ~AlignRate_scaled + timepoint, "Alignment Corrected")

# 3. SVA-corrected
vst_mat <- assay(vst(gdds, blind=FALSE))
mod <- model.matrix(~timepoint, data = treatmentinfo)
mod0 <- model.matrix(~1, data = treatmentinfo)
svobj <- sva::svaseq(vst_mat, mod, mod0)
treatmentinfo$SV1 <- svobj$sv[,1]
treatmentinfo$SV2 <- svobj$sv[,2]
results$sva <- run_correction_pca(gcount_filt, treatmentinfo, ~SV1 + SV2 + timepoint, "SVA Corrected")

# 4. TIN corrected
treatmentinfo$TIN_scaled <- scale(treatmentinfo$TIN)
results$tin <- run_correction_pca(gcount_filt, treatmentinfo, ~TIN_scaled + timepoint, "TIN Corrected")

# 5. ComBat-Seq corrected
combat_counts <- ComBat_seq(as.matrix(gcount), batch = treatmentinfo$spat_status)
vst_combat <- vst(DESeqDataSetFromMatrix(countData = combat_counts, colData = treatmentinfo, design = ~timepoint), blind = TRUE)
pca_data_combat <- plotPCA(vst_combat, intgroup = "timepoint", returnData = TRUE, ntop = 25741)
percent_var_combat <- round(100 * attr(pca_data_combat, "percentVar"))
spread_combat <- range(pca_data_combat$PC1)
plot_combat <- ggplot(pca_data_combat, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percent_var_combat[1], "%")) +
  ylab(paste0("PC2: ", percent_var_combat[2], "%")) +
  coord_fixed() +
  theme_classic() +
  ggtitle("ComBat-Seq Corrected")
results$combat <- list(
  plot = plot_combat,
  variance = percent_var_combat[1],
  pc1_range = spread_combat,
  pca_data = pca_data_combat
)

# Summary Table
summary_table <- tibble(
  Method = names(results),
  PC1_Variance = sapply(results, function(x) x$variance),
  PC1_Range = sapply(results, function(x) diff(x$pc1_range))
)
print(summary_table)

#Combine plots
library(patchwork)
results$raw$plot + results$align$plot + results$sva$plot + results$tin$plot + results$combat$plot + plot_layout(ncol = 2)


#See plots individually

# Show each plot separately
print(results$raw$plot)
print(results$align$plot)
print(results$sva$plot)
print(results$tin$plot)
print(results$combat$plot)

#Note,all but ComBat-Seq are very similar because they're just being used as a covariate in the VST model. Will now try to residualize the data. 

# Build DESeq2 dataset and compute VST
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = treatmentinfo,
                               design = ~ timepoint)

gdds <- DESeq(gdds)
gvst <- vst(gdds, blind = FALSE)


# Function: Residualize expression matrix and return DESeqTransform object for plotPCA
residualize_vst <- function(vst_obj, coldata, vars_to_remove) {
  vst_mat <- assay(vst_obj)
  formula <- as.formula(paste("~", paste(vars_to_remove, collapse = " + ")))
  covariate_matrix <- model.matrix(formula, data = coldata)
  
  residuals_mat <- matrix(NA, nrow = nrow(vst_mat), ncol = ncol(vst_mat))
  rownames(residuals_mat) <- rownames(vst_mat)
  colnames(residuals_mat) <- colnames(vst_mat)
  
  for (i in 1:nrow(vst_mat)) {
    fit <- lm(vst_mat[i, ] ~ covariate_matrix)
    residuals_mat[i, ] <- residuals(fit)
  }
  
  vst_corrected <- vst_obj
  assay(vst_corrected) <- residuals_mat
  return(vst_corrected)
}

# Scale numeric covariates
library(sva)
treatmentinfo$AlignRate_scaled <- scale(treatmentinfo$AlignRate)
treatmentinfo$TIN_scaled <- scale(treatmentinfo$TIN)

# Calculate surrogate variables 
mod <- model.matrix(~ timepoint, data = treatmentinfo)
mod0 <- model.matrix(~ 1, data = treatmentinfo)
svobj <- svaseq(assay(gvst), mod, mod0)
treatmentinfo$SV1 <- svobj$sv[,1]
treatmentinfo$SV2 <- svobj$sv[,2]

# Apply individual corrections
vst_TIN <- residualize_vst(gvst, treatmentinfo, c("TIN_scaled"))
vst_Align <- residualize_vst(gvst, treatmentinfo, c("AlignRate_scaled"))
vst_SVA <- residualize_vst(gvst, treatmentinfo, c("SV1", "SV2"))
vst_All <- residualize_vst(gvst, treatmentinfo, c("TIN_scaled", "AlignRate_scaled", "SV1", "SV2"))

# Plotting example (DESeq2-style)
library(ggplot2)
pca_df <- plotPCA(vst_SVA, intgroup = "timepoint", returnData = TRUE, ntop = 25741)
percent_var <- round(100 * attr(pca_df, "percentVar"))

p <- ggplot(pca_df, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
  geom_point(size = 6) +
  xlab(paste0("PC1: ", percent_var[1], "%")) +
  ylab(paste0("PC2: ", percent_var[2], "%")) +
  coord_fixed() +
  theme_classic()
print(p)


library(ggplot2)

plot_custom_pca <- function(vst_object, title) {
  pca_df <- plotPCA(vst_object, intgroup = "timepoint", returnData = TRUE, ntop = 25741)
  percent_var <- round(100 * attr(pca_df, "percentVar"))
  
  ggplot(pca_df, aes(PC1, PC2, color = timepoint, shape = timepoint)) +
    geom_point(size = 6) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    coord_fixed() +
    theme_classic() +
    ggtitle(title) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
}

# Plot for each correction
ggplot_TIN <- plot_custom_pca(vst_TIN, "TIN Corrected")
ggplot_Align <- plot_custom_pca(vst_Align, "Alignment Rate Corrected")
ggplot_SVA <- plot_custom_pca(vst_SVA, "SVA Corrected")

# Display the plots
print(ggplot_TIN)
print(ggplot_Align)
print(ggplot_SVA)
print(ggplot_All)

#Testing with SVA corrected some WGCNA network parameters

# Prepare WGCNA input using raw and SVA-corrected matrices
library(WGCNA)
options(stringsAsFactors = FALSE)

# Raw expression matrix
expr_raw <- as.data.frame(t(assay(gvst)))

# SVA-corrected expression matrix
expr_sva <- as.data.frame(t(assay(vst_SVA)))

# Use the same top 5000 most variable genes based on the raw data
var_genes_raw <- apply(expr_raw, 2, var)
top_genes <- names(sort(var_genes_raw, decreasing = TRUE))[1:5000]

# Subset both matrices to the same genes
expr_raw_top <- expr_raw[, top_genes]
expr_sva_top <- expr_sva[, top_genes]

# Soft-thresholding power selection for raw
powers <- c(seq(from = 1, to=19, by=2), c(21:30))
sft_raw <- pickSoftThreshold(expr_raw_top, powerVector = powers, verbose = 5)

# Soft-thresholding power selection for SVA
sft_sva <- pickSoftThreshold(expr_sva_top, powerVector = powers, verbose = 5)

# Plot the results (optional: export to file)
par(mfrow = c(1,2))
plot(sft_raw$fitIndices[,1], -sign(sft_raw$fitIndices[,3])*sft_raw$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main="Raw - Scale independence")
text(sft_raw$fitIndices[,1], -sign(sft_raw$fitIndices[,3])*sft_raw$fitIndices[,2],
     labels=powers, cex=0.8, col="red")
abline(h=0.8, col="blue")

plot(sft_sva$fitIndices[,1], -sign(sft_sva$fitIndices[,3])*sft_sva$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main="SVA - Scale independence")
text(sft_sva$fitIndices[,1], -sign(sft_sva$fitIndices[,3])*sft_sva$fitIndices[,2],
     labels=powers, cex=0.8, col="red")
abline(h=0.8, col="blue")

hist(as.vector(assay(gvst)), breaks = 100, main = "Distribution of VST-transformed values (gvst)", xlab = "Expression")
summary(as.vector(assay(gvst)))



# Create DESeq2 object
dds_sva <- DESeqDataSetFromMatrix(gcount_filt, treatmentinfo, design = ~SV1 + SV2 + timepoint)

# Run DESeq
dds_sva <- DESeq(dds_sva)

# Get VST matrix including SVA covariates
vst_sva_cov <- vst(dds_sva, blind = FALSE)

# Raw VST matrix (no SVA)
expr_raw <- assay(gvst)  # from earlier step

# SVA-corrected VST matrix (via covariates)
expr_sva_cov <- assay(vst_sva_cov)


# Raw VST matrix (no SVA)
expr_raw <- assay(gvst)  # from earlier step

# SVA-corrected VST matrix (via covariates)
expr_sva_cov <- assay(vst_sva_cov)

datExpr_raw <- t(expr_raw)
datExpr_sva_cov <- t(expr_sva_cov)


# Match sample IDs
datExpr_raw <- datExpr_raw[match(treatmentinfo$sampleID, rownames(datExpr_raw)), ]
datExpr_sva_cov <- datExpr_sva_cov[match(treatmentinfo$sampleID, rownames(datExpr_sva_cov)), ]

# Check for NAs
datExpr_raw <- datExpr_raw[complete.cases(datExpr_raw), ]
datExpr_sva_cov <- datExpr_sva_cov[complete.cases(datExpr_sva_cov), ]

powers <- c(1:20)

# Raw
sft_raw <- pickSoftThreshold(datExpr_raw, powerVector = powers, verbose = 5)

# SVA as covariate
sft_sva_cov <- pickSoftThreshold(datExpr_sva_cov, powerVector = powers, verbose = 5)

par(mfrow = c(1, 2))  # 1 row, 2 columns for side-by-side plots

# --- Plot 1: Raw VST ---
plot(sft_raw$fitIndices[,1], 
     -sign(sft_raw$fitIndices[,3]) * sft_raw$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale-Free Topology Fit (R²)",
     type = "n", main = "Raw VST - Scale Independence")
text(sft_raw$fitIndices[,1], 
     -sign(sft_raw$fitIndices[,3]) * sft_raw$fitIndices[,2], 
     labels = powers, col = "red")
abline(h = 0.8, col = "blue", lty = 2)

# --- Plot 2: SVA VST ---
plot(sft_sva_cov$fitIndices[,1], 
     -sign(sft_sva_cov$fitIndices[,3]) * sft_sva_cov$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale-Free Topology Fit (R²)",
     type = "n", main = "SVA-Corrected VST - Scale Independence")
text(sft_sva_cov$fitIndices[,1], 
     -sign(sft_sva_cov$fitIndices[,3]) * sft_sva_cov$fitIndices[,2], 
     labels = powers, col = "red")
abline(h = 0.8, col = "blue", lty = 2)

#Mean connectivity side by side
par(mfrow = c(1, 2))  # Side-by-side 

# --- Plot 1: Raw VST ---
plot(sft_raw$fitIndices[,1], sft_raw$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Raw VST - Mean Connectivity")
text(sft_raw$fitIndices[,1], sft_raw$fitIndices[,5], 
     labels = powers, col = "red")

# --- Plot 2: SVA VST ---
plot(sft_sva_cov$fitIndices[,1], sft_sva_cov$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "SVA VST - Mean Connectivity")
text(sft_sva_cov$fitIndices[,1], sft_sva_cov$fitIndices[,5], 
     labels = powers, col = "red")

power_raw <- 15      
power_sva <- 15      

# RAW VST network
net_raw <- blockwiseModules(datExpr_raw, power = power_raw,
                            TOMType = "unsigned", minModuleSize = 30,
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, saveTOMFileBase = "TOM_raw",
                            verbose = 3)

# SVA-corrected VST network
net_sva <- blockwiseModules(datExpr_sva_cov, power = power_sva,
                            TOMType = "unsigned", minModuleSize = 30,
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, saveTOMFileBase = "TOM_sva",
                            verbose = 3)

table(net_raw$colors)
table(net_sva$colors)

# Convert numeric labels to color names
moduleColors_raw <- labels2colors(net_raw$colors)
moduleColors_sva <- labels2colors(net_sva$colors)

table(moduleColors_raw)
table(moduleColors_sva)

# === Prepare multi-trait matrix ===
datTraits <- data.frame(
  timepoint = as.numeric(treatmentinfo$timepoint),
  larvae_released = c(1,1,1,0,0,0,0,0,0),
  larvae_compressed = c(0,0,0,1,1,1,0,0,0),
  spat = c(0,0,0,0,0,0,1,1,1)
)

# Make sure sample order matches
rownames(datTraits) <- treatmentinfo$sampleID
rownames(MEs_raw) <- rownames(datTraits)  # same for MEs_sva

# === Correlation and p-values ===
moduleTraitCor_raw <- cor(MEs_raw, datTraits, use = "p")
moduleTraitP_raw <- corPvalueStudent(moduleTraitCor_raw, nrow(datTraits))

moduleTraitCor_sva <- cor(MEs_sva, datTraits, use = "p")
moduleTraitP_sva <- corPvalueStudent(moduleTraitCor_sva, nrow(datTraits))

# Create text matrix
textMatrix <- matrix(paste(signif(moduleTraitCor_raw, 2), "\n(",
                           signif(moduleTraitP_raw, 1), ")"),
                     nrow = nrow(moduleTraitCor_raw),
                     ncol = ncol(moduleTraitCor_raw),
                     dimnames = dimnames(moduleTraitCor_raw))

# Adjust font and PDF dimensions
pdf("Raw_ModuleTraitRelationships_fixed.pdf", height = 12, width = 8)

labeledHeatmap(Matrix = moduleTraitCor_raw,
               xLabels = colnames(datTraits),
               yLabels = rownames(moduleTraitCor_raw),
               ySymbols = rownames(moduleTraitCor_raw),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.6,
               cex.lab.x = 1.0,
               cex.lab.y = 0.5,
               zlim = c(-1, 1),
               main = "Raw Module-Trait Relationships")

dev.off()


#Complex heatmaps
# Assume MEs_raw, MEs_sva, and datTraits are already prepared

# Load required packages
# Load required packages
library(ComplexHeatmap)
library(circlize)
library(grid)

# Assume MEs_raw, MEs_sva, datTraits, and module color vectors are already prepared

# Get module colors directly from WGCNA output if available
# For example: moduleColors_raw <- net_raw$colors and moduleColors_sva <- net_sva$colors

# Convert module colors to unique names matching MEs
moduleNames_raw <- names(MEs_raw)
moduleColors_named_raw <- gsub("^ME", "", moduleNames_raw)
names(moduleColors_named_raw) <- moduleNames_raw
rownames(moduleTraitCor_raw) <- moduleColors_named_raw
rownames(moduleTraitP_raw) <- moduleColors_named_raw

moduleNames_sva <- names(MEs_sva)
moduleColors_named_sva <- gsub("^ME", "", moduleNames_sva)
names(moduleColors_named_sva) <- moduleNames_sva
rownames(moduleTraitCor_sva) <- moduleColors_named_sva
rownames(moduleTraitP_sva) <- moduleColors_named_sva

# 3. Generate heatmap text matrices
heatmappval_raw <- signif(moduleTraitP_raw, 1)
rownames(heatmappval_raw) <- moduleColors_named_raw

heatmappval_sva <- signif(moduleTraitP_sva, 1)
rownames(heatmappval_sva) <- moduleColors_named_sva

# 4. Dendrograms
row_dend_raw <- hclust(dist(moduleTraitCor_raw))
col_dend_raw <- hclust(dist(t(moduleTraitCor_raw)))

row_dend_sva <- hclust(dist(moduleTraitCor_sva))
col_dend_sva <- hclust(dist(t(moduleTraitCor_sva)))

# 5. Raw Module-Trait Heatmap
pdf("Raw_ModuleTrait_Heatmap.pdf", width = 8, height = 11)
Heatmap(moduleTraitCor_raw,
        name = "Correlation",
        column_title = "Raw Module–Trait Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_rows = row_dend_raw,
        cluster_columns = col_dend_raw,
        show_column_names = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          label <- formatC(heatmappval_raw[i, j], format = "e", digits = 1)
          grid.text(label, x, y,
                    gp = gpar(fontsize = 8,
                              fontface = ifelse(heatmappval_raw[i, j] <= 0.05, "bold", "plain")))
        })
dev.off()

# 6. SVA-Corrected Module-Trait Heatmap
pdf("SVA_ModuleTrait_Heatmap.pdf", width = 8, height = 11)
Heatmap(moduleTraitCor_sva,
        name = "Correlation",
        column_title = "SVA-Corrected Module–Trait Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_rows = row_dend_sva,
        cluster_columns = col_dend_sva,
        show_column_names = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10),
        cell_fun = function(j, i, x, y, width, height, fill) {
          label <- formatC(heatmappval_sva[i, j], format = "e", digits = 1)
          grid.text(label, x, y,
                    gp = gpar(fontsize = 8,
                              fontface = ifelse(heatmappval_sva[i, j] <= 0.05, "bold", "plain")))
        })
dev.off()


# Convert numeric labels to color names if needed
moduleColors_raw <- labels2colors(net_raw$colors)
moduleColors_sva <- labels2colors(net_sva$colors)

# Count module sizes
table(moduleColors_raw)
table(moduleColors_sva)


# Convert color labels for clarity
moduleColors_raw <- labels2colors(net_raw$colors)
moduleColors_sva <- labels2colors(net_sva$colors)

# Count module sizes
raw_sizes <- table(moduleColors_raw)
sva_sizes <- table(moduleColors_sva)

# Combine into a data frame
module_df <- merge(as.data.frame(raw_sizes), as.data.frame(sva_sizes),
                   by = "Var1", all = TRUE)
names(module_df) <- c("Module", "Raw", "SVA")
module_df[is.na(module_df)] <- 0

# Barplot of module sizes comparison
pdf("Module_Size_Comparison.pdf", width = 10, height = 12)
ggplot(module_df, aes(x = reorder(Module, -Raw), y = Raw, fill = "Raw")) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_bar(aes(y = SVA, fill = "SVA"), stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Raw" = "steelblue", "SVA" = "tomato")) +
  labs(x = "Module Color", y = "Number of Genes", fill = "Dataset") +
  theme_minimal(base_size = 12) +
  coord_flip()
dev.off()

# Load required packages
library(ComplexHeatmap)
library(circlize)
library(grid)
library(WGCNA)
library(flashClust)
library(ggplot2)

# Assume MEs_raw, MEs_sva, datTraits, and module color vectors are already prepared

# Convert color labels for clarity
moduleColors_raw <- labels2colors(net_raw$colors)
moduleColors_sva <- labels2colors(net_sva$colors)

# Count module sizes
raw_sizes <- table(moduleColors_raw)
sva_sizes <- table(moduleColors_sva)

# Convert to data frames with consistent structure
df_raw <- data.frame(Module = names(raw_sizes), Count = as.vector(raw_sizes), Dataset = "Raw")
df_sva <- data.frame(Module = names(sva_sizes), Count = as.vector(sva_sizes), Dataset = "SVA")



# Count module sizes
raw_sizes <- table(moduleColors_raw)
sva_sizes <- table(moduleColors_sva)

# Ensure all color names are included in both
all_colors <- union(names(raw_sizes), names(sva_sizes))
raw_sizes_full <- setNames(as.numeric(raw_sizes[all_colors]), all_colors)
sva_sizes_full <- setNames(as.numeric(sva_sizes[all_colors]), all_colors)
raw_sizes_full[is.na(raw_sizes_full)] <- 0
sva_sizes_full[is.na(sva_sizes_full)] <- 0

# Build tidy data frame
combined_df <- data.frame(
  Module = rep(all_colors, 2),
  Count = c(raw_sizes_full, sva_sizes_full),
  Dataset = rep(c("Raw", "SVA"), each = length(all_colors))
)

# Plot 
library(ggplot2)
ggplot(combined_df, aes(x = reorder(Module, Count), y = Count, fill = Module)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Dataset, ncol = 1, scales = "free_y") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = "Module Color", y = "Gene Count", title = "Module Sizes by Dataset")

#  Which number corresponds to "turquoise"
turquoise_label <- which(labels2colors(0:50) == "turquoise")
turquoise_label
#1, associated with larval stages

which(labels2colors(0:50) == "blue")
#3, associated with spat

# Print mapping of module numbers to color names
module_label_map <- data.frame(
  ModuleNumber = 0:50,
  ModuleColor = labels2colors(0:50)
)
# print(module_label_map)

# Example output:
#   ModuleNumber     ModuleColor
# 1             0            grey
# 2             1       turquoise
# 3             2            blue
# 4             3           brown
# 5             4          yellow
# 6             5           green
# 7             6             red
# 8             7           black
# 9             8            pink
# 10            9         magenta
# 11           10          purple
# 12           11     greenyellow
# 13           12             tan
# 14           13          salmon
# 15           14            cyan
# 16           15    midnightblue
# 17           16       lightcyan
# 18           17          grey60
# 19           18      lightgreen
# 20           19     lightyellow
# 21           20       royalblue
# 22           21         darkred
# 23           22       darkgreen
# 24           23   darkturquoise
# 25           24        darkgrey
# 26           25          orange
# 27           26      darkorange
# 28           27           white
# 29           28         skyblue
# 30           29     saddlebrown
# 31           30       steelblue
# 32           31   paleturquoise
# 33           32          violet
# 34           33  darkolivegreen
# 35           34     darkmagenta
# 36           35         sienna3
# 37           36     yellowgreen
# 38           37        skyblue3
# 39           38           plum1
# 40           39      orangered4
# 41           40   mediumpurple3
# 42           41 lightsteelblue1
# 43           42      lightcyan1
# 44           43           ivory
# 45           44     floralwhite
# 46           45     darkorange2
# 47           46          brown4
# 48           47         bisque4
# 49           48   darkslateblue
# 50           49           plum2
# 51           50        thistle2


