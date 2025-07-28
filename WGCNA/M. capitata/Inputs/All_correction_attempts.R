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
  left_join(alignment_table, by = c("sampleID" = "SampleID"))

#Correction based on TINS
# Load TIN data
tin_data <- read.table("multiqc_tin.txt", header = TRUE)

# Rename column for clarity
colnames(tin_data) <- c("SRR", "TIN")

tin_data <- left_join(tin_data, alignment_table, by = "SRR")  # Adds SampleID column

treatmentinfo <- left_join(treatmentinfo, tin_data, by = c("sampleID" = "SampleID"))

treatmentinfo$TIN_scaled <- scale(treatmentinfo$TIN)

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


