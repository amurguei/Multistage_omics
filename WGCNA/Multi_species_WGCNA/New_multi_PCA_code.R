#Orthogroups WGCNA

#Packages-------------------------------------------------------------------------------
if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse')
if ("RColorBrewer" %in% rownames(installed.packages()) == 'FALSE') install.packages('RColorBrewer') 
if ("WGCNA" %in% rownames(installed.packages()) == 'FALSE') install.packages('WGCNA') 
if ("flashClust" %in% rownames(installed.packages()) == 'FALSE') install.packages('flashClust') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("pheatmap" %in% rownames(installed.packages()) == 'FALSE') install.packages('pheatmap') 
if ("magrittr" %in% rownames(installed.packages()) == 'FALSE') install.packages('magrittr') 
if ("vegan" %in% rownames(installed.packages()) == 'FALSE') install.packages('vegan') 
if ("factoextra" %in% rownames(installed.packages()) == 'FALSE') install.packages('factoextra') 


if ("preprocessCore" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("preprocessCore") 
if ("impute" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("impute") 
if ("genefilter" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("genefilter") 
if ("DESeq2" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install('DESeq2')
if ("ComplexHeatmap" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install('ComplexHeatmap') 
if ("goseq" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install('goseq')
if ("clusterProfiler" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install('clusterProfiler') 


if ("tidyverse" %in% rownames(installed.packages()) == 'FALSE') install.packages('tidyverse')
if ("RColorBrewer" %in% rownames(installed.packages()) == 'FALSE') install.packages('RColorBrewer')
if ("WGCNA" %in% rownames(installed.packages()) == 'FALSE') install.packages('WGCNA') 
if ("flashClust" %in% rownames(installed.packages()) == 'FALSE') install.packages('flashClust') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("pheatmap" %in% rownames(installed.packages()) == 'FALSE') install.packages('pheatmap') 
if ("magrittr" %in% rownames(installed.packages()) == 'FALSE') install.packages('magrittr') 
if ("vegan" %in% rownames(installed.packages()) == 'FALSE') install.packages('vegan') 
if ("factoextra" %in% rownames(installed.packages()) == 'FALSE') install.packages('factoextra') 
if ("dendsort" %in% rownames(installed.packages()) == 'FALSE') install.packages('dendsort') 
if ("NbClust" %in% rownames(installed.packages()) == 'FALSE') install.packages('NbClust')
if ("VenDiagram"%in% rownames(installed.packages()) == 'FALSE') install.packages('VennDiagram')
if ("patchwork"%in% rownames(installed.packages()) == 'FALSE') install.packages('patchwork')


if ("genefilter" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("genefilter") 
if ("DESeq2" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install('DESeq2') 
if ("ComplexHeatmap" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install('ComplexHeatmap') 
if ("goseq" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install('goseq')
if ("simplifyEnrichment" %in% rownames(installed.packages()) == 'FALSE')BiocManager::install("simplifyEnrichment")
if ("clusterProfiler" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install('clusterProfiler') 
if ("sva" %in% rownames(installed.packages()) == 'FALSE') BiocManager::install("sva")



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
#Loading the packages

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
library("dendsort")
library("ggplot2")
library("dplyr")
library("sva")

setwd("C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA")
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")


# Load necessary library
library(readr)

# Read the CSV file with specified column types
gcount_raw <- read_csv("merged_data_acropora_montipora_pocillopora_spis.csv", 
                       col_types = cols(
                         DRR318288 = col_integer(), 
                         DRR318289 = col_integer(), 
                         DRR318291 = col_integer(), 
                         DRR318293 = col_integer(), 
                         DRR318294 = col_integer(), 
                         DRR318295 = col_integer(), 
                         DRR318297 = col_integer(), 
                         DRR318298 = col_integer(), 
                         DRR318299 = col_integer(), 
                         AH1 = col_integer(), 
                         AH2 = col_integer(), 
                         SRR14333320 = col_integer(), 
                         SRR14333321 = col_integer(), 
                         SRR14333322 = col_integer(), 
                         SRR14333323 = col_integer(), 
                         SRR14333324 = col_integer(), 
                         SRR14333325 = col_integer(), 
                         SRR14333326 = col_integer(), 
                         SRR14333327 = col_integer()
                       ))

treatmentinfo <- read_csv("treatmentinfo_merged.csv")

#Renaming columns so they're more manageable

# Identify columns that need to be renamed (7 onward)
old_names <- colnames(gcount_raw)[7:ncol(gcount_raw)]

# Match old names with new_names in treatmentinfo
new_names <- treatmentinfo$new_names[match(old_names, treatmentinfo$sampleID)]

# Replace only the column names from column 7 onward
colnames(gcount_raw)[7:ncol(gcount_raw)] <- ifelse(is.na(new_names), old_names, new_names)

#Switching the order for Spis spat

# Define column ranges
spis_spat_cols <- 34:36  # Current positions of Spis spat
spis_larvae_cols <- 40:42  # Current positions of Spis larvae

# Swap the columns by reordering
# Define column ranges
spis_spat_cols <- 34:36  # Current positions of Spis spat
spis_larvae_cols <- 40:42  # Current positions of Spis larvae

# Swap the columns
gcount_raw <- gcount_raw[, c(1:33, spis_larvae_cols, 37:39, spis_spat_cols)]

#Now I'm making a dataframe and calling Orthogroup rownames
gcount_raw <- as.data.frame(gcount_raw)
rownames(gcount_raw) <- gcount_raw$Orthogroup  # If gene IDs are in a column
# Extract only the count data (columns 7 to 42)
gcount <- gcount_raw[, 7:42]

# Set filter values for PoverA
filt <- filterfun(pOverA(0.083, 10))

# Apply filter to count data only
gfilt <- genefilter(gcount, filt)

# Identify genes to keep by count filter
gkeep <- gcount[gfilt, ]

# Identify gene names to keep
gn.keep <- rownames(gkeep)

# Subset the original dataset to keep metadata and filtered count data
gcount_filt <- gcount[rownames(gcount) %in% gn.keep, ]

# Check how many rows before and after filtering
nrow(gcount)  # Before filtering
#7129
nrow(gcount_filt) # After filtering
#7118

#Now I'm going to try to remove batch effects because it's looking a little tetric right now. I'm trying Combat-Seq

# Load necessary library
library(sva)
library(dplyr)  # Ensure dplyr is loaded for data manipulation

# Convert treatmentinfo to a dataframe (avoid tibble issues)
treatmentinfo <- as.data.frame(treatmentinfo)

# Ensure gcount_filt is a matrix
gcount_filt <- as.matrix(gcount_filt)

# Extract batch information correctly (align column names with new_names)
batch_info <- treatmentinfo$Species[match(colnames(gcount_filt), treatmentinfo$new_names)]

# Ensure batch_info is a factor
batch_info <- as.factor(batch_info)

# Apply ComBat-seq batch correction
gcount_corrected <- ComBat_seq(counts = gcount_filt, batch = batch_info, group = NULL)

# Check the adjusted count matrix
head(gcount_corrected)


# Normalize our read counts using VST-normalization in DESeq2
# Construct the DESeq2 dataset
#Merge the timepoint columns into a new column, group. Set group as a factor.
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I","II","III"))
treatmentinfo$Species <- factor(treatmentinfo$Species, levels = c("Acropora_tenuis", "Montipora_capitata", "Pocillopora_acuta","Stylophora_pistillata"))


gdds <- DESeqDataSetFromMatrix(countData = gcount_corrected,
                               colData = treatmentinfo,
                               design = ~ Species * timepoint)
# Log-transform the count data using a variance stabilizing transforamtion (vst). 

SF.gdds <- estimateSizeFactors( gdds ) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than for to use vst
print(sizeFactors(SF.gdds)) #View size factors
#Larva1_Atenuis Larva2_Atenuis Larva3_Atenuis  Meta1_Atenuis  Meta2_Atenuis  Meta3_Atenuis  Spat1_Atenuis  Spat2_Atenuis 
#1.8633992      2.1859570      2.1389978      1.2525420      1.3537586      1.1058002      1.7772996      1.9022468 
#Spat3_Atenuis    Larva1_Mcap    Larva2_Mcap    Larva3_Mcap     Meta1_Mcap     Meta2_Mcap     Meta3_Mcap     Spat1_Mcap 
#1.4566066      2.4887247      3.0634007      3.1470177      3.0737340      2.6238828      3.1272286      0.7344009 
#Spat2_Mcap     Spat3_Mcap    Larva1_Pacu    Larva2_Pacu    Larva3_Pacu     Meta1_Pacu     Meta2_Pacu     Meta3_Pacu 
#0.9863250      1.5349595      0.6585072      0.8927192      1.0516380      0.6744375      0.6988360      0.7895888 
#Spat1_Pacu     Spat2_Pacu     Spat3_Pacu    Larva1_Spis    Larva2_Spis    Larva3_Spis     Meta1_Spis     Meta2_Spis 
#0.9387776      0.9030795      0.8710691      0.3277546      0.2142616      0.2788691      0.3028476      0.4435857 
#Meta3_Spis     Spat1_Spis     Spat2_Spis     Spat3_Spis 
#0.5622346      0.4970006      0.4693681      0.4603734 
gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size

# Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.

datExpr <- as.data.frame(t(assay(gvst))) #transpose to output to a new data frame with the column names as row names. And make all data numeric



#Centering, no ComBat, No Z-scoring. 
  
## From raw counts (no ComBat-seq)
dds1 <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = treatmentinfo,
                               design = ~ Species * timepoint)
vst1 <- vst(dds1, blind = FALSE)

# Transpose and center
datExpr1 <- as.data.frame(t(assay(vst1)))
datExpr1_centered <- scale(datExpr1, center = TRUE, scale = FALSE)

# PCA
pca1 <- prcomp(datExpr1_centered, center = FALSE, scale. = FALSE)
plot (pca1)


# Define Z-score function per species
zscore_within_species <- function(x, groups) {
  tapply(x, groups, function(vals) scale(vals, center = TRUE, scale = TRUE)) |>
    unlist(use.names = FALSE)
}

# Helper to convert prcomp result to tidy PCA dataframe
get_pca_df <- function(pca_obj, metadata, sample_names) {
  df <- as.data.frame(pca_obj$x)
  df$Sample <- sample_names
  df$Species <- metadata$Species[match(df$Sample, metadata$new_names)]
  df$timepoint <- metadata$timepoint[match(df$Sample, metadata$new_names)]
  return(df)
}

# Helper to create ggplot from PCA dataframe
plot_pca <- function(pca_df, pca_obj, title) {
  percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))[1:2]
  ggplot(pca_df, aes(x = PC1, y = PC2, color = timepoint, shape = Species)) +
    geom_point(size = 5) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_classic() +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    ggtitle(title)
}

# Helper to create ggplot from PCA dataframe
plot_pca <- function(pca_df, pca_obj, title) {
  percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))[1:2]
  ggplot(pca_df, aes(x = PC1, y = PC2, color = timepoint, shape = Species)) +
    geom_point(size = 5) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_classic() +
    scale_shape_manual(values = c(15, 16, 17, 18)) +
    ggtitle(title)
}

#Treatment 1. Scenario 1: VST only (species-aware), with centering

dds1 <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = treatmentinfo,
                               design = ~ Species * timepoint)
vst1 <- vst(dds1, blind = FALSE)
datExpr1 <- as.data.frame(t(assay(vst1)))
datExpr1_centered <- scale(datExpr1, center = TRUE, scale = FALSE)
pca1 <- prcomp(datExpr1_centered, center = FALSE, scale. = FALSE)
pca_df1 <- get_pca_df(pca1, treatmentinfo, rownames(datExpr1))
p1 <- plot_pca(pca_df1, pca1, "1. VST + Centering")

#Treatment 2. vst + Z-scoring per species
datExpr2 <- datExpr1
datExpr2$Species <- treatmentinfo$Species[match(rownames(datExpr2), treatmentinfo$new_names)]
datExpr2_zscored <- as.data.frame(apply(datExpr2[, -ncol(datExpr2)], 2, zscore_within_species, groups = datExpr2$Species))
rownames(datExpr2_zscored) <- rownames(datExpr2)
pca2 <- prcomp(datExpr2_zscored, center = FALSE, scale. = FALSE)
pca_df2 <- get_pca_df(pca2, treatmentinfo, rownames(datExpr2))
p2 <- plot_pca(pca_df2, pca2, "2. VST + Z-scoring per Species")

# Z-score per species
datExpr2$Species <- treatmentinfo$Species[match(rownames(datExpr2), treatmentinfo$new_names)]
datExpr2_zscored <- as.data.frame(apply(datExpr2[, -ncol(datExpr2)], 2, zscore_within_species, groups = datExpr2$Species))

# Remove any problematic genes
datExpr2_zscored <- datExpr2_zscored[, apply(datExpr2_zscored, 2, function(x) all(is.finite(x)))]
rownames(datExpr2_zscored) <- rownames(datExpr2)

# Now PCA will work
pca2 <- prcomp(datExpr2_zscored, center = FALSE, scale. = FALSE)
pca_df2 <- get_pca_df(pca2, treatmentinfo, rownames(datExpr2_zscored))
p2 <- plot_pca(pca_df2, pca2, "2. VST + Z-scoring per Species")

#How many genes were removed? 
cat("Genes before filtering:", ncol(datExpr2) - 1, "\n")  # subtract Species col
cat("Genes after Z-scoring cleanup:", ncol(datExpr2_zscored), "\n")


# Scenario 3: ComBat-seq + VST + Centering

gdds3 <- DESeqDataSetFromMatrix(countData = gcount_corrected,
                                colData = treatmentinfo,
                                design = ~ Species * timepoint)
vst3 <- vst(gdds3, blind = FALSE)
datExpr3 <- as.data.frame(t(assay(vst3)))
datExpr3_centered <- scale(datExpr3, center = TRUE, scale = FALSE)
pca3 <- prcomp(datExpr3_centered, center = FALSE, scale. = FALSE)
pca_df3 <- get_pca_df(pca3, treatmentinfo, rownames(datExpr3))
p3 <- plot_pca(pca_df3, pca3, "3. ComBat-seq + VST + Centering")

#Scenario 4: ComBat-Seq + VST + Centering + Z scoring per species

datExpr4 <- datExpr3
datExpr4$Species <- treatmentinfo$Species[match(rownames(datExpr4), treatmentinfo$new_names)]
datExpr4_zscored <- as.data.frame(apply(datExpr4[, -ncol(datExpr4)], 2, zscore_within_species, groups = datExpr4$Species))
rownames(datExpr4_zscored) <- rownames(datExpr4)
pca4 <- prcomp(datExpr4_zscored, center = FALSE, scale. = FALSE)
pca_df4 <- get_pca_df(pca4, treatmentinfo, rownames(datExpr4))
p4 <- plot_pca(pca_df4, pca4, "4. ComBat-seq + VST + Z-scoring")

# Z-score per species
datExpr4$Species <- treatmentinfo$Species[match(rownames(datExpr4), treatmentinfo$new_names)]
datExpr4_zscored <- as.data.frame(apply(datExpr4[, -ncol(datExpr4)], 2, zscore_within_species, groups = datExpr4$Species))

# Remove any problematic genes (columns with NaN/Inf/NA)
datExpr4_zscored <- datExpr4_zscored[, apply(datExpr4_zscored, 2, function(x) all(is.finite(x)))]

# Set rownames
rownames(datExpr4_zscored) <- rownames(datExpr4)

# Now PCA will work
pca4 <- prcomp(datExpr4_zscored, center = FALSE, scale. = FALSE)
pca_df4 <- get_pca_df(pca4, treatmentinfo, rownames(datExpr4_zscored))
p4 <- plot_pca(pca_df4, pca4, "4. ComBat-seq + VST + Z-scoring")

#How many genes were removed? 
cat("Genes before filtering:", ncol(datExpr4) - 1, "\n")  # subtract Species col
cat("Genes after Z-scoring cleanup:", ncol(datExpr4_zscored), "\n")


library(patchwork)
(p1 | p2) / (p3 | p4)

plot(p1)
plot(p2)
plot(p3)
plot(p4)

# Mean expression per sample (after VST, no centering)
mean_expr_5 <- rowMeans(datExpr3)  # datExpr3 = from ComBat-seq + VST

# PC1 scores for each sample
pc1_scores_5 <- pca3$x[, 1]  # from Scenario 5 (ComBat-seq + VST + Centering)

# Plot mean expression vs PC1
plot(mean_expr_5, pc1_scores_5,
     xlab = "Mean expression (VST)", ylab = "PC1 score",
     main = "Mean Expression vs PC1 Score",
     pch = 19)
abline(lm(pc1_scores_5 ~ mean_expr_5), col = "red")

cor(mean_expr_5, pc1_scores_5)
summary(lm(pc1_scores_5 ~ mean_expr_5))

# Mean expression per sample after centering
mean_expr_centered <- rowMeans(datExpr1_centered)

# PC1 scores after centering
pc1_scores_centered <- pca1$x[, 1]

plot(mean_expr_centered, pc1_scores_centered,
     xlab = "Mean expression (centered VST)",
     ylab = "PC1 score",
     main = "Mean Expression vs PC1 Score (Centered)",
     pch = 19)
abline(lm(pc1_scores_centered ~ mean_expr_centered), col = "red")

cor(mean_expr_centered, pc1_scores_centered)
summary(lm(pc1_scores_centered ~ mean_expr_centered))

# Mean expression per sample after Z-scoring
mean_expr_z4 <- rowMeans(datExpr4_zscored)

# PC1 scores after PCA on Z-scored data
pc1_scores_z4 <- pca4$x[, 1]


# Mean expression per sample after Z-scoring
mean_expr_z4 <- rowMeans(datExpr4_zscored)

# PC1 scores after PCA on Z-scored data
pc1_scores_z4 <- pca4$x[, 1]


plot(mean_expr_z4, pc1_scores_z4,
     xlab = "Mean expression (Z-scored, Scenario 4)",
     ylab = "PC1 score",
     main = "Mean Expression vs PC1 Score (Z-scored)",
     pch = 19)
abline(lm(pc1_scores_z4 ~ mean_expr_z4), col = "red")


#Now getting my final PCA plot. I'm doing: ComBat-Seq + VST (species-aware) + centered. 

# Color palette by timepoint
# High-contrast color-blind palette for timepoints
timepoint_palette <- c(
  "I"   = "#e9a3c9",  # pink-purple
  "II"  = "#f7f7f7",  # soft white
  "III" = "#a1d76a"   # green
)

plot_pca_cb <- function(pca_df, pca_obj, title) {
  percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))[1:2]
  
  ggplot(pca_df, aes(x = PC1, y = PC2, fill = timepoint, shape = Species)) +
    geom_point(size = 5, color = "black", stroke = 0.8) +  # black outline
    scale_fill_manual(values = timepoint_palette, name = "Timepoint") +
    scale_shape_manual(
      values = c(
        "Acropora_tenuis" = 21,
        "Montipora_capitata" = 22,
        "Pocillopora_acuta" = 24,
        "Stylophora_pistillata" = 23
      ),
      name = "Species"
    ) +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, color = "black")),  # fixes black-only issue
      shape = guide_legend(override.aes = list(fill = "grey80", color = "black"))  # default fill for clarity
    ) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    ) +
    ggtitle(title)
}

p3_cb <- plot_pca_cb(pca_df3, pca3, "Scenario 3: Final Aesthetic Version ✨")
p3_cb



#Testing with centroids-----------------------------------------

library(dplyr)
library(ggplot2)

# Color-blind-friendly palette
timepoint_palette <- c(
  "I"   = "#e9a3c9",  # pink
  "II"  = "lightgoldenrod2",  # yellow
  "III" = "#a1d76a"   # green
)

# Ensure consistent group structure for ellipses
pca_df3$timepoint <- factor(pca_df3$timepoint, levels = c("I", "II", "III"))


plot_pca_cb <- function(pca_df, pca_obj, title = "Global Gene Expression By Life Stage and Species") {
  percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))[1:2]
  
  ggplot(pca_df, aes(x = PC1, y = PC2)) +
    # Add confidence ellipses by timepoint
    stat_ellipse(
      data = pca_df,
      mapping = aes(x = PC1, y = PC2, group = timepoint, fill = timepoint),
      geom = "polygon",
      alpha = 0.3,
      color = NA,
      level = 0.95,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    # Sample points with shape and fill
    geom_point(aes(fill = timepoint, shape = Species), size = 5, color = "black", stroke = 0.8) +
    
    # Color and shape scales
    scale_fill_manual(
      values = timepoint_palette,
      name = "Timepoint"
    ) +
    scale_shape_manual(
      values = c(
        "Acropora_tenuis" = 21,
        "Montipora_capitata" = 22,
        "Pocillopora_acuta" = 24,
        "Stylophora_pistillata" = 23
      ),
      labels = c(
        "Acropora_tenuis" = expression(italic("A. tenuis")),
        "Montipora_capitata" = expression(italic("M. capitata")),
        "Pocillopora_acuta" = expression(italic("P. acuta")),
        "Stylophora_pistillata" = expression(italic("S. pistillata"))
      ),
      name = "Species"
    ) +
    
    # Tidy legends
    guides(
      fill = guide_legend(override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(override.aes = list(fill = "grey80", color = "black"))
    ) +
    
    # Labels + Theme
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(title) +
    coord_fixed() +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}

p3_cb <- plot_pca_cb(pca_df3, pca3)
p3_cb


ggsave("PCA_Global_Expression.png", plot = p3_cb, width = 7, height = 5.5, dpi = 600)


#Getting PCA loadings

# Get loadings matrix (genes x PCs)
loadings <- pca3$rotation

# Convert to data frame for easy manipulation
loadings_df <- as.data.frame(loadings)
loadings_df$gene <- rownames(loadings_df)

# Rank genes by absolute loading for PC1 and PC2
top_pc1 <- loadings_df %>%
  arrange(desc(abs(PC1))) %>%
  slice(1:20)

top_pc2 <- loadings_df %>%
  arrange(desc(abs(PC2))) %>%
  slice(1:20)

top_both <- loadings_df %>%
  mutate(total_contrib = abs(PC1) + abs(PC2)) %>%
  arrange(desc(total_contrib)) %>%
  slice(1:20)

head(top_pc1)
head(top_pc2)


# Get PCA loadings (genes x PCs)
loadings <- pca3$rotation

# Convert to data frame
loadings_df <- as.data.frame(loadings)
loadings_df$Orthogroup <- rownames(loadings_df)

# Top 50 genes by absolute loading on PC1
top_pc1 <- loadings_df %>%
  arrange(desc(abs(PC1))) %>%
  slice(1:50)

# Top 50 genes by absolute loading on PC2
top_pc2 <- loadings_df %>%
  arrange(desc(abs(PC2))) %>%
  slice(1:50)

top_pc1_pc2 <- loadings_df %>%
  mutate(total_loading = abs(PC1) + abs(PC2)) %>%
  arrange(desc(total_loading)) %>%
  slice(1:50)

write.csv(top_pc1, "top50_PC1_loadings.csv", row.names = FALSE)
write.csv(top_pc2, "top50_PC2_loadings.csv", row.names = FALSE)
write.csv(top_pc1_pc2, "top50_PC1_PC2_combined.csv", row.names = FALSE)


