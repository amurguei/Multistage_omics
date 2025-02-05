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



library(installr)
updateR()
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
library("edgeR")

setwd("C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA")
#setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")

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

# Normalize our read counts using VST-normalization in DESeq2
# Construct the DESeq2 dataset

#Merge the timepoint columns into a new column, group. Set group as a factor.
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I","II","III"))
treatmentinfo$Species <- factor(treatmentinfo$Species, levels = c("Acropora_tenuis", "Montipora_capitata", "Pocillopora_acuta","Stylophora_pistillata"))


gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = treatmentinfo,
                               design = ~ Species * timepoint)
# Log-transform the count data using a variance stabilizing transforamtion (vst). 

#SF.gdds <- estimateSizeFactors( gdds ) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than for to use vst
#print(sizeFactors(SF.gdds)) #View size factors
#print(sizeFactors(SF.gdds)) #View size factors
#Larva1_Atenuis Larva2_Atenuis Larva3_Atenuis  Meta1_Atenuis  Meta2_Atenuis  Meta3_Atenuis  Spat1_Atenuis  Spat2_Atenuis 
#2.1484437      2.5045719      2.4650281      1.6930436      1.8298762      1.4530130      2.2826295      2.4530962 
#Spat3_Atenuis    Larva1_Mcap    Larva2_Mcap    Larva3_Mcap     Meta1_Mcap     Meta2_Mcap     Meta3_Mcap     Spat1_Mcap 
#1.9067000      2.5355610      3.0515890      3.0873104      2.9781568      2.5851077      3.0724323      0.9420058 
#Spat2_Mcap     Spat3_Mcap    Larva1_Pacu    Larva2_Pacu    Larva3_Pacu     Meta1_Pacu     Meta2_Pacu     Meta3_Pacu 
#1.2404501      1.6396740      0.6839530      0.9259933      1.0532392      0.6829982      0.7023775      0.7466266 
#Spat1_Pacu     Spat2_Pacu     Spat3_Pacu    Larva1_Spis    Larva2_Spis    Larva3_Spis     Meta1_Spis     Meta2_Spis 
#0.9111311      0.8727706      0.8633586      0.2405015      0.1598825      0.2085578      0.2252045      0.3610082 
#Meta3_Spis     Spat1_Spis     Spat2_Spis     Spat3_Spis 
#0.4935508      0.4431366      0.4204062      0.4138748 

#gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size

# Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.

datExpr <- as.data.frame(t(assay(gdds))) #transpose to output to a new data frame with the column names as row names. And make all data numeric

#Check for genes and samples with too many missing values with goodSamplesGenes. There shouldn't be any because we performed pre-filtering
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #Should return TRUE if not, the R chunk below will take care of flagged data
# [1] allOK is TRUE

#Remove flagged samples if the allOK is FALSE, not used here.  

#ncol(datExpr) #number genes before
#if (!gsg$allOK) #If the allOK is FALSE...
#{
# Optionally, print the gene and sample names that are flagged:
#if (sum(!gsg$goodGenes)>0)
#printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
#if (sum(!gsg$goodSamples)>0)
#printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
#datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
#}
#ncol(datExpr) #number genes after

#Cluster the samples to look for obvious outliers
#Look for outliers by examining the sample tree:

sampleTree = hclust(dist(datExpr), method = "average")

#Plot the sample tree
pdf(paste0('sampleTree_noTransform','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()


#PCA--------------------------------------------------------------------------------------------------------------
counts_matrix <- assay(gdds) # Extract raw counts
pca <- prcomp(t(counts_matrix)) # PCA
plot(pca$x[,1], pca$x[,2]) # Quick scatterplot


gPCAdata <- plotPCA(gdds, intgroup = c("timepoint", "Species"), returnData=TRUE, ntop=5000)

percentVar <- round(100 * attr(gPCAdata, "percentVar"))


# Load necessary libraries
library(DESeq2)
library(ggplot2)

log_counts <- log2(counts(gdds, normalized=FALSE) + 1) # Log-transform raw counts
pca_res <- prcomp(t(log_counts), center=TRUE, scale.=FALSE)


# Compute percent variance explained
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)

# Create PCA dataframe
pca_data <- as.data.frame(pca_res$x)
pca_data$timepoint <- colData(gdds)$timepoint
pca_data$Species <- colData(gdds)$Species

# Ensure timepoint and Species are factors with correct levels
pca_data$timepoint <- factor(pca_data$timepoint, levels = c("I", "II", "III"))
pca_data$Species <- factor(pca_data$Species, levels = c("Acropora_tenuis", "Montipora_capitata", "Pocillopora_acuta", "Stylophora_pistillata"))

# Create PCA visualization with ggplot (same aesthetics as before)
allgenesfilt_PCA_visual <- ggplot(pca_data, aes(x = PC1, y = PC2, color = timepoint, shape = Species)) + 
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(-50, 50) +
  coord_fixed() +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank()
  ) +
  scale_shape_manual(values = c(15, 16, 17, 18)) # Different shapes for each species

# Print the PCA plot
print(allgenesfilt_PCA_visual)

#Combat_seq + untransformed data-------------------------------------------------
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
#print(sizeFactors(SF.gdds)) #View size factors
#Larva1_Atenuis Larva2_Atenuis Larva3_Atenuis  Meta1_Atenuis  Meta2_Atenuis  Meta3_Atenuis  Spat1_Atenuis  Spat2_Atenuis 
#2.1484437      2.5045719      2.4650281      1.6930436      1.8298762      1.4530130      2.2826295      2.4530962 
#Spat3_Atenuis    Larva1_Mcap    Larva2_Mcap    Larva3_Mcap     Meta1_Mcap     Meta2_Mcap     Meta3_Mcap     Spat1_Mcap 
#1.9067000      2.5355610      3.0515890      3.0873104      2.9781568      2.5851077      3.0724323      0.9420058 
#Spat2_Mcap     Spat3_Mcap    Larva1_Pacu    Larva2_Pacu    Larva3_Pacu     Meta1_Pacu     Meta2_Pacu     Meta3_Pacu 
#1.2404501      1.6396740      0.6839530      0.9259933      1.0532392      0.6829982      0.7023775      0.7466266 
#Spat1_Pacu     Spat2_Pacu     Spat3_Pacu    Larva1_Spis    Larva2_Spis    Larva3_Spis     Meta1_Spis     Meta2_Spis 
#0.9111311      0.8727706      0.8633586      0.2405015      0.1598825      0.2085578      0.2252045      0.3610082 
#Meta3_Spis     Spat1_Spis     Spat2_Spis     Spat3_Spis 
#0.4935508      0.4431366      0.4204062      0.4138748 

#gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size

# Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.

datExpr <- as.data.frame(t(assay(gdds))) #transpose to output to a new data frame with the column names as row names. And make all data numeric

#Check for genes and samples with too many missing values with goodSamplesGenes. There shouldn't be any because we performed pre-filtering
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK #Should return TRUE if not, the R chunk below will take care of flagged data
# [1] allOK is TRUE

#Remove flagged samples if the allOK is FALSE, not used here.  

#ncol(datExpr) #number genes before
#if (!gsg$allOK) #If the allOK is FALSE...
#{
# Optionally, print the gene and sample names that are flagged:
#if (sum(!gsg$goodGenes)>0)
#printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
#if (sum(!gsg$goodSamples)>0)
#printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
#datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
#}
#ncol(datExpr) #number genes after

#Cluster the samples to look for obvious outliers
#Look for outliers by examining the sample tree:

sampleTree = hclust(dist(datExpr), method = "average")

#Plot the sample tree
pdf(paste0('sampleTree_ComBat_noTrans','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()


#PCA--------------------------------------------------------------------------------------------------------------
gPCAdata <- plotPCA(gvst, intgroup = c("timepoint", "Species"), returnData=TRUE, ntop=5000)

percentVar <- round(100 * attr(gPCAdata, "percentVar"))



#CombatSeq + VST-----------------------------------------------------------------------------------

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

# Load edgeR package for TMM normalization
library(edgeR)

# Merge the timepoint columns into a new column, group. Set group as a factor.
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I", "II", "III"))
treatmentinfo$Species <- factor(treatmentinfo$Species, levels = c("Acropora_tenuis", "Montipora_capitata", "Pocillopora_acuta", "Stylophora_pistillata"))

# Create the DESeqDataSet object (no transformation yet)
gdds <- DESeqDataSetFromMatrix(countData = gcount_corrected,
                               colData = treatmentinfo,
                               design = ~ Species * timepoint)

# Use edgeR to calculate TMM normalization factors
dgelist <- DGEList(counts = gcount_corrected)  # Create DGEList object
dgelist <- calcNormFactors(dgelist, method = "TMM")  # Calculate TMM normalization factors

# You can check the normalization factors to see if they make sense
print(dgelist$samples$norm.factors)

# Apply TMM normalization and transform the counts (no VST, just normalizing by TMM)
norm_counts <- cpm(dgelist, normalized.lib.sizes = TRUE)  # CPM with TMM normalization

# Transpose the normalized gene count matrix so that the gene IDs are rows and the sample IDs are columns
datExpr <- as.data.frame(t(norm_counts))  # Transpose to output to a new data frame with column names as row names, and make all data numeric

# Check for genes and samples with too many missing values with goodSamplesGenes. There shouldn't be any because we performed pre-filtering
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK  # Should return TRUE if not, the R chunk below will take care of flagged data

variances <- apply(norm_counts, 1, var)
hist(variances)  # Check the distribution of variances
pca_result <- prcomp(t(norm_counts_top))  # Ensure samples are columns, genes are rows


# Select top 5000 most variable genes based on the normalized data
# You can use the row variances of the normalized counts to select the top genes
variances <- apply(norm_counts, 1, var)  # Calculate the variance for each gene
top_genes <- order(variances, decreasing = TRUE)[1:500]  # Select top 5000 genes by variance
norm_counts_top <- norm_counts[top_genes, ]  # Subset the top genes

# Perform PCA
pca_result <- prcomp(t(norm_counts_top))  # Use transposed data for PCA (samples as columns)

# Get the percentage of variance explained by each principal component
percentVar <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)))

# Prepare PCA data for ggplot
gPCAdata <- data.frame(pca_result$x)  # Get the PCA results (coordinates of samples in PCA space)
gPCAdata$timepoint <- treatmentinfo$timepoint  # Add the 'timepoint' and 'Species' columns
gPCAdata$Species <- treatmentinfo$Species

# Ensure timepoint and species are factors with correct levels
gPCAdata$timepoint <- factor(gPCAdata$timepoint, levels = c("I", "II", "III"))
gPCAdata$Species <- factor(gPCAdata$Species, levels = c("Acropora_tenuis", "Montipora_capitata", "Pocillopora_acuta", "Stylophora_pistillata"))

# Create PCA visualization
allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(x = PC1, y = PC2, color = timepoint, shape = Species)) + 
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(-50, 50) +
  coord_fixed() +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank()
  ) +
  scale_shape_manual(values = c(15, 16, 17, 18))  # Different shapes for each species

# Print the PCA plot
print(allgenesfilt_PCA_visual)

variances <- apply(norm_counts, 1, var)
hist(variances)  # Check the distribution of variances



gPCAdata <- data.frame(pca_result$x)  # Use the PCA result to create the data for plotting
gPCAdata$timepoint <- treatmentinfo$timepoint
gPCAdata$Species <- treatmentinfo$Species

ggplot(data = gPCAdata, aes(x = PC1, y = PC2, color = timepoint, shape = Species)) + 
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(-50, 50) +
  coord_fixed() +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank()
  ) +
  scale_shape_manual(values = c(15, 16, 17, 18))  # Different shapes for each species

sum(is.na(gPCAdata$PC1))  # Check for NAs in PC1
sum(is.na(gPCAdata$PC2))  # Check for NAs in PC2
sum(is.na(gPCAdata$timepoint))  # Check for NAs in timepoint
sum(is.na(gPCAdata$Species))  # Check for NAs in Species


