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

#Set working directory
setwd("C:/Users/amurg/OneDrive/Desktop/PhD Haifa/Multistage Omics/R scripts/A_tenuis")
library(readr)
#Load data
treatmentinfo <- read_csv("Sample_info_A_tenuis.csv")
gcount <- read_csv("gene_count_matrix_Acropora.csv")
#Some formatting
colnames(gcount)[colnames(gcount) == "Geneid"] <- "gene_id"
gcount <- gcount[, -c(2:6)]
gcount <- as.data.frame(gcount)
rownames(gcount) <- gcount$gene_id
gcount <- gcount %>% mutate(across(-gene_id, as.integer))
gcount <- gcount[, -1] # Remove `gene_id` column, as it's now in the row names


# Set filter values for PoverA
# Smallest sample size per treatment is 4, so 4/12 (12 samples) is 0.33 (rounded to 2 decimal places)
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
nrow(gcount) #Before 30327
nrow(gcount_filt) #After 18917

# Normalize our read counts using VST-normalization in DESeq2
# Construct the DESeq2 dataset

#Merge Set group as a factor.
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I","II","III"))

#Create a DESeqDataSet design from gene count matrix and labels. Here we set the design to look at 
#any differences in gene expression across samples attributed to depth.

#Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = treatmentinfo,
                               design = ~timepoint)

# Log-transform the count data using a variance stabilizing transforamtion (vst). 

SF.gdds <- estimateSizeFactors( gdds ) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than four to use vst
print(sizeFactors(SF.gdds)) #View size factors
#DRR318288 DRR318289 DRR318290 DRR318291 DRR318292 DRR318293 DRR318294 DRR318295 DRR318296 DRR318297 DRR318298 DRR318299 
#1.1187900 1.2877618 1.2765853 1.2701737 0.8221285 0.8614881 0.9337686 0.7456923 1.0930993 1.0393193 1.1166223 0.8606887
gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size

#Compile WGCNA dataset--------------------------------------------------------------------------
# Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.

datExpr <- as.data.frame(t(assay(gvst))) #transpose to output to a new data frame with the column names as row names. And make all data numeric

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

gPCAdata <- plotPCA(gvst, intgroup = c("timepoint"), returnData=TRUE, ntop=18917) #use ntop to specify all genes

percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data

# Explicitly set the levels of the timepoint variable
gPCAdata$timepoint <- factor(gPCAdata$timepoint, levels = c("I", "II", "III"))

allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(PC1, PC2)) + 
  geom_point(aes(shape = timepoint, colour = timepoint), size = 6) +  # Increase point size
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(-50, 50) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),            # Increase axis text size
    axis.title = element_text(size = 18, face = "bold"),  # Increase axis label size and make it bold
    legend.text = element_text(size = 14),          # Increase legend text size
    legend.title = element_text(size = 16),         # Increase legend title size
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank()
  )

print(allgenesfilt_PCA_visual)


library(ggrepel) # Optional but helps avoid overlapping text

allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(PC1, PC2)) + 
  geom_point(aes(shape = timepoint, colour = timepoint), size = 6) +  # Increase point size
  geom_text_repel(aes(label = rownames(gPCAdata)), size = 4) +       # Add sample names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(-50, 50) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),            # Increase axis text size
    axis.title = element_text(size = 18, face = "bold"),  # Increase axis label size and make it bold
    legend.text = element_text(size = 14),          # Increase legend text size
    legend.title = element_text(size = 16),         # Increase legend title size
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank()
  )

print(allgenesfilt_PCA_visual)

# Cluster the samples to look for obvious outliers
#Look for outliers by examining the sample tree:

sampleTree = hclust(dist(datExpr), method = "average")

# Plot the sample tree
pdf(paste0('sampleTree','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()


#Filtered data: 

# Samples to remove
remove_samples <- c("DRR318292", "DRR318296", "DRR318290")

# Filter metadata
treatmentinfo <- treatmentinfo[!treatmentinfo$sampleID %in% remove_samples, ]

# Filter expression data
gvst_filtered <- gvst[, !colnames(gvst) %in% remove_samples]

# Transpose expression matrix for WGCNA and PERMANOVA
datExpr <- as.data.frame(t(assay(gvst_filtered)))

# Make sure sample names match
treatmentinfo <- treatmentinfo[match(rownames(datExpr), treatmentinfo$sampleID), ]
stopifnot(all(rownames(datExpr) == treatmentinfo$sampleID))

#PCA with filtered data
gPCAdata <- plotPCA(gvst_filtered, intgroup = c("timepoint"), returnData = TRUE, ntop = ncol(datExpr))


percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data

# Explicitly set the levels of the timepoint variable
gPCAdata$timepoint <- factor(gPCAdata$timepoint, levels = c("I", "II", "III"))

allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(PC1, PC2)) + 
  geom_point(aes(shape = timepoint, colour = timepoint), size = 6) +  # Increase point size
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(-50, 50) +
  coord_fixed() +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16),            # Increase axis text size
    axis.title = element_text(size = 18, face = "bold"),  # Increase axis label size and make it bold
    legend.text = element_text(size = 14),          # Increase legend text size
    legend.title = element_text(size = 16),         # Increase legend title size
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank()
  )

print(allgenesfilt_PCA_visual)

library(vegan)

# Bray-Curtis distance matrix
bray_dist <- vegdist(datExpr, method = "bray")

# Make sure timepoint is a factor
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I", "II", "III"))

# Run PERMANOVA
permanova_result <- adonis2(datExpr ~ timepoint,
                            data = treatmentinfo,
                            method = "bray",
                            permutations = 9999,
                            by = "term")

print(permanova_result)

