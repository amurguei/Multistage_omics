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



#library(installr)
#updateR()
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
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


# Load necessary library
library(readr)
library(WGCNA)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(magrittr)

# Set working directory
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")

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

# Renaming columns to make them more manageable
old_names <- colnames(gcount_raw)[7:ncol(gcount_raw)]
new_names <- treatmentinfo$new_names[match(old_names, treatmentinfo$sampleID)]
colnames(gcount_raw)[7:ncol(gcount_raw)] <- ifelse(is.na(new_names), old_names, new_names)

# Switching the order for Spis spat
spis_spat_cols <- 34:36
spis_larvae_cols <- 40:42
gcount_raw <- gcount_raw[, c(1:33, spis_larvae_cols, 37:39, spis_spat_cols)]

# Create a dataframe and call Orthogroup rownames
gcount_raw <- as.data.frame(gcount_raw)
rownames(gcount_raw) <- gcount_raw$Orthogroup

# Extract count data
gcount <- gcount_raw[, 7:42]

# Set filter values for PoverA
filt <- filterfun(pOverA(0.083, 10))

# Apply filter to count data
gfilt <- genefilter(gcount, filt)

# Identify genes to keep by count filter
gkeep <- gcount[gfilt, ]
gn.keep <- rownames(gkeep)

# Subset the original dataset to keep metadata and filtered count data
gcount_filt <- gcount[rownames(gcount) %in% gn.keep, ]

# Remove batch effects with ComBat-Seq
library(sva)
batch_info <- treatmentinfo$Species[match(colnames(gcount_filt), treatmentinfo$new_names)]
batch_info <- as.factor(batch_info)

# Apply ComBat-seq batch correction
gcount_corrected <- ComBat_seq(counts = as.matrix(gcount_filt), batch = batch_info, group = NULL)

# Normalize read counts using VST normalization in DESeq2
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I", "II", "III"))
treatmentinfo$Species <- factor(treatmentinfo$Species, levels = c("Acropora_tenuis", "Montipora_capitata", "Pocillopora_acuta", "Stylophora_pistillata"))

gdds <- DESeqDataSetFromMatrix(countData = gcount_corrected, colData = treatmentinfo, design = ~ Species * timepoint)

SF.gdds <- estimateSizeFactors(gdds)  # Estimate size factors
print(sizeFactors(SF.gdds))  # View size factors
print(sizeFactors(SF.gdds))  # View size factors
#Larva1_Atenuis Larva2_Atenuis Larva3_Atenuis  Meta1_Atenuis  Meta2_Atenuis  Meta3_Atenuis 
#1.8633992      2.1859570      2.1389978      1.2525420      1.3537586      1.1058002 
#Spat1_Atenuis  Spat2_Atenuis  Spat3_Atenuis    Larva1_Mcap    Larva2_Mcap    Larva3_Mcap 
#1.7772996      1.9022468      1.4566066      2.4887247      3.0634007      3.1470177 
#Meta1_Mcap     Meta2_Mcap     Meta3_Mcap     Spat1_Mcap     Spat2_Mcap     Spat3_Mcap 
#3.0737340      2.6238828      3.1272286      0.7344009      0.9863250      1.5349595 
#Larva1_Pacu    Larva2_Pacu    Larva3_Pacu     Meta1_Pacu     Meta2_Pacu     Meta3_Pacu 
#0.6585072      0.8927192      1.0516380      0.6744375      0.6988360      0.7895888 
#Spat1_Pacu     Spat2_Pacu     Spat3_Pacu    Larva1_Spis    Larva2_Spis    Larva3_Spis 
#0.9387776      0.9030795      0.8710691      0.3277546      0.2142616      0.2788691 
#Meta1_Spis     Meta2_Spis     Meta3_Spis     Spat1_Spis     Spat2_Spis     Spat3_Spis 
#0.3028476      0.4435857      0.5622346      0.4970006      0.4693681      0.4603734 
 

gvst <- vst(gdds, blind = FALSE)  # Apply variance stabilizing transformation

# Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.
datExpr1 <- as.data.frame(t(assay(gvst)))  # Transpose

# Center the dataset (mean = 0 for each gene) but do not scale (keep variance as it is)
datExpr1_centered <- scale(datExpr1, center = TRUE, scale = FALSE)

# Cluster the samples to look for obvious outliers
sampleTree <- hclust(dist(datExpr1_centered), method = "average")

# Plot the sample tree
pdf(paste0('sampleTree_ComBat_centered','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# PCA
gPCAdata <- plotPCA(gvst, intgroup = c("timepoint", "Species"), returnData=TRUE)
percentVar <- round(100 * attr(gPCAdata, "percentVar"))

# Ensure timepoint and species are factors with correct levels
gPCAdata$timepoint <- factor(gPCAdata$timepoint, levels = c("I", "II", "III"))
gPCAdata$Species <- factor(gPCAdata$Species, levels = c("Acropora_tenuis", "Montipora_capitata", "Pocillopora_acuta", "Stylophora_pistillata"))

# Create PCA visualization
allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(x = PC1, y = PC2, color = timepoint, shape = Species)) + 
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
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


#PERMANOVA
library(vegan)

# Make sure sample order in metadata matches expression matrix
treatmentinfo <- treatmentinfo[match(rownames(datExpr1), treatmentinfo$new_names), ]

# Confirm alignment
stopifnot(all(rownames(datExpr1) == treatmentinfo$new_names))

# PERMANOVA with interaction term
permanova_result <- adonis2(datExpr1 ~ Species * timepoint, 
                            data = treatmentinfo, 
                            method = "bray", 
                            permutations = 9999,
                            by = "term")


print(permanova_result)

bray_dist <- vegdist(datExpr1, method = "bray")

# By Species
bd_species <- betadisper(bray_dist, treatmentinfo$Species)
anova(bd_species)
#Analysis of Variance Table

#Response: Distances
#Df     Sum Sq    Mean Sq F value Pr(>F)
#Groups     3 0.00001408 4.6930e-06  0.0693 0.9759
#Residuals 32 0.00216610 6.7691e-05               

# By Timepoint
bd_timepoint <- betadisper(bray_dist, treatmentinfo$timepoint)
anova(bd_timepoint)
#Analysis of Variance Table

#Response: Distances
#Df     Sum Sq   Mean Sq F value Pr(>F)
#Groups     2 0.00004896 2.448e-05  1.0511  0.361
#Residuals 33 0.00076858 2.329e-05 

# By interaction
treatmentinfo$group <- interaction(treatmentinfo$Species, treatmentinfo$timepoint)
bd_interaction <- betadisper(bray_dist, treatmentinfo$group)
anova(bd_interaction)
#anova(bd_interaction)
#Analysis of Variance Table

#Response: Distances
#Df     Sum Sq    Mean Sq F value  Pr(>F)  
#Groups    11 0.00035584 3.2349e-05  3.0731 0.01038 *
#  Residuals 24 0.00025264 1.0527e-05                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Boxplot of distances to centroids for Species × Timepoint

library(pairwiseAdonis)

# Use same distance matrix (Bray-Curtis on uncentered VST)
bray_dist <- vegdist(datExpr1, method = "bray")

# Make sure interaction group exists
treatmentinfo$group <- interaction(treatmentinfo$Species, treatmentinfo$timepoint)

# Run pairwise PERMANOVA across all group pairs
pairwise_results <- pairwise.adonis2(bray_dist ~ group, data = treatmentinfo, 
                                     permutations = 9999, p.adjust.m = "BH")  # BH = Benjamini-Hochberg correction
print(pairwise_results)

timepoints <- levels(treatmentinfo$timepoint)
pairwise_species_by_stage <- lapply(timepoints, function(stage) {
  subset_data <- datExpr1[treatmentinfo$timepoint == stage, ]
  subset_meta <- droplevels(treatmentinfo[treatmentinfo$timepoint == stage, ])
  dist_stage <- vegdist(subset_data, method = "bray")
  pairwise.adonis2(dist_stage ~ Species, data = subset_meta, permutations = 9999, p.adjust.m = "BH")
})
names(pairwise_species_by_stage) <- timepoints

print(pairwise_species_by_stage)

boxplot(bd_interaction,
        main = "Distances to Group Centroids (Species × Timepoint)",
        xlab = "Species × Timepoint",
        ylab = "Distance to Centroid",
        las = 2,             # Rotate x-axis labels
        cex.axis = 0.8,      # Shrink axis label size if too crowded
        col = "lightblue")


# Ordination plot of dispersions
plot(bd_interaction,
     main = "Dispersion Ordination (PCoA of Bray-Curtis Distances)",
     hull = FALSE,     # Don’t draw convex hulls
     ellipse = TRUE,   # Show ellipses for each group
     label = TRUE)     # Label group centroids


# Choose a set of soft-thresholding powers
powers <- c(seq(from = 1, to = 19, by = 1), c(21:30)) 
powerVector <- c(seq(1, 10, by = 1), seq(12, 20, by = 1))

# Load the doParallel package
library(doParallel)

# Register a parallel backend
registerDoParallel(cores = 4)  # Set the number of cores

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr1_centered, powerVector = powers, verbose = 5)

# Stop the parallel backend when done
# stopImplicitCluster()

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(paste0('network+combat+centering', '.pdf'))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

#The lowest scale-free topology fit index R^2 recommended by Langfelder and Horvath is 0.8. 
#From the graph, it appears that our soft thresholding power is 19 because it is the lowest 
#power before the R^2=0.8 threshold that maximizes with model fit (s17 is the number right above the red line).

### Network construction and module detection:
# Co-expression adjacency and topological overlap matrix similarity
# Co-expression similarity and adjacency, using the soft thresholding power 19 and translate the adjacency into topological overlap matrix to calculate 
# the corresponding dissimilarity. I will use a signed network because we have a relatively high softPower, according 
# to >12 (https://peterlangfelder.com/2018/11/25/__trashed/). 
#Moreover, in expression data where you are interested in when expression on one gene increases or decreases with expression level of another you would use a signed network (when you are interested in the direction of change, correlation and anti-correlation, you use a signed network).

options(stringsAsFactors = FALSE)
enableWGCNAThreads() #Allow multi-threading within WGCNA


#Version two of this graph (with dots)
dev.new()


a1 <- ggplot(sft$fitIndices, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft$fitIndices, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)


softPower=7 #Set softPower to 7
adjacency=adjacency(datExpr1_centered, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(adjacency, file = "C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/adjacency_centered7sft.RData")
save(TOM, file ="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/TOM_7sft.RData")
save(dissTOM, file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/disstomsft7_centered.RData") 
load(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/disstomsft7.RData")

#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Multisp_dissTOMClustering_centered7sft.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()

# Module identification is cutting the branches off the tree in the dendrogram above. We want large modules, so we set the minimum module size 
# relatively high (minimum size = 30).

minModuleSize = 30 #default value used most often
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods) #list modules and respective sizes
#dynamicMods
#dynamicMods
#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24 
#1440  540  457  316  304  290  289  273  230  212  210  195  195  169  160  157  157  119  114  112  104  101   99   95 
#25   26   27   28   29   30   31   32   33   34   35   36   37   38 
#93   78   64   62   60   56   52   51   51   47   46   45   39   36 

dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
#table(dynamicColors)
#dynamicColors
#black           blue          brown           cyan      darkgreen       darkgrey    darkmagenta darkolivegreen 
#289            540            457            169            101             95             47             51 
#darkorange        darkred  darkturquoise          green    greenyellow         grey60      lightcyan     lightgreen 
#78            104             99            304            210            157            157            119 
#lightyellow        magenta   midnightblue         orange  paleturquoise           pink          plum1         purple 
#114            230            160             93             52            273             36            212 
#red      royalblue    saddlebrown         salmon        sienna3        skyblue       skyblue3      steelblue 
#290            112             60            195             46             62             39             56 
#tan      turquoise         violet          white         yellow    yellowgreen 
#195           1440             51             64            316             45 

pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/dissTOMColorClustering_centered7sft.pdf", width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar or choose not to merge
# Plot module similarity based on eigengene value

#Calculate eigengenes
MEList = moduleEigengenes(datExpr1_centered, colors = dynamicColors, softPower = 7)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster again and plot the results
METree = flashClust(as.dist(MEDiss), method = "average")

pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/eigengeneClustering1_WGCNA.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >85% eigengene similarity (most studies use80-90% similarity)

MEDissThres= 0.15 #merge modules that are 85% similar

pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/sft7eigengeneClustering2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr1_centered, dynamicColors, cutHeight= MEDissThres, verbose =3)

# Check the structure of the merge object
print(names(merge))
#[1] "colors"    "dendro"    "oldDendro" "cutHeight" "oldMEs"    "newMEs"    "allOK"    

# Assign mergedColors
mergedColors = merge$colors

# Check if mergedColors is assigned correctly
print(head(mergedColors))
#[1] "darkgrey"   "darkred"    "skyblue"    "lightgreen" "turquoise"  "turquoise" 

# Plot dendrogram and colors
pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/sft7mergedClusters.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
dev.off()

#Save new colors

moduleColors = mergedColors # Rename to moduleColors
colorOrder = c("grey", standardColors(50)); # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;

# Assign merged module eigengenes
MEs = merge$newMEs

# Check if MEs is assigned correctly
print(head(MEs))

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/sft7eigengeneClustering3.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# Relating modules to group, quantifying moduletrait associations
#Prepare trait data. Data has to be numeric, so I replaced the group for numeric values
colnames(treatmentinfo)
str(treatmentinfo)
head(treatmentinfo$group)

allTraits <- names(treatmentinfo$group)
allTraits$Larvae_Atenuis <- c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
allTraits$Meta_Atenuis   <- c(0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
allTraits$Spat_Atenuis   <- c(0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
allTraits$Larvae_Mcap    <- c(0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
allTraits$Meta_Mcap      <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
allTraits$Spat_Mcap      <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
allTraits$Larvae_Pacu    <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
allTraits$Meta_Pacu      <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
allTraits$Spat_Pacu      <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0)
allTraits$Larvae_Spis    <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0)
allTraits$Meta_Spis      <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0)
allTraits$Spat_Spis      <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1)

datTraits <- as.data.frame(allTraits)
dim(datTraits)
#[1] 36 12

rownames(datTraits) <- treatmentinfo$sample_id
print(datTraits)
#This is creating an issue

#Define numbers of genes and samples
nGenes = ncol(datExpr1_centered)
nSamples = nrow(datExpr1_centered)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr1_centered, moduleColors,softPower=17)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)
#[1] "MEsteelblue"      "MEturquoise"      "MEdarkred"        "MEpurple"         "MElightgreen"    
#[6] "MEdarkgrey"       "MEskyblue"        "MEmidnightblue"   "MEmagenta"        "MEtan"           
#[11] "MEdarkgreen"      "MEroyalblue"      "MEblack"          "MEgrey60"         "MEdarkmagenta"   
#[16] "MEsienna3"        "MEyellow"         "MEcyan"           "MElightcyan"      "MEdarkolivegreen"
#[21] "MEskyblue3"       "MEbrown"          "MEpink"           "MEred"            "MEsaddlebrown"   
#[26] "MEdarkturquoise"  "MEorange"         "MEwhite"          "MEyellowgreen"    "MEblue"          
#[31] "MEdarkorange"     "MElightyellow"    "MEsalmon"         "MEviolet"         "MEpaleturquoise" 
#[36] "MEplum1"


moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/sft7_Life-stage clustering based on module-trait correlation.pdf", width = 20)
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr1_centered)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)
#head(textMatrix)
#[,1]             [,2]             [,3]             [,4]            [,5]           
#[1,] "-0.52\n(0.001)" "-0.52\n(0.001)" "-0.52\n(0.001)" "0.14\n(0.4)"   "0.14\n(0.4)"  
#[2,] "0.17\n(0.3)"    "0.12\n(0.5)"    "-0.3\n(0.07)"   "-0.25\n(0.1)"  "-0.3\n(0.08)" 
#[3,] "0.1\n(0.6)"     "0.045\n(0.8)"   "-0.19\n(0.3)"   "-0.29\n(0.09)" "-0.33\n(0.05)"
#[4,] "-0.021\n(0.9)"  "0.084\n(0.6)"   "-0.11\n(0.5)"   "-0.27\n(0.1)"  "-0.3\n(0.07)" 
#[5,] "0.13\n(0.5)"    "0.1\n(0.6)"     "-0.26\n(0.1)"   "-0.34\n(0.04)" "-0.29\n(0.09)"
#[6,] "-0.22\n(0.2)"   "0.44\n(0.007)"  "-0.22\n(0.2)"   "-0.35\n(0.03)" "-0.18\n(0.3)" 
#[,6]            [,7]            [,8]            [,9]            [,10]          [,11]         
#[1,] "0.19\n(0.3)"   "0.21\n(0.2)"   "0.17\n(0.3)"   "0.19\n(0.3)"   "0.21\n(0.2)"  "0.15\n(0.4)" 
#[2,] "0.63\n(3e-05)" "0.13\n(0.4)"   "-0.36\n(0.03)" "0.25\n(0.1)"   "-0.17\n(0.3)" "-0.2\n(0.2)" 
#[3,] "0.69\n(3e-06)" "0.36\n(0.03)"  "-0.028\n(0.9)" "-0.32\n(0.05)" "0.037\n(0.8)" "-0.2\n(0.2)" 
#[4,] "0.59\n(2e-04)" "0.37\n(0.02)"  "-0.2\n(0.2)"   "-0.18\n(0.3)"  "0.43\n(0.01)" "-0.11\n(0.5)"
#[5,] "0.68\n(5e-06)" "-0.28\n(0.09)" "0.39\n(0.02)"  "-0.1\n(0.5)"   "-0.19\n(0.3)" "0.12\n(0.5)" 
#[6,] "0.57\n(3e-04)" "-0.34\n(0.04)" "0.039\n(0.8)"  "0.29\n(0.08)"  "-0.24\n(0.2)" "0.13\n(0.4)" 
#[,12]          
#[1,] "0.15\n(0.4)"  
#[2,] "0.27\n(0.1)"  
#[3,] "0.12\n(0.5)"  
#[4,] "-0.28\n(0.09)"
#[5,] "0.053\n(0.8)" 
#[6,] "0.087\n(0.6)" 

pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/sft7_Module-trait-relationships.pdf")
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), 
               cex.lab.y= 0.55, cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , 
               zlim = c(-1,1), main = paste("Module-trait relationships"))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), cex.lab.y= 0.55, 
               cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , zlim = c(-1,1), 
               main = paste("Module-trait relationships"))
dev.off() 

#Clustering Life stage first


# Step 1: Create trait matrix (life stage first, then species)
datTraits_mod <- data.frame(
  Larvae_Atenuis = c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  Larvae_Mcap    = c(0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  Larvae_Pacu    = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  Larvae_Spis    = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0),
  Meta_Atenuis   = c(0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  Meta_Mcap      = c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  Meta_Pacu      = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),
  Meta_Spis      = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0),
  Spat_Atenuis   = c(0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  Spat_Mcap      = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  Spat_Pacu      = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0),
  Spat_Spis      = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1)
)

# Step 2: Assign row names from your metadata
rownames(datTraits_mod) <- treatmentinfo$sample_id  # or treatmentinfo$new_names if renamed

# Step 3: Compute module eigengenes and correlations
nSamples <- nrow(datExpr1_centered)

MEs0 <- moduleEigengenes(datExpr1_centered, moduleColors, softPower = 10)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, datTraits_mod, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Step 4: Create labeled text matrix
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Step 5: Plot the heatmap
pdf(file = "sft7module_trait_heatmap_stage_first.pdf", width = 14, height = 10)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits_mod),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  cex.lab.y = 0.55,
  cex.lab.x = 0.55,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = TRUE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Module–Trait Relationships by Life Stage-Species Groups"
)

dev.off()

# Plot as clustered Heatmap
#add bold sigignificant p-values, dendrogram with WGCNA MEtree cut-off, module clusters

#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)

htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

#library(dendsort)
row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

pdf(file = "sft7v1Module-life stage-relationship-heatmap3.pdf", height = 11.5, width = 8)

ht = Heatmap(moduleTraitCor, 
             name = "Eigengene", 
             column_title = "Module-Life stage Eigengene Correlation", 
             col = blueWhiteRed(50), 
             row_names_side = "left", 
             row_dend_side = "left",
             width = unit(4, "in"), 
             height = unit(8.5, "in"), 
             column_order = 1:ncol(moduleTraitCor),  # Adjusting for 12 columns
             column_dend_reorder = FALSE, 
             cluster_columns = hclust(dist(t(moduleTraitCor)), method = "average"), 
             column_split = 3, 
             column_dend_height = unit(0.5, "in"),
             cluster_rows = METree, 
             row_split = 10, 
             row_gap = unit(2.5, "mm"), 
             border = TRUE,
             cell_fun = function(j, i, x, y, w, h, col) {
               if(heatmappval[i, j] <= 0.05) {
                 grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
               }
               else {
                 grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
               }},
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10, alpha = 0.75, border = TRUE, fill = htmap.colors))
draw(ht)
dev.off()

#Version without column clustering: 

# Prepare p-values
heatmappval <- signif(moduleTraitPvalue, 1)

# Clean up module names (row colors)
htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

# Use pretty italic labels
pretty_labels <- parse(text = c(
  "Larvae~italic('A. tenuis')", "Larvae~italic('M. capitata')", "Larvae~italic('P. acuta')", "Larvae~italic('S. pistillata')",
  "Meta~italic('A. tenuis')",   "Meta~italic('M. capitata')",   "Meta~italic('P. acuta')",   "Meta~italic('S. pistillata')",
  "Spat~italic('A. tenuis')",   "Spat~italic('M. capitata')",   "Spat~italic('P. acuta')",   "Spat~italic('S. pistillata')")
)

# PDF output
pdf(file = "module_trait_heatmap_cleaned.pdf", height = 11.5, width = 8)

ht = Heatmap(
  moduleTraitCor,
  name = "Eigengene",
  column_title = "Module–Life Stage Eigengene Correlation",
  col = blueWhiteRed(50),
  row_names_side = "left",
  row_dend_side = "left",
  width = unit(4, "in"),
  height = unit(8.5, "in"),
  column_order = 1:ncol(moduleTraitCor),
  cluster_columns = FALSE,  #  Disable column clustering
  column_labels = pretty_labels,
  cluster_rows = METree,
  row_split = 10,
  row_gap = unit(2.5, "mm"),
  border = TRUE,
  cell_fun = function(j, i, x, y, w, h, col) {
    if (heatmappval[i, j] <= 0.05) {
      grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
    } else {
      grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
    }
  },
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10, alpha = 0.75, border = TRUE, fill = htmap.colors)  # his is what gives row color boxes
)

draw(ht)
dev.off()

#By Life-stage - species groups

#Same but let´s try now with Dattraits grouped per life stage rather than life stage & species

# Reorder the columns to have larvae first, then meta, and then spat for every species

allTraits_justLifeStage <- list(
  Larvae         = c(1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0),
  Meta           = c(0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0),
  Spat           = c(0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1)
)

# Combine the individual lists into a data frame
datTraits_justLifeStage <- as.data.frame(allTraits_justLifeStage)

print(datTraits_justLifeStage)

#Define numbers of genes and samples
nGenes = ncol(datExpr1_centered)
nSamples = nrow(datExpr1_centered)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr1_centered, moduleColors,softPower=10)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)

moduleTraitCor = cor(MEs, datTraits_justLifeStage, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/sftt7_Life stages clustering based on module-trait correlation v3.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# Required packages
library(dendextend)
library(WGCNA)

# Assuming you have already defined `moduleTraitCor` as:
# moduleTraitCor <- cor(MEs, datTraits_justLifeStage, use = "p")

# Step 1: Cluster as dendrogram
d <- as.dendrogram(hclust(dist(t(moduleTraitCor)), method = "average"))

# Step 2: Rotate to desired order
desired_order <- c("Larvae", "Meta", "Spat")  # Your custom trait column order
d_rotated <- rotate(d, order = desired_order)

# Step 3: Plot the dendrogram
pdf("sft7_Life_stage_dendrogram_rotated.pdf", width = 6, height = 6)
plot(d_rotated,
     main = "Group clustering based on module–trait correlation",
     cex.main = 1.8,
     cex.lab = 1.5,
     cex.axis = 1.2,
     horiz = FALSE)  # Set to TRUE for horizontal layout if preferred
dev.off()


#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr1_centered)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)
head(textMatrix)
#=[,1]             [,2]             [,3]           
#[1,] "0.026\n(0.9)"   "-0.036\n(0.8)"  "0.01\n(1)"    
#[2,] "-0.065\n(0.7)"  "-0.43\n(0.008)" "0.5\n(0.002)" 
#[3,] "0.12\n(0.5)"    "-0.3\n(0.07)"   "0.18\n(0.3)"  
#[4,] "0.3\n(0.08)"    "-0.31\n(0.07)"  "0.0085\n(1)"  
#[5,] "-0.4\n(0.01)"   "0.19\n(0.3)"    "0.22\n(0.2)"  
#[6,] "-0.68\n(5e-06)" "0.25\n(0.1)"    "0.43\n(0.009)"pdf(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Module-life stage relationships no species.pdf", width = 14, height = 10)

labeledHeatmap(
  Matrix = moduleTraitCor, 
  xLabels = names(datTraits_justLifeStage),  
  yLabels = names(MEs), 
  ySymbols = names(MEs), 
  cex.lab.y = 0.55,  # Increase y-axis label size
  cex.lab.x = 0.55,  # Increase x-axis label size
  colors = blueWhiteRed(50), 
  textMatrix = textMatrix, 
  setStdMargins = TRUE, 
  cex.text = 0.5,  # Increase text size inside the heatmap
  zlim = c(-1,1), 
  main = "Module-life stage relationships"
)

dev.off()


# Plot as clustered Heatmap
#add bold sigignificant p-values, dendrogram with WGCNA MEtree cut-off, module clusters

#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)

htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

#library(dendsort)
row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

pdf(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/sft7nosp_Module-life stage-relationship-heatmap3.pdf", height = 11.5, width = 8)

ht = Heatmap(moduleTraitCor, 
             name = "Eigengene", 
             column_title = "Module-Life stage Eigengene Correlation", 
             col = blueWhiteRed(50), 
             row_names_side = "left", 
             row_dend_side = "left",
             width = unit(4, "in"), 
             height = unit(8.5, "in"), 
             column_order = 1:ncol(moduleTraitCor),  # Adjusting for 12 columns
             column_dend_reorder = FALSE, 
             cluster_columns = hclust(dist(t(moduleTraitCor)), method = "average"), 
             column_split = 3, 
             column_dend_height = unit(0.5, "in"),
             cluster_rows = METree, 
             row_split = 10, 
             row_gap = unit(2.5, "mm"), 
             border = TRUE,
             cell_fun = function(j, i, x, y, w, h, col) {
               if(heatmappval[i, j] <= 0.05) {
                 grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
               }
               else {
                 grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
               }},
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 10, alpha = 0.75, border = TRUE, fill = htmap.colors))
draw(ht)
dev.off()


library(ComplexHeatmap)
library(circlize)
library(grid)
library(dendextend)



# Step 1: Rotated top dendrogram in desired order
d <- as.dendrogram(hclust(dist(t(moduleTraitCor)), method = "average"))
d_rotated <- rotate(d, order = c("Larvae", "Metamorphosed", "Spat"))
column_order_indices <- order.dendrogram(d_rotated)

# Step 2: Simplified column labels (parse for better formatting)
pretty_labels <- parse(text = c("Larvae", "Metamorphosed", "Spat"))

# Step 3: Create the heatmap object
ht <- Heatmap(
  matrix = moduleTraitCor,
  name = "Eigengene",
  col = blueWhiteRed(50),
  column_order = column_order_indices,
  column_labels = pretty_labels[column_order_indices],
  cluster_columns = d_rotated,
  column_dend_reorder = FALSE,
  column_dend_height = unit(0.5, "in"),
  row_split = 8,
  cluster_rows = METree,
  column_names_rot = 45,
  row_gap = unit(2.5, "mm"),
  row_names_side = "left",
  row_dend_side = "left",
  border = TRUE,
  width = unit(5, "in"),
  height = unit(9, "in"),
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 10, alpha = 0.75, border = TRUE, fill = htmap.colors),
  cell_fun = function(j, i, x, y, w, h, col) {
    pval <- heatmappval[i, j]
    label <- sprintf("%s", pval)
    font <- if (pval <= 0.05) "bold" else "plain"
    grid.text(label, x, y, gp = gpar(fontsize = 8, fontface = font))
  }
)

# Step 4: Save with proper layout
pdf("sft7_module_trait_heatmap_final_simple_labels.pdf", width = 8, height = 12)
draw(
  ht,
  heatmap_legend_side = "right",
  column_title = "All Species Module–Life Stage Eigengene Correlation",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  padding = unit(c(6, 6, 10, 6), "mm")  # Top, right, bottom, left
)
dev.off()

png("module_trait_heatmap_final_simple_labels.png",
    width = 3000, height = 3300, res = 300)  # 300 DPI, large canvas

draw(
  ht,
  heatmap_legend_side = "right",
  column_title = "All Species Module–Life Stage Eigengene Correlation",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  padding = unit(c(10, 10, 14, 10), "mm")  # More breathing room
)

dev.off()

# Create dataframe that associates module colors with clusters based on the new heatmap
MEcluster1 <- data.frame(moduleColor = c("steelblue"), moduleCluster = c(1))
MEcluster2 <- data.frame(moduleColor = c("turquoise","darkred","purple","lightgreen","darkgrey","skyblue","midnightblue","magenta","tan","darkgreen","royalblue"), moduleCluster = c(2))
MEcluster3 <- data.frame(moduleColor = c("black","grey60", "darkmagenta", "sienna3","yellow"), moduleCluster = c(3))
MEcluster4 <- data.frame(moduleColor = c("cyan","lightcyan","darkolivegreen","skyblue3"), moduleCluster = c(4))
MEcluster5 <- data.frame(moduleColor = c("brown", "pink", "red"), moduleCluster = c(5))
MEcluster6 <- data.frame(moduleColor = c("saddlebrown","darkturquoise","orange","white"), moduleCluster = c(6))
MEcluster7 <- data.frame(moduleColor = c("yellowgreen","blue","darkorange"), moduleCluster = c(7))
MEcluster8 <- data.frame(moduleColor = c("lightyellow", "salmon","violet","paleturquoise", "plum1"), moduleCluster = c(8))

moduleCluster = bind_rows(MEcluster1, MEcluster2, MEcluster3, MEcluster4, MEcluster5, MEcluster6, MEcluster7, MEcluster8)

chooseTopHubInEachModule(
  datExpr1_centered, 
  Colors, 
  power = 7, 
  type = "signed")

head(moduleCluster) 
head(moduleCluster) 
#"moduleColor moduleCluster
#1   steelblue             1
#2   turquoise             2
#3     darkred             2
#4      purple             2
#5  lightgreen             2
#6    darkgrey             2


# View module eigengene data
head(MEs)
#Looks OK

names(MEs)
#Looks good

Strader_MEs <- MEs
Strader_MEs$group <- treatmentinfo$group
Strader_MEs$timepoint <- treatmentinfo$timepoint
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)

head(Strader_MEs$sample_id)
#[1] "Larva1_Atenuis" "Larva2_Atenuis" "Larva3_Atenuis" "Meta1_Atenuis"  "Meta2_Atenuis"  "Meta3_Atenuis" 
head(Strader_MEs$timepoint)
#[1] I  I  I  II II II
#Levels: I II III
head(Strader_MEs$group)
#[1] Acropora_tenuis.I  Acropora_tenuis.I  Acropora_tenuis.I  Acropora_tenuis.II Acropora_tenuis.II
#[6] Acropora_tenuis.II
#12 Levels: Acropora_tenuis.I Montipora_capitata.I Pocillopora_acuta.I ... Stylophora_pistillata.III


# Cluster 1
C1_Strader_MEs <- select(Strader_MEs, MEsteelblue)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)

# Cluster 2
C2_Strader_MEs <- select(Strader_MEs, MEturquoise, MEdarkred, MEpurple, MElightgreen, 
                         MEdarkgrey, MEskyblue, MEmidnightblue, MEmagenta, MEtan, 
                         MEdarkgreen, MEroyalblue)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)

# Cluster 3
C3_Strader_MEs <- select(Strader_MEs, MEblack, MEgrey60, MEdarkmagenta, MEsienna3, MEyellow)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)

# Cluster 4
C4_Strader_MEs <- select(Strader_MEs, MEcyan, MElightcyan, MEdarkolivegreen, MEskyblue3)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)

# Cluster 5
C5_Strader_MEs <- select(Strader_MEs, MEbrown, MEpink, MEred)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)

# Cluster 6
C6_Strader_MEs <- select(Strader_MEs, MEsaddlebrown, MEdarkturquoise, MEorange, MEwhite)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)

# Cluster 7
C7_Strader_MEs <- select(Strader_MEs, MEyellowgreen, MEblue, MEdarkorange)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)

# Cluster 8
C8_Strader_MEs <- select(Strader_MEs, MElightyellow, MEsalmon, MEviolet, MEpaleturquoise, MEplum1)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)

#Create a expression profile dataframe
expressionProfile_data_group <- data.frame(
  group = Strader_MEs$group, 
  cluster1 = C1_Strader_MEs$Mean, 
  cluster2 = C2_Strader_MEs$Mean, 
  cluster3 = C3_Strader_MEs$Mean, 
  cluster4 = C4_Strader_MEs$Mean,
  cluster5 = C5_Strader_MEs$Mean,
  cluster6 = C6_Strader_MEs$Mean,
  cluster7 = C7_Strader_MEs$Mean,
  cluster8 = C8_Strader_MEs$Mean
)

# Create a lookup table
group_map <- c(
  "Acropora_tenuis.I"   = "Larvae_Atenuis",
  "Montipora_capitata.I"= "Larvae_Mcap",
  "Pocillopora_acuta.I" = "Larvae_Pacu",
  "Stylophora_pistillata.I" = "Larvae_Spis",
  
  "Acropora_tenuis.II"   = "Meta_Atenuis",
  "Montipora_capitata.II"= "Meta_Mcap",
  "Pocillopora_acuta.II" = "Meta_Pacu",
  "Stylophora_pistillata.II" = "Meta_Spis",
  
  "Acropora_tenuis.III"   = "Spat_Atenuis",
  "Montipora_capitata.III"= "Spat_Mcap",
  "Pocillopora_acuta.III" = "Spat_Pacu",
  "Stylophora_pistillata.III" = "Spat_Spis"
)

# Apply the mapping
expressionProfile_data_group$group <- group_map[as.character(expressionProfile_data_group$group)]

# Now make it a factor in the desired order
expressionProfile_data_group$group <- factor(
  expressionProfile_data_group$group,
  levels = colnames(datTraits_mod)
)

# Extract the mean columns from each cluster data frame
cluster_means <- cbind(
  C1_Strader_MEs$Mean,
  C2_Strader_MEs$Mean,
  C3_Strader_MEs$Mean,
  C4_Strader_MEs$Mean,
  C5_Strader_MEs$Mean,
  C6_Strader_MEs$Mean,
  C7_Strader_MEs$Mean,
  C8_Strader_MEs$Mean
)

expressionProfile_data_group$group <- factor(
  expressionProfile_data_group$group,
  levels = c(
    "Larvae_Atenuis", "Larvae_Mcap", "Larvae_Pacu", "Larvae_Spis",
    "Meta_Atenuis", "Meta_Mcap", "Meta_Pacu", "Meta_Spis",
    "Spat_Atenuis", "Spat_Mcap", "Spat_Pacu", "Spat_Spis"
  )
)

# Step 2: Create boxplots for each cluster
plots <- list()
for (i in 1:8) {
  cluster_col <- paste0("cluster", i)
  
  p <- expressionProfile_data_group %>%
    select(group, all_of(cluster_col)) %>%
    ggplot(aes(x = group, y = .data[[cluster_col]], fill = group)) +
    geom_boxplot(width = .5, outlier.shape = NA, alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 20, size = 5, color = "red") +
    geom_point(pch = 21, size = 3, position = position_jitter(width = 0.2)) +
    scale_fill_manual(values = rep(RColorBrewer::brewer.pal(12, "Paired"), length.out = 12)) +
    scale_x_discrete(labels = gsub("_", "\n", levels(expressionProfile_data_group$group))) +
    labs(x = "Species & Life Stage", y = "Mean Module Eigengene", title = paste("Cluster", i)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray")
  
  plots[[i]] <- p
}

# Step 3: Combine and export
combined_plot <- wrap_plots(plots, ncol = 2)

ggsave(
  filename = "expression_eigengene_profiles_species_lifestage_ordered.png",
  plot = combined_plot,
  width = 20,
  height = 24,
  units = "in",
  dpi = 300
)

