#WGCNA for proteomics

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
if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
  #BiocManager::install("geneLenDataBase")
  
  
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
library("goseq")


# Load necessary libraries
library(dplyr)
library(DESeq2)
library(readr)

setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Proteomics/Inputs")
setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics")


# Upload data--------------------------------------------------------------------------

treatmentinfo <- read_csv("Sample_info.csv")

# Load and rename the 'Gene' column to 'gene_id'
proteome_data <- read.csv("Protein_count_matrix_Spis.csv", header=TRUE)

# Create a unique identifier column
proteome_data <- proteome_data %>%
  mutate(Unique_ID = paste0("Spis_Protein", row_number()))

# Remove unnecessary columns for gcount function to work 
# Remove 'Gene' and 'Protein.Full.Name' columns
proteome_data <- proteome_data %>%
  select(-Gene, -Protein.Full.Name)

# Check for duplicate and missing values in 'Protein.Accession'
sum(is.na(proteome_data$Protein.Accession))  # Check for missing values
sum(duplicated(proteome_data$Protein.Accession))  # Check for duplicates

# Crop the dataset to keep only 5221 rows
proteome_data <- proteome_data %>%
  slice_head(n = 5221)

# Move 'Unique_ID' to the first column
proteome_data <- proteome_data %>%
  select(Unique_ID, everything())

write_csv(proteome_data, "Cleaned_Protein_count_matrix_Spis.csv")

# Load the cleaned data with Unique_ID as row names
gcount <- as.data.frame(read.csv("Cleaned_Protein_count_matrix_Spis.csv", row.names = "Unique_ID", header = TRUE))

# Remove columns for the WGCNA pipeline
gcount <- gcount %>%
  select(-Protein.Accession)

# Troubleshooting for integer values issue
# Check if there are any non-integer values in gcount
any(!apply(gcount, 2, function(x) all(x == floor(x))))
# [1] TRUE... Check with Jeana if next steps are acceptable 

# Preserve row names
row_names <- rownames(gcount)

# Convert all values in gcount to integers
gcount <- as.data.frame(lapply(gcount, function(x) as.integer(round(x))))

# Reapply row names
rownames(gcount) <- row_names

# Check again for non-integer values
any(!apply(gcount, 2, function(x) all(x == floor(x))))
# [1] FALSE

# Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = gcount,
                               colData = treatmentinfo,
                               design = ~timepoint)
# Log-transform the count data using a variance stabilizing transformation (vst). 
SF.gdds <- estimateSizeFactors(gdds) # estimate size factors to determine if we can use vst to transform our data. Size factors should be less than 4 to use vst
print(sizeFactors(SF.gdds)) # View size factors

gvst <- vst(gdds, blind = FALSE) # apply a variance stabilizing transformation to minimize effects of small counts and normalize with respect to library size

# Compile WGCNA dataset--------------------------------------------------------------------------
# Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.
datExpr <- as.data.frame(t(assay(gvst))) # transpose to output to a new data frame with the column names as row names and make all data numeric

# Set the row names of datExpr to the Unique_IDs from gcount
rownames(datExpr) <- gcount$Unique_ID


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
pdf(paste0('sampleTree','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#PCA-----------------------------------------------------------------------------

gPCAdata <- plotPCA(gvst, intgroup = c("timepoint"), returnData=TRUE, ntop=1000) #use ntop to specify all genes
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
#Some troubleshooting steps so all samples are included
range(gPCAdata$PC1)
range(gPCAdata$PC2)
ylim_range <- range(gPCAdata$PC2, na.rm = TRUE)

ylim_range <- range(gPCAdata$PC2, na.rm = TRUE)

# Create the PCA plot with updated limits if needed
allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(PC1, PC2)) + 
  geom_point(aes(shape = timepoint, colour = timepoint), size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(ylim_range) +  # Use updated y-axis limits
  coord_fixed() +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank()
  )

# Print the PCA plot
print(allgenesfilt_PCA_visual)


# Explicitly set the levels of the timepoint variable
gPCAdata$timepoint <- factor(gPCAdata$timepoint, levels = c("I", "II", "III"))

allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(PC1, PC2)) + 
  geom_point(aes(shape = timepoint, colour = timepoint), size = 6) +  # Increase point size
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ylim(ylim_range) +  # Use updated y-axis limits
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

ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/PCA_timepointntop=1000.png", allgenesfilt_PCA_visual, width = 11, height = 8)

# Load necessary libraries
library(vegan)
library(ggplot2)


# Assuming gvst is your DESeqTransform object, extract the assay (expression) data
gvst_matrix <- assay(gvst)

# Check the structure of the extracted matrix (optional)
str(gvst_matrix)

# Transpose the matrix and calculate the Bray-Curtis dissimilarity
# Transposing so that samples are in rows and genes/features are in columns
dissimilarity_matrix <- vegdist(t(gvst_matrix), method = "bray")

# Optionally, inspect the dissimilarity matrix
print(dissimilarity_matrix)

# If you plan to use the dissimilarity matrix for downstream analysis (e.g., NMDS or clustering), you can proceed accordingly.
# For example, you could run NMDS with this matrix:
nmds <- metaMDS(dissimilarity_matrix)

# Visualize the NMDS (optional)
plot(nmds, type = "t")


#Choose a set of soft-thresholding powers
powers <- c(seq(from = 1, to=19, by=2), c(21:30)) #Create a string of numbers from 1 through 10, and even numbers from 10 through 20
powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2))

#Call the network topology analysis function
library(doParallel)
sft <-pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(paste0('network', '.pdf'))
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

# The lowest scale-free topology fit index R^2 recommended by Langfelder and Horvath is 0.8. 
# From the graph, it appears that our soft thresholding power is 19 because it is the lowest 
# power before the R^2=0.8 threshold that maximizes with model fit (s17 is the number right above the red line).


### Network construction and module detection:
# Co-expression adjacency and topological overlap matrix similarity
# Co-expression similarity and adjacency, using the soft thresholding power 19 and translate the adjacency into topological overlap matrix to calculate 
# the corresponding dissimilarity. I will use a signed network because we have a relatively high softPower, according 
# to >12 (https://peterlangfelder.com/2018/11/25/__trashed/). 
# Moreover, in expression data where you are interested in when expression on one gene increases or decreases with expression level of another you would use a signed network (when you are interested in the direction of change, correlation and anti-correlation, you use a signed network).

options(stringsAsFactors = FALSE)
enableWGCNAThreads() #Allow multi-threading within WGCNA


#Version two of this graph (with dots)

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

#Run analysis. Tengo que correr esto en un compu más fancy. 

softPower=19 #Set softPower to 19
adjacency=adjacency(datExpr, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
save(TOM, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/TOMspisprot.RData")
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(dissTOM, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/dissTOMspisprot.RData")
load("~/Multistage_omics/R scripts/S. pistillata/Proteomics/dissTOMspisprot.RData")

# Clustering using TOM
#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/dissTOMClusteringProt.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()
dev.off()

# Module identification is cutting the branches off the tree in the dendrogram above. We want large modules, so we set the minimum module size 
# relatively high (minimum size = 30).

minModuleSize = 30 #default value used most often
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods) #list modules and respective sizes
save(dynamicMods, geneTree, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/dyMod_geneTree.RData")

dyMod_geneTree <- load(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/dyMod_geneTree.RData")
table(dynamicMods) #list modules and respective sizes
dynamicMods
#1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
#105  98  97  95  91  88  87  84  79  76  74  74  72  71  70  69  69  69  69  69  67  65  64  63  61 
#26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50 
#61  61  60  59  56  56  56  55  55  55  55  54  54  54  53  52  52  52  52  52  51  50  50  49  49 
#51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75 
#49  48  48  48  47  47  46  46  46  45  44  44  44  44  43  43  42  42  41  41  41  41  41  40  39 
#76  77  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 
#39  39  39  39  38  38  37  37  37  37  37  36  36  35  34  34  33  32  32  32  31  30  30  30 

dyMod_geneTree
#[1] "dynamicMods" "geneTree"   

# Plot the module assignment under the gene dendrogram
dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
#antiquewhite2   antiquewhite4         bisque4           black            blue           blue2 
#34              44              50              87              98              39 
#brown          brown2          brown4           coral          coral1          coral2 
#97              39              51              35              45              44 
#coral3            cyan       darkgreen        darkgrey     darkmagenta  darkolivegreen 
#34              71              65              63              55              55 
#darkolivegreen4      darkorange     darkorange2         darkred   darkseagreen3   darkseagreen4 
#40              61              52              67              36              46 
#darkslateblue   darkturquoise      darkviolet      firebrick4     floralwhite           green 
#50              64              39              41              52              91 
#greenyellow          grey60        honeydew       honeydew1      indianred4           ivory 
#74              69              36              46              41              52 
#lavenderblush2  lavenderblush3      lightcoral       lightcyan      lightcyan1      lightgreen 
#37              46              41              69              52              69 
#lightpink3      lightpink4  lightslateblue  lightsteelblue lightsteelblue1     lightyellow 
#37              47              30              41              52              69 
#magenta        magenta4          maroon    mediumorchid   mediumpurple1   mediumpurple2 
#79              37              47              44              30              41 
#mediumpurple3   mediumpurple4    midnightblue    navajowhite1    navajowhite2          orange 
#53              33              70              37              48              61 
#orangered1      orangered3      orangered4   paleturquoise  palevioletred2  palevioletred3 
#30              42              54              56              37              48 
#pink           pink4            plum           plum1           plum2           plum3 
#84              31              42              54              49              39 
#purple             red       royalblue     saddlebrown          salmon         salmon2 
#76              88              69              59              72              38 
#salmon4         sienna3         sienna4         skyblue        skyblue1        skyblue2 
#48              55              32              60              43              44 
#skyblue3        skyblue4       steelblue             tan         thistle        thistle1 
#54              32              56              74              38              49 
#thistle2        thistle3       turquoise          violet           white          yellow 
#49              39             105              56              61              95 
#yellow3         yellow4     yellowgreen 
#32              43              55 


pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/dissTOMColorClustering.pdf", width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar or choose not to merge
# Plot module similarity based on eigengene value

#Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower = 19)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster again and plot the results
METree = flashClust(as.dist(MEDiss), method = "average")

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/eigengeneClustering1.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >85% eigengene similarity (most studies use80-90% similarity)

MEDissThres= 0.15 #merge modules that are 85% similar

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/eigengeneClusteringV2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
#mergeCloseModules: Merging modules whose distance is less than 0.15
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 99 module eigengenes in given set.
#Calculating new MEs...
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 99 module eigengenes in given set.

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/mergedClusters.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

#Save new colors

moduleColors = mergedColors # Rename to moduleColors
colorOrder = c("grey", standardColors(50)); # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
ncol(MEs) 
#99

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/eigengeneClustering3.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# Relating modules to group, quantifying module–trait associations
#Prepare trait data. Data has to be numeric, so I replaced the group for numeric values
#Check order with Treatment info

allTraits <- names(treatmentinfo$timepoint)
allTraits$larvae_released <- c(1,1,1,1,1,0,0,0,0,0)
allTraits$spat <- c(0,0,0,0,0,1,1,1,1,1)

datTraits <- as.data.frame(allTraits)
dim(datTraits)
#10 2

rownames(datTraits) <- treatmentinfo$sampleID
print(datTraits)
#larvae_released spat
#Larva_4               1    0
#Larva_5               1    0
#Larva_6               1    0
#Larva_8               1    0
#Larva_9               1    0
#Spat_2                0    1
#Spat_3                0    1
#Spat_4                0    1
#Spat_5                0    1
#Spat_6                0    1

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors, softPower = 19)$eigengenes
MEs = orderMEs(MEs0)

# Calculate the correlation matrix of module eigengenes with themselves
moduleCor = cor(MEs, use = "p")

# Convert the correlation matrix to a dissimilarity matrix
dissimilarity = 1 - abs(moduleCor)

# Ensure the dissimilarity matrix is square
print(dim(dissimilarity))

# Perform hierarchical clustering on the dissimilarity matrix
moduleTraitTree = hclust(as.dist(dissimilarity), method = "average")

# Check the class of moduleTraitTree to ensure it's "hclust"
print(class(moduleTraitTree))

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/Life-stage clustering based on module-trait correlation2.pdf", width = 12, height = 8)
plot(moduleTraitTree, 
     main = "Group clustering based on module-trait correlation", 
     sub = "", 
     xlab = "", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)
#[1,] "-0.28\n(0.4)"  "0.28\n(0.4)" 
#[2,] "-0.13\n(0.7)"  "0.13\n(0.7)" 
#[3,] "0.4\n(0.3)"    "-0.4\n(0.3)" 
#[4,] "-0.063\n(0.9)" "0.063\n(0.9)"
#[5,] "-0.34\n(0.3)"  "0.34\n(0.3)" 
#[6,] "0.18\n(0.6)"   "-0.18\n(0.6)"

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. Pistillata/Proteomics/Module-trait-relationships2.pdf")
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), 
               cex.lab.y= 0.55, cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , 
               zlim = c(-1,1), main = paste("Module-trait relationships"))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), cex.lab.y= 0.55, 
               cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , zlim = c(-1,1), 
               main = paste("Module-trait relationships"))
dev.off()

#In two pages

# Define the number of rows for each page
rowsPerPage = 50

# Calculate the number of pages required
numPages = ceiling(nrow(moduleTraitCor) / rowsPerPage)

# Open the PDF file
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. Pistillata/Proteomics/Module-trait-relationships2pages.pdf", width = 12, height = 8)

for (page in 1:numPages) {
  # Determine the row indices for the current page
  startRow = (page - 1) * rowsPerPage + 1
  endRow = min(page * rowsPerPage, nrow(moduleTraitCor))
  
  # Extract the subset of the data for the current page
  subsetTraitCor = moduleTraitCor[startRow:endRow, ]
  subsetTextMatrix = textMatrix[startRow:endRow, ]
  
  # Plot the heatmap for the current page
  labeledHeatmap(
    Matrix = subsetTraitCor,
    xLabels = names(datTraits),
    yLabels = rownames(subsetTraitCor),
    ySymbols = rownames(subsetTraitCor),
    cex.lab.y = 0.55,
    cex.lab.x = 0.55,
    colors = blueWhiteRed(50),
    textMatrix = subsetTextMatrix,
    setStdMargins = TRUE,
    cex.text = 0.25,
    textAdj = c(0.5, 0.5),  # Adjust text alignment if needed
    zlim = c(-1, 1),
    main = paste("Module-trait relationships (Page", page, "of", numPages, ")")
  )
}



#Create list of p-values for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)

# Define color scheme for heatmap
htmap.colors <- gsub("ME", "", names(MEs))

# Ensure dendrograms are correctly sorted
row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

# Create heatmap
pdf(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. Pistillata/Proteomics/Module-trait-relationship-heatmap3.pdf", height = 11.5, width = 8)
ht = Heatmap(
  moduleTraitCor,
  name = "Eigengene",
  column_title = "Module-Trait Eigengene Correlation",
  col = blueWhiteRed(50),
  row_names_side = "left",
  row_dend_side = "left",
  width = unit(4, "in"),
  height = unit(8.5, "in"),
  column_dend_reorder = FALSE,
  cluster_columns = col_dend,
  column_split = 3,  # Adjust as needed
  column_dend_height = unit(0.5, "in"),
  cluster_rows = row_dend,
  row_split = 10,  # Adjust as needed
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
  row_names_gp = gpar(fontsize = 5, alpha = 0.75, border = TRUE, fill = htmap.colors)
)
draw(ht)
dev.off()

# Create dataframe that associates module colors with clusters based on the above heatmap and the MEtree dendrogram.I did this manually looking at the clusters in the previous plot.
MEcluster1 <- data.frame(moduleColor = c("darkorange","palevioletred3","ivory","coral1","lightpink4","yellowgreen","bisque4","grey60","magenta4","darkorange2","blue2","antiquewhite","darkseagreen4"), moduleCluster = c(1))
MEcluster2 <- data.frame(moduleColor = c("white","cyan","firebrick4","palevioletred2","lightslateblue","coral","coral3","darkviolet","honeydew1","lavenderblush3"), moduleCluster = c(2))
MEcluster3 <- data.frame(moduleColor = c("darkolivegreen4","darkslateblue","skyblue3","magenta","coral2","purple","sienna3","orangered4","salmon","plum2","lightcyan1","plum1","orange","antiquewhite2","paleturquoise","thistle","violet","brown2","plum"), moduleCluster = c(3))
MEcluster4 <- data.frame(moduleColor = c("thistle2","red","salmon4","plum3","navajowhite2","darkred","tan"), moduleCluster = c(4))
MEcluster5 <- data.frame(moduleColor = c("yellow","mediumpurple4","darkmagenta","lightsteelblue1","lightpink3","lightcoral","pink","thistle3","lightgreen","indianred4","skyblue","thistle1","yellow4","mediumorchid","black"), moduleCluster = c(5))
MEcluster6 <- data.frame(moduleColor = c("greenyellow","maroon","brown4","turquoise"), moduleCluster = c(6)) 
MEcluster7 <- data.frame(moduleColor = c("salmon2","skyblue1","mediumpurple1","orangered1","mediumpurple2","brown","floralwhite","yellow3","skyblue2","blue","darkseagreen3","lavenderblush2","darkolivegreen","darkgrey","pink4","honeydew","lightyellow","lightcyan","mediumpurple3","saddlebrown","skyblue4"), moduleCluster = c(7)) 
MEcluster8 <- data.frame(moduleColor = c("orangered3","darkgreen","navajowhite1","royalblue"), moduleCluster = c(8)) 
MEcluster9 <- data.frame(moduleColor = c("darkturquoise","steelblue","green","midnightblue","sienna4"), moduleCluster = c(9))
MEcluster10 <- data.frame(moduleColor = c("lightsteelblue"), moduleCluster = c(10)) 

moduleCluster = bind_rows(MEcluster1, MEcluster2, MEcluster3, MEcluster4, MEcluster5, MEcluster6, MEcluster7, MEcluster8, MEcluster9, MEcluster10)
head(moduleCluster) 
# moduleColor moduleCluster
#1     darkorange             1
#2 palevioletred3             1
#3          ivory             1
#4         coral1             1
#5     lightpink4             1
#6    yellowgreen             1
head(MEs)
#MEsaddlebrown MEantiquewhite2   MEyellow  MEsienna3 MEdarkolivegreen MElightslateblue    MEplum2    MEorange
#Larva_4   -0.57668924      -0.2851064 -0.3410496  0.3387905      -0.51521862      -0.30110563 -0.3378313 -0.37857721
#Larva_5    0.13396779       0.2026294  0.2473703  0.3456226       0.33226954       0.36527654  0.4407567  0.02224021
#Larva_6    0.04690997      -0.2812836  0.4299356 -0.1390453       0.19166981       0.23921869  0.1657482  0.26163203
#Larva_8    0.26060292       0.4895164  0.4729596 -0.1167882      -0.04743883       0.09663953 -0.2543979  0.37797339
#Larva_9   -0.31461528      -0.3376395 -0.1832719 -0.5288867      -0.49192266      -0.12085670 -0.1387295 -0.45586552
#Spat_2     0.17812094       0.4009117  0.1761350  0.1854952       0.34037030       0.37615062  0.1520776  0.34135156

names(MEs)
#[1] "MEsaddlebrown"     "MEantiquewhite2"   "MEyellow"          "MEsienna3"         "MEdarkolivegreen"  "MElightslateblue" 
#[7] "MEplum2"           "MEorange"          "MEcoral2"          "MEfirebrick4"      "MEdarkred"         "MEturquoise"      
#[13] "MEgreenyellow"     "MEtan"             "MEplum3"           "MEdarkmagenta"     "MEhoneydew1"       "MEantiquewhite4"  
#[19] "MEindianred4"      "MEnavajowhite2"    "MElavenderblush3"  "MElightpink3"      "MEbrown2"          "MEcoral3"         
#[25] "MElightcyan"       "MEbrown"           "MElightcyan1"      "MEdarkorange2"     "MElightsteelblue1" "MEmagenta"        
#[31] "MEhoneydew"        "MEskyblue2"        "MEmagenta4"        "MEpink4"           "MEdarkolivegreen4" "MEdarkslateblue"  
#[37] "MEdarkseagreen4"   "MEthistle"         "MEplum1"           "MEpurple"          "MEdarkseagreen3"   "MEplum"           
#[43] "MEskyblue4"        "MEviolet"          "MEthistle1"        "MEblue2"           "MEthistle3"        "MEmediumorchid"   
#[49] "MEthistle2"        "MEcoral"           "MEpink"            "MElightyellow"     "MEblue"            "MEdarkorange"     
#[55] "MEdarkgrey"        "MEsalmon"          "MEmediumpurple1"   "MEorangered3"      "MElavenderblush2"  "MEsalmon2"        
#[61] "MEskyblue1"        "MEroyalblue"       "MEfloralwhite"     "MEyellow3"         "MEmediumpurple3"   "MEmidnightblue"   
#[67] "MEmediumpurple2"   "MEnavajowhite1"    "MEorangered1"      "MEsienna4"         "MEgreen"           "MElightsteelblue" 
#[73] "MEdarkturquoise"   "MEsteelblue"       "MEdarkgreen"       "MElightpink4"      "MEmediumpurple4"   "MEyellowgreen"    
#[79] "MEwhite"           "MElightgreen"      "MEpalevioletred3"  "MEcyan"            "MEpaleturquoise"   "MEred"            
#[85] "MEskyblue"         "MEskyblue3"        "MEblack"           "MEivory"           "MEyellow4"         "MEcoral1"         
#[91] "MEorangered4"      "MEbrown4"          "MEmaroon"          "MElightcoral"      "MEsalmon4"         "MEbisque4"        
#[97] "MEpalevioletred2"  "MEdarkviolet"      "MEgrey60"         

Strader_MEs <- MEs
Strader_MEs$timepoint <- treatmentinfo$timepoint
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)
#        MEsalmon4   MEbisque4 MEpalevioletred2 MEdarkviolet    MEgrey60 timepoint sample_id
#Larva_4  0.3308714 -0.38271770       -0.3408834   -0.2407381 -0.35434518         I   Larva_4
#Larva_5  0.1758668  0.05012932        0.3402429   -0.1451983 -0.06530827         I   Larva_5
#Larva_6 -0.3418352 -0.31700728       -0.3922246   -0.2119718 -0.26032574         I   Larva_6
#Larva_8  0.3063943  0.34760685        0.4075936    0.4793531  0.35767232         I   Larva_8
#Larva_9  0.3464471  0.34218641        0.3251685    0.5288320  0.36100128         I   Larva_9
#Spat_2  -0.3826482 -0.28193127       -0.1340736    0.3941905  0.02022685       III    Spat_2
#note: it's longer
head(Strader_MEs$sample_id)
#[1] "Larva_4" "Larva_5" "Larva_6" "Larva_8" "Larva_9" "Spat_2" 
head(Strader_MEs$timepoint)
#[1] I   I   I   I   I   III
#Levels: I III

# Calculate 10 over-arching expression patterns using mean eigengene for each module in a cluster. Use MEcluster as models for first lines
# Create a column to the Strader_MEs data frame containing lifestage groups


MEcluster1 <- data.frame(moduleColor = c("darkorange","palevioletred3","ivory","coral1","lightpink4","yellowgreen","bisque4","grey60","magenta4","darkorange2","blue2","antiquewhite","darkseagreen4"), moduleCluster = c(1))
MEcluster2 <- data.frame(moduleColor = c("white","cyan","firebrick4","palevioletred2","lightslateblue","coral","coral3","darkviolet","honeydew1","lavenderblush3"), moduleCluster = c(2))
MEcluster3 <- data.frame(moduleColor = c("darkolivegreen4","darkslateblue","skyblue3","magenta","coral2","purple","sienna3","orangered4","salmon","plum2","lightcyan1","plum1","orange","antiquewhite2","paleturquoise","thistle","violet","brown2","plum"), moduleCluster = c(3))
MEcluster4 <- data.frame(moduleColor = c("thistle2","red","salmon4","plum3","navajowhite2","darkred","tan"), moduleCluster = c(4))
MEcluster5 <- data.frame(moduleColor = c("yellow","mediumpurple4","darkmagenta","lightsteelblue1","lightpink3","lightcoral","pink","thistle3","lightgreen","indianred4","skyblue","thistle1","yellow4","mediumorchid","black"), moduleCluster = c(5))
MEcluster6 <- data.frame(moduleColor = c("greenyellow","maroon","brown4","turquoise"), moduleCluster = c(6)) 
MEcluster7 <- data.frame(moduleColor = c("salmon2","skyblue1","mediumpurple1","orangered1","mediumpurple2","brown","floralwhite","yellow3","skyblue2","blue","darkseagreen3","lavenderblush2","darkolivegreen","darkgrey","pink4","honeydew","lightyellow","lightcyan","mediumpurple3","saddlebrown","skyblue4"), moduleCluster = c(7)) 
MEcluster8 <- data.frame(moduleColor = c("orangered3","darkgreen","navajowhite1","royalblue"), moduleCluster = c(8)) 
MEcluster9 <- data.frame(moduleColor = c("darkturquoise","steelblue","green","midnightblue","sienna4"), moduleCluster = c(9))
MEcluster10 <- data.frame(moduleColor = c("lightsteelblue"), moduleCluster = c(10)) 


C1_Strader_MEs <- select(Strader_MEs, MEdarkorange:MEdarkseagreen4)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)
C2_Strader_MEs <- select(Strader_MEs, MEwhite:MElavenderblush3)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)
C3_Strader_MEs <- select(Strader_MEs, MEdarkolivegreen4:MEplum)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)
C4_Strader_MEs <- select(Strader_MEs, MEthistle2:MEtan)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)
C5_Strader_MEs <- select(Strader_MEs, MEyellow:MEblack)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)
C6_Strader_MEs <- select(Strader_MEs, MEgreenyellow:MEturquoise)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)
C7_Strader_MEs <- select(Strader_MEs, MEsalmon2:MEskyblue4)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)
C8_Strader_MEs <- select(Strader_MEs, MEorangered3:MEroyalblue)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)
C9_Strader_MEs <- select(Strader_MEs, MEdarkturquoise:MEsienna4)
C9_Strader_MEs$Mean <- rowMeans(C9_Strader_MEs)
C10_Strader_MEs <- select(Strader_MEs, MElightsteelblue)
C10_Strader_MEs$Mean <- rowMeans(C10_Strader_MEs)

expressionProfile_data <- as.data.frame(cbind(group = Strader_MEs$timepoint, cluster1= C1_Strader_MEs$Mean, cluster2 = C2_Strader_MEs$Mean, 
                                              cluster3 = C3_Strader_MEs$Mean, cluster4 = C4_Strader_MEs$Mean,cluster5 = C5_Strader_MEs$Mean,cluster6 = C6_Strader_MEs$Mean,
                                              cluster7 = C7_Strader_MEs$Mean,cluster8 = C8_Strader_MEs$Mean,cluster9 = C9_Strader_MEs$Mean,cluster10 = C10_Strader_MEs$Mean))

head(expressionProfile_data)
head(expressionProfile_data)
#group    cluster1      cluster2    cluster3    cluster4     cluster5
#1     1  0.10320896 -5.048773e-05 -0.12323730 -0.02507255  0.007017195
#2     1 -0.01383511 -1.154881e-01 -0.21077120 -0.07585809 -0.022924872
#3     1  0.19623224  9.740973e-03  0.29907236  0.31178897  0.043720273
#4     1 -0.20113447 -1.500497e-01 -0.02446714 -0.08480652 -0.048995771
#5     1 -0.04922112 -2.390850e-02 -0.10546081 -0.01227533 -0.044635575
#6     2 -0.19230672 -1.413776e-03 -0.18855158 -0.05282130 -0.014666968
#cluster6     cluster7     cluster8   cluster9  cluster10
#1  0.33992461  0.141714813  0.003559169 -0.2232839 -0.2589222
#2  0.14428724  0.010618508 -0.104755289 -0.3018798 -0.2950282
#3  0.41461859 -0.009817634 -0.308445626 -0.3337249 -0.2765966
#4  0.06611453 -0.182262383 -0.162729829 -0.1250049 -0.3607362
#5  0.34357056 -0.044409215 -0.151645459 -0.3259148 -0.3286812
#6 -0.18938191 -0.218022573 -0.155001587  0.3770063  0.2995247

# Assuming you have a list or data frame 'expressionProfile_data' with cluster data
# And your cluster-specific data frames are named like C1_Strader_MEs, C2_Strader_MEs, etc.This is a loop, run it together. 
#These ones came out without a group, check later if that's needed
# Assuming you have a list or data frame 'expressionProfile_data' with cluster data
# And your cluster-specific data frames are named like C1_Strader_MEs, C2_Strader_MEs, etc.This is a loop, run it together. 
#These ones came out without a group, check later if that's needed
for (i in 1:10) {
  cluster_name = paste0("cluster", i)
  strader_ME_name = paste0("C", i, "_Strader_MEs")
  
  # Extract mean eigengene values for the current cluster
  meanEigenClust <- expressionProfile_data[[cluster_name]]
  
  # Save the mean eigengene values to a CSV file
  meanEigenFilename = paste0("meanEigenClust", i, ".csv")
  write.csv(meanEigenClust, file = meanEigenFilename)
  
  # Save the Strader ME data for the current cluster to a CSV file
  straderMEFilename = paste0("C", i, "_Strader_MEs.csv")
  write.csv(get(strader_ME_name), file = straderMEFilename)
}
 
expressionProfile_data$group
#> expressionProfile_data$group
#[1] 1 1 1 1 1 2 2 2 2 2


#cols.num <- c(2:11)
#expressionProfile_data[cols.num] <- sapply(expressionProfile_data[cols.num],as.numeric)

expressionProfile_data$group <- factor(expressionProfile_data$group)
sapply(expressionProfile_data, class)
#group  cluster1  cluster2  cluster3  cluster4  cluster5  cluster6  cluster7 
#"factor" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" 
#cluster8  cluster9 cluster10 
#"numeric" "numeric" "numeric"

dim(expressionProfile_data)
#[1] 10 11

head(expressionProfile_data)
#group    cluster1      cluster2    cluster3    cluster4     cluster5
#1     1  0.10320896 -5.048773e-05 -0.12323730 -0.02507255  0.007017195
#2     1 -0.01383511 -1.154881e-01 -0.21077120 -0.07585809 -0.022924872
#3     1  0.19623224  9.740973e-03  0.29907236  0.31178897  0.043720273
#4     1 -0.20113447 -1.500497e-01 -0.02446714 -0.08480652 -0.048995771
#5     1 -0.04922112 -2.390850e-02 -0.10546081 -0.01227533 -0.044635575
#6     2 -0.19230672 -1.413776e-03 -0.18855158 -0.05282130 -0.014666968
#cluster6     cluster7     cluster8   cluster9  cluster10
#1  0.33992461  0.141714813  0.003559169 -0.2232839 -0.2589222
#2  0.14428724  0.010618508 -0.104755289 -0.3018798 -0.2950282
#3  0.41461859 -0.009817634 -0.308445626 -0.3337249 -0.2765966
#4  0.06611453 -0.182262383 -0.162729829 -0.1250049 -0.3607362
#5  0.34357056 -0.044409215 -0.151645459 -0.3259148 -0.3286812
#6 -0.18938191 -0.218022573 -0.155001587  0.3770063  0.2995247

library(ggplot2)

group_order = c("1","2","3")

Cluster1Plot <- expressionProfile_data %>%
  select(group, cluster1) %>% 
  #group_by(group) %>% 
  ggplot(aes(x=group, y=cluster1, fill=group)) +
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5), alpha = 0.7) +
  #stat_summary(fun=mean, geom="line", aes(group=group, color = group), position = position_dodge(width = 0.5))  + 
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch = 21, size=5, position = position_dodge(width = 1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + #Axis titles
  #ylim(-0.5,1) +
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster1") +
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, color = "black"))+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") 

Cluster2Plot <- expressionProfile_data %>%
  select(group, cluster2) %>% 
  #group_by(group) %>% 
  ggplot(aes(x=group, y=cluster2, fill=group)) +
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5), alpha = 0.7) +
  #stat_summary(fun=mean, geom="line", aes(group=group, color = group), position = position_dodge(width = 0.5))  + 
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch = 21, size=5, position = position_dodge(width = 1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + #Axis titles
  #ylim(-0.5,1) +
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster2") +
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, color = "black"))+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") 
plot(Cluster2Plot)

Cluster3Plot <- expressionProfile_data %>%
  select(group, cluster3) %>% 
  ggplot(aes(x=group, y=cluster3, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigengene") +
  ggtitle("Cluster3") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank(),
        axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster4Plot <- expressionProfile_data %>%
  select(group, cluster4) %>% 
  ggplot(aes(x=group, y=cluster4, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigengene") +
  ggtitle("Cluster4") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank(),
        axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")


Cluster5Plot <- expressionProfile_data %>%
  select(group, cluster5) %>% 
  ggplot(aes(x=group, y=cluster5, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigengene") +
  ggtitle("Cluster5") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank(),
        axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

plot(Cluster5Plot)

Cluster6Plot <- expressionProfile_data %>%
  select(group, cluster6) %>% 
  ggplot(aes(x=group, y=cluster6, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigengene") +
  ggtitle("Cluster6") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank(),
        axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster7Plot <- expressionProfile_data %>%
  select(group, cluster7) %>% 
  ggplot(aes(x=group, y=cluster7, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigengene") +
  ggtitle("Cluster7") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank(),
        axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")
Cluster8Plot <- expressionProfile_data %>%
  select(group, cluster8) %>% 
  ggplot(aes(x=group, y=cluster8, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigengene") +
  ggtitle("Cluster8") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank(),
        axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster9Plot <- expressionProfile_data %>%
  select(group, cluster9) %>% 
  ggplot(aes(x=group, y=cluster9, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigengene") +
  ggtitle("Cluster9") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank(),
        axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster10Plot <- expressionProfile_data %>%
  select(group, cluster10) %>% 
  ggplot(aes(x=group, y=cluster10, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigengene") +
  ggtitle("Cluster10") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank(),
        axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")


expressionProfiles <- (Cluster1Plot + Cluster2Plot + Cluster3Plot) / (Cluster4Plot + Cluster5Plot + Cluster6Plot) / (Cluster7Plot + Cluster8Plot + Cluster9Plot)/ (Cluster10Plot+Cluster10Plot+Cluster10Plot)
ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Proteomics/expression_eigengene_Profiles_withMean.pdf", expressionProfiles, height = 25, width = 28, units = "in")

#  Gene relationship to trait and important modules: Gene Significance and Module Membership

#We quantify associations of individual genes with life stage by defining Gene Significance GS as the absolute value of the correlation between the gene and the time_point. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 

#Define variable weight containing the weight column of datTrait

treatmentinfo <- filter(treatmentinfo, sampleID %in% rownames(datExpr))
datExpr <- datExpr[treatmentinfo$sampleID,]
group <- data.frame(group = as.numeric(as.factor(treatmentinfo$timepoint)))
MEs <- MEs[treatmentinfo$sampleID,]

group_spis_proteome <- datTraits
dim(group_spis_proteome)
#[1] 10  2
# Colors of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Use group_p_acuta for geneTraitSignificance
geneTraitSignificance_spis_proteome = as.data.frame(cor(datExpr, group_spis_proteome, use = "p"))
GSPvalue_spis_proteome = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_spis_proteome), nSamples))

names(geneTraitSignificance_spis_proteome) = paste("GS.", names(group_spis_proteome), sep="")
names(GSPvalue_spis_proteome) = paste("p.GS.", names(group_spis_proteome), sep="")

#Use group for geneTraitSignificance for groups
geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(group), sep="")
names(GSPvalue) = paste("p.GS.", names(group), sep="")

### Summary output of network analysis results
library(readr)
#Annotation file, needs a few modifications

#Annot <- read_csv("Annot_Spis_Proteins.csv")

#Add the Unique_ID
#Annot<- Annot %>%
#  mutate(Unique_ID = paste0("Spis_Protein", row_number()))
#Crop the dataset to keep only 5221 rows, because a lot of rows were coming empty
#Annot <- Annot %>%
#  slice_head(n = 5221)
#  Move 'Unique_ID' to the first column
#Annot <- Annot %>%
#  select(Unique_ID, everything())

#Write CSV and re-load
#write_csv(Annot, "Annot_Spis_Proteins_mod.csv")

#Load Annotation file
annot <- read_csv("Annot_Spis_Proteins_mod.csv")

probes = names(datExpr)
probes2annot = match(probes, annot$Unique_ID)
nrow(annot)#5221
head(annot)
#A tibble: 6 x 14
#Unique_ID     Protein.Accession     Gene  Protein.Full.Name Entry
#<chr>         <chr>                 <chr> <chr>             <chr>
#  1 Spis_Protein1 A0A059U7T2            NA    Beta-tubulin (Fr~ A0A0~
#                                                              2 Spis_Protein2 A0A059UGC1|A0A059UGU~ ||CY~ Cytochrome b (Fr~ A0A0~
#                                                                                                                          3 Spis_Protein3 A0A0G2SJD9            SLC4~ Anion exchange p~ A0A0~
#                                                                                                                          4 Spis_Protein4 A0A0G2SJL4|A0A2B4SEX6 SLC4~ Anion exchange p~ A0A0~
#                                                                                                                          5 Spis_Protein5 A0A0G2SJL5            SLC2~ Solute carrier f~ A0A0~
#             6 Spis_Protein6 A0A1B0Y2D9            CA11  Alpha carbonic a~ A0A1~
sum(is.na(probes2annot))#0, this number must be 0.

# Match probes to annot Unique_ID
probes2annot <- match(probes, annot$Unique_ID)



# Create the starting data frame
protInfo0 = data.frame(
  Unique_ID = annot$Unique_ID,
  Protein.Accession = annot$Protein.Accession[probes2annot],
  Protein.Full.Name = annot$Protein.Full.Name[probes2annot],
  Gene.Ontology.IDs = annot$Gene.ontology.IDs[probes2annot],
  Gene.ontology.GO = annot$Gene.ontology.GO[probes2annot],
  Gene = annot$Gene[probes2annot],
  Entry = annot$Entry[probes2annot],
  Entry.name = annot$Entry.name[probes2annot],
  Protein.names =annot$Protein.names[probes2annot],
  moduleColor = dynamicColors,
  geneTraitSignificance,
  GSPvalue,
  geneTraitSignificance_spis_proteome,
  GSPvalue_spis_proteome
)

# Order modules by their significance for time_point
modOrder = order(-abs(cor(MEs, group, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(protInfo0)
  protInfo0 = data.frame(protInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  names(protInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(protInfo0$moduleColor, -abs(protInfo0$GS.group))
protInfo = protInfo0[geneOrder, ]

#Add module cluster information
protInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")
dim(protInfo)
#[1] 5221  214

head(protInfo)
#Looks good

write.csv(protInfo, file = "ProtInfo.csv")

