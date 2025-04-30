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
load(file = "C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/disstomsft7.RData")

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


