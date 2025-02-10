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
pdf(paste0('sampleTree_ComBat','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()


#PCA--------------------------------------------------------------------------------------------------------------
gPCAdata <- plotPCA(gvst, intgroup = c("timepoint", "Species"), returnData=TRUE, ntop=500)

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
  scale_shape_manual(values = c(15, 16, 17, 18)) # Different shapes for each species

# Print the PCA plot
print(allgenesfilt_PCA_visual)


#Choose a set of soft-thresholding powers
powers <- c(seq(from = 1, to=19, by=1), c(21:30)) #Create a string of numbers from 1 through 10, and even numbers from 10 through 20
powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 1))

# Load the doParallel package
library(doParallel)

# Register a parallel backend with the desired number of cores
# Adjust the number of cores as needed
registerDoParallel(cores = 4)  # Set the number of cores
#Call the network topology analysis function
sft <-pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Don't forget to stop the parallel backend when you're done
#stopImplicitCluster()

# This will unregister the parallel backend and prevent further messages

#Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9

#Scale-free topology fit index as a function of the soft-thresholding power
pdf(paste0('network+combat', '.pdf'))
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


softPower=10 #Set softPower to 17# let's start with 10
adjacency=adjacency(datExpr, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(adjacency, file = "C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/TOM_10sft.RData")
save(TOM, file ="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/TOM_10sft.RData")
save(dissTOM, file = "C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/disstomsft10.RData") 
load(file = "C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/disstomsft10.RData")

#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/Multisp_dissTOMClustering.pdf", width=20, height=20)
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
#1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
#1020  476  439  331  293  267  245  234  230  205  189  170  151  150  144  143  130  123  120  120  119  114  105   95   92   92   91   84   78 
#30   31   32   33   34   35   36   37   38   39   40   41   42   43   44   45   46   47   48   49 
#77   67   65   63   62   60   58   57   57   57   56   55   55   49   46   44   41   35   32   32

dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
#table(dynamicColors)
dynamicColors
#bisque4           black            blue           brown          brown4            cyan       darkgreen        darkgrey     darkmagenta 
#35             245             476             439              41             150             114              95              62 
#darkolivegreen      darkorange     darkorange2         darkred   darkslateblue   darkturquoise     floralwhite           green     greenyellow 
#63              92              44             119              32             105              46             293             189 
#grey60           ivory       lightcyan      lightcyan1      lightgreen lightsteelblue1     lightyellow         magenta   mediumpurple3 
#130              49             143              55             123              55             120             230              56 
#midnightblue          orange      orangered4   paleturquoise            pink           plum1           plum2          purple             red 
#144              92              57              67             234              57              32             205             267 
#royalblue     saddlebrown          salmon         sienna3         skyblue        skyblue3       steelblue             tan       turquoise 
#120              78             151              60              84              57              77             170            1020 
#violet           white          yellow     yellowgreen 
#65              91             331              58 

pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/dissTOMColorClustering.pdf", width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar or choose not to merge
# Plot module similarity based on eigengene value

#Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower = 17)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster again and plot the results
METree = flashClust(as.dist(MEDiss), method = "average")

pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/eigengeneClustering1.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >85% eigengene similarity (most studies use80-90% similarity)

MEDissThres= 0.15 #merge modules that are 85% similar

pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/eigengeneClustering2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
#mergeCloseModules: Merging modules whose distance is less than 0.15
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 49 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 43 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 42 module eigengenes in given set.
#Calculating new MEs...
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 42 module eigengenes in given set.

# Check the structure of the merge object
print(names(merge))
#[1] "colors"    "dendro"    "oldDendro" "cutHeight" "oldMEs"    "newMEs"    "allOK"    

# Assign mergedColors
mergedColors = merge$colors


# Check if mergedColors is assigned correctly
print(head(mergedColors))
#[1] "darkmagenta" "salmon"      "darkmagenta" "pink"        "blue"        "blue"     

# Plot dendrogram and colors
pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/mergedClusters.pdf", width=20, height=20)
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
#42 modules, looking good

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/eigengeneClustering3.pdf")
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
#36 12

rownames(datTraits) <- treatmentinfo$sample_id
print(datTraits)

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=10)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)

# [1] "MEdarkolivegreen"  "MEdarkslateblue"   "MElightyellow"     "MEdarkred"         "MEturquoise"       "MEroyalblue"       "MEsalmon"         
#[8] "MEbisque4"         "MEmidnightblue"    "MEpaleturquoise"   "MEblack"           "MEpink"            "MEpurple"          "MEdarkmagenta"    
#[15] "MEblue"            "MEskyblue3"        "MEorangered4"      "MEskyblue"         "MEdarkorange"      "MEmediumpurple3"   "MEbrown4"         
#[22] "MEivory"           "MEyellowgreen"     "MEgrey60"          "MElightcyan"       "MEcyan"            "MEgreenyellow"     "MEgreen"          
#[29] "MEdarkturquoise"   "MEorange"          "MEbrown"           "MEsaddlebrown"     "MEdarkgrey"        "MEmagenta"         "MElightgreen"     
#[36] "MElightsteelblue1" "MEred"             "MEplum2"           "MEsteelblue"       "MEfloralwhite"     "MEdarkgreen"       "MEdarkorange2"    


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
#36 12

rownames(datTraits) <- treatmentinfo$sample_id
print(datTraits)

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=17)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)

# [1] "MEdarkolivegreen"  "MEdarkslateblue"   "MElightyellow"     "MEdarkred"         "MEturquoise"       "MEroyalblue"       "MEsalmon"         
#[8] "MEbisque4"         "MEmidnightblue"    "MEpaleturquoise"   "MEblack"           "MEpink"            "MEpurple"          "MEdarkmagenta"    
#[15] "MEblue"            "MEskyblue3"        "MEorangered4"      "MEskyblue"         "MEdarkorange"      "MEmediumpurple3"   "MEbrown4"         
#[22] "MEivory"           "MEyellowgreen"     "MEgrey60"          "MElightcyan"       "MEcyan"            "MEgreenyellow"     "MEgreen"          
#[29] "MEdarkturquoise"   "MEorange"          "MEbrown"           "MEsaddlebrown"     "MEdarkgrey"        "MEmagenta"         "MElightgreen"     
#[36] "MElightsteelblue1" "MEred"             "MEplum2"           "MEsteelblue"       "MEfloralwhite"     "MEdarkgreen"       "MEdarkorange2"    


moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/Life-stage clustering based on module-trait correlation.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)
#[,1]             [,2]             [,3]             [,4]            [,5]            [,6]            [,7]            [,8]           
#[1,] "-0.52\n(0.001)" "-0.52\n(0.001)" "-0.52\n(0.001)" "0.14\n(0.4)"   "0.13\n(0.4)"   "0.22\n(0.2)"   "0.22\n(0.2)"   "0.18\n(0.3)"  
#[2,] "0.44\n(0.007)"  "0.56\n(4e-04)"  "0.3\n(0.08)"    "0.04\n(0.8)"   "0.033\n(0.9)"  "0.19\n(0.3)"   "-0.12\n(0.5)"  "-0.22\n(0.2)" 
#[3,] "-0.071\n(0.7)"  "-0.27\n(0.1)"   "0.31\n(0.06)"   "-0.29\n(0.08)" "-0.33\n(0.05)" "0.67\n(7e-06)" "0.22\n(0.2)"   "-0.33\n(0.05)"
#[4,] "0.1\n(0.6)"     "0.11\n(0.5)"    "-0.26\n(0.1)"   "-0.27\n(0.1)"  "-0.31\n(0.06)" "0.62\n(5e-05)" "0.43\n(0.01)"  "-0.076\n(0.7)"
#[5,] "0.16\n(0.4)"    "0.13\n(0.5)"    "-0.3\n(0.07)"   "-0.25\n(0.1)"  "-0.3\n(0.07)"  "0.64\n(3e-05)" "0.28\n(0.1)"   "-0.39\n(0.02)"
#[6,] "-0.057\n(0.7)"  "-0.0055\n(1)"   "0.012\n(0.9)"   "-0.26\n(0.1)"  "-0.3\n(0.08)"  "0.57\n(2e-04)" "-0.067\n(0.7)" "0.29\n(0.09)" 
#[,9]            [,10]           [,11]           [,12]          
#[1,] "0.21\n(0.2)"   "0.19\n(0.3)"   "0.14\n(0.4)"   "0.15\n(0.4)"  
#[2,] "-0.18\n(0.3)"  "-0.35\n(0.03)" "-0.35\n(0.04)" "-0.33\n(0.05)"
#[3,] "0.11\n(0.5)"   "-0.08\n(0.6)"  "-0.1\n(0.5)"   "0.16\n(0.3)"  
#[4,] "-0.34\n(0.04)" "0.023\n(0.9)"  "-0.17\n(0.3)"  "0.15\n(0.4)"  
#[5,] "0.13\n(0.4)"   "-0.1\n(0.6)"   "-0.21\n(0.2)"  "0.23\n(0.2)"  
#[6,] "-0.23\n(0.2)"  "0.39\n(0.02)"  "-0.066\n(0.7)" "-0.28\n(0.09)"

pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/Module-trait-relationships.pdf")
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), 
               cex.lab.y= 0.55, cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , 
               zlim = c(-1,1), main = paste("Module-trait relationships"))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), cex.lab.y= 0.55, 
               cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , zlim = c(-1,1), 
               main = paste("Module-trait relationships"))
dev.off() 


#Testing by showing first from each particular group

# Reorder the columns to have larvae first, then meta, and then spat for every species

allTraits_mod <- list(
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
  
  
# Combine the individual lists into a data frame
datTraits_mod <- as.data.frame(allTraits_mod)

print(datTraits_mod)

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=10)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)

moduleTraitCor = cor(MEs, datTraits_mod, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/Life stage-species clustering based on module-trait correlation v2.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

pdf(file="C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/Life_stage_species_clustering.pdf", width = 12, height = 10)

par(mar = c(10, 5, 5, 2))  # Adjust margins to prevent label cutoff

plot(moduleTraitTree, 
     main = "Group clustering based on module-trait correlation", 
     sub="", 
     xlab="", 
     cex.lab = 1.8, 
     cex.axis = 1.5, 
     cex.main = 2, 
     las = 2)  # Rotate x-axis labels

dev.off()


moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Life-stage clustering based on module-trait correlation.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)
head(textMatrix)
#[,1]             [,2]            [,3]            [,4]            [,5]             [,6]           
#[1,] "-0.52\n(0.001)" "0.14\n(0.4)"   "0.22\n(0.2)"   "0.19\n(0.3)"   "-0.52\n(0.001)" "0.13\n(0.4)"  
#[2,] "0.44\n(0.007)"  "0.04\n(0.8)"   "-0.12\n(0.5)"  "-0.35\n(0.03)" "0.56\n(4e-04)"  "0.033\n(0.9)" 
#[3,] "-0.071\n(0.7)"  "-0.29\n(0.08)" "0.22\n(0.2)"   "-0.08\n(0.6)"  "-0.27\n(0.1)"   "-0.33\n(0.05)"
#[4,] "0.1\n(0.6)"     "-0.27\n(0.1)"  "0.43\n(0.01)"  "0.023\n(0.9)"  "0.11\n(0.5)"    "-0.31\n(0.06)"
#[5,] "0.16\n(0.4)"    "-0.25\n(0.1)"  "0.28\n(0.1)"   "-0.1\n(0.6)"   "0.13\n(0.5)"    "-0.3\n(0.07)" 
#[6,] "-0.057\n(0.7)"  "-0.26\n(0.1)"  "-0.067\n(0.7)" "0.39\n(0.02)"  "-0.0055\n(1)"   "-0.3\n(0.08)" 
#[,7]            [,8]            [,9]             [,10]           [,11]           [,12]          
#[1,] "0.18\n(0.3)"   "0.14\n(0.4)"   "-0.52\n(0.001)" "0.22\n(0.2)"   "0.21\n(0.2)"   "0.15\n(0.4)"  
#[2,] "-0.22\n(0.2)"  "-0.35\n(0.04)" "0.3\n(0.08)"    "0.19\n(0.3)"   "-0.18\n(0.3)"  "-0.33\n(0.05)"
#[3,] "-0.33\n(0.05)" "-0.1\n(0.5)"   "0.31\n(0.06)"   "0.67\n(7e-06)" "0.11\n(0.5)"   "0.16\n(0.3)"  
#[4,] "-0.076\n(0.7)" "-0.17\n(0.3)"  "-0.26\n(0.1)"   "0.62\n(5e-05)" "-0.34\n(0.04)" "0.15\n(0.4)"  
#[5,] "-0.39\n(0.02)" "-0.21\n(0.2)"  "-0.3\n(0.07)"   "0.64\n(3e-05)" "0.13\n(0.4)"   "0.23\n(0.2)"  
#[6,] "0.29\n(0.09)"  "-0.066\n(0.7)" "0.012\n(0.9)"   "0.57\n(2e-04)" "-0.23\n(0.2)"  "-0.28\n(0.09)"

pdf(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Module-trait-relationships_v2.pdf", width = 14, height = 10)

labeledHeatmap(
  Matrix = moduleTraitCor, 
  xLabels = names(datTraits_mod),  
  yLabels = names(MEs), 
  ySymbols = names(MEs), 
  cex.lab.y = 0.55,  # Increase y-axis label size
  cex.lab.x = 0.55,  # Increase x-axis label size
  colors = blueWhiteRed(50), 
  textMatrix = textMatrix, 
  setStdMargins = TRUE, 
  cex.text = 0.5,  # Increase text size inside the heatmap
  zlim = c(-1,1), 
  main = "Module-trait relationships"
)

dev.off()

pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Life_stage_species_clustering.pdf", width = 14  , height = 12)

par(mar = c(10, 5, 5, 2))  # Adjust margins to prevent label cutoff

plot(moduleTraitTree, 
     main = "Group clustering based on module-trait correlation", 
     sub="", 
     xlab="", 
     cex.lab = 1.8, 
     cex.axis = 1.5, 
     cex.main = 2, 
     las = 2)  # Rotate x-axis labels

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

pdf(file = "C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/Module-trait-relationship-heatmap3.pdf", height = 11.5, width = 8)

ht = Heatmap(moduleTraitCor, 
             name = "Eigengene", 
             column_title = "Module-Trait Eigengene Correlation", 
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


# Open the PDF device to save the plot
pdf(file = "C:/Users/amurgueitio/Documents/Multistage_Omics/Multi_species_WGCNA/Module-trait-relationship-heatmap3.pdf", 
    height = 11.5, width = 8)

# Generate the heatmap object
ht = Heatmap(moduleTraitCor, 
             name = "Eigengene", 
             column_title = "Module-Trait Eigengene Correlation", 
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

# Draw the heatmap to the PDF device
draw(ht)

# Close the PDF device to save the plot
dev.off()


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
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=10)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)

moduleTraitCor = cor(MEs, datTraits_justLifeStage, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Life stages clustering based on module-trait correlation v3.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)
head(textMatrix)
#[,1]           [,2]             [,3]           
#[1,] "0.011\n(0.9)" "-0.045\n(0.8)"  "0.034\n(0.8)" 
#[2,] "0.0035\n(1)"  "0.013\n(0.9)"   "-0.017\n(0.9)"
#[3,] "-0.13\n(0.4)" "-0.61\n(9e-05)" "0.74\n(3e-07)"
#[4,] "0.16\n(0.3)"  "-0.26\n(0.1)"   "0.098\n(0.6)" 
#[5,] "0.05\n(0.8)"  "-0.45\n(0.005)" "0.41\n(0.01)" 
#[6,] "0.0043\n(1)"  "-0.047\n(0.8)"  "0.043\n(0.8)" 
pdf(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Module-life stage relationships no species.pdf", width = 14, height = 10)

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

pdf(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Module-life stage-relationship-heatmap3.pdf", height = 11.5, width = 8)

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


#Re-adjusting columns

# Create list of p-values for eigengene correlation with specific life stages 
heatmappval <- signif(moduleTraitPvalue, 1)

htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

# library(dendsort)
row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

# Identify and set the column order as 'Larvae', 'Meta', 'Spat'
# Assuming your columns are named something like "Larvae", "Meta", "Spat"
desired_column_order <- c("Larvae", "Meta", "Spat")  # Adjust these names based on actual column names in your data

# Make sure the column names are in the correct order
column_order_indices <- match(desired_column_order, colnames(moduleTraitCor))

# Save heatmap to PDF
pdf(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Module-life_stage-relationship-heatmap3.pdf", height = 11.5, width = 8)

ht = Heatmap(moduleTraitCor, 
             name = "Eigengene", 
             column_title = "Module-Life stage Eigengene Correlation", 
             col = blueWhiteRed(50), 
             row_names_side = "left", 
             row_dend_side = "left",
             width = unit(4, "in"), 
             height = unit(8.5, "in"), 
             column_order = column_order_indices,  # Apply the new order here
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

# Create list of p-values for eigengene correlation with specific life stages 
heatmappval <- signif(moduleTraitPvalue, 1)

htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

# library(dendsort)
row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

# Identify and set the column order as 'Larvae', 'Meta', 'Spat'
# Assuming your columns are named something like "Larvae", "Meta", "Spat"
desired_column_order <- c("Larvae", "Meta", "Spat")  # Adjust these names based on actual column names in your data

# Make sure the column names are in the correct order
column_order_indices <- match(desired_column_order, colnames(moduleTraitCor))

# Save heatmap to PDF
pdf(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Module-life_stage-relationship-heatmap3.pdf", height = 11.5, width = 8)

ht = Heatmap(moduleTraitCor, 
             name = "Eigengene", 
             column_title = "Module-Life stage Eigengene Correlation", 
             col = blueWhiteRed(50), 
             row_names_side = "left", 
             row_dend_side = "left",
             width = unit(4, "in"), 
             height = unit(8.5, "in"), 
             column_order = column_order_indices,  # Apply the new order here
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

# Create list of p-values for eigengene correlation with specific life stages 
heatmappval <- signif(moduleTraitPvalue, 1)

htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

# library(dendsort)
row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

# Identify and set the column order as 'Larvae', 'Meta', 'Spat'
# Assuming your columns are named something like "Larvae", "Meta", "Spat"
desired_column_order <- c("Larvae", "Meta", "Spat")  # Adjust these names based on actual column names in your data

# Make sure the column names are in the correct order
column_order_indices <- match(desired_column_order, colnames(moduleTraitCor))

# Save heatmap to PDF
pdf(file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/Module-life_stage-relationship-heatmap3 v2.pdf", height = 11.5, width = 8)

ht = Heatmap(moduleTraitCor, 
             name = "Eigengene", 
             column_title = "Module-Life stage Eigengene Correlation", 
             col = blueWhiteRed(50), 
             row_names_side = "left", 
             row_dend_side = "left",
             width = unit(4, "in"), 
             height = unit(8.5, "in"), 
             column_order = column_order_indices,  # Apply the new order here
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

colnames(moduleTraitCor)

# Create dataframe that associates module colors with clusters based on the new heatmap
MEcluster1 <- data.frame(moduleColor = c("darkolivegreen"), moduleCluster = c(1))
MEcluster2 <- data.frame(moduleColor = c("darkslateblue"), moduleCluster = c(2))
MEcluster3 <- data.frame(moduleColor = c("lightyellow","darkred", "turquoise", "royalblue","salmon","bisque4", "midnightblue","paleturquoise","black","pink","purple","darkmagenta","blue","skyblue3"), moduleCluster = c(3))
MEcluster4 <- data.frame(moduleColor = c("orangered4"), moduleCluster = c(4))
MEcluster5 <- data.frame(moduleColor = c("skyblue", "darkorange", "mediumpurple3", "brown4","ivory"), moduleCluster = c(5))
MEcluster6 <- data.frame(moduleColor = c("yellowgreen"), moduleCluster = c(6))
MEcluster7 <- data.frame(moduleColor = c("grey60","lightcyan","cyan","greenyellow"), moduleCluster = c(7))
MEcluster8 <- data.frame(moduleColor = c("green", "darkturquoise","orange","brown", "saddlebrown", "darkgrey","magenta"), moduleCluster = c(8))
MEcluster9 <- data.frame(moduleColor = c("lightgreen"), moduleCluster = c(9))
MEcluster10 <- data.frame(moduleColor = c("lightsteelblue1", "plum2", "red", "steelblue", "floralwhite", "darkgreen", "darkorange2"), moduleCluster = c(10))

moduleCluster = bind_rows(MEcluster1, MEcluster2, MEcluster3, MEcluster4, MEcluster5, MEcluster6, MEcluster7, MEcluster8, MEcluster9, MEcluster10)



chooseTopHubInEachModule(
  datExpr, 
  Colors, 
  power = 10, 
  type = "signed")

#moduleCluster = bind_rows(MEcluster1, MEcluster2, MEcluster3, MEcluster4, MEcluster5, MEcluster6, MEcluster7, MEcluster8, MEcluster9, MEcluster10)
#bisque4           black            blue           brown          brown4            cyan       darkgreen 
#"OG0009145"     "OG0008065"     "OG0008976"     "OG0009965"     "OG0009675"     "OG0010874"     "OG0010523" 
#darkgrey     darkmagenta  darkolivegreen      darkorange     darkorange2         darkred   darkslateblue 
#"OG0005819"     "OG0012697"     "OG0012759"     "OG0010837"     "OG0009824"     "OG0011576"     "OG0005793" 
#darkturquoise     floralwhite           green     greenyellow          grey60           ivory       lightcyan 
#"OG0008175"     "OG0006584"     "OG0009384"     "OG0005481"     "OG0011233"     "OG0007705"     "OG0009822" 
#lightgreen lightsteelblue1     lightyellow         magenta   mediumpurple3    midnightblue          orange 
#"OG0006154"     "OG0011559"     "OG0011846"     "OG0007431"     "OG0009291"     "OG0011109"     "OG0012502" 
#orangered4   paleturquoise            pink           plum2          purple             red       royalblue 
#"OG0009830"     "OG0009622"     "OG0006854"     "OG0009122"     "OG0011257"     "OG0007346"     "OG0009401" 
#saddlebrown          salmon         skyblue        skyblue3       steelblue       turquoise     yellowgreen 
#"OG0012210"     "OG0007290"     "OG0008849"     "OG0007104"     "OG0010950"     "OG0005704"     "OG0007997" 

head(moduleCluster) 

#   moduleColor moduleCluster
#1 darkolivegreen             1
#2  darkslateblue             2
#3    lightyellow             3
#4        darkred             3
#5      turquoise             3
#6      royalblue             3

# View module eigengene data
head(MEs)
#MEdarkolivegreen MEdarkslateblue MElightyellow  MEdarkred MEturquoise  MEroyalblue   MEsalmon    MEbisque4
#Larva1_Atenuis       -0.2855066       0.2478781   -0.04414259 0.04843330  0.08172188  0.006446895 0.05743817  0.018910620
#Larva2_Atenuis       -0.2867903       0.2280316   -0.04242803 0.05257476  0.10051992 -0.036746615 0.02412668 -0.006716871
#Larva3_Atenuis       -0.2873017       0.2602053   -0.03040436 0.06831929  0.07643965 -0.064433935 0.02746105  0.091172681
#Meta1_Atenuis        -0.2871411       0.2765646   -0.17086580 0.05680712  0.07752099  0.012155348 0.05933017 -0.008878983
#Meta2_Atenuis        -0.2870194       0.3205526   -0.13337134 0.05789281  0.07304967 -0.013096854 0.06274494  0.020449007
#Meta3_Atenuis        -0.2862190       0.3250330   -0.14274554 0.07365532  0.05757292 -0.008196814 0.05654361  0.025607501

names(MEs)
# [1] "MEdarkolivegreen"  "MEdarkslateblue"   "MElightyellow"     "MEdarkred"         "MEturquoise"      
#[6] "MEroyalblue"       "MEsalmon"          "MEbisque4"         "MEmidnightblue"    "MEpaleturquoise"  
#[11] "MEblack"           "MEpink"            "MEpurple"          "MEdarkmagenta"     "MEblue"           
#[16] "MEskyblue3"        "MEorangered4"      "MEskyblue"         "MEdarkorange"      "MEmediumpurple3"  
#[21] "MEbrown4"          "MEivory"           "MEyellowgreen"     "MEgrey60"          "MElightcyan"      
#[26] "MEcyan"            "MEgreenyellow"     "MEgreen"           "MEdarkturquoise"   "MEorange"         
#[31] "MEbrown"           "MEsaddlebrown"     "MEdarkgrey"        "MEmagenta"         "MElightgreen"     
#[36] "MElightsteelblue1" "MEred"             "MEplum2"           "MEsteelblue"       "MEfloralwhite"    
#[41] "MEdarkgreen"       "MEdarkorange2"    

Strader_MEs <- MEs
Strader_MEs$timepoint <- treatmentinfo$timepoint
Strader_MEs$sample_id <- rownames(Strader_MEs)
Strader_MEs$group <- treatmentinfo$group
head(Strader_MEs)
#Looks good
head(Strader_MEs$sample_id)
#[1] "Larva1_Atenuis" "Larva2_Atenuis" "Larva3_Atenuis" "Meta1_Atenuis"  "Meta2_Atenuis"  "Meta3_Atenuis" 
head(Strader_MEs$timepoint)
#[1] I  I  I  II II II
#Levels: I II III
head(Strader_MEs$group)
#[1] "Larvae_Atenuis" "Larvae_Atenuis" "Larvae_Atenuis" "Meta_Atenuis"   "Meta_Atenuis"   "Meta_Atenuis"  

# Select modules based on the new cluster assignments
C1_Strader_MEs <- select(Strader_MEs, MEdarkolivegreen)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)

C2_Strader_MEs <- select(Strader_MEs, MEdarkslateblue)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)

C3_Strader_MEs <- select(Strader_MEs, MElightyellow, MEdarkred, MEturquoise, MEroyalblue, MEsalmon, MEbisque4, 
                         MEmidnightblue, MEpaleturquoise, MEblack, MEpink, MEpurple, MEdarkmagenta, MEblue, MEskyblue3)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)

C4_Strader_MEs <- select(Strader_MEs, MEorangered4)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)

C5_Strader_MEs <- select(Strader_MEs, MEskyblue, MEdarkorange, MEmediumpurple3, MEbrown4, MEivory)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)

C6_Strader_MEs <- select(Strader_MEs, MEyellowgreen)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)

C7_Strader_MEs <- select(Strader_MEs, MEgrey60, MElightcyan, MEcyan, MEgreenyellow)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)

C8_Strader_MEs <- select(Strader_MEs, MEgreen, MEdarkturquoise, MEorange, MEbrown, MEsaddlebrown, MEdarkgrey, MEmagenta)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)

C9_Strader_MEs <- select(Strader_MEs, MElightgreen)
C9_Strader_MEs$Mean <- rowMeans(C9_Strader_MEs)

C10_Strader_MEs <- select(Strader_MEs, MElightsteelblue1, MEplum2, MEred, MEsteelblue, MEfloralwhite, MEdarkgreen, MEdarkorange2)
C10_Strader_MEs$Mean <- rowMeans(C10_Strader_MEs)

# Create the final dataframe with new clusters and timepoint
expressionProfile_data_timepoint <- as.data.frame(cbind(
  group = Strader_MEs$timepoint, 
  cluster1 = C1_Strader_MEs$Mean, 
  cluster2 = C2_Strader_MEs$Mean, 
  cluster3 = C3_Strader_MEs$Mean, 
  cluster4 = C4_Strader_MEs$Mean,
  cluster5 = C5_Strader_MEs$Mean,
  cluster6 = C6_Strader_MEs$Mean,
  cluster7 = C7_Strader_MEs$Mean,
  cluster8 = C8_Strader_MEs$Mean,
  cluster9 = C9_Strader_MEs$Mean,
  cluster10 = C10_Strader_MEs$Mean
))



# Create the final dataframe with group rather than timepoint
expressionProfile_data_group <- as.data.frame(cbind(
  group = Strader_MEs$group, 
  cluster1 = C1_Strader_MEs$Mean, 
  cluster2 = C2_Strader_MEs$Mean, 
  cluster3 = C3_Strader_MEs$Mean, 
  cluster4 = C4_Strader_MEs$Mean,
  cluster5 = C5_Strader_MEs$Mean,
  cluster6 = C6_Strader_MEs$Mean,
  cluster7 = C7_Strader_MEs$Mean,
  cluster8 = C8_Strader_MEs$Mean,
  cluster9 = C9_Strader_MEs$Mean,
  cluster10 = C10_Strader_MEs$Mean
))

# Extract the mean columns from each cluster data frame
cluster_means <- cbind(
  C1_Strader_MEs$Mean,
  C2_Strader_MEs$Mean,
  C3_Strader_MEs$Mean,
  C4_Strader_MEs$Mean,
  C5_Strader_MEs$Mean,
  C6_Strader_MEs$Mean,
  C7_Strader_MEs$Mean,
  C8_Strader_MEs$Mean,
  C9_Strader_MEs$Mean,
  C10_Strader_MEs$Mean
)

# Assuming you have a list or data frame 'expressionProfile_data' with cluster data
# And your cluster-specific data frames are named like C1_Strader_MEs, C2_Strader_MEs, etc.This is a loop, run it together. 

for (i in 1:10) {
  cluster_name = paste0("cluster", i)
  strader_ME_name = paste0("C", i, "_Strader_MEs")
  
  # Extract mean eigengene values for the current cluster
  meanEigenClust <- expressionProfile_data_group[[cluster_name]]
  
  # Save the mean eigengene values to a CSV file
  meanEigenFilename = paste0("meanEigenClustgroup", i, ".csv")
  write.csv(meanEigenClust, file = meanEigenFilename)
  
  # Save the Strader ME data for the current cluster to a CSV file
  straderMEFilename = paste0("C", i, "_Strader_MEsgroup.csv")
  write.csv(get(strader_ME_name), file = straderMEFilename)
}


expressionProfile_data_timepoint$group <- as.factor(expressionProfile_data_timepoint$group)

group_order = c("1","2","3")

Cluster1Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster1) %>% 
  #group_by(group) %>% 
  ggplot(aes(x=group, y=cluster1, fill=group)) +
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5), alpha = 0.7) +
  #stat_summary(fun=mean, geom="line", aes(group=group, color = group), position = position_dodge(width = 0.5))  + 
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch = 21, size=5, position = position_dodge(width = 1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + #Axis titles
  #ylim(-0.5,1) +
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster1") +
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, color = "black"))+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") 

Cluster2Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster2) %>% 
  ggplot(aes(x=group, y=cluster2, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster2") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")  


Cluster3Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster3) %>% 
  ggplot(aes(x=group, y=cluster3, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster3") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")  

Cluster4Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster4) %>% 
  ggplot(aes(x=group, y=cluster4, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster4") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")  

Cluster5Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster5) %>% 
  ggplot(aes(x=group, y=cluster5, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster5") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")  

Cluster6Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster6) %>% 
  ggplot(aes(x=group, y=cluster6, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster6") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster7Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster7) %>% 
  ggplot(aes(x=group, y=cluster7, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster7") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")  

Cluster8Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster8) %>% 
  ggplot(aes(x=group, y=cluster8, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster8") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")  

Cluster9Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster9) %>% 
  ggplot(aes(x=group, y=cluster9, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster9") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")  

Cluster10Plot <- expressionProfile_data_timepoint %>%
  select(group, cluster10) %>% 
  ggplot(aes(x=group, y=cluster10, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + 
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster10") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")


library(patchwork)

expressionProfiles <- (Cluster1Plot + Cluster2Plot + Cluster3Plot) / (Cluster4Plot + Cluster5Plot + Cluster6Plot) / (Cluster7Plot + Cluster8Plot + Cluster9Plot)/ (Cluster10Plot+Cluster10Plot+Cluster10Plot)
ggsave("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA/expression_eigengene_Profiles_withMean.pdf", expressionProfiles, height = 25, width = 28, units = "in")

