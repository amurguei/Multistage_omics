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


setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_GFF")

#Upload data--------------------------------------------------------------------------

treatmentinfo <- read_csv("~/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/5-Pacu-SampleInfo.csv")

#gene count matrix
gcount <- as.data.frame(read.csv("Pacu_gene_count_matrix_newGFF.csv", row.names="gene_id"), colClasses = double, header=TRUE)

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
nrow(gcount) #Before 33705                                        

nrow(gcount_filt) #After 20891

# Normalize our read counts using VST-normalization in DESeq2
# Construct the DESeq2 dataset

#Merge the timepoint columns into a new column, group. Set group as a factor.
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I","II","III"))

#Create a DESeqDataSet design from gene count matrix and labels. Here we set the design to look at 
#any differences in gene expression across samples attributed to depth.

#Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = treatmentinfo,
                               design = ~timepoint)

# Log-transform the count data using a variance stabilizing transforamtion (vst). 

SF.gdds <- estimateSizeFactors( gdds ) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than for to use vst
print(sizeFactors(SF.gdds)) #View size factors
#SRR3051863 SRR3051864 SRR3051865 SRR3051866 SRR3051867 SRR3051868 SRR3051869 SRR3051870 SRR3051871 
#0.8683862  1.1838339  1.3386539  0.7629672  0.7573425  0.8809547  1.1700139  1.1345417  1.1324302 

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


#Cluster the samples to look for obvious outliers
#Look for outliers by examining the sample tree:

sampleTree = hclust(dist(datExpr), method = "average")

#Plot the sample tree
pdf(paste0('sampleTree','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()


#PCA--------------------------------------------------------------------------------------------------------------
gPCAdata <- plotPCA(gvst, intgroup = c("timepoint"), returnData=TRUE)

percentVar <- round(100 * attr(gPCAdata, "percentVar"))

allgenesfilt_PCA <- ggplot(gPCAdata, aes(PC1, PC2, shape = timepoint)) + 
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_shape_manual(values = c("I" = 1, "II" = 2, "III" = 3)) +
  xlim(-40, 40) + 
  ylim(-40, 40) +
  coord_fixed() +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.background = element_blank()
  ) +
  theme(legend.position = "none")

allgenesfilt_PCA

ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_GFF/PCA_timepointntop=1000.png", allgenesfilt_PCA, width = 11, height = 8)

#Other PCA

levels(gPCAdata$timepoint)

allgenesfilt_PCA_visual <- 
  ggplot(data=gPCAdata, aes(PC1, PC2, shape=timepoint, colour=timepoint)) + 
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  
  scale_shape_manual(values = c("I" = 1, "II" = 2, "III" = 3)) +
  
  scale_colour_manual(values = c("I" = 1, "II" = 2, "III" = 3)) +
  ylim(-50, 50)+
  coord_fixed()+
  theme_classic() + #Set background color
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines 
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) # + #Set the plot background
#theme(legend.position = ("none")) #set title attributes
allgenesfilt_PCA_visual


ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_GFF/PCA_timepointntop=1000.png", allgenesfilt_PCA_visual, width = 11, height = 8)


#some AI modification
library(ggplot2)

# Check the levels of the timepoint variable
levels(gPCAdata$timepoint)

# Explicitly set the levels of the timepoint variable
gPCAdata$timepoint <- factor(gPCAdata$timepoint, levels = c("I", "II", "III"))

# Create the PCA plot
allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(PC1, PC2, shape = timepoint, colour = timepoint)) + 
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_shape_manual(values = c("I" = 1, "II" = 2, "III" = 3)) +
  scale_colour_manual(values = c("I" = 1, "II" = 2, "III" = 3)) +
  ylim(-50, 50) +
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

allgenesfilt_PCA_visual <- 
  ggplot(data = gPCAdata, aes(PC1, PC2)) + 
  geom_point(aes(shape = timepoint, colour = timepoint), size = 5) +
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
  )

print(allgenesfilt_PCA_visual)

ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_GFF/PCA_timepointntop=1000.png", allgenesfilt_PCA_visual, width = 11, height = 8)

#Choose a set of soft-thresholding powers
powers <- c(seq(from = 1, to=19, by=2), c(21:30)) #Create a string of numbers from 1 through 10, and even numbers from 10 through 20
powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2))

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
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")# dev.off()


#The lowest scale-free topology fit index R^2 recommended by Langfelder and Horvath is 0.8. 
#From the graph, it appears that our soft thresholding power is 19 because it is the lowest 
#power before the R^2=0.8 threshold that maximizes with model fit (s17 is the number right above the red line).

### Network construction and module detection:
# Co-expression adjacency and topological overlap matrix similarity
# Co-expression similarity and adjacency, using the soft thresholding power 19 and translate the adjacency into topological overlap matrix to calculate 
# the corresponding dissimilarity. I will use a signed network because we have a relatively high softPower, according 
# to >12 (https://peterlangfelder.com/2018/11/25/__trashed/). 
# Moreover, in expression data where you are interested in when expression on one gene increases or decreases with expression level of another you would use a signed network (when you are interested in the direction of change, correlation and anti-correlation, you use a signed network).

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


#Run analysis. 

softPower=17 #Set softPower to 17
adjacency=adjacency(datExpr, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(adjacency, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/New_genome_P_acutaadjTOM.RData")
save(TOM, file ="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/New_genome_P_acutaTOM.RData")
save(dissTOM, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/New_genome_P_acutadisstom.RData") 
load(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/New_genome_P_acutadisstom.RData")


#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/dissTOMClustering.pdf", width=20, height=20)
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
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19 
#1 2472 1401 1130 1044 1034 1020  843  732  730  729  646  604  598  570  549  494  462  456  450 
#20   21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39 
#442  436  407  391  388  361  309  303  289  283  275  205  179  164  120  116   86   67   57   48 

save(dynamicMods, geneTree, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/dyMod_geneTree.RData")

dyMod_geneTree <- load(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/dyMod_geneTree.RData")

dyMod_geneTree
## [1] "dynamicMods" "geneTree"

# Plot the module assignment under the gene dendrogram
dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)

#black           blue          brown           cyan      darkgreen       darkgrey 
#843           1401           1130            570            407            388 
#darkmagenta darkolivegreen     darkorange        darkred  darkturquoise          green 
#120            164            309            436            391           1034 
#greenyellow           grey         grey60      lightcyan     lightgreen    lightyellow 
#646              1            462            494            456            450 
#magenta   midnightblue         orange     orangered4  paleturquoise           pink 
#730            549            361             48            205            732 
#plum1         purple            red      royalblue    saddlebrown         salmon 
#57            729           1020            442            283            598 
#sienna3        skyblue       skyblue3      steelblue            tan      turquoise 
#116            289             67            275            604           2472 
#violet          white         yellow    yellowgreen 
#179            303           1044             86 
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/dissTOMColorClustering.pdf", width=20, height=20)
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

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/eigengeneClustering1.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >85% eigengene similarity (most studies use80-90% similarity)

MEDissThres= 0.15 #merge modules that are 85% similar

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/eigengeneClustering2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

#Merging modules whose distance is less than 0.15
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 40 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 19 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 16 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 15 module eigengenes in given set.
#Calculating new MEs...
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 15 module eigengenes in given set.
#mergedColors= merge$colors
#mergedMEs= merge$newMEs

#Troubleshooting

# Run mergeCloseModules
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight=MEDissThres, verbose=3)

# Check the structure of the merge object
print(names(merge))
#[1] "colors"    "dendro"    "oldDendro" "cutHeight" "oldMEs"    "newMEs"    "allOK"    

# Assign mergedColors
mergedColors = merge$colors


# Check if mergedColors is assigned correctly
print(head(mergedColors))
#[1] "salmon"         "saddlebrown"    "darkolivegreen" "lightyellow"    "lightyellow"   
#[6] "green

# Plot dendrogram and colors
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/mergedClusters.pdf", width=20, height=20)
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
#print(head(MEs))
#MEblack MEgreenyellow MElightyellow MEpaleturquoise    MEgreen MEdarkorange MEsaddlebrown
#SRR3051863 -0.4407852   -0.46848430    0.47066636      0.29508062 -0.2016300  0.224244257    0.11576492
#SRR3051864 -0.4531115   -0.44811825    0.46303303      0.31473068 -0.2437415  0.119888251    0.06494867
#SRR3051865 -0.4330213   -0.42969863    0.42788466      0.35581076 -0.2339743  0.142629798    0.05759451
#SRR3051866  0.3882864    0.10940587   -0.11908833      0.08892962  0.4865318  0.315576854    0.36732763
#SRR3051867  0.4090003    0.05656749   -0.04051124     -0.13308778  0.5974472  0.494796738    0.22929293
#SRR3051868  0.2666735    0.14118551   -0.25319891      0.34735544  0.2958598  0.007137438    0.47062872
#MEorangered4 MEyellowgreen   MEsalmon MEdarkmagenta MEdarkolivegreen   MEsienna3  MEskyblue3
#SRR3051863  -0.01762650    0.11362612 -0.1810208    0.25304733        0.3085635  0.16442661  0.30144461
#SRR3051864   0.06695769    0.11924942 -0.1247619    0.21550504        0.4033185  0.19678147  0.16288701
#SRR3051865   0.24075081    0.18008516 -0.1090534    0.25624627        0.3642847  0.17803709  0.07480525
#SRR3051866  -0.31464770    0.07383415 -0.3342176   -0.41601491       -0.4701605 -0.91771482 -0.61829246
#SRR3051867  -0.60576308   -0.93804010 -0.3667384   -0.05277386       -0.5001653 -0.06980354 -0.60281274
#SRR3051868   0.63921226    0.13075065 -0.2562510   -0.75071064       -0.3442587  0.05273902  0.26574807
#MEgrey
#SRR3051863 -0.15039310
#SRR3051864 -0.09770823
#SRR3051865  0.24429827
#SRR3051866 -0.33071748
#SRR3051867  0.12804038
#SRR3051868  0.19401132
ncol(MEs) 
# [1] 15
#Instead of 40 modules, we now have 15

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/eigengeneClustering3.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# Relating modules to group, quantifying module–trait associations
#Prepare trait data. Data has to be numeric, so I replaced the group for numeric values

allTraits <- names(treatmentinfo$timepoint)
allTraits$larvae_released <- c(1,1,1,0,0,0,0,0,0)
allTraits$lavae_compressed <- c(0,0,0,1,1,1,0,0,0)
allTraits$spat <- c(0,0,0,0,0,0,1,1,1)


datTraits <- as.data.frame(allTraits)
dim(datTraits)
#[1] 9 3

rownames(datTraits) <- treatmentinfo$sampleID
print(datTraits)
#larvae_released lavae_compressed spat
#1               1                0    0
#2               1                0    0
#3               1                0    0
#4               0                1    0
#5               0                1    0
#6               0                1    0
#7               0                0    1
#8               0                0    1
#9               0                0    1

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=17)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)

#[1] "MEblack"          "MEgreenyellow"    "MElightyellow"    "MEpaleturquoise"  "MEgreen"         
#[6] "MEdarkorange"     "MEsaddlebrown"    "MEorangered4"     "MEyellowgreen"    "MEsalmon"        
#[11] "MEdarkmagenta"    "MEdarkolivegreen" "MEsienna3"        "MEskyblue3"       "MEgrey" 

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/Life-stage clustering based on module-trait correlation.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)

#[,1]             [,2]            [,3]            
#[1,] "-0.94\n(2e-04)" "0.75\n(0.02)"  "0.19\n(0.6)"   
#[2,] "-0.95\n(8e-05)" "0.22\n(0.6)"   "0.73\n(0.02)"  
#[3,] "0.96\n(3e-05)"  "-0.29\n(0.4)"  "-0.67\n(0.05)" 
#[4,] "0.68\n(0.04)"   "0.21\n(0.6)"   "-0.9\n(0.001)" 
#[5,] "-0.48\n(0.2)"   "0.98\n(7e-06)" "-0.5\n(0.2)"   
#[6,] "0.34\n(0.4)"    "0.58\n(0.1)"   "-0.92\n(4e-04)"

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/Module-trait-relationships.pdf")
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), 
               cex.lab.y= 0.55, cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , 
               zlim = c(-1,1), main = paste("Module-trait relationships"))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),  yLabels = names(MEs), ySymbols = names(MEs), cex.lab.y= 0.55, 
               cex.lab.x= 0.55, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.25, textAdj = , zlim = c(-1,1), 
               main = paste("Module-trait relationships"))
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

pdf(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/Module-trait-relationship-heatmap3.pdf", height = 11.5, width = 8)
ht=Heatmap(moduleTraitCor, name = "Eigengene", column_title = "Module-Trait Eigengene Correlation", 
           col = blueWhiteRed(50), 
           row_names_side = "left", row_dend_side = "left",
           width = unit(4, "in"), height = unit(8.5, "in"), 
           column_order = 1:3, column_dend_reorder = FALSE, cluster_columns = hclust(dist(t(moduleTraitCor)), method = "average"), column_split = 3, column_dend_height = unit(0.5, "in"),
           cluster_rows = METree, row_split = 10, row_gap = unit(2.5, "mm"), border = TRUE,
           cell_fun = function(j, i, x, y, w, h, col) {
             if(heatmappval[i, j] <= 0.05) {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
             }
             else {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
             }},
           column_names_gp =  gpar(fontsize = 10),
           row_names_gp = gpar(fontsize = 10, alpha = 0.75, border = TRUE, fill = htmap.colors))
draw(ht)
dev.off()


# Create dataframe that associates module colors with clusters based on the above heatmap and the MEtree dendrogram.I did this manually looking at the clusters in the previous plot.
MEcluster1 <- data.frame(moduleColor = c("black", "greenyellow"), moduleCluster = c(1))
MEcluster2 <- data.frame(moduleColor = c("lightyellow", "paleturquoise"), moduleCluster = c(2))
MEcluster3 <- data.frame(moduleColor = c("green"), moduleCluster = c(3))
MEcluster4 <- data.frame(moduleColor = c("darkorange","saddlebrown"), moduleCluster = c(4))
MEcluster5 <- data.frame(moduleColor = c("grey"), moduleCluster = c(5))
MEcluster6 <- data.frame(moduleColor = c("orangered4"), moduleCluster = c(6)) 
MEcluster7 <- data.frame(moduleColor = c("yellowgreen"), moduleCluster = c(7)) 
MEcluster8 <- data.frame(moduleColor = c("salmon"), moduleCluster = c(8)) 
MEcluster9 <- data.frame(moduleColor = c("darkmagenta","darkolivegreen"), moduleCluster = c(9)) 
MEcluster10 <- data.frame(moduleColor = c("sienna3", "skyblue3"), moduleCluster = c(10)) 

moduleCluster = bind_rows(MEcluster1, MEcluster2, MEcluster3, MEcluster4, MEcluster5, MEcluster6, MEcluster7, MEcluster8, MEcluster9, MEcluster10)

head(moduleCluster) 
#   moduleColor moduleCluster
#1         black             1
#2   greenyellow             1
#3   lightyellow             2
#4 paleturquoise             2
#5         green             3
#6    darkorange             4

# View module eigengene data
head(MEs)

#MEblack MEgreenyellow MElightyellow MEpaleturquoise    MEgreen MEdarkorange
#SRR3051863 -0.4407852   -0.46848430    0.47066636      0.29508062 -0.2016300  0.224244257
#SRR3051864 -0.4531115   -0.44811825    0.46303303      0.31473068 -0.2437415  0.119888251
#SRR3051865 -0.4330213   -0.42969863    0.42788466      0.35581076 -0.2339743  0.142629798
#SRR3051866  0.3882864    0.10940587   -0.11908833      0.08892962  0.4865318  0.315576854
#SRR3051867  0.4090003    0.05656749   -0.04051124     -0.13308778  0.5974472  0.494796738
#SRR3051868  0.2666735    0.14118551   -0.25319891      0.34735544  0.2958598  0.007137438
#MEsaddlebrown MEorangered4 MEyellowgreen   MEsalmon MEdarkmagenta MEdarkolivegreen
#SRR3051863    0.11576492  -0.01762650    0.11362612 -0.1810208    0.25304733        0.3085635
#SRR3051864    0.06494867   0.06695769    0.11924942 -0.1247619    0.21550504        0.4033185
#SRR3051865    0.05759451   0.24075081    0.18008516 -0.1090534    0.25624627        0.3642847
#SRR3051866    0.36732763  -0.31464770    0.07383415 -0.3342176   -0.41601491       -0.4701605
#SRR3051867    0.22929293  -0.60576308   -0.93804010 -0.3667384   -0.05277386       -0.5001653
#SRR3051868    0.47062872   0.63921226    0.13075065 -0.2562510   -0.75071064       -0.3442587
#MEsienna3  MEskyblue3      MEgrey
#SRR3051863  0.16442661  0.30144461 -0.15039310
#SRR3051864  0.19678147  0.16288701 -0.09770823
#SRR3051865  0.17803709  0.07480525  0.24429827
#SRR3051866 -0.91771482 -0.61829246 -0.33071748
#SRR3051867 -0.06980354 -0.60281274  0.12804038
#SRR3051868  0.05273902  0.26574807  0.19401132
names(MEs)
#[1] "MEblack"          "MEgreenyellow"    "MElightyellow"    "MEpaleturquoise" 
#[5] "MEgreen"          "MEdarkorange"     "MEsaddlebrown"    "MEorangered4"    
#[9] "MEyellowgreen"    "MEsalmon"         "MEdarkmagenta"    "MEdarkolivegreen"
#[13] "MEsienna3"        "MEskyblue3"       "MEgrey" 

Strader_MEs <- MEs
Strader_MEs$timepoint <- treatmentinfo$timepoint
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)

            #MEblack MEgreenyellow MElightyellow MEpaleturquoise    MEgreen MEdarkorange
#SRR3051863 -0.4407852   -0.46848430    0.47066636      0.29508062 -0.2016300  0.224244257
#SRR3051864 -0.4531115   -0.44811825    0.46303303      0.31473068 -0.2437415  0.119888251
#SRR3051865 -0.4330213   -0.42969863    0.42788466      0.35581076 -0.2339743  0.142629798
#SRR3051866  0.3882864    0.10940587   -0.11908833      0.08892962  0.4865318  0.315576854
#SRR3051867  0.4090003    0.05656749   -0.04051124     -0.13308778  0.5974472  0.494796738
#SRR3051868  0.2666735    0.14118551   -0.25319891      0.34735544  0.2958598  0.007137438
#MEsaddlebrown MEorangered4 MEyellowgreen   MEsalmon MEdarkmagenta MEdarkolivegreen
#SRR3051863    0.11576492  -0.01762650    0.11362612 -0.1810208    0.25304733        0.3085635
#SRR3051864    0.06494867   0.06695769    0.11924942 -0.1247619    0.21550504        0.4033185
#SRR3051865    0.05759451   0.24075081    0.18008516 -0.1090534    0.25624627        0.3642847
#SRR3051866    0.36732763  -0.31464770    0.07383415 -0.3342176   -0.41601491       -0.4701605
#SRR3051867    0.22929293  -0.60576308   -0.93804010 -0.3667384   -0.05277386       -0.5001653
#SRR3051868    0.47062872   0.63921226    0.13075065 -0.2562510   -0.75071064       -0.3442587
#MEsienna3  MEskyblue3      MEgrey timepoint  sample_id
#SRR3051863  0.16442661  0.30144461 -0.15039310         I SRR3051863
#SRR3051864  0.19678147  0.16288701 -0.09770823         I SRR3051864
#SRR3051865  0.17803709  0.07480525  0.24429827         I SRR3051865
#SRR3051866 -0.91771482 -0.61829246 -0.33071748        II SRR3051866
#SRR3051867 -0.06980354 -0.60281274  0.12804038        II SRR3051867
#SRR3051868  0.05273902  0.26574807  0.19401132        II SRR3051868

head(Strader_MEs$sample_id)
#[1] "SRR3051863" "SRR3051864" "SRR3051865" "SRR3051866" "SRR3051867" "SRR3051868"

head(Strader_MEs$timepoint)
#[1] I  I  I  II II II
#Levels: I II III
# Calculate 10 over-arching expression patterns using mean eigengene for each module in a cluster. Use MEcluster as models for first lines
# Create a column to the Strader_MEs data frame containing age-depth groups

C1_Strader_MEs <- select(Strader_MEs, MEblack: MEgreenyellow)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)
C2_Strader_MEs <- select(Strader_MEs, MElightyellow: MEpaleturquoise)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)
C3_Strader_MEs <- select(Strader_MEs, MEgreen)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)
C4_Strader_MEs <- select(Strader_MEs, MEdarkorange: MEsaddlebrown)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)
C5_Strader_MEs <- select(Strader_MEs, MEgrey)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)
C6_Strader_MEs <- select(Strader_MEs,MEorangered4)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)
C7_Strader_MEs <- select(Strader_MEs, MEyellowgreen)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)
C8_Strader_MEs <- select(Strader_MEs, MEsalmon)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)
C9_Strader_MEs <- select(Strader_MEs, MEdarkmagenta, MEdarkolivegreen)
C9_Strader_MEs$Mean <- rowMeans(C9_Strader_MEs)
C10_Strader_MEs <- select(Strader_MEs, MEsienna3: MEskyblue3)
C10_Strader_MEs$Mean <- rowMeans(C10_Strader_MEs)

expressionProfile_data <- as.data.frame(cbind(group = Strader_MEs$timepoint, cluster1= C1_Strader_MEs$Mean, cluster2 = C2_Strader_MEs$Mean, 
                                              cluster3 = C3_Strader_MEs$Mean, cluster4 = C4_Strader_MEs$Mean,cluster5 = C5_Strader_MEs$Mean,cluster6 = C6_Strader_MEs$Mean,
                                              cluster7 = C7_Strader_MEs$Mean,cluster8 = C8_Strader_MEs$Mean,cluster9 = C9_Strader_MEs$Mean,cluster10 = C10_Strader_MEs$Mean))



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
  meanEigenClust <- expressionProfile_data[[cluster_name]]
  
  # Save the mean eigengene values to a CSV file
  meanEigenFilename = paste0("meanEigenClust", i, ".csv")
  write.csv(meanEigenClust, file = meanEigenFilename)
  
  # Save the Strader ME data for the current cluster to a CSV file
  straderMEFilename = paste0("C", i, "_Strader_MEs.csv")
  write.csv(get(strader_ME_name), file = straderMEFilename)
}

library(ggplot2)

expressionProfile_data$group <- as.factor(expressionProfile_data$group)

group_order = c("1","2","3")

Cluster1Plot <- expressionProfile_data %>%
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

Cluster2Plot <- expressionProfile_data %>%
  select(group, cluster2) %>% 
  #group_by(group) %>% 
  ggplot(aes(x=group, y=cluster2, fill=group)) +
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
  ggtitle("Cluster 2") +
  theme_bw() + theme(panel.border = element_rect(color="black", fill=NA, size=0.75), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, color = "black"))+
  geom_hline(yintercept = 0, linetype="dashed", color = "grey") 

Cluster3Plot <- expressionProfile_data %>%
  select(group, cluster3) %>% 
  ggplot(aes(x=group, y=cluster3, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + # Axis titles
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster 3") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster4Plot <- expressionProfile_data %>%
  select(group, cluster4) %>% 
  ggplot(aes(x=group, y=cluster4, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + # Axis titles
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster 4") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")


Cluster5Plot <- expressionProfile_data %>%
  select(group, cluster5) %>% 
  ggplot(aes(x=group, y=cluster5, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + # Axis titles
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster 5") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster6Plot <- expressionProfile_data %>%
  select(group, cluster6) %>% 
  ggplot(aes(x=group, y=cluster6, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + # Axis titles
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster 6") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster7Plot <- expressionProfile_data %>%
  select(group, cluster7) %>% 
  ggplot(aes(x=group, y=cluster7, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + # Axis titles
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster 7") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster8Plot <- expressionProfile_data %>%
  select(group, cluster8) %>% 
  ggplot(aes(x=group, y=cluster8, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + # Axis titles
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster 8") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster9Plot <- expressionProfile_data %>%
  select(group, cluster9) %>% 
  ggplot(aes(x=group, y=cluster9, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + # Axis titles
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster 9") +
  theme_bw() + 
  theme(panel.border=element_rect(color="black", fill=NA, size=0.75), 
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        axis.line=element_blank()) +
  theme(axis.text=element_text(size=15, color="black"),
        axis.title=element_text(size=18, color="black")) +
  geom_hline(yintercept=0, linetype="dashed", color="grey")

Cluster10Plot <- expressionProfile_data %>%
  select(group, cluster10) %>% 
  ggplot(aes(x=group, y=cluster10, fill=group)) +
  geom_boxplot(width=.5, outlier.shape=NA, position=position_dodge(width=0.5), alpha=0.7) +
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red") + 
  geom_point(pch=21, size=5, position=position_dodge(width=1)) +
  scale_fill_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) +
  scale_color_manual(values=c("1"="deeppink1", "2"="darkturquoise", "3"="lightpink3")) + 
  scale_x_discrete(labels=c("1"="swimming", "2"="round", "3"="spat")) +
  xlab("Group") + # Axis titles
  ylab("Mean Module Eigenegene") +
  ggtitle("Cluster 10") +
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
ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/expression_eigengene_Profiles_withMean.pdf", expressionProfiles, height = 25, width = 28, units = "in")

#  Gene relationship to trait and important modules: Gene Significance and Module Membership

#We quantify associations of individual genes with life stage by defining Gene Significance GS as the absolute value of the correlation between the gene and the time_point. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 

#Define variable weight containing the weight column of datTrait

treatmentinfo <- filter(treatmentinfo, sampleID %in% rownames(datExpr))
datExpr <- datExpr[treatmentinfo$sampleID,]
group <- data.frame(group = as.numeric(as.factor(treatmentinfo$timepoint)))
MEs <- MEs[treatmentinfo$sampleID,]

# Creating the group data frame for P. acuta so we can get GS per lifestage
group_p_acuta <- data.frame(
  I = c(1, 1, 1, 0, 0, 0, 0, 0, 0),
  II = c(0, 0, 0, 1, 1, 1, 0, 0, 0),
  III = c(0, 0, 0, 0, 0, 0, 1, 1, 1),
  row.names = c("SRR3051863", "SRR3051864", "SRR3051865", "SRR3051866", "SRR3051867", "SRR3051868", "SRR3051869", "SRR3051870", "SRR3051871")
)

# Display the data frame
print(group_p_acuta)
dim(group_p_acuta)
# [1] 9 3

# Colors of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Use group_p_acuta for geneTraitSignificance
geneTraitSignificance_p_acuta = as.data.frame(cor(datExpr, group_p_acuta, use = "p"))
GSPvalue_p_acuta = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_p_acuta), nSamples))

names(geneTraitSignificance_p_acuta) = paste("GS.", names(group_p_acuta), sep="")
names(GSPvalue_p_acuta) = paste("p.GS.", names(group_p_acuta), sep="")

#Use group for geneTraitSignificance for groups
geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(group), sep="")
names(GSPvalue) = paste("p.GS.", names(group), sep="")

### Summary output of network analysis results
library(readr)

annot <- read_delim(
  "~/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/Pocillopora_acuta_HIv2.genes.EggNog_results.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
annot <- rename(annot, "gene_id"=`#query`)
probes = names(datExpr)
probes2annot = match(probes, annot$gene_id)
nrow(annot)#16784
head(annot)
#A tibble: 6 x 21
#gene_id     seed_ortholog    evalue score eggNOG_OGs max_annot_lvl COG_category Description
#<chr>       <chr>             <dbl> <dbl> <chr>      <chr>         <chr>        <chr>      
#  1 Pocillopor~ 45351.EDO487~ 2.16e-120 364   COG0620@1~ 33208|Metazoa E            Cobalamin-~
#  2 Pocillopor~ 45351.EDO386~ 3.18e-123 355   COG0450@1~ 33208|Metazoa O            negative r~
#  3 Pocillopor~ 192875.XP_00~ 1.70e-183 526   COG3239@1~ 33154|Opisth~ I            Fatty acid~
#  4 Pocillopor~ 45351.EDO287~ 2.94e- 48 172   2ED36@1|r~ 33208|Metazoa -            -          
#  5 Pocillopor~ 10224.XP_006~ 3.19e- 20  92.4 COG2801@1~ 33208|Metazoa OU           K02A2.6-li~
#  6 Pocillopor~ 106582.XP_00~ 3.7 e- 14  68.2 2CSTD@1|r~ 33208|Metazoa S            Phosphatid~
sum(is.na(probes2annot))#7808, this number must be 0.
# Sums to annot was 7808 because the eggnog files don't have annotations for all the genes,
# To sort this I created a new dataframe including also the genes without annotation
# with empty annotation values. 

#Troubleshooting steps

length(probes)   # Number of probes 20891
length(annot$gene_id)   # Number of gene IDs 16784 

# Prepare the full list of probes (using original case)
probes <- names(datExpr)
all_probes <- data.frame(gene_id = probes, stringsAsFactors = FALSE)  # Rename to match annotation

# Ensure column names are consistent
annot$gene_id <- trimws(annot$gene_id)
all_probes$gene_id <- trimws(all_probes$gene_id)

# Merge the dataframes
annotated_probes <- merge(all_probes, annot, by = "gene_id", all.x = TRUE)

# Check the number of unmatched probes
missing_annot_count <- sum(is.na(annotated_probes$Description))
cat("Number of probes without annotation:", missing_annot_count, "\n")#7808

# View the head of the merged dataframe
head(annotated_probes)
#                                  gene_id         seed_ortholog    evalue score
#1 Pocillopora_acuta_HIv2___RNAseq.10002_t        45351.EDO27354  2.41e-93   317
#2 Pocillopora_acuta_HIv2___RNAseq.10010_t   6087.XP_002166004.2  1.28e-38   164
#3  Pocillopora_acuta_HIv2___RNAseq.1016_t                  <NA>        NA    NA
#4 Pocillopora_acuta_HIv2___RNAseq.10171_t 106582.XP_004539380.1  1.68e-57   225
#5 Pocillopora_acuta_HIv2___RNAseq.10263_t 106582.XP_004546958.1 3.40e-183   571
#6 Pocillopora_acuta_HIv2___RNAseq.10431_t  10224.XP_002733515.1 1.46e-118   374

# write csv for annotated probes
write.csv(annotated_probes, file = "Annotated_non_annotated_genes.csv")

probes2annot = match(probes,annotated_probes$gene_id)
nrow(annotated_probes)#20891
sum(is.na(probes2annot))#0 

#Since now the sum is 0, I can continue to the following steps. 

# Create the starting data frame
geneInfo0 = data.frame(
  gene_id = annotated_probes$gene_id,
  Accession = annot$seed_ortholog[probes2annot],
  Score = annot$score[probes2annot],
  eValue = annot$evalue[probes2annot],
  Description = annot$Description[probes2annot],
  KEGG = annot$KEGG_ko[probes2annot],
  Annotation.GO.ID = annot$GOs[probes2annot],
  moduleColor = dynamicColors,
  geneTraitSignificance,
  GSPvalue,
  geneTraitSignificance_p_acuta,
  GSPvalue_p_acuta
)

# Order modules by their significance for time_point
modOrder = order(-abs(cor(MEs, group, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.group))
geneInfo = geneInfo0[geneOrder, ]
head(geneInfo)
#gene_id.gene_id
#Pocillopora_acuta_HIv2___TS.g12975.t1         Pocillopora_acuta_HIv2___TS.g24751.t1
#Pocillopora_acuta_HIv2___TS.g21138.t3     Pocillopora_acuta_HIv2___RNAseq.g19228.t1
#Pocillopora_acuta_HIv2___RNAseq.g18587.t1  Pocillopora_acuta_HIv2___RNAseq.g1930.t1
#Pocillopora_acuta_HIv2___RNAseq.g17030.t1     Pocillopora_acuta_HIv2___TS.g11376.t1
#Pocillopora_acuta_HIv2___TS.g26697.t1a     Pocillopora_acuta_HIv2___RNAseq.g5083.t1
#Pocillopora_acuta_HIv2___RNAseq.g29921.t1 Pocillopora_acuta_HIv2___RNAseq.g17537.t1
#gene_id.seed_ortholog gene_id.evalue gene_id.score
#Pocillopora_acuta_HIv2___TS.g12975.t1                      <NA>             NA            NA
#Pocillopora_acuta_HIv2___TS.g21138.t3                      <NA>             NA            NA

#Add module cluster information
geneInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")
dim(geneInfo)
#[1] 20891    67
head(geneInfo)
#looks good

write.csv(geneInfo, file = "Gene_info_P.acuta_final.csv")


#The rest of the lines are older versions of this code I'm just keeping in case
#there's need to revise something later. 

#GS & MM old versions

treatmentinfo <- filter(treatmentinfo, sampleID %in% rownames(datExpr))
datExpr <- datExpr[treatmentinfo$sampleID,]
group <- data.frame(group = as.numeric(as.factor(treatmentinfo$timepoint)))
MEs <- MEs[treatmentinfo$sampleID,]

#Colors of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(group), sep="")
names(GSPvalue) = paste("p.GS.", names(group), sep="")

### Summary output of network analysis results

#### Make a dataframe that connects traits, genes, and gene annotation
library(readr)

annot <- read_delim(
  "~/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_gff/Pocillopora_acuta_HIv2.genes.EggNog_results.txt",
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

annot <- rename(annot, "gene_id"=`#query`)
probes = names(datExpr)
probes2annot = match(probes, annot$gene_id)
nrow(annot)#16784
head(annot)#
sum(is.na(probes2annot))#7808, time for troubleshooting madame

#troubleshooting steps
#length(probes)   # Number of probes 20891
#length(annot$gene_id)   # Number of gene IDs 16784

#Sums to annot was 7808 because the eggnog files don't have annotations for all the genes,
#To sort this I created a new dataframe including also the genes without annotation
#with empty annotation values. 

# Rename the column to "gene_id" for consistency
annot <- rename(annot, gene_id = `#query`)

# Prepare the full list of probes (using original case)
probes <- names(datExpr)
all_probes <- data.frame(gene_id = probes, stringsAsFactors = FALSE)  # Rename to match annotation

# Ensure column names are consistent
annot$gene_id <- trimws(annot$gene_id)
all_probes$gene_id <- trimws(all_probes$gene_id)

# Merge the dataframes
annotated_probes <- merge(all_probes, annot, by = "gene_id", all.x = TRUE)

# Check the number of unmatched probes
missing_annot_count <- sum(is.na(annotated_probes$Description))
cat("Number of probes without annotation:", missing_annot_count, "\n")#7808

# View the head of the merged dataframe
head(annotated_probes)
#write csv for annotated probes
write.csv(annotated_probes, file = "Annotated_non_annotated_genes.csv")


probes2annot = match(probes,annotated_probes$gene_id)
nrow(annotated_probes)#20891
sum(is.na(probes2annot))#0 

#Create the starting data frame

geneInfo0 = data.frame(
  gene_id = annotated_probes,
  Accession = annot$seed_ortholog[probes2annot],
  Score = annot$score[probes2annot],
  eValue = annot$evalue[probes2annot],
  Description = annot$Description[probes2annot],
  KEGG = annot$KEGG_ko[probes2annot],
  Annotation.GO.ID = annot$GOs[probes2annot],
  moduleColor = dynamicColors,
  geneTraitSignificance,
  GSPvalue
)

#Order modules by their significance for time_point
modOrder = order(-abs(cor(MEs, group, use = "p")))

#Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

#Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.group));
geneInfo = geneInfo0[geneOrder, ]
head(geneInfo)

#See and save geneInfo as a CSV
dim(geneInfo)
#20891 60
head(geneInfo)
#looking good

#Add module cluster information
geneInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")
dim(geneInfo)
head(geneInfo)

write.csv(geneInfo, file = "Gene_info_P.acuta_nolifestage.csv")
