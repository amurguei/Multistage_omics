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

#Sorting issues with WGCNA
#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("impute")
#BiocManager::install("preprocessCore")
#install.packages("WGCNA")


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

setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff")


#Data--------------------------------------------------------------------------------------
treatmentinfo <- read_csv("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/5-Mcap-SampleInfo.csv")

gcount <- as.data.frame(read.csv("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/Mcap_transcript_count_matrix.csv", row.names="gene_id"), colClasses = double, header=TRUE)


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
#any differences in gene expression across samples attributed to depth.

#Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = treatmentinfo,
                               design = ~timepoint)

# Log-transform the count data using a variance stabilizing transforamtion (vst). 

SF.gdds <- estimateSizeFactors( gdds ) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than for to use vst
print(sizeFactors(SF.gdds)) #View size factors
#AH1       AH2       AH3       AH4       AH5       AH6       AH7       AH8       AH9 
#1.0805764 1.3101459 1.3130424 1.2556892 1.1001138 1.3045822 0.5273869 0.7045234 0.9337767  

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

#PCA-----------------------------------------------------------------------------

gPCAdata <- plotPCA(gvst, intgroup = c("timepoint"), returnData=TRUE, ntop=1000) #use ntop to specify all genes

percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data

#allgenesfilt_PCA <- ggplot(gPCAdata, aes(PC1, PC2, shape=timepoint)) + 
#geom_point(size=3) +
#xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#scale_shape_manual(values = c("I"=1, "II"=2, "III"=3)) +
#xlim(-40,40)+ 
#ylim(-40,40)+
#coord_fixed()+
#theme_bw() + #Set background color
#theme(panel.border = element_blank(), # Set border
#panel.grid.major = element_blank(), #Set major gridlines 
#panel.grid.minor = element_blank(), #Set minor gridlines
#axis.line = element_line(colour = "black"), #Set axes color
#plot.background=element_blank()) # + #Set the plot background
#theme(legend.position = ("none")) #set title attributes
#allgenesfilt_PCA
#ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/PCA_timepoint.png", allgenesfilt_PCA, width=11, height=8)


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

#another attempt

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

ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/PCA_timepointntop=1000.png", allgenesfilt_PCA_visual, width = 11, height = 8)


# Cluster the samples to look for obvious outliers
#Look for outliers by examining the sample tree:

sampleTree = hclust(dist(datExpr), method = "average")

# Plot the sample tree
pdf(paste0('sampleTree','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

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

#Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
#Scale-free topology fit index as a function of the soft-thresholding power
pdf(paste0('network','.pdf'))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# # # this line corresponds to using an R^2 cut-off
abline(h=0.8,col="red")

# # # Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#The lowest scale-free topology fit index R^2 recommended by Langfelder and Horvath is 0.8. 
#From the graph, it appears that our soft thresholding power is 15 because it is the lowest 
#power before the R^2=0.8 threshold that maximizes with model fit (s17 is the number right above the red line).

### Network construction and module detection:
# Co-expression adjacency and topological overlap matrix similarity
# Co-expression similarity and adjacency, using the soft thresholding power 15 and translate the adjacency into topological overlap matrix to calculate 
# the corresponding dissimilarity. I will use a signed network because we have a relatively high softPower, according 
# to >12 (https://peterlangfelder.com/2018/11/25/__trashed/). 
# Moreover, in expression data where you are interested in when expression on one gene increases or decreases with expression level of another you would use a signed network (when you are interested in the direction of change, correlation and anti-correlation, you use a signed network).

options(stringsAsFactors = FALSE)
enableWGCNAThreads() #Allow multi-threading within WGCNA

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

softPower=15 #Set softPower to 15
adjacency=adjacency(datExpr, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(adjacency, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/adj.RData")
adjacency <- load(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/adj.RData")
save(dissTOM, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/disstomTOM.RData") 
dissTOM <- load (file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/disstomTOM.RData")
# Clustering using TOM
#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")
save(geneTree, file = "geneTree.RData")
geneTree <- load (file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/geneTree.RData")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/dissTOMClusteringfixedgff.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()

# Module identification is cutting the branches off the tree in the dendrogram above. We want large modules, so we set the minimum module size 
# relatively high (minimum size = 30).

minModuleSize = 30 #default value used most often
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

table(dynamicMods) #list modules and respective sizes
#1     2     3     4     5     6     7     8     9    10    11    12    13 
#11757 11257   714   506   238   235   233   197   175   153   145    76    55 
save(dynamicMods, geneTree, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/dyMod_geneTree.RData")

dyMod_geneTree <- load(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/dyMod_geneTree.RData")

dyMod_geneTree
#[1] "dynamicMods" "geneTree"   

# Plot the module assignment under the gene dendrogram
dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
#dynamicColors
#black        blue       brown       green greenyellow     magenta        pink      purple         red      salmon         tan   turquoise 
#233       11257         714         238         145         175         197         153         235          55          76       11757 
#yellow 
#506 
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/dissTOMColorClusteringgff.pdf", width=20, height=20)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar or choose not to merge
# Plot module similarity based on eigengene value

#Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower = 15)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster again and plot the results
METree = flashClust(as.dist(MEDiss), method = "average")

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/eigengeneClustering1.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >85% eigengene similarity (most studies use80-90% similarity)

MEDissThres= 0.15 #merge modules that are 85% similar

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/eigengeneClustering2V2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
#mergeCloseModules: Merging modules whose distance is less than 0.15
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 13 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 12 module eigengenes in given set.
#Calculating new MEs...
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 12 module eigengenes in given set.

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/Mcap_mergedClusters.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

#Save new colors

moduleColors = mergedColors # Rename to moduleColors
colorOrder = c("grey", standardColors(50)); # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
ncol(MEs) 
#[1] 12

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/eigengeneClustering3.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# Relating modules to group, quantifying module–trait associations
#Prepare trait data. Data has to be numeric, so I replaced the group for numeric values
#Check order with Treatment info

allTraits <- names(treatmentinfo$timepoint)
allTraits$larvae_released <- c(1,1,1,0,0,0,0,0,0)
allTraits$lavae_compressed <- c(0,0,0,1,1,1,0,0,0)
allTraits$spat <- c(0,0,0,0,0,0,1,1,1)


datTraits <- as.data.frame(allTraits)
dim(datTraits)
#[1] 9 3

#rownames(datTraits) <- timepoint$sample_id
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
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=15)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)
#[1] "MEblue"        "MEgreenyellow" "MEyellow"      "MEbrown"       "MEpink"        "MEturquoise"   "MEred"         "MEblack"       "MEgreen"      
#[10] "MEsalmon"      "MEpurple"      "MEtan"      


#Correlations of traits and eigengenes

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))


moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/Life-stage clustering based on module-trait correlation.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()


#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)
#[,1]           [,2]            [,3]           
#[1,] "0.49\n(0.2)"  "0.51\n(0.2)"   "-1\n(1e-08)"  
#[2,] "0.76\n(0.02)" "-0.67\n(0.05)" "-0.093\n(0.8)"
#[3,] "0.12\n(0.7)"  "0.00058\n(1)"  "-0.13\n(0.7)" 
#[4,] "-0.18\n(0.7)" "-0.38\n(0.3)"  "0.55\n(0.1)"  
#[5,] "0.096\n(0.8)" "0.11\n(0.8)"   "-0.21\n(0.6)" 
#[6,] "-0.49\n(0.2)" "-0.5\n(0.2)"   "0.99\n(8e-08)"

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/Module-trait-relationships.pdf")
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

pdf(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/Module-trait-relationship-heatmap3.pdf", height = 11.5, width = 8)
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
MEcluster1 <- data.frame(moduleColor = c("blue"), moduleCluster = c(1))
MEcluster2 <- data.frame(moduleColor = c("greenyellow"), moduleCluster = c(2))
MEcluster3 <- data.frame(moduleColor = c("yellow"), moduleCluster = c(3))
MEcluster4 <- data.frame(moduleColor = c("brown"), moduleCluster = c(4))
MEcluster5 <- data.frame(moduleColor = c("pink"), moduleCluster = c(5))
MEcluster6 <- data.frame(moduleColor = c("turquoise"), moduleCluster = c(6)) 
MEcluster7 <- data.frame(moduleColor = c("red"), moduleCluster = c(7)) 
MEcluster8 <- data.frame(moduleColor = c("black","green"), moduleCluster = c(8)) 
MEcluster9 <- data.frame(moduleColor = c("salmon"), moduleCluster = c(9)) 
MEcluster10 <- data.frame(moduleColor = c("purple", "tan"), moduleCluster = c(10)) 

moduleCluster = bind_rows(MEcluster1, MEcluster2, MEcluster3, MEcluster4, MEcluster5, MEcluster6, MEcluster7, MEcluster8, MEcluster9, MEcluster10)

head(moduleCluster) 

#  moduleColor moduleCluster
#moduleColor moduleCluster
#1        blue             1
#2 greenyellow             2
#3      yellow             3
#4       brown             4
#5        pink             5
#6   turquoise             6

head(MEs)
#MEblue MEgreenyellow     MEyellow      MEbrown       MEpink MEturquoise
#AH1 0.2347669     0.5357287  0.099081038 -0.002959745  0.105036951  -0.2373928
#AH2 0.2271438     0.3601512  0.064788493 -0.083632463  0.035507273  -0.2325282
#AH3 0.2313268     0.1758096  0.012340035 -0.162665612 -0.005171406  -0.2237745
#AH4 0.2434107    -0.3207421 -0.023175105 -0.237635473 -0.007581367  -0.2446720
#AH5 0.2386577    -0.3275927  0.022223785 -0.131134356  0.092289145  -0.2361615
#AH6 0.2332680    -0.2923570  0.001771422 -0.162599680  0.070302880  -0.2301644
#MEred     MEblack       MEgreen   MEsalmon   MEpurple      MEtan
#AH1 -0.3526648 -0.05598583 -0.0002329925 -0.3610804 -0.4907164 -0.3412926
#AH2 -0.2398858  0.02882405  0.1028206114 -0.2264035 -0.4728992 -0.3064592
#AH3 -0.1137728  0.08015971  0.2097099756 -0.1842593 -0.4024340 -0.3433712
#AH4  0.2268384 -0.05984661  0.1616428387  0.3091669  0.2945990  0.4622720
#AH5  0.2258356 -0.18663474  0.0869670891  0.4387961  0.3139318  0.4094073
#AH6  0.2221884 -0.09199443  0.1321559244  0.2378661  0.2825892  0.4659501

names(MEs)
#[1] "MEblue"        "MEgreenyellow" "MEyellow"      "MEbrown"       "MEpink"       
#[6] "MEturquoise"   "MEred"         "MEblack"       "MEgreen"       "MEsalmon"     
#[11] "MEpurple"      "MEtan"

Strader_MEs <- MEs
Strader_MEs$timepoint <- treatmentinfo$timepoint
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)
#MEblue MEgreenyellow     MEyellow      MEbrown       MEpink MEturquoise      MEred     MEblack       MEgreen   MEsalmon   MEpurple      MEtan
#AH1 0.2347669     0.5357287  0.099081038 -0.002959745  0.105036951  -0.2373928 -0.3526648 -0.05598583 -0.0002329925 -0.3610804 -0.4907164 -0.3412926
#AH2 0.2271438     0.3601512  0.064788493 -0.083632463  0.035507273  -0.2325282 -0.2398858  0.02882405  0.1028206114 -0.2264035 -0.4728992 -0.3064592
#AH3 0.2313268     0.1758096  0.012340035 -0.162665612 -0.005171406  -0.2237745 -0.1137728  0.08015971  0.2097099756 -0.1842593 -0.4024340 -0.3433712
#AH4 0.2434107    -0.3207421 -0.023175105 -0.237635473 -0.007581367  -0.2446720  0.2268384 -0.05984661  0.1616428387  0.3091669  0.2945990  0.4622720
#AH5 0.2386577    -0.3275927  0.022223785 -0.131134356  0.092289145  -0.2361615  0.2258356 -0.18663474  0.0869670891  0.4387961  0.3139318  0.4094073
#AH6 0.2332680    -0.2923570  0.001771422 -0.162599680  0.070302880  -0.2301644  0.2221884 -0.09199443  0.1321559244  0.2378661  0.2825892  0.4659501
#timepoint sample_id
#AH1         I       AH1
#AH2         I       AH2
#AH3         I       AH3
#AH4        II       AH4
#AH5        II       AH5
#AH6        II       AH6

head(Strader_MEs$sample_id)
#[1] "AH1" "AH2" "AH3" "AH4" "AH5" "AH6"

head(Strader_MEs$timepoint)
#[1] I  I  I  II II II
#Levels: I II III

# Calculate 10 over-arching expression patterns using mean eigengene for each module in a cluster. Use MEcluster as models for first lines
# Create a column to the Strader_MEs data frame containing age-depth groups

C1_Strader_MEs <- select(Strader_MEs, MEblue)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)
C2_Strader_MEs <- select(Strader_MEs, MEgreenyellow)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)
C3_Strader_MEs <- select(Strader_MEs, MEyellow)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)
C4_Strader_MEs <- select(Strader_MEs, MEbrown)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)
C5_Strader_MEs <- select(Strader_MEs, MEpink)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)
C6_Strader_MEs <- select(Strader_MEs,MEturquoise)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)
C7_Strader_MEs <- select(Strader_MEs, MEred)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)
C8_Strader_MEs <- select(Strader_MEs, MEblack:MEgreen)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)
C9_Strader_MEs <- select(Strader_MEs, MEsalmon)
C9_Strader_MEs$Mean <- rowMeans(C9_Strader_MEs)
C10_Strader_MEs <- select(Strader_MEs, MEpurple:MEtan)
C10_Strader_MEs$Mean <- rowMeans(C10_Strader_MEs)

expressionProfile_data <- as.data.frame(cbind(group = Strader_MEs$timepoint, cluster1= C1_Strader_MEs$Mean, cluster2 = C2_Strader_MEs$Mean, 
                                              cluster3 = C3_Strader_MEs$Mean, cluster4 = C4_Strader_MEs$Mean,cluster5 = C5_Strader_MEs$Mean,cluster6 = C6_Strader_MEs$Mean,
                                              cluster7 = C7_Strader_MEs$Mean,cluster8 = C8_Strader_MEs$Mean,cluster9 = C9_Strader_MEs$Mean,cluster10 = C10_Strader_MEs$Mean))
head(expressionProfile_data)

#cluster1   cluster2     cluster3     cluster4     cluster5   cluster6   cluster7
#1 0.2347669  0.5357287  0.099081038 -0.002959745  0.105036951 -0.2373928 -0.3526648
#2 0.2271438  0.3601512  0.064788493 -0.083632463  0.035507273 -0.2325282 -0.2398858
#3 0.2313268  0.1758096  0.012340035 -0.162665612 -0.005171406 -0.2237745 -0.1137728
#4 0.2434107 -0.3207421 -0.023175105 -0.237635473 -0.007581367 -0.2446720  0.2268384
#5 0.2386577 -0.3275927  0.022223785 -0.131134356  0.092289145 -0.2361615  0.2258356
#6 0.2332680 -0.2923570  0.001771422 -0.162599680  0.070302880 -0.2301644  0.2221884
#cluster8   cluster9  cluster10
#1 -0.02810941 -0.3610804 -0.4160045
#2  0.06582233 -0.2264035 -0.3896792
#3  0.14493484 -0.1842593 -0.3729026
#4  0.05089811  0.3091669  0.3784355
#5 -0.04983382  0.4387961  0.3616695
#6  0.02008075  0.2378661  0.3742697

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

expressionProfile_data$group
#> expressionProfile_data$group
#[1] "III" "III" "III" "II"  "II"  "II"  "I"   "I"   "I" 

#cols.num <- c(2:11)
#expressionProfile_data[cols.num] <- sapply(expressionProfile_data[cols.num],as.numeric)

expressionProfile_data$group <- factor(expressionProfile_data$group)

sapply(expressionProfile_data, class)
# [1] group  cluster1  cluster2  cluster3  cluster4  cluster5  cluster6  cluster7  cluster8  cluster9 cluster10 
#"factor" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" 

dim(expressionProfile_data)
## [1]  9 11

head(expressionProfile_data)
# head(expressionProfile_data)
#group  cluster1   cluster2     cluster3     cluster4     cluster5   cluster6   cluster7    cluster8   cluster9  cluster10
#1     1 0.2347669  0.5357287  0.099081038 -0.002959745  0.105036951 -0.2373928 -0.3526648 -0.02810941 -0.3610804 -0.4160045
#2     1 0.2271438  0.3601512  0.064788493 -0.083632463  0.035507273 -0.2325282 -0.2398858  0.06582233 -0.2264035 -0.3896792
#3     1 0.2313268  0.1758096  0.012340035 -0.162665612 -0.005171406 -0.2237745 -0.1137728  0.14493484 -0.1842593 -0.3729026
#4     2 0.2434107 -0.3207421 -0.023175105 -0.237635473 -0.007581367 -0.2446720  0.2268384  0.05089811  0.3091669  0.3784355
#5     2 0.2386577 -0.3275927  0.022223785 -0.131134356  0.092289145 -0.2361615  0.2258356 -0.04983382  0.4387961  0.3616695
#6     2 0.2332680 -0.2923570  0.001771422 -0.162599680  0.070302880 -0.2301644  0.2221884  0.02008075  0.2378661  0.3742697
 


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
ggsave("E:\\Users\\amurgueitio\\Documents\\Multistage_omics\\R scripts\\M. capitata\\New_genome\\fixed_gff\\expression_eigengene_Profiles_withMean.pdf", expressionProfiles, height = 25, width = 28, units = "in")

#  Gene relationship to trait and important modules: Gene Significance and Module Membership

#We quantify associations of individual genes with life stage by defining Gene Significance GS as the absolute value of the correlation between the gene and the time_point. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 

#Define variable weight containing the weight column of datTrait

treatmentinfo <- filter(treatmentinfo, sampleID %in% rownames(datExpr))
datExpr <- datExpr[treatmentinfo$sampleID,]
group <- data.frame(group = as.numeric(as.factor(treatmentinfo$timepoint)))
MEs <- MEs[treatmentinfo$sampleID,]

# Creating the group data frame for P. acuta so we can get GS per lifestage
# Create the data
group_m_capitata <- data.frame(
  I = c(1, 1, 1, 0, 0, 0, 0, 0, 0),
  II = c(0, 0, 0, 1, 1, 1, 0, 0, 0),
  III = c(0, 0, 0, 0, 0, 0, 1, 1, 1),
  row.names = c("AH1", "AH2", "AH3", "AH4", "AH5", "AH6", "AH7", "AH8", "AH9")
)
# Display the data frame
print(group_m_capitata)
#   I II III
#AH1 1  0   0
#AH2 1  0   0
#AH3 1  0   0
#AH4 0  1   0
#AH5 0  1   0
#AH6 0  1   0
#AH7 0  0   1
#AH8 0  0   1
#AH9 0  0   1

dim(group_m_capitata)
# [1] 9 3

# Colors of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Use group_p_acuta for geneTraitSignificance
geneTraitSignificance_m_capitata = as.data.frame(cor(datExpr, group_m_capitata, use = "p"))
GSPvalue_m_capitata = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_m_capitata), nSamples))

names(geneTraitSignificance_m_capitata) = paste("GS.", names(group_m_capitata), sep="")
names(GSPvalue_m_capitata) = paste("p.GS.", names(group_m_capitata), sep="")

#Use group for geneTraitSignificance for groups
geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(group), sep="")
names(GSPvalue) = paste("p.GS.", names(group), sep="")

### Summary output of network analysis results
library(readr)

annot <- read_delim("Montipora_capitata_HIv3.genes.EggNog_results.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
#24074 genes are annotated 

annot <- rename(annot, "gene_id"=`#query`)
probes = names(datExpr)
probes2annot = match(probes, annot$gene_id)
nrow(annot)#24074
head(annot)
sum(is.na(probes2annot))#10383, this number must be 0.
# Sums to annot was 10383 because the eggnog files don't have annotations for all the genes,
# To sort this I created a new dataframe including also the genes without annotation
# with empty annotation values. 

#Troubleshooting steps

length(probes)   # Number of probes 25741
length(annot$gene_id)   # Number of gene IDs 24072 

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
cat("Number of probes without annotation:", missing_annot_count, "\n")#10383

# View the head of the merged dataframe
head(annotated_probes)
#gene_id      seed_ortholog    evalue score
#1 Montipora_capitata_HIv3___RNAseq.10136_t               <NA>        NA    NA
#2 Montipora_capitata_HIv3___RNAseq.10187_t     45351.EDO41826 1.73e-147   441
#3 Montipora_capitata_HIv3___RNAseq.10207_t 7668.SPU_012679-tr  0.00e+00   914
#4 Montipora_capitata_HIv3___RNAseq.10214_t     45351.EDO43644  4.92e-61   205
#5 Montipora_capitata_HIv3___RNAseq.10304_t     45351.EDO43884  3.51e-37   138
#6 Montipora_capitata_HIv3___RNAseq.10384_t     45351.EDO33972  2.80e-37   142

# write csv for annotated probes
write.csv(annotated_probes, file = "Annotated_non_annotated_genesMontipora.csv")

probes2annot = match(probes,annotated_probes$gene_id)
nrow(annotated_probes)#25741
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
  geneTraitSignificance_m_capitata,
  GSPvalue_m_capitata
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
#head(geneInfo)
#Montipora_capitata_HIv3___RNAseq.g5451.t1      Montipora_capitata_HIv3___TS.g44009.t1       45351.EDO37896  85.9
#Montipora_capitata_HIv3___RNAseq.g23296.t1 Montipora_capitata_HIv3___RNAseq.g12326.t1 27679.XP_010336303.1 671.0
#Montipora_capitata_HIv3___RNAseq.g46628.t1 Montipora_capitata_HIv3___RNAseq.g41190.t1       45351.EDO48351 225.0
#Montipora_capitata_HIv3___RNAseq.g22907.t1 Montipora_capitata_HIv3___RNAseq.g45900.t1  400682.PAC_15720747 248.0
#Montipora_capitata_HIv3___RNAseq.g12630.t1 Montipora_capitata_HIv3___RNAseq.g44118.t1  6087.XP_004211416.1  81.3
#Montipora_capitata_HIv3___RNAseq.g32207.t1 Montipora_capitata_HIv3___RNAseq.g43801.t1       45351.EDO36281 333.0

#Add module cluster information
geneInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")
dim(geneInfo)
#[1] 20891    41
head(geneInfo)
#looks good

write.csv(geneInfo, file = "Gene_info_M.capitata_final.csv")


#The rest of the lines are older versions of this code I'm just keeping in case
#there's need to revise something later. 

#GS & MM old versions

#Prepare trait data. Data has to be numeric, so I will substitute time_points and type for numeric values. Not sure if the character thing is correct...

# Create the data
Mcap <- data.frame(
  I = c(1, 1, 1, 0, 0, 0, 0, 0, 0),
  II = c(0, 0, 0, 1, 1, 1, 0, 0, 0),
  III = c(0, 0, 0, 0, 0, 0, 1, 1, 1),
  row.names = c("AH1", "AH2", "AH3", "AH4", "AH5", "AH6", "AH7", "AH8", "AH9")
)

# Convert columns to characters
Mcap$I <- as.character(Mcap$I)
Mcap$II <- as.character(Mcap$II)
Mcap$III <- as.character(Mcap$III)

# Display the data frame
print(Mcap)


#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Correlations of traits with eigengenes
moduleTraitCor = cor(MEs, Mcap, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

##  Gene relationship to trait and important modules: Gene Significance and Module Membership

#We quantify associations of individual genes with life stage by defining Gene Significance GS as the absolute value of the correlation between the gene and the time_point. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 

#Define variable weight containing the weight column of datTrait

treatmentinfo <- filter(treatmentinfo, sampleID %in% rownames(datExpr))
datExpr <- datExpr[treatmentinfo$sampleID,]
group <- data.frame(group=as.numeric(as.factor(treatmentinfo$timepoint)))
MEs <- MEs[treatmentinfo$sampleID,]

#Colours of the modules

modNames = substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, group, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- paste("GS.", names(group), sep = "")
names(GSPvalue) <- paste("p.GS.", names(group), sep = "")

# Summary output of network analysis results

# Make a dataframe that connects traits, genes, and gene annotation

#Import annotation file

library(readr)

annot <- read_delim("Montipora_capitata_HIv3.genes.EggNog_results.txt",
                    delim = "\t",
                    escape_double = FALSE,
                    trim_ws = TRUE)


annot <- rename(annot, "gene_id"=`#query`)

#Match up genes in datExpr to annotation file

names(annot)
#[1] "gene_id"        "seed_ortholog"  "evalue"         "score"         
#[5] "eggNOG_OGs"     "max_annot_lvl"  "COG_category"   "Description"   
#[9] "Preferred_name" "GOs"            "EC"             "KEGG_ko"       
#[13] "KEGG_Pathway"   "KEGG_Module"    "KEGG_Reaction"  "KEGG_rclass"   
#[17] "BRITE"          "KEGG_TC"        "CAZy"           "BiGG_Reaction" 
#[21] "PFAMs"   

#This is not working
probes = names(datExpr)
probes2annot = match(probes, annot$gene_id)
nrow(annot)
#24072
head(annot)

# The following is the number of probes without annotation... Should return 0.
sum(is.na(probes2annot))
#10385

#This Chat GPT solution is working???
# Ensure both probes and Annot$gene_id are character vectors
probes <- as.character(probes)
annot$gene_id <- as.character(annot$gene_id)

# Merge based on gene_id
merged_data <- merge(data.frame(probes), Annot, by.x = "probes", by.y = "gene_id", all.x = TRUE)

# Check for probes without annotation
sum(is.na(merged_data$gene_id))
#0


#Create the starting dataset

# Create the starting data frame
geneInfo0 <- data.frame(
  gene_id = probes,
  Accession = merged_data$seed_ortholog,
  Score = merged_data$score,
  eValue = merged_data$evalue,
  Description = merged_data$Description,
  KEGG = merged_data$KEGG_ko,
  Annotation.GO.ID = merged_data$GOs,
  moduleColor = dynamicColors,  # Replace dynamicColors with your actual module color data
  geneTraitSignificance,
  GSPvalue
)

# Order modules by their significance for time_point
modOrder <- order(-abs(cor(MEs, group, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(
    geneInfo0,
    geneModuleMembership[, modOrder[mod]],
    MMPvalue[, modOrder[mod]]
  )
  names(geneInfo0) <- c(
    oldNames,
    paste("MM.", modNames[modOrder[mod]], sep = ""),
    paste("p.MM.", modNames[modOrder[mod]], sep = "")
  )
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.group))
geneInfo <- geneInfo0[geneOrder, ]

# View and save geneInfo as a CSV
dim(geneInfo)  # Check dimensions
#[1] 25741    34
head(geneInfo)  # View the first few rows

# Write geneInfo to a CSV file
write.csv(geneInfo, file = "geneinfo_didntwork.csv", row.names = FALSE)


#Second version (Erin's original)
geneInfo0 = data.frame(gene_id = probes,
                       Accession = annot$seed_ortholog[probes2annot],
                       Score = annot$score[probes2annot],
                       eValue = annot$evalue[probes2annot],
                       Description = annot$Description[probes2annot],
                       KEGG = annot$KEGG_ko[probes2annot],
                       Annotation.GO.ID = annot$GOs[probes2annot],
                       moduleColor = dynamicColors,
                       geneTraitSignificance,
                       GSPvalue)
#Order modules by their significance for time_point
modOrder = order(-abs(cor(MEs, group, use = "p")))

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

dim(geneInfo)
head(geneInfo)

