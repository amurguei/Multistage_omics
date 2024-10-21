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


setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata")

setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Transcriptomics/Inputs")

#Upload data--------------------------------------------------------------------------

treatmentinfo <- read_csv("5-Spis-SampleInfo.csv")

#gene count matrix
gcount <- as.data.frame(read.csv("4-Spis-GeneCountMatrix.csv", row.names="gene_id"), colClasses = double, header=TRUE)

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
nrow(gcount) #Before 25769                                        

nrow(gcount_filt) #After 15998

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

gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transforamtion to minimize effects of small counts and normalize wrt library size

#Compile WGCNA dataset--------------------------------------------------------------------------
# Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.

datExpr <- as.data.frame(t(assay(gvst))) #transpose to output to a new data frame with the column names as row names. And make all data numeric

#Check for genes and samples with too many missing values with goodSamplesGenes. There shouldn't be any because we performed pre-filtering
#gsg = goodSamplesGenes(datExpr, verbose = 3)
#gsg$allOK #Should return TRUE if not, the R chunk below will take care of flagged data
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

#sampleTree = hclust(dist(datExpr), method = "average")

#Plot the sample tree
#pdf(paste0('sampleTree','.pdf'))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
#dev.off()

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

ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/PCA_timepointntop=1000.png", allgenesfilt_PCA_visual, width = 11, height = 8)

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



#PERMANOVA
# Load necessary library
library(vegan)

# Check the data frames
str(datExpr)         # Should be a matrix or dataframe with genes as columns and samples as rows
str(treatmentinfo)   # Should have a sampleID and timepoint (life stages)

# Ensure the row names of datExpr match the sampleID in treatmentinfo
rownames(datExpr) <- treatmentinfo$sampleID

# Convert timepoint (life stage) to factor if not already
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I", "II", "III"))

# Running PERMANOVA with Bray-Curtis dissimilarity
permanova_result <- adonis2(datExpr ~ timepoint, data = treatmentinfo, method = "bray", permutations = 9999)
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = datExpr ~ timepoint, data = treatmentinfo, permutations = 999, method = "bray")
#Df  SumOfSqs      R2      F Pr(>F)   
#timepoint  2 0.0069751 0.65401 5.6708  0.003 **
#  Residual   6 0.0036900 0.34599                 
#Total      8 0.0106651 1.00000                 
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# View results
print(permanova_result)

# Load the vegan package
library(vegan)
library(cluster)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

pairwise.adonis(datExpr ~ timepoint, sim.method = "bray")

# Assuming datExpr is a data frame with sample IDs as row names
datExpr <- as.data.frame(datExpr)
datExpr$sampleID <- rownames(datExpr)

# Merge the expression data with the treatment info
combined_data <- merge(datExpr, treatmentinfo, by = "sampleID")

# Now, the first few columns will be your expression data, and the last one will be the timepoint
library(pairwiseAdonis)

# Prepare the expression data (excluding the sampleID and timepoint columns)
expr_data <- combined_data[, -which(names(combined_data) %in% c("sampleID", "timepoint"))]

# Run pairwise.adonis
result <- pairwise.adonis(expr_data, combined_data$timepoint, sim.method = "bray", perm=9999)

# View the results
print(result)

dim(expr_data) # should match with the number of samples in treatmentinfo
head(combined_data) # check combined data to see if timepoints align with samples


#Choose a set of soft-thresholding powers
powers <- c(seq(from = 1, to=19, by=2), c(21:30)) #Create a string of numbers from 1 through 10, and even numbers from 10 through 20
powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2))

#Call the network topology analysis function
sft <-pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#Plot the results
# sizeGrWindow(9, 5)
# par(mfrow = c(1, 2))
# cex1 = 0.9
# 
# # Scale-free topology fit index as a function of the soft-thresholding power
# pdf(paste0('network', '.pdf'))
# plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
#      xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
#      main = paste("Scale independence"))
# text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
#      labels = powers, cex = cex1, col = "red")
# abline(h = 0.8, col = "red")
# 
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
#      xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
# dev.off()


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
save(TOM, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/TOMspis.RData")
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(dissTOM, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/disTOMSpis.RData")

# Clustering using TOM
#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/dissTOMClusteringV2.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()

# Module identification is cutting the branches off the tree in the dendrogram above. We want large modules, so we set the minimum module size 
# relatively high (minimum size = 30).

minModuleSize = 30 #default value used most often
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods) #list modules and respective sizes
save(dynamicMods, geneTree, file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/dyMod_geneTree.RData")

dyMod_geneTree <- load(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/dyMod_geneTree.RData")
table(dynamicMods) #list modules and respective sizes
dynamicMods
#1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23 
#878 369 356 327 283 236 228 227 211 205 204 185 179 179 173 172 168 150 147 146 145 139 135 
#24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46 
#132 132 132 131 127 127 126 126 125 124 123 123 123 122 122 122 121 120 119 115 114 113 109 
#47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69 
#109 107 107 106 106 106 102 100  99  98  98  97  97  95  93  91  91  91  90  90  85  85  80 
#70  71  72  73  74  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92 
#80  80  79  79  79  78  77  76  76  75  75  74  73  72  72  71  71  71  70  70  69  69  69 
#93  94  95  96  97  98  99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 
#66  66  66  66  65  65  65  65  64  64  63  62  62  62  61  61  60  60  59  59  58  58  57 
#116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 
#57  56  55  55  53  53  53  53  52  52  52  51  50  50  49  49  48  48  47  47  46  46  46 
#139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 
#45  45  45  44  44  43  43  43  43  42  42  42  41  41  40  40  39  39  39  38  38  38  38 
#162 163 164 165 166 167 168 
#37  36  36  35  35  34  33 


dyMod_geneTree <- load(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/dyMod_geneTree.RData")

dyMod_geneTree
## [1] "dynamicMods" "geneTree"

# Plot the module assignment under the gene dendrogram
dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
# [1] dynamicColors
# dynamicColors
#aliceblue   antiquewhite1   antiquewhite2   antiquewhite4         bisque4 
#42              55              69              93             109 
#black            blue           blue1           blue2           blue3 
#228             369              37              77              48 
#blue4      blueviolet           brown          brown1          brown2 
#62              62             356              47              78 
#brown3          brown4      chocolate3      chocolate4           coral 
#36             109              43              55              70 
#coral1          coral2          coral3          coral4  cornflowerblue 
#95              91              69              53              42 
#cyan  darkgoldenrod4       darkgreen        darkgrey     darkmagenta 
#179              38             139             132             123 
#darkolivegreen darkolivegreen1 darkolivegreen2 darkolivegreen4      darkorange 
#124              48              63              79             132 
#darkorange2         darkred   darkseagreen1   darkseagreen2   darkseagreen3 
#113             145              43              56              70 
#darkseagreen4   darkslateblue   darkturquoise      darkviolet        deeppink 
#97             107             135              76              62 
#deeppink1       deeppink2       firebrick      firebrick2      firebrick3 
#47              36              38              49              64 
#firebrick4     floralwhite           green          green3          green4 
#79             114             283              43              57 
#greenyellow          grey60        honeydew       honeydew1      indianred1 
#204             168              71              97              38 
#indianred2      indianred3      indianred4           ivory   lavenderblush 
#49              64              79             115              43 
#lavenderblush1  lavenderblush2  lavenderblush3      lightblue2      lightblue3 
#57              71              98              38              50 
#lightblue4      lightcoral       lightcyan      lightcyan1      lightgreen 
#65              80             172             119             150 
#lightpink1      lightpink2      lightpink3      lightpink4   lightskyblue3 
#44              58              71              98              39 
#lightskyblue4  lightslateblue  lightsteelblue lightsteelblue1     lightyellow 
#50              65              80             120             147 
#magenta        magenta2        magenta3        magenta4          maroon 
#211              44              58              72              99 
#mediumorchid   mediumorchid4    mediumpurple   mediumpurple1   mediumpurple2 
#91              39              51              65              80 
#mediumpurple3   mediumpurple4    midnightblue       mistyrose        moccasin 
#121              69             173              53              45 
#navajowhite    navajowhite1    navajowhite2    navajowhite3          orange 
#59              72             100              42             132 
#orange4       orangered      orangered1      orangered3      orangered4 
#39              52              65              85             122 
#paleturquoise   palevioletred  palevioletred1  palevioletred2  palevioletred3 
#126              45              59              73             102 
#pink           pink2           pink3           pink4            plum 
#227              40              52              66              85 
#plum1           plum2           plum3           plum4      powderblue 
#122             107              76              61              46 
#purple         purple2             red       royalblue      royalblue2 
#205              35             236             146              33 
#royalblue3     saddlebrown          salmon         salmon1         salmon2 
#45             127             179              60              74 
#salmon4         sienna1         sienna2         sienna3         sienna4 
#106              40              52             123              66 
#skyblue        skyblue1        skyblue2        skyblue3        skyblue4 
#127              90              91             122              66 
#slateblue      slateblue1       steelblue             tan            tan2 
#53              41             126             185              34 
#tan3            tan4         thistle        thistle1        thistle2 
#46              60              75             106             106 
#thistle3        thistle4          tomato         tomato2       turquoise 
#75              61              46              35             878 
#violet           white      whitesmoke          yellow         yellow2 
#125             131              41             327              53 
#yellow3         yellow4     yellowgreen 
#66              90             123 

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/dissTOMColorClusteringV2.pdf", width=20, height=20)
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

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/eigengeneClustering1V2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >85% eigengene similarity (most studies use80-90% similarity)

MEDissThres= 0.15 #merge modules that are 85% similar

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/eigengeneClustering2V2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
# [1] mergeCloseModules: Merging modules whose distance is less than 0.15
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 168 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 96 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 83 module eigengenes in given set.
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 81 module eigengenes in given set.
#Calculating new MEs...
#multiSetMEs: Calculating module MEs.
#Working on set 1 ...
#moduleEigengenes: Calculating 81 module eigengenes in given set.
 

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/mergedClustersV2.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

#Save new colors

moduleColors = mergedColors # Rename to moduleColors
colorOrder = c("grey", standardColors(50)); # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
ncol(MEs) 
# [1] 81

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. Pistillata/eigengeneClustering3V2.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()


# Relating modules to group, quantifying module–trait associations
#Prepare trait data. Data has to be numeric, so I replaced the group for numeric values
#Check order with Treatment info

allTraits <- names(treatmentinfo$group)
allTraits$larvae_released <- c(0,0,0,0,0,0,1,1,1)
allTraits$lavae_compressed <- c(0,0,0,1,1,1,0,0,0)
allTraits$spat <- c(1,1,1,0,0,0,0,0,0)


datTraits <- as.data.frame(allTraits)
dim(datTraits)
#9, 3

rownames(datTraits) <- treatmentinfo$sample_id
print(datTraits)

#larvae_released lavae_compressed spat
#1               0                0    1
#2               0                0    1
#3               0                0    1
#4               0                1    0
#5               0                1    0
#6               0                1    0
#7               1                0    0
#8               1                0    0
#9               1                0    0

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=19)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)
#[1] "MEcornflowerblue"  "MEyellow2"         "MEblack"           "MEyellow3"        
#[5] "MEgreen3"          "MEindianred3"      "MEsalmon2"         "MEdeeppink1"      
#[9] "MEnavajowhite3"    "MElavenderblush2"  "MEmediumpurple4"   "MEdarkolivegreen4"
#[13] "MEcoral2"          "MEroyalblue2"      "MEblue4"           "MEdarkgrey"       
#[17] "MEgrey60"          "MElightcyan1"      "MEdarkseagreen2"   "MEcyan"           
#[21] "MEdarkorange"      "MEblue"            "MEindianred4"      "MEbisque4"        
#[25] "MEmediumpurple"    "MEfirebrick4"      "MEpink4"           "MEcoral"          
#[29] "MEdarkolivegreen1" "MEdarkseagreen1"   "MEsalmon4"         "MEdarkmagenta"    
#[33] "MEskyblue"         "MElightpink2"      "MElightslateblue"  "MEtan3"           
#[37] "MEbrown4"          "MEblue3"           "MEchocolate4"      "MEantiquewhite4"  
#[41] "MEmediumorchid"    "MEpink"            "MEfirebrick3"      "MEaliceblue"      
#[45] "MElavenderblush"   "MEbrown"           "MElightblue2"      "MEpalevioletred"  
#[49] "MEdeeppink2"       "MEfirebrick"       "MEhoneydew"        "MEindianred1"     
#[53] "MEwhitesmoke"      "MEbrown3"          "MElightpink1"      "MEmagenta3"       
#[57] "MEsaddlebrown"     "MEbrown1"          "MEskyblue2"        "MEantiquewhite1"  
#[61] "MEmoccasin"        "MEblue1"           "MEdarkgoldenrod4"  "MEivory"          
#[65] "MEpowderblue"      "MElightcyan"       "MEmagenta2"        "MEorangered3"     
#[69] "MEtomato2"         "MEcoral3"          "MEcoral4"          "MEdarkturquoise"  
#[73] "MElightsteelblue"  "MEbrown2"          "MEindianred2"      "MEslateblue1"     
#[77] "MEdeeppink"        "MEmediumorchid4"   "MEorangered"       "MEchocolate3"     
#[81] "MElightskyblue3" 

#Correlations of traits and eigengenes

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))


moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. Pistillata/Life-stage clustering based on module-trait correlation.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#Plot module trait correlations as a heatmap

textMatrix = paste(signif(moduleTraitCor, 2), "\n(" ,signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
head(textMatrix)

pdf(file="E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. Pistillata/Module-trait-relationships.pdf")
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

pdf(file = "E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. Pistillata/Module-trait-relationship-heatmap3.pdf", height = 11.5, width = 8)
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
MEcluster1 <- data.frame(moduleColor = c("cornflowerblue","yellow2","black","yellow3","green3","indianred3"), moduleCluster = c(1))
MEcluster2 <- data.frame(moduleColor = c("salmon2","deeppink1","navajowhite3","lavenderblush2","mediumpurple4","darkolivegreen4","coral2","royalblue2"), moduleCluster = c(2))
MEcluster3 <- data.frame(moduleColor = c("blue4","darkgrey","grey60","lightcyan1","darkseagreen2","cyan","darkorange","blue","indianred4","bisque4","mediumpurple","firebrick4","pink4"), moduleCluster = c(3))
MEcluster4 <- data.frame(moduleColor = c("coral","darkolivegreen1","darkseagreen1","salmon4","darkmagenta","skyblue","lightpink2","lightslateblue","tan3"), moduleCluster = c(4))
MEcluster5 <- data.frame(moduleColor = c("brown4","blue3","chocolate4","antiquewhite4","mediumorchid","pink","firebrick3","aliceblue","lavenderblush","brown","lightblue2","palevioletred"), moduleCluster = c(5))
MEcluster6 <- data.frame(moduleColor = c("deeppink2","firebrick","honeydew","indianred1","whitesmoke"), moduleCluster = c(6)) 
MEcluster7 <- data.frame(moduleColor = c("brown3","lightpink1","magenta3","saddlebrown","brown1","skyblue2","antiquewhite1","moccasin","blue1","darkgoldenrod4"), moduleCluster = c(7)) 
MEcluster8 <- data.frame(moduleColor = c("ivory","powderblue","lightcyan","magenta2","orangered3","tomato2","coral3","coral4","darkturquoise","lightsteelblue"), moduleCluster = c(8)) 
MEcluster9 <- data.frame(moduleColor = c("brown2","indianred2","slateblue1"), moduleCluster = c(9))
MEcluster10 <- data.frame(moduleColor = c("deeppink", "mediumorchid4","orangered","chocolate3","lightskyblue3"), moduleCluster = c(10)) 

moduleCluster = bind_rows(MEcluster1, MEcluster2, MEcluster3, MEcluster4, MEcluster5, MEcluster6, MEcluster7, MEcluster8, MEcluster9, MEcluster10)
head(moduleCluster) 
#   moduleColor moduleCluster
#1 cornflowerblue             1
#2        yellow2             1
#3          black             1
#4        yellow3             1
#5         green3             1
#6     indianred3             1

head(MEs)

#MEcornflowerblue   MEyellow2    MEblack  MEyellow3      MEgreen3 MEindianred3 MEsalmon2  MEdeeppink1 MEnavajowhite3 MElavenderblush2
#SRR14333319       -0.1539588  0.03527640  0.2393682  0.1739224  0.0001629248  -0.08626927 0.1265013  0.007753072     0.16200765        0.2699916
#SRR14333320       -0.1593908  0.11680039  0.3802084  0.1779685  0.0420750591   0.21402703 0.1477487 -0.004279204     0.25706623        0.2823233
#SRR14333321       -0.1315325  0.18434902  0.3481822  0.2427760 -0.0170750525   0.07349463 0.1332277 -0.007477901     0.28814489        0.3037046
#SRR14333322       -0.1281852 -0.05125104 -0.4038646 -0.4142822 -0.3450500283  -0.37505329 0.1328867  0.090654597     0.05868169       -0.4123771
#SRR14333323        0.1380989  0.02630253  0.1525994  0.2315720  0.5306128991   0.08192228 0.1659785  0.329909077     0.15221631        0.2660582
#SRR14333324        0.7909674  0.53604264  0.3410017  0.4578577  0.6017282607   0.77670371 0.3183431  0.472889576     0.04975688       -0.2109061
#            MEmediumpurple4 MEdarkolivegreen4   MEcoral2 MEroyalblue2    MEblue4 MEdarkgrey   MEgrey60 MElightcyan1 MEdarkseagreen2     MEcyan
#SRR14333319      0.08515150        0.04279293  0.2169719   0.07460750 0.22199562  0.1984324  0.2565821    0.2680967       0.3053653  0.3405317
#SRR14333320      0.03272885        0.09148874  0.3461497   0.14433872 0.29519651  0.2697398  0.3088139    0.3669669       0.3966324  0.4203874
#SRR14333321     -0.01986668        0.14388249  0.3783050   0.24833076 0.31456287  0.3054259  0.3722781    0.3193472       0.3889432  0.3701000
#SRR14333322     -0.26463990       -0.20124412 -0.3486717   0.04971767 0.10918099  0.1000814  0.3571311    0.4822638      -0.4093580 -0.4900082
#SRR14333323      0.28286473        0.77920968  0.4906035   0.67867684 0.06777992  0.2849949  0.1708722   -0.2961981      -0.1460393 -0.2838789
#SRR14333324     -0.25423285       -0.25930004 -0.1123693  -0.33144648 0.17926018  0.1901359 -0.3271556   -0.3515044       0.1816994 -0.1218740
#            MEdarkorange     MEblue MEindianred4    MEbisque4 MEmediumpurple MEfirebrick4    MEpink4     MEcoral MEdarkolivegreen1 MEdarkseagreen1
#SRR14333319   0.30133183  0.4141724    0.3288598  0.077000850     0.17336266   0.36556120  0.2676408  0.23887562        0.34585066       0.2345640
#SRR14333320   0.38227547  0.5022334    0.3837632  0.110240642     0.20943562   0.38004352  0.3038318  0.21697257        0.35635711       0.2212203
#SRR14333321   0.39933982  0.4828019    0.3501371  0.132470939     0.22275578   0.32450550  0.3489944  0.11047590        0.34253222       0.1561683
#SRR14333322  -0.15413973 -0.2443654   -0.4350062  0.194484516     0.02836625   0.10414580  0.3020141 -0.87109082       -0.76747218      -0.6403511
#SRR14333323  -0.05145070 -0.1931623   -0.1456246 -0.007493273    -0.63002089  -0.47645412 -0.3031070 -0.09316642       -0.06695815      -0.3037390
#SRR14333324  -0.02940979 -0.1395911    0.1812338  0.231458181     0.25294830   0.02358894  0.1845928  0.30930960       -0.20136230      -0.2877859
#             MEsalmon4 MEdarkmagenta  MEskyblue MElightpink2 MElightslateblue      MEtan3    MEbrown4     MEblue3 MEchocolate4 MEantiquewhite4  MEmoccasin
#SRR14333319  0.3845950     0.3828713  0.2069695   0.17759994        0.1866408  0.29488395 -0.04976654  0.07703347  -0.12747860      -0.3026429 -0.10693542
#SRR14333320  0.3537090     0.3922119  0.2051167   0.35007310        0.3551318  0.34674053 -0.09283600  0.15015729  -0.05823484      -0.3853304 -0.05501081
#SRR14333321  0.2947369     0.3763757  0.1460994   0.28687815        0.3304538  0.29128433 -0.14556760  0.06504874  -0.10916165      -0.4066914 -0.04417096
#SRR14333322 -0.4577735    -0.3960708 -0.4879148  -0.09894001       -0.2706137 -0.19712717 -0.23562422 -0.03493395   0.30875294       0.2932076  0.85993477
#SRR14333323 -0.2884531    -0.5676556 -0.6450686   0.04603941       -0.3248785  0.07433987 -0.52460027 -0.57282876  -0.48674560      -0.1383312 -0.38348936
#SRR14333324 -0.1766409    -0.1707170 -0.1870647  -0.43725651       -0.6732900 -0.73561406 -0.24593885 -


names(MEs)
# [1] "MEcornflowerblue"  "MEyellow2"         "MEblack"          
# [4] "MEyellow3"         "MEgreen3"          "MEindianred3"     
# [7] "MEsalmon2"         "MEdeeppink1"       "MEnavajowhite3"   
# [10] "MElavenderblush2"  "MEmediumpurple4"   "MEdarkolivegreen4"
# [13] "MEcoral2"          "MEroyalblue2"      "MEblue4"          
# [16] "MEdarkgrey"        "MEgrey60"          "MElightcyan1"     
# [19] "MEdarkseagreen2"   "MEcyan"            "MEdarkorange"     
# [22] "MEblue"            "MEindianred4"      "MEbisque4"        
# [25] "MEmediumpurple"    "MEfirebrick4"      "MEpink4"          
# [28] "MEcoral"           "MEdarkolivegreen1" "MEdarkseagreen1"  
# [31] "MEsalmon4"         "MEdarkmagenta"     "MEskyblue"        
# [34] "MElightpink2"      "MElightslateblue"  "MEtan3"           
# [37] "MEbrown4"          "MEblue3"           "MEchocolate4"     
# [40] "MEantiquewhite4"   "MEmediumorchid"    "MEpink"           
# [43] "MEfirebrick3"      "MEaliceblue"       "MElavenderblush"  
# [46] "MEbrown"           "MElightblue2"      "MEpalevioletred"  
# [49] "MEdeeppink2"       "MEfirebrick"       "MEhoneydew"       
# [52] "MEindianred1"      "MEwhitesmoke"      "MEbrown3"         
# [55] "MElightpink1"      "MEmagenta3"        "MEsaddlebrown"    
# [58] "MEbrown1"          "MEskyblue2"        "MEantiquewhite1"  
# [61] "MEmoccasin"        "MEblue1"           "MEdarkgoldenrod4" 
# [64] "MEivory"           "MEpowderblue"      "MElightcyan"      
# [67] "MEmagenta2"        "MEorangered3"      "MEtomato2"        
# [70] "MEcoral3"          "MEcoral4"          "MEdarkturquoise"  
# [73] "MElightsteelblue"  "MEbrown2"          "MEindianred2"     
# [76] "MEslateblue1"      "MEdeeppink"        "MEmediumorchid4"  
# [79] "MEorangered"       "MEchocolate3"      "MElightskyblue3"

Strader_MEs <- MEs
Strader_MEs$timepoint <- treatmentinfo$timepoint
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)

# Muting the dataset for display
# MEcornflowerblue   MEyellow2    MEblack  MEyellow3      MEgreen3 MEindianred3 MEsalmon2  MEdeeppink1 MEnavajowhite3 MElavenderblush2
# SRR14333319       -0.1539588  0.03527640  0.2393682  0.1739224  0.0001629248  -0.08626927 0.1265013  0.007753072     0.16200765        0.2699916
# SRR14333320       -0.1593908  0.11680039  0.3802084  0.1779685  0.0420750591   0.21402703 0.1477487 -0.004279204     0.25706623        0.2823233
# SRR14333321       -0.1315325  0.18434902  0.3481822  0.2427760 -0.0170750525   0.07349463 0.1332277 -0.007477901     0.28814489        0.3037046
# SRR14333322       -0.1281852 -0.05125104 -0.4038646 -0.4142822 -0.3450500283  -0.37505329 0.1328867  0.090654597     0.05868169       -0.4123771
# SRR14333323        0.1380989  0.02630253  0.1525994  0.2315720  0.5306128991   0.08192228 0.1659785  0.329909077     0.15221631        0.2660582
# SRR14333324        0.7909674  0.53604264  0.3410017  0.4578577  0.6017282607   0.77670371 0.3183431  0.472889576     0.04975688       -0.2109061
# MEmediumpurple4 MEdarkolivegreen4   MEcoral2 MEroyalblue2    MEblue4 MEdarkgrey   MEgrey60 MElightcyan1 MEdarkseagreen2     MEcyan
# SRR14333319      0.08515150        0.04279293  0.2169719   0.07460750 0.22199562  0.1984324  0.2565821    0.2680967       0.3053653  0.3405317
# SRR14333320      0.03272885        0.09148874  0.3461497   0.14433872 0.29519651  0.2697398  0.3088139    0.3669669       0.3966324  0.4203874
# SRR14333321     -0.01986668        0.14388249  0.3783050   0.24833076 0.31456287  0.3054259  0.3722781    0.3193472       0.3889432  0.3701000
# SRR14333322     -0.26463990       -0.20124412 -0.3486717   0.04971767 0.10918099  0.1000814  0.3571311    0.4822638      -0.4093580 -0.4900082
# SRR14333323      0.28286473        0.77920968  0.4906035   0.67867684 0.06777992  0.2849949  0.1708722   -0.2961981      -0.1460393 -0.2838789
# SRR14333324     -0.25423285       -0.25930004 -0.1123693  -0.33144648 0.17926018  0.1901359 -0.3271556   -0.3515044       0.1816994 -0.1218740
# MEdarkorange     MEblue MEindianred4    MEbisque4 MEmediumpurple MEfirebrick4    MEpink4     MEcoral MEdarkolivegreen1 MEdarkseagreen1
# SRR14333319   0.30133183  0.4141724    0.3288598  0.077000850     0.17336266   0.36556120  0.2676408  0.23887562        0.34585066       0.2345640
# SRR14333320   0.38227547  0.5022334    0.3837632  0.110240642     0.20943562   0.38004352  0.3038318  0.21697257        0.35635711       0.2212203
# SRR14333321   0.39933982  0.4828019    0.3501371  0.132470939     0.22275578   0.32450550  0.11047590  0.34253222       0.1561683
# SRR14333322  -0.15413973 -0.2443654   -0.4350062  0.194484516     0.02836625   0.10414580  0.3020141 -0.87109082       -0.76747218      -0.6403511
# SRR14333323  -0.05145070 -0.1931623   -0.1456246 -0.007493273    -0.63002089  -0.47645412 -0.3031070 -0.09316642       -0.06695815      -0.3037390
# SRR14333324  -0.02940979 -0.1395911    0.1812338  0.231458181     0.25294830   0.02358894  0.1845928  0.30930960       -0.20136230      -0.2877859
# MEsalmon4 MEdarkmagenta  MEskyblue MElightpink2 MElightslateblue      MEtan3    MEbrown4     MEblue3 MEchocolate4 MEantiquewhite4
# SRR14333319  0.3845950     0.3828713  0.2069695   0.17759994        0.1866408  0.29488395 -0.04976654  0.07703347  -0.12747860      -0.3026429
# SRR14333320  0.3537090     0.3922119  0.2051167   0.35007310        0.3551318  0.34674053 -0.09283600  0.15015729  -0.05823484      -0.3853304
# SRR14333321  0.2947369     0.3763757  0.1460994   0.28687815        0.3304538  0.29128433 -0.14556760  0.06504874  -0.10916165      -0.4066914
# SRR14333322 -0.4577735    -0.3960708 -0.4879148  -0.09894001       -0.2706137 -0.19712717 -0.23562422 -0.03493395   0.30875294       0.2932076
# SRR14333323 -0.2884531    -0.5676556 -0.6450686   0.04603941       -0.3248785  0.07433987 -0.52460027 -0.57282876  -0.48674560      -0.1383312
# SRR14333324 -0.1766409    -0.1707170 -0.1870647  -0.43725651       -0.6732900 -0.73561406 -0.24593885 -0.53284538  -0.52623947      -0.1377024
#MEmoccasin    MEblue1 MEdarkgoldenrod4    MEivory MEpowderblue MElightcyan  MEmagenta2 MEorangered3  MEtomato2   MEcoral3    MEcoral4 MEdarkturquoise
#SRR14333319 -0.10693542 -0.05281997       0.07302638 -0.2622696  -0.11018395 -0.36804607 -0.36340118  -0.13525207 0.03567272 -0.1865693 -0.17560233     -0.15445101
#SRR14333320 -0.05501081  0.09465338       0.06790039 -0.2626980  -0.12273246 -0.43790652 -0.39643256  -0.09393749 0.04808416 -0.2180750  0.01284896     -0.16081975
#SRR14333321 -0.04417096  0.12552565       0.13100719 -0.2057022  -0.06298302 -0.40982138 -0.37701275  -0.04949148 0.09960000 -0.1831822  0.06866431     -0.05972167
#SRR14333322  0.85993477  0.32610861       0.48109601  0.7993526   0.72572868  0.52002487  0.68651718   0.73520160 0.39720895  0.5454977  0.34396532      0.65630926
#SRR14333323 -0.38348936 -0.45139412      -0.37519268  0.2049975   0.22578654  0.27957859  0.03996922   0.09589084 0.13279525  0.1507223  0.02506655      0.53343744
#SRR14333324 -0.26647548 -0.21427254      -0.47703936 -0.3536608  -0.30872568 -0.04959461  0.06240819  -0.04578688 0.17477302  0.5489901  0.71673726      0.01101735
#            MElightsteelblue    MEbrown2 MEindianred2 MEslateblue1  MEdeeppink MEmediumorchid4 MEorangered MEchocolate3 MElightskyblue3 timepoint
#SRR14333319     -0.050355185 -0.09514398   -0.2909084  -0.27662380  0.05507943      -0.2660329 -0.28785215   -0.1519834      0.01368040       III
#SRR14333320     -0.066222234 -0.17832487   -0.3511874  -0.42591215 -0.01992492      -0.3910929 -0.38209651   -0.1264370     -0.02660291       III
#SRR14333321      0.006919412 -0.26106483   -0.2989501  -0.47068807 -0.11558449      -0.3247613 -0.37079548   -0.1017485     -0.09405620       III
#SRR14333322      0.173182136 -0.60666449    0.2123955  -0.04242742 -0.50533401       0.1596205  0.18466439   -0.2071656     -0.18001566        II
#SRR14333323      0.616909680  0.44990667    0.6005451   0.56841160  0.23088164      -0.3989967 -0.06368611   -0.3943578     -0.30446239        II
#SRR14333324      0.407776374  0.48097587    0.4687780   0.23900674  0.34635842       0.4595326  0.73724638    0.3820427      0.84293622        II
#              sample_id
#SRR14333319 SRR14333319
#SRR14333320 SRR14333320
#SRR14333321 SRR14333321
#SRR14333322 SRR14333322
#SRR14333323 SRR14333323
#SRR14333324 SRR14333324

head(Strader_MEs$sample_id)
#head(Strader_MEs$sample_id)
#[1] "SRR14333319" "SRR14333320" "SRR14333321" "SRR14333322" "SRR14333323" "SRR14333324"

head(Strader_MEs$timepoint)

#[1] III III III II  II  II 
#Levels: I II III

# Calculate 10 over-arching expression patterns using mean eigengene for each module in a cluster. Use MEcluster as models for first lines
# Create a column to the Strader_MEs data frame containing lifestage groups

C1_Strader_MEs <- select(Strader_MEs, MEcornflowerblue:MEindianred3)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)
C2_Strader_MEs <- select(Strader_MEs, MEsalmon2:MEroyalblue2)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)
C3_Strader_MEs <- select(Strader_MEs, MEblue4:MEpink4)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)
C4_Strader_MEs <- select(Strader_MEs, MEcoral:MEtan3)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)
C5_Strader_MEs <- select(Strader_MEs, MEbrown4:MEpalevioletred)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)
C6_Strader_MEs <- select(Strader_MEs, MEdeeppink2:MEwhitesmoke)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)
C7_Strader_MEs <- select(Strader_MEs, MEbrown3:MEdarkgoldenrod4)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)
C8_Strader_MEs <- select(Strader_MEs, MEivory:MElightsteelblue)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)
C9_Strader_MEs <- select(Strader_MEs, MEbrown2:MEslateblue1)
C9_Strader_MEs$Mean <- rowMeans(C9_Strader_MEs)
C10_Strader_MEs <- select(Strader_MEs, MEdeeppink:MElightskyblue3)
C10_Strader_MEs$Mean <- rowMeans(C10_Strader_MEs)

expressionProfile_data <- as.data.frame(cbind(group = Strader_MEs$timepoint, cluster1= C1_Strader_MEs$Mean, cluster2 = C2_Strader_MEs$Mean, 
                                              cluster3 = C3_Strader_MEs$Mean, cluster4 = C4_Strader_MEs$Mean,cluster5 = C5_Strader_MEs$Mean,cluster6 = C6_Strader_MEs$Mean,
                                              cluster7 = C7_Strader_MEs$Mean,cluster8 = C8_Strader_MEs$Mean,cluster9 = C9_Strader_MEs$Mean,cluster10 = C10_Strader_MEs$Mean))

head(expressionProfile_data)
# group           cluster1            cluster2             cluster3           cluster4           cluster5            cluster6
# 1   III 0.0347503219340546   0.123222186039259    0.270687176989076  0.272538961872467 -0.107472623284725   0.171613300130117
# 2   III  0.128614763032405   0.162195631904521    0.333043128229262  0.310837000237138 -0.161955351931919   0.157414075776753
# 3   III  0.116699044378824   0.183531355242665    0.333204809408833  0.259444960467366 -0.225079015210858   0.162653241364556
# 4    II -0.286281057965992  -0.111874017082919 -0.00424689338392289 -0.465261565935405 -0.126839019550724    0.11630501896857
# 5    II  0.193518011389986   0.393189595546525   -0.154598640345414 -0.241059990482721 -0.297093191480431  -0.695201587099655
# 6    II  0.584050228911511 -0.0409081462228039   0.0350294352655437 -0.284491296438426 -0.139882974571766 -0.0391768519235037
#                cluster7           cluster8           cluster9          cluster10
# 1 -0.000750212792476973 -0.177045791163233 -0.220892062967631 -0.127421725932282
# 2    0.0386003385746464 -0.169789095359217 -0.318474807381929 -0.189230857906726
# 3    0.0266636289466788 -0.117273096232804 -0.343567676884239 -0.201389205401437
# 4     0.266186930158272  0.558298824762024 -0.145565482479009 -0.109646072952204
# 5    -0.477598077155284  0.230515389848139  0.539621137729709 -0.186124265731688
# 6    -0.182729735523429  0.116393429498723  0.396253539299589  0.553623273155967

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
#[1] "III" "III" "III" "II"  "II"  "II"  "I"   "I"   "I" 

#cols.num <- c(2:11)
#expressionProfile_data[cols.num] <- sapply(expressionProfile_data[cols.num],as.numeric)

expressionProfile_data$group <- factor(expressionProfile_data$group)


#I want clusters to be numeric, not factors, so I'm using the following custom-made function. Simply using 'as.numeric' I loose the real numbers
#as.double.factor <- function(x) {as.numeric(levels(x))[x]} 
#these steps are from Federica's code, but didn't work for me. I think I reached the same solution though. 
#expressionProfile_data$cluster1 <- as.double.factor(expressionProfile_data$cluster1)
#expressionProfile_data$cluster2 <- as.double.factor(expressionProfile_data$cluster2)
#expressionProfile_data$cluster3 <- as.double.factor(expressionProfile_data$cluster3)
#expressionProfile_data$cluster4 <- as.double.factor(expressionProfile_data$cluster4)
#expressionProfile_data$cluster5 <- as.double.factor(expressionProfile_data$cluster5)
#expressionProfile_data$cluster6 <- as.double.factor(expressionProfile_data$cluster6)
#expressionProfile_data$cluster7 <- as.double.factor(expressionProfile_data$cluster7)
#expressionProfile_data$cluster8 <- as.double.factor(expressionProfile_data$cluster8)
#expressionProfile_data$cluster9 <- as.double.factor(expressionProfile_data$cluster9)
#expressionProfile_data$cluster10 <- as.double.factor(expressionProfile_data$cluster10)

sapply(expressionProfile_data, class)
# [1] group  cluster1  cluster2  cluster3  cluster4  cluster5  cluster6  cluster7
#  "factor" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric"
#  cluster8  cluster9 cluster10
# "numeric" "numeric" "numeric"

dim(expressionProfile_data)
## [1]  9 11

head(expressionProfile_data)


library(ggplot2)

group_order = c("I","II","III")

Cluster1Plot <- expressionProfile_data %>%
  select(group, cluster1) %>% 
  #group_by(group) %>% 
  ggplot(aes(x=group, y=cluster1, fill=group)) +
  geom_boxplot(width=.5, outlier.shape= NA, position = position_dodge(width = 0.5), alpha = 0.7) +
  #stat_summary(fun=mean, geom="line", aes(group=group, color = group), position = position_dodge(width = 0.5))  + 
  stat_summary(fun=mean, geom="point", shape=20, size=14, color="red", fill="red")  + 
  geom_point(pch = 21, size=5, position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values = c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round","III"="spat")) +
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
  scale_fill_manual(values = c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values = c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round","III"="spat")) +
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
  scale_fill_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round", "III"="spat")) +
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
  scale_fill_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round", "III"="spat")) +
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
  scale_fill_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round", "III"="spat")) +
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
  scale_fill_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round", "III"="spat")) +
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
  scale_fill_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round", "III"="spat")) +
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
  scale_fill_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round", "III"="spat")) +
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
  scale_fill_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round", "III"="spat")) +
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
  scale_fill_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) +
  scale_color_manual(values=c("I"="deeppink1", "II"="darkturquoise", "III"="lightpink3")) + 
  scale_x_discrete(labels=c("I"="swimming", "II"="round", "III"="spat")) +
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
ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/expression_eigengene_Profiles_withMean.pdf", expressionProfiles, height = 25, width = 28, units = "in")

#Second version with plot sparcerer

expressionProfilesv2 <- (Cluster1Plot + Cluster2Plot + Cluster3Plot + Cluster4Plot) /
  (Cluster5Plot + Cluster6Plot + Cluster7Plot + Cluster8Plot) /
  (Cluster9Plot + Cluster10Plot + plot_spacer() + plot_spacer())

ggsave("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/expression_eigengene_Profiles_withMeanv2.pdf", expressionProfilesv2, height = 25, width = 28, units = "in")
save.image()         

#  Gene relationship to trait and important modules: Gene Significance and Module Membership

#We quantify associations of individual genes with life stage by defining Gene Significance GS as the absolute value of the correlation between the gene and the time_point. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 

#Define variable weight containing the weight column of datTrait
treatmentinfo <- filter(treatmentinfo, sampleID %in% rownames(datExpr))

# Ensure 'treatmentinfo$timepoint' is a factor
treatmentinfo$timepoint <- as.factor(treatmentinfo$timepoint)

# Create a mapping from factor levels to numeric values
timepoint_levels <- levels(treatmentinfo$timepoint)
level_mapping <- setNames(seq_along(timepoint_levels), timepoint_levels)

# Apply the mapping to create the numeric group values
group <- data.frame(group = level_mapping[as.character(treatmentinfo$timepoint)])

# Display the group dataframe
print(group)
#print(group)
group
#1     3
#2     3
#3     3
#4     2
#5     2
#6     2
#7     1
#8     1
#9     1


datExpr <- datExpr[treatmentinfo$sampleID,]
MEs <- MEs[treatmentinfo$sampleID,]

group_s_pistillata <- data.frame(
  I = c(0, 0, 0, 0, 0, 0, 1, 1, 1),
  II = c(0, 0, 0, 1, 1, 1, 0, 0, 0),
  III = c(1, 1, 1, 0, 0, 0, 0, 0, 0),
  row.names = c("SRR14333319", "SRR14333320", "SRR14333321", "SRR14333322", "SRR14333323", "SRR14333324", "SRR14333325", "SRR14333326", "SRR14333327")
)

# Display the data frame
print(group_s_pistillata)
print(group_s_pistillata)
#I II III
#SRR14333319 0  0   1
#SRR14333320 0  0   1
#SRR14333321 0  0   1
#SRR14333322 0  1   0
#SRR14333323 0  1   0
#SRR14333324 0  1   0
#SRR14333325 1  0   0
#SRR14333326 1  0   0
#SRR14333327 1  0   0
dim(group_s_pistillata)
#[1] 9 3
# Colors of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Use group_s_pistillata for geneTraitSignificance
geneTraitSignificance_s_pistillata = as.data.frame(cor(datExpr, group_s_pistillata, use = "p"))
GSPvalue_s_pistillata = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_s_pistillata), nSamples))


names(geneTraitSignificance_s_pistillata) = paste("GS.", names(group_s_pistillata), sep="")
names(GSPvalue_s_pistillata) = paste("p.GS.", names(group_s_pistillata), sep="")

#Use group for geneTraitSignificance for groups
geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(group), sep="")
names(GSPvalue) = paste("p.GS.", names(group), sep="")


#Import annotation file

annot <- read_csv("Spis_annot.csv")

#Match up genes in datExpr to annotation file
names(annot)
# [1] "gene_id"         "Source"          "Hit accession"  
#[4] "Hit description" "Query length"    "Hit length"     
#[7] "Max bit score"   "Total bit score" "Identity"       
#[10] "Identity %"      "Coverage %"      "Expect"         
#[13] "GO_IDs"  
probes = names(datExpr)
probes2annot = match(probes, annot$gene_id)
nrow(annot)
head(annot)

# The following is the number of probes without annotation... Should return 0.
sum(is.na(probes2annot))
#[1] 0

geneInfo0 = data.frame(
  gene_id = probes,
  `Hit accession` = annot$`Hit accession`[probes2annot],
  `Hit description` = annot$`Hit description`[probes2annot],
  `Query length` = annot$`Query length`[probes2annot],
  `Hit length` = annot$`Hit length`[probes2annot],
  `Max bit score` = annot$`Max bit score`[probes2annot],
  `Total bit score` = annot$`Total bit score`[probes2annot],
  Identity = annot$Identity[probes2annot],
  `Identity %` = annot$`Identity %`[probes2annot],
  `Coverage %` = annot$`Coverage %`[probes2annot],
  Expect = annot$Expect[probes2annot],
  GO_IDs = annot$GO_IDs[probes2annot],
  moduleColor = dynamicColors,
  geneTraitSignificance = geneTraitSignificance,
  GSPvalue = GSPvalue,
  geneTraitSignificance_s_pistillata = geneTraitSignificance_s_pistillata,
  GSPvalue_s_pistillata = GSPvalue_s_pistillata
)


#Order modules by their significance for time_point
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

#            Query.length Hit.length Max.bit.score Total.bit.score
#SpisGene7898           374        361        103.22          103.22
#SpisGene20122          740        567        300.06          300.06
#SpisGene5722           474        474        304.68          304.68
#SpisGene1436           450        366        252.68          252.68
#SpisGene11499          195        197         68.17           68.17
#Identity Identity.. Coverage..   Expect
#SpisGene7898    95/340     0.2794     0.8770 2.02e-23
#SpisGene20122  156/375     0.4160     0.4878 3.43e-90
#SpisGene5722   186/485     0.3835     0.9937 4.69e-96
#SpisGene1436   151/370     0.4081     0.8022 1.22e-77
#SpisGene11499   47/138     0.3406     0.7026 3.26e-11

#Add module cluster information
geneInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")
dim(geneInfo)
#[1] 15998   184
head(geneInfo)
#looks like it has far too many MM, check 
write.csv(geneInfo, file = "GeneInfo_S_pistillata.csv")


#Old code, keeping just to double-check
# Define group containing columns for each lifestage I think this is not at use, revise
group <- data.frame(
  Swimming = c(0, 0, 0, 0, 0, 0, 1, 1, 1),
  Round = c(0, 0, 0, 1, 1, 1, 0, 0, 0),
  Spat = c(1, 1, 1, 0, 0, 0, 0, 0, 0)
)

# Calculate gene significance for each lifestage
geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS", colnames(group), sep = ".")
names(GSPvalue) = paste("p.GS", colnames(group), sep = ".")


# We quantify associations of individual genes with lifestage groups by defining Gene Significance GS as the absolute value of the correlation between the gene 
# and the group. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene 
# expression profile. THese steps have some issues. Try with this code (from 326 & then from 425 https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/809f516ba78b0fcb162e52d2ba84b75ed204fe13/3-WGCNA/WGCNA.Rmd#L346)

# Define variable weight containing the weight column of datTrait... Not sure if I interpreted this correctly
#group <- as.data.frame(c(3,3,3,
#                         2,2,2,
#                         1,1,1))


group <- read.table(text = "
SRR14333319 0 0 1
SRR14333320 0 0 1
SRR14333321 0 0 1
SRR14333322 0 1 0
SRR14333323 0 1 0
SRR14333324 0 1 0
SRR14333325 1 0 0
SRR14333326 1 0 0
SRR14333327 1 0 0
", header = FALSE)

# Assign column names manually, leaving the first column unnamed
colnames(group) <- c("", "Swimming", "Round", "Spat")

# View the data
print(group)

names(group) = "group"
dim(group)
#[1] 9 4

#Colors of the modules

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(group), sep="");
names(GSPvalue) = paste("p.GS.", names(group), sep="")

## Summary output of network analysis results
# Make a dataframe that connects traits, genes, and gene annotation

GO.annot <- read_csv("Spis_annot.csv")

#Match up genes in datExpr to annotation file

names(GO.annot)
#[1] "gene_id"         "Source"          "Hit accession"  
#[4] "Hit description" "Query length"    "Hit length"     
#[7] "Max bit score"   "Total bit score" "Identity"       
#[10] "Identity %"      "Coverage %"      "Expect"         
#[13] "GO_IDs"         

probes = names(datExpr)
probes2annot = match(probes, GO.annot$gene_id)

# The following is the number of probes without annotation... Should return 0.
sum(is.na(probes2annot))
#0

#Create the starting data frame
geneInfo0 = data.frame(gene_id = probes,
                       Annotation.GO.ID = GO.annot$GO_IDs[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

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


#Add module cluster in geneInfo
geneInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")
dim(geneInfo)
head(geneInfo)
#[1] 15998    12

#Save geneInfo as a CSV
geneInfo$Annotation.GO.ID <- gsub(";NA", "", geneInfo$Annotation.GO.ID) #Remove NAs
geneInfo$Annotation.GO.ID <- gsub("NA", "", geneInfo$Annotation.GO.ID) #Remove NAs
write.csv(geneInfo, file = "geneInfoLifestage.csv")

#Old code

# Define variable weight containing the weight column of datTrait... Not sure if I interpreted this correctly
group <- as.data.frame(c(3,3,3,
                         2,2,2,
                         1,1,1))

names(group) = "group"
dim(group)
#[1] 9 1

#Colors of the modules

modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(group), sep="");
names(GSPvalue) = paste("p.GS.", names(group), sep="")

## Summary output of network analysis results
# Make a dataframe that connects traits, genes, and gene annotation

GO.annot <- read_csv("Spis_annot.csv")

#Match up genes in datExpr to annotation file

names(GO.annot)
#[1] "gene_id"         "Source"          "Hit accession"  
#[4] "Hit description" "Query length"    "Hit length"     
#[7] "Max bit score"   "Total bit score" "Identity"       
#[10] "Identity %"      "Coverage %"      "Expect"         
#[13] "GO_IDs"         

probes = names(datExpr)
probes2annot = match(probes, GO.annot$gene_id)

# The following is the number of probes without annotation... Should return 0.
sum(is.na(probes2annot))
#0

#Create the starting data frame
geneInfo0 = data.frame(gene_id = probes,
                       Annotation.GO.ID = GO.annot$GO_IDs[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

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
# gene_id
# SpisGene7898   SpisGene7898
# SpisGene20122  SpisGene20122
# SpisGene5722   SpisGene5722
# SpisGene1436   SpisGene1436
# SpisGene11499  SpisGene11499

# Annotation.GO.ID
# SpisGene7898  GO:0005575,GO:0016020,GO:0016021,GO:0044425
# SpisGene20122  GO:0003674,GO:0005488,GO:0005515,GO:0005575,GO:0005634,GO:0005886,GO:0007165,GO:0008150,GO:0009987,GO:0016020,GO:0019899,GO:0019902,GO:0019903,GO:0031344,GO:0043226,GO:0043227,GO:0043229,GO:0043231,GO:0044087,GO:0044424,GO:0044464,GO:0050789,GO:0050794,GO:0050896,GO:0051128,GO:0051489,GO:0051716,GO:0060491,GO:0065007
# SpisGene5722  GO:0003674,GO:0003824,GO:0005575,GO:0006629,GO:0006644,GO:0006793,GO:0006796,GO:0008150,GO:0008152,GO:0008610,GO:0008654,GO:0009058,GO:0009987,GO:0016020,GO:0016021,GO:0016740,GO:0016746,GO:0019637,GO:0044237,GO:0044238,GO:0044249,GO:0044255,GO:0044425,GO:0044699,GO:0044710,GO:0044711,GO:0044763,GO:0071704,GO:0090407,GO:1901576
# SpisGene1436  GO:0003674,GO:0003824,GO:0006576,GO:0006595,GO:0006596,GO:0006597,GO:0006807,GO:0008150,GO:0008152,GO:0008215,GO:0009058,GO:0009308,GO:0009309,GO:0009987,GO:0016740,GO:0016765,GO:0016768,GO:0034641,GO:0042401,GO:0044106,GO:0044237,GO:0044249,GO:0044271,GO:0071704,GO:1901564,GO:1901566,GO:1901576
# SpisGene11499  <NA>

# moduleColor   GS.group p.GS.group   MM.pink  p.MM.pink    MM.blue
# SpisGene7898  aliceblue -0.7863153 0.01196459 0.7089349 0.03249139 -0.5631009
# SpisGene20122 aliceblue -0.7839731 0.01239931 0.7372489 0.02341340 -0.5502762
# SpisGene5722  aliceblue -0.7774503 0.01366558 0.7344025 0.02423911 -0.4870729
# SpisGene1436  aliceblue -0.7498714 0.01997531 0.7077358 0.03291948 -0.4842566
# SpisGene11499 aliceblue -0.7215222 0.02821506 0.6338781 0.06677818 -0.5263933

# p.MM.blue MM.antiquewhite4 p.MM.antiquewhite4 MM.darkgrey p.MM.darkgrey
# SpisGene7898  0.1144121        0.7087439         0.03255934  -0.8520581  0.0035373623
# SpisGene20122 0.1247483        0.6265572         0.07098622  -0.8414390  0.0044596973
# SpisGene5722  0.1835757        0.6639469         0.05115226  -0.9014457  0.0008975073
# SpisGene1436  0.1865042        0.6455778         0.06038478  -0.8750175  0.0020068841
# SpisGene11499 0.1454296        0.6073069         0.08282843  -0.7852512  0.0121608009

# MM.grey60  p.MM.grey60   MM.black p.MM.black   MM.pink4 p.MM.pink4
# SpisGene7898  -0.8808253 1.709068e-03 -0.5114635  0.1593117 -0.6334306 0.06703080
# SpisGene20122 -0.9532465 6.949904e-05 -0.4650361  0.2071883 -0.6815650 0.04319386
# SpisGene5722  -0.8501398 3.693238e-03 -0.5137190  0.1571673 -0.7643147 0.01647251
# SpisGene1436  -0.8753950 1.986516e-03 -0.4832590  0.1875477 -0.6905702 0.03945425
# SpisGene11499 -0.8615708 2.830678e-03 -0.4288995  0.2493545 -0.6518292 0.05713348

# MM.brown  p.MM.brown   MM.blue4 p.MM.blue4 MM.darkorange
# SpisGene7898  0.7822617 0.012723580 -0.6904390 0.03950721 -0.6866943
# SpisGene20122 0.8393775 0.004656044 -0.6385664 0.06416756 -0.5927274
# SpisGene5722  0.8102901 0.008089304 -0.7160547 0.03002419 -0.6352877
# SpisGene1436  0.8085659 0.008334558 -0.5890347 0.09513156 -0.6033596
# SpisGene11499 0.8179330 0.007060922 -0.4722916 0.19923723 -0.5916078

# p.MM.darkorange  MM.coral2 p.MM.coral2 MM.darkseagreen2
# SpisGene7898  0.04103707 -0.7442401 0.02146462 -0.4404821
# SpisGene20122 0.09256061 -0.6789106 0.04433805 -0.3951494
# SpisGene5722  0.06598642 -0.6030783 0.08558276 -0.4347247
# SpisGene1436  0.08539777 -0.6625653 0.05181306 -0.4560672
# SpisGene11499 0.09333555 -0.6423455 0.06211068 -0.4904708

# p.MM.darkseagreen2 MM.palevioletred p.MM.palevioletred MM.firebrick4
# SpisGene7898  0.2353765 0.7597064 0.0175415434 -0.4256893
# SpisGene20122 0.2925323 0.9091287 0.0006808072 -0.3894806
# SpisGene5722  0.2422704 0.8148298 0.0074670465 -0.4567079
# SpisGene1436  0.2172561 0.8439351 0.0042296355 -0.3131862
# SpisGene11499 0.1800776 0.8006171 0.0095305019 -0.2829920

# p.MM.firebrick4 MM.mediumorchid4 p.MM.mediumorchid4 MM.lightcyan1
# SpisGene7898  0.2533052 0.7173241 0.02959759 -0.6559750
# SpisGene20122 0.3001366 0.7347777 0

#Add module cluster in geneInfo
geneInfo <- left_join(geneInfo, moduleCluster, by = "moduleColor")
dim(geneInfo)
head(geneInfo)
#[1] 15998   168

#Save geneInfo as a CSV
geneInfo$Annotation.GO.ID <- gsub(";NA", "", geneInfo$Annotation.GO.ID) #Remove NAs
geneInfo$Annotation.GO.ID <- gsub("NA", "", geneInfo$Annotation.GO.ID) #Remove NAs
write.csv(geneInfo, file = "geneInfo.csv")

