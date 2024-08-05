# GO Enrichment WGCNA data
setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff")
load("Workspaces_Mcap_fixedGFF") # Load all WCGNA data

# Prepare the dataframe
genes.GO <- as.data.frame(t(datExpr))
genes.GO <- cbind(gene_id = rownames(genes.GO), genes.GO)
rownames(genes.GO) <- NULL

head(genes.GO)
#gene_id      AH1
#1    Montipora_capitata_HIv3___TS.g49315.t1a 7.442739
#2  Montipora_capitata_HIv3___RNAseq.g8204.t1 9.229137
#3  Montipora_capitata_HIv3___RNAseq.g9577.t1 8.566735
#4 Montipora_capitata_HIv3___RNAseq.g23056.t1 9.426291
#5 Montipora_capitata_HIv3___RNAseq.g37262.t1 8.643819
#6  Montipora_capitata_HIv3___RNAseq.g5416.t1 8.134063
#AH2      AH3      AH4      AH5      AH6      AH7
#1 7.676553 7.927178 7.442739 7.754977 7.442739 7.442739
#2 9.417425 9.398208 9.749351 9.899797 9.728030 8.076629
#3 8.695685 8.792582 8.576484 8.459184 8.440691 7.916886
#4 9.453078 9.442767 9.264178 9.467357 9.541120 9.533989
#5 8.503989 8.627807 8.663587 8.696978 8.655319 8.126522
#6 8.071614 8.105080 8.159240 8.143120 8.051393 7.916886
#AH8      AH9
#1 8.007147 8.151469
#2 8.340371 8.168615
#3 8.035718 8.124876
#4 9.477473 9.524112
#5 8.215343 8.552708
#6 7.832487 7.958688
# Select interesting Clusters (the ones with a pattern by lifestage)
# Here I looked at the difference in mean module eigengene value (boxplot profiles)
# and significance in heatmap and decided to analyze cluster. I'll start with cluster 6

## Cluster 6 modules

genes_clust6.GO <- genes.GO[,c(1:10)]
geneColor <- geneInfo %>% select(gene_id, moduleColor) # Make a dataframe containing just gene_id and moduleColor

module_cluster6 <- geneColor %>% filter(moduleColor == "turquoise")

geneColor$gene_id <- as.factor(geneColor$gene_id) # Make factor for merge

genes_clust6.GO$gene_id <- as.factor(genes_clust6.GO$gene_id) # Make factor for merge
genes_clust6.GO <- merge(module_cluster6, genes_clust6.GO)

# I need a vector of gene lengths not present in the annotation file, can get it from the GFF 
## Generate a vector of gene lengths 

# Import file
gff <- rtracklayer::import("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff/Montipora_capitata_HIv3.genes.gff3") # If this doesn't work, restart R and try again 

transcripts <- subset(gff, type == "transcript") # Keep only transcripts 

transcripts_gr <- makeGRangesFromDataFrame(transcripts, keep.extra.columns=TRUE) # Extract length information

transcript_lengths <- width(transcripts_gr) # Isolate length of each transcript

seqnames <- transcripts_gr$ID # Extract list of gene ID 

lengths <- cbind(seqnames, transcript_lengths)

lengths <- as.data.frame(lengths) # Convert to data frame

# Add in length to the annotation file    
annotated_probes$length <- lengths$transcript_lengths[match(annotated_probes$gene, lengths$seqnames)]

# Back to GOseq, requires a vector of all genes and all differentially expressed genes

# Build a data frame that links the gene IDs, modules, and counts of expressed genes and the gene lengths

GO.annot <- subset(annotated_probes, select= c(gene_id, GOs, length)) # Select only relevant information

GOref <- merge(genes.GO, GO.annot, by.x="gene_id")

####### Cluster 6 ######
cluster6_gene.vector <- as.vector(genes_clust6.GO$gene_id)
cluster6_gene.vector=as.integer(GOref$gene_id %in% cluster6_gene.vector)
names(cluster6_gene.vector)=GOref$gene_id
head(cluster6_gene.vector)
#Montipora_capitata_HIv3___RNAseq.10136_t Montipora_capitata_HIv3___RNAseq.10187_t 
#0                                        0 
#Montipora_capitata_HIv3___RNAseq.10207_t Montipora_capitata_HIv3___RNAseq.10214_t 
#0                                        1 
#Montipora_capitata_HIv3___RNAseq.10304_t Montipora_capitata_HIv3___RNAseq.10384_t 
#0                                        0

cluster6_ID.vector <- as.character(GOref$gene_id) # Make ID vector
head(cluster6_ID.vector)
#head(cluster6_ID.vector)
#[1] "Montipora_capitata_HIv3___RNAseq.10136_t" "Montipora_capitata_HIv3___RNAseq.10187_t"
#[3] "Montipora_capitata_HIv3___RNAseq.10207_t" "Montipora_capitata_HIv3___RNAseq.10214_t"
#[5] "Montipora_capitata_HIv3___RNAseq.10304_t" "Montipora_capitata_HIv3___RNAseq.10384_t"
dim(cluster6_ID.vector)
#NULL 

# Make length vector for cluster 6
cluster6_Length.vector <- GOref$length # Make length vector
cluster6_Length.vector <- as.numeric(cluster6_Length.vector)
head(cluster6_Length.vector)
# [1] 50267 12774  9051 10360  6889  2773
 

# Calculate Probability Weighting Function for cluster 6
pwf.cluster6 <- nullp(cluster6_gene.vector, cluster6_ID.vector, bias.data=cluster6_Length.vector) # Weight vector by length of gene

#Prepare GO term dataframe
GO.annot <- select(geneInfo, gene_id, Annotation.GO.ID)
splitted <- strsplit(as.character(GO.annot$Annotation.GO.ID), ";") #split into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GO.annot$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")
head(GO.terms) #GO.ID
#1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   -
#  2 GO:0000902,GO:0000904,GO:0001558,GO:0003674,GO:0005488,GO:0005515,GO:0005575,GO:0005622,GO:0005623,GO:0005634,GO:0005654,GO:0005737,GO:0005768,GO:0005829,GO:0005856,GO:0005884,GO:0005886,GO:0006082,GO:0006725,GO:0006807,GO:0006810,GO:0006897,GO:0006907,GO:0006928,GO:0006935,GO:0006950,GO:0006996,GO:0007010,GO:0007015,GO:0007163,GO:0007275,GO:0007399,GO:0007409,GO:0007411,GO:0007596,GO:0007599,GO:0008064,GO:0008150,GO:0008152,GO:0008361,GO:0009605,GO:0009611,GO:0009653,GO:0009966,GO:0009968,GO:0009987,GO:0010638,GO:0010646,GO:0010648,GO:0010720,GO:0010769,GO:0010770,GO:0010810,GO:0010811,GO:0010975,GO:0012505,GO:0015629,GO:0016020,GO:0016043,GO:0016192,GO:0016477,GO:0016604,GO:0016607,GO:0022008,GO:0022603,GO:0022604,GO:0022607,GO:0023051,GO:0023057,GO:0030027,GO:0030029,GO:0030030,GO:0030031,GO:0030032,GO:0030036,GO:0030154,GO:0030155,GO:0030182,GO:0030334,GO:0030335,GO:0030516,GO:0030832,GO:0030833,GO:0030834,GO:0030836,GO:0030838,GO:0031175,GO:0031252,GO:0031334,GO:0031344,GO:0031346,GO:0031410,GO:0031529,GO:0031941,GO:0031974,GO:0031981,GO:0031982,GO:0032101,GO:0032231,GO:0032233,GO:0032271,GO:0032273,GO:0032501,GO:0032502,GO:0032535,GO:0032794,GO:0032879,GO:0032956,GO:0032970,GO:0032989,GO:0032990,GO:0032991,GO:0033043,GO:0034641,GO:0035020,GO:0035021,GO:0040008,GO:0040011,GO:0040012,GO:0040017,GO:0042060,GO:0042221,GO:0042330,GO:0042995,GO:0043226,GO:0043227,GO:0043228,GO:0043229,GO:0043231,GO:0043232,GO:0043233,GO:0043243,GO:0043244,GO:0043254,GO:0044085,GO:0044087,GO:0044089,GO:0044237,GO:0044281,GO:0044351,GO:0044352,GO:0044354,GO:0044422,GO:0044424,GO:0044428,GO:0044430,GO:0044444,GO:0044446,GO:0044451,GO:0044464,GO:0044877,GO:0045111,GO:0045595,GO:0045597,GO:0045664,GO:0045785,GO:0046415,GO:0046483,GO:0046578,GO:0046580,GO:0048468,GO:0048518,GO:0048519,GO:0048522,GO:0048523,GO:0048583,GO:0048585,GO:0048638,GO:0048666,GO:0048667,GO:0048699,GO:0048731,GO:0048812,GO:0048841,GO:0048856,GO:0048858,GO:0048869,GO:0048870,GO:0050767,GO:0050770,GO:0050789,GO:0050793,GO:0050794,GO:0050817,GO:0050878,GO:0050896,GO:0050920,GO:0051056,GO:0051058,GO:0051094,GO:0051128,GO:0051129,GO:0051130,GO:0051179,GO:0051234,GO:0051239,GO:0051270,GO:0051272,GO:0051492,GO:0051493,GO:0051495,GO:0051496,GO:0051638,GO:0051639,GO:0051674,GO:0051695,GO:0051960,GO:0060284,GO:0061387,GO:0061564,GO:0065007,GO:0065008,GO:0070013,GO:0071704,GO:0071840,GO:0071944,GO:0072521,GO:0090066,GO:0097435,GO:0097485,GO:0097581,GO:0097708,GO:0098657,GO:0099080,GO:0099081,GO:0099512,GO:0099513,GO:0110020,GO:0110053,GO:0120025,GO:0120031,GO:0120035,GO:0120036,GO:0120039,GO:1900024,GO:1900026,GO:1901360,GO:1901564,GO:1901879,GO:1901881,GO:1902531,GO:1902532,GO:1902667,GO:1902743,GO:1902745,GO:1902903,GO:1902905,GO:2000026,GO:2000145,GO:2000147,GO:2000812,GO:2000813
tail(GO.terms)
#25738  GO:0000070,GO:0000086,GO:0000226,GO:0000228,GO:0000278,GO:0000280,GO:0000775,GO:0000776,GO:0000793,GO:0000794,GO:0000819,GO:0000922,GO:0001701,GO:0001824,GO:0003674,GO:0003824,GO:0004672,GO:0004674,GO:0005488,GO:0005515,GO:0005575,GO:0005622,GO:0005623,GO:0005634,GO:0005654,GO:0005694,GO:0005737,GO:0005813,GO:0005815,GO:0005819,GO:0005829,GO:0005856,GO:0006464,GO:0006468,GO:0006793,GO:0006796,GO:0006807,GO:0006996,GO:0007010,GO:0007017,GO:0007049,GO:0007051,GO:0007052,GO:0007059,GO:0007088,GO:0007098,GO:0007275,GO:0007346,GO:0008150,GO:0008152,GO:0009790,GO:0009792,GO:0009889,GO:0009891,GO:0009893,GO:0009987,GO:0010389,GO:0010556,GO:0010557,GO:0010564,GO:0010604,GO:0010638,GO:0010639,GO:0010948,GO:0015630,GO:0016043,GO:0016301,GO:0016310,GO:0016740,GO:0016772,GO:0016773,GO:0019219,GO:0019222,GO:0019538,GO:0019899,GO:0019902,GO:0019903,GO:0022402,GO:0022406,GO:0022607,GO:0030030,GO:0030031,GO:0030496,GO:0030997,GO:0031023,GO:0031323,GO:0031325,GO:0031326,GO:0031328,GO:0031974,GO:0031981,GO:0032204,GO:0032206,GO:0032210,GO:0032212,GO:0032501,GO:0032502,GO:0032886,GO:0032991,GO:0033043,GO:0033044,GO:0036211,GO:0043009,GO:0043085,GO:0043170,GO:0043226,GO:0043227,GO:0043228,GO:0043229,GO:0043231,GO:0043232,GO:0043233,GO:0043392,GO:0043412,GO:0044085,GO:0044092,GO:0044093,GO:0044237,GO:0044238,GO:0044260,GO:0044267,GO:0044422,GO:0044424,GO:0044427,GO:0044428,GO:0044430,GO:0044444,GO:0044446,GO:0044464,GO:0044770,GO:0044772,GO:0044782,GO:0044839,GO:0045786,GO:0045935,GO:0046602,GO:0046605,GO:0046606,GO:0046777,GO:0048285,GO:0048518,GO:0048519,GO:0048522,GO:0048523,GO:0048856,GO:0050789,GO:0050790,GO:0050794,GO:0051052,GO:0051054,GO:0051098,GO:0051100,GO:0051101,GO:0051128,GO:0051129,GO:0051130,GO:0051171,GO:0051173,GO:0051179,GO:0051225,GO:0051276,GO:0051299,GO:0051338,GO:0051347,GO:0051493,GO:0051494,GO:0051640,GO:0051641,GO:0051726,GO:0051783,GO:0051972,GO:0051973,GO:0051983,GO:0051988,GO:0060255,GO:0060271,GO:0065007,GO:0065008,GO:0065009,GO:0070013,GO:0070507,GO:0070925,GO:0071704,GO:0071840,GO:0080090,GO:0090307,GO:0097711,GO:0098687,GO:0098813,GO:0120031,GO:0120036,GO:0140014,GO:0140056,GO:0140096,GO:1901564,GO:1901987,GO:1901990,GO:1902749,GO:1902850,GO:1903047,GO:1903126,GO:1904353,GO:1904355,GO:1904356,GO:1904358,GO:2000112,GO:2000278,GO:2000573,GO:2001252                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            GO:0000070,GO:0000086,GO:0000226,GO:0000228,GO:0000278,GO:0000280,GO:0000775,GO:0000776,GO:0000793,GO:0000794,GO:0000819,GO:0000922,GO:0001701,GO:0001824,GO:0003674,GO:0003824,GO:0004672,GO:0004674,GO:0005488,GO:0005515,GO:0005575,GO:0005622,GO:0005623,GO:0005634,GO:0005654,GO:0005694,GO:0005737,GO:0005813,GO:0005815,GO:0005819,GO:0005829,GO:0005856,GO:0006464,GO:0006468,GO:0006793,GO:0006796,GO:0006807,GO:0006996,GO:0007010,GO:0007017,GO:0007049,GO:0007051,GO:0007052,GO:0007059,GO:0007088,GO:0007098,GO:0007275,GO:0007346,GO:0008150,GO:0008152,GO:0009790,GO:0009792,GO:0009889,GO:0009891,GO:0009893,GO:0009987,GO:0010389,GO:0010556,GO:0010557,GO:0010564,GO:0010604,GO:0010638,GO:0010639,GO:0010948,GO:0015630,GO:0016043,GO:0016301,GO:0016310,GO:0016740,GO:0016772,GO:0016773,GO:0019219,GO:0019222,GO:0019538,GO:0019899,GO:0019902,GO:0019903,GO:0022402,GO:0022406,GO:0022607,GO:0030030,GO:0030031,GO:0030496,GO:0030997,GO:0031023,GO:0031323,GO:0031325,GO:0031326,GO:0031328,GO:0031974,GO:0031981,GO:0032204,GO:0032206,GO:0032210,GO:0032212,GO:0032501,GO:0032502,GO:0032886,GO:0032991,GO:0033043,GO:0033044,GO:0036211,GO:0043009,GO:0043085,GO:0043170,GO:0043226,GO:0043227,GO:0043228,GO:0043229,GO:0043231,GO:0043232,GO:0043233,GO:0043392,GO:0043412,GO:0044085,GO:0044092,GO:0044093,GO:0044237,GO:0044238,GO:0044260,GO:0044267,GO:0044422,GO:0044424,GO:0044427,GO:0044428,GO:0044430,GO:0044444,GO:0044446,GO:0044464,GO:0044770,GO:0044772,GO:0044782,GO:0044839,GO:0045786,GO:0045935,GO:0046602,GO:0046605,GO:0046606,GO:0046777,GO:0048285,GO:0048518,GO:0048519,GO:0048522,GO:0048523,GO:0048856,GO:0050789,GO:0050790,GO:0050794,GO:0051052,GO:0051054,GO:0051098,GO:0051100,GO:0051101,GO:0051128,GO:0051129,GO:0051130,GO:0051171,GO:0051173,GO:0051179,GO:0051225,GO:0051276,GO:0051299,GO:0051338,GO:0051347,GO:0051493,GO:0051494,GO:0051640,GO:0051641,GO:0051726,GO:0051783,GO:0051972,GO:0051973,GO:0051983,GO:0051988,GO:0060255,GO:0060271,GO:0065007,GO:0065008,GO:0065009,GO:0070013,GO:0070507,GO:0070925,GO:0071704,GO:0071840,GO:0080090,GO:0090307,GO:0097711,GO:0098687,GO:0098813,GO:0120031,GO:0120036,GO:0140014,GO:0140056,GO:0140096,GO:1901564,GO:1901987,GO:1901990,GO:1902749,GO:1902850,GO:1903047,GO:1903126,GO:1904353,GO:1904355,GO:1904356,GO:1904358,GO:2000112,GO:2000278,GO:2000573,GO:2001252
#25739                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      -

GO.terms$GO.ID<- as.character(GO.terms$GO.ID)
GO.terms$GO.ID <- replace_na(GO.terms$GO.ID, "unknown")
GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)
GO.terms$gene_id <- as.factor(GO.terms$gene_id)
GO.terms$GO.ID <- gsub(" ", "", GO.terms$GO.ID)
dim(GO.terms)
#25741     2

#Troubleshooting once again for the next step
# Remove rows with invalid GO IDs
GO.terms <- GO.terms[!GO.terms$GO.ID %in% c("unknown", "-"), ]

# Split GO terms if they are in a single string separated by commas
splitted <- strsplit(as.character(GO.terms$GO.ID), ",")

# Reformat the dataframe
GO.terms <- data.frame(
  gene_id = rep(GO.terms$gene_id, sapply(splitted, length)),
  GO.ID = unlist(splitted)
)

# Convert GO.ID back to factor
GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)


##Find enriched GO terms, "selection-unbiased testing for category enrichment amongst significantly expressed genes for RNA-seq data"
GOwall.cluster6<- goseq(pwf.cluster6, GOref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
dim(GOwall.cluster6)
#[1] 20369     7

head(GOwall.cluster6)
# category over_represented_pvalue under_represented_pvalue numDEInCat numInCat
#5708  GO:0016572            8.243368e-05                0.9999862         23       28
#3903  GO:0008361            2.052767e-04                0.9998739        146      257
#6714  GO:0022029            2.775207e-04                0.9999095         36       51
#6645  GO:0021885            3.266174e-04                0.9998894         37       53
#10727 GO:0043525            4.166623e-04                0.9998398         43       64
#6601  GO:0021795            5.790695e-04                0.9998088         33       47
#term ontology
#5708                                             <NA>     <NA>
#  3903                          regulation of cell size       BP
#6714                     telencephalon cell migration       BP
#6645                         forebrain cell migration       BP
#10727 positive regulation of neuron apoptotic process       BP
#6601                   cerebral cortex cell migration       BP

# Filter GO terms with p-values < 0.05
cluster6.GO.05 <- GOwall.cluster6$category[GOwall.cluster6$over_represented_pvalue < 0.05]

# Convert to dataframe
cluster6.GO.05 <- data.frame(category = cluster6.GO.05)

# Merge with the original GO results to get detailed information
cluster6.GO.05 <- merge(cluster6.GO.05, GOwall.cluster6, by = "category")

# Order by ontology, p-value, and number of DE genes
cluster6.GO.05 <- cluster6.GO.05[order(cluster6.GO.05$ontology, cluster6.GO.05$over_represented_pvalue, -cluster6.GO.05$numDEInCat),]

# Convert term column to factor
cluster6.GO.05$term <- as.factor(cluster6.GO.05$term)

# Display the dimensions of the filtered dataframe
dim(cluster6.GO.05)  # Number of significant GO terms
#587   7

# Save significant terms for cluster 6
write.csv(cluster6.GO.05, file = "GO_WGCNA_cluster6.csv", row.names = FALSE)
