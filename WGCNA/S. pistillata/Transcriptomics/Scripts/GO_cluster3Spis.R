#GO Enrichment WGCNA data
setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata")
load("Workspace S. pistillata") #load all WCGNA data

#Prepare the dataframe
genes.GO <- as.data.frame(t(datExpr))
genes.GO <- cbind(gene_id = rownames(genes.GO), genes.GO)
rownames(genes.GO) <- NULL

head(genes.GO)
#gene_id SRR14333319 SRR14333320 SRR14333321 SRR14333322
#1 SpisGene3668    7.744240    7.568799    7.804169    6.307419
#2 SpisGene3666   10.980818   10.829867   10.828419   13.191497
#3 SpisGene3665    9.357614    9.455700    9.419425    9.446066
#4 SpisGene3664    9.791704    9.865025   10.249185   12.512649
#5 SpisGene3663    6.064986    6.437467    6.511859    6.447295
#6 SpisGene3662    7.523644    7.574262    7.542644    7.124676


#Getting the length from the GFF for the next steps. 
#Import gff
gff <- rtracklayer::import("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata/Spis.genome.annotation.gff3") # If this doesn't work, restart R and try again 
transcripts <- subset(gff, type == "mRNA")
transcripts_gr <- makeGRangesFromDataFrame(transcripts, keep.extra.columns=TRUE) # Extract length information
transcript_lengths <- width(transcripts_gr) # Isolate length of each transcript
seqnames <- transcripts_gr$ID # Extract list of gene ID 

lengths <- cbind(seqnames, transcript_lengths)

lengths <- as.data.frame(lengths) # Convert to data frame

write.csv(lengths, "lengths_to_modify.csv", row.names = FALSE)#I made a CSV to change the names...
#I changed SpisXXXX for SpisGeneXXXX, now realizing there's another issue with the
#lengths per isoform... So re-loading and trying to modify code by taking an average
#length per gene of the 3 isoforms... 

lengths <- read_csv("lenghts_to_modify.csv", col_types = cols(transcript_lengths = col_number()))

# Rename the column 'seqnames' to 'gene_id'
lengths <- lengths %>% rename(gene_id = seqnames)

# Remove transcript identifiers to aggregate at the gene level
lengths$gene_id <- sub("\\.t\\d+", "", lengths$gene_id)  # Remove transcript variants (.t1, .t2, etc.)

# Averaging the lengths
gene_lengths <- lengths %>%
  group_by(gene_id) %>%
  summarise(average_length = mean(transcript_lengths))

# Save the aggregated data
write_csv(gene_lengths, "averaged_gene_lengths.csv")

# Add in length to the annotation file
# Example conversion, assuming probes is a character vector of genes
probes <- data.frame(gene = probes, stringsAsFactors = FALSE)

# Join probes with gene_lengths based on gene_id
probes <- dplyr::left_join(probes, gene_lengths, by = c("gene" = "gene_id"))

# Select interesting Clusters (the ones with a pattern by lifestage)
# Here I looked at the difference in mean module eigengene value (boxplot profiles)
# and significance in heatmap and decided to analyze cluster. I'll start with cluster 3
## Cluster 3 modules

genes_clust3.GO <- genes.GO[,c(1:10)]
geneColor <- geneInfo %>% select(gene_id, moduleColor) # Make a dataframe containing just gene_id and moduleColor

# Specify the colors for cluster3
cluster3_colors <- c("blue4", "darkgrey", "grey60", "lightcyan1", "darkseagreen2", 
                     "cyan", "darkorange", "blue", "indianred4", "bisque4", 
                     "mediumpurple", "firebrick4", "pink4")

# Filter the dataframe for the specified module colors
module_cluster3 <- geneColor %>% filter(moduleColor %in% cluster3_colors)

# Make a dataframe containing just gene_id and moduleColor of cluster3
module_cluster3 <- module_cluster3 %>% select(gene_id, moduleColor)

# Display the first few rows of the filtered dataframe
head(module_cluster3)

#gene_id moduleColor
#1   SpisGene109     bisque4
#2  SpisGene5636     bisque4
#3  SpisGene3973     bisque4
#4  SpisGene6220     bisque4
#5 SpisGene25678     bisque4
#6  SpisGene7710     bisque4

geneColor$gene_id <- as.factor(geneColor$gene_id) # Make factor for merge

genes_clust3.GO$gene_id <- as.factor(genes_clust3.GO$gene_id) # Make factor for merge
genes_clust3.GO <- merge(module_cluster3, genes_clust3.GO)

#Build a data frame that links the gene IDs, modules, and counts of expressed genes and the gene lengths.
probes <- probes %>%
  rename(gene_id = gene)#Intermediate step

probes <- probes %>%
  left_join(annot %>% select(gene_id, GO_IDs), by = "gene_id")

GO.annot <- subset(probes, select= c(gene_id, GO_IDs, average_length)) # Select only relevant information

GOref <- merge(genes.GO, GO.annot, by ="gene_id")
head(GOref)
# head(GOref)
#gene_id SRR14333319 SRR14333320 SRR14333321 SRR14333322
#1     SpisGene1   12.025179   10.368011   10.195286   11.487078
#2    SpisGene10    6.005917    6.589392    6.283911    5.746702
#3   SpisGene100    7.355780    7.495694    7.619552    8.061834
#4  SpisGene1000   12.116459   11.099984   10.695502   13.718149
#5 SpisGene10003    5.797574    5.959827    6.071229    5.033130
#6 SpisGene10005    5.712400    5.546757    5.760242    5.033130
dim(GOref)
#[1] 15998    12

####### Cluster 3 ######
cluster3_gene.vector <- as.vector(genes_clust3.GO$gene_id)
cluster3_gene.vector=as.integer(GOref$gene_id %in% cluster3_gene.vector)
names(cluster3_gene.vector)=GOref$gene_id
head(cluster3_gene.vector)
#SpisGene1    SpisGene10   SpisGene100  SpisGene1000 
#0             0             0             0 
#SpisGene10003 SpisGene10005 
#0             0 

cluster3_ID.vector <- as.character(GOref$gene_id) #Make ID vector
head(cluster3_ID.vector)
#[1] "SpisGene1"     "SpisGene10"    "SpisGene100"   "SpisGene1000" 
#[5] "SpisGene10003" "SpisGene10005"
dim(cluster3_ID.vector)
#NULL


cluster3_Length.vector <- GOref$average_length #Make length vector
head(cluster3_Length.vector)
#[1]  7845.333  6462.000  7573.000  2436.000  3892.000 10256.000

# Recalculate the Probability Weighting Function (PWF)
pwf.cluster3 <- nullp(cluster3_gene.vector, cluster3_ID.vector, bias.data=cluster3_Length.vector)
#Warning message:
#In pcls(G) : initial point very close to some inequality constraints

#Prepare GO term dataframe
GO.annot <- select(geneInfo, gene_id, GO_IDs)
splitted <- strsplit(as.character(GO.annot$GO_IDs), ";") #split into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GO.annot$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")
head(GO.terms)
#gene_id
#1  SpisGene7898
#2 SpisGene20122
#3  SpisGene5722
#4  SpisGene1436
#5 SpisGene11499
#6 SpisGene16839
#GO.ID
#1                                                                                                                                                                                                                                                                                                                                                      GO:0005575,GO:0016020,GO:0016021,GO:0044425
#2                                                                   GO:0003674,GO:0005488,GO:0005515,GO:0005575,GO:0005634,GO:0005886,GO:0007165,GO:0008150,GO:0009987,GO:0016020,GO:0019899,GO:0019902,GO:0019903,GO:0031344,GO:0043226,GO:0043227,GO:0043229,GO:0043231,GO:0044087,GO:0044424,GO:0044464,GO:0050789,GO:0050794,GO:0050896,GO:0051128,GO:0051489,GO:0051716,GO:0060491,GO:0065007
#3                                                        GO:0003674,GO:0003824,GO:0005575,GO:0006629,GO:0006644,GO:0006793,GO:0006796,GO:0008150,GO:0008152,GO:0008610,GO:0008654,GO:0009058,GO:0009987,GO:0016020,GO:0016021,GO:0016740,GO:0016746,GO:0019637,GO:0044237,GO:0044238,GO:0044249,GO:0044255,GO:0044425,GO:0044699,GO:0044710,GO:0044711,GO:0044763,GO:0071704,GO:0090407,GO:1901576
#4                                                                                         GO:0003674,GO:0003824,GO:0006576,GO:0006595,GO:0006596,GO:0006597,GO:0006807,GO:0008150,GO:0008152,GO:0008215,GO:0009058,GO:0009308,GO:0009309,GO:0009987,GO:0016740,GO:0016765,GO:0016768,GO:0034641,GO:0042401,GO:0044106,GO:0044237,GO:0044249,GO:0044271,GO:0071704,GO:1901564,GO:1901566,GO:1901576
#5                        
tail(GO.terms)
gene_id
#15993 SpisGene7328
#15994  SpisGene299
#15995 SpisGene8282
#15996 SpisGene4110
#15997 SpisGene1964
#15998 SpisGene5363
#GO.ID
#15993                                                                                                                                                                                                                                                                                                                                                                                                                                                                          <NA>
# 15994                                                                                                                                                                                                                                        GO:0003674,GO:0005488,GO:0005515,GO:0005575,GO:0005634,GO:0005730,GO:0005737,GO:0005794,GO:0008150,GO:0043226,GO:0043227,GO:0043228,GO:0043229,GO:0043231,GO:0043232,GO:0044422,GO:0044424,GO:0044428,GO:0044444,GO:0044446,GO:0044464
#15995 <NA>

GO.terms$GO.ID<- as.character(GO.terms$GO.ID)
GO.terms$GO.ID <- replace_na(GO.terms$GO.ID, "unknown")
GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)
GO.terms$gene_id <- as.factor(GO.terms$gene_id)
GO.terms$GO.ID <- gsub(" ", "", GO.terms$GO.ID)
dim(GO.terms)
#15998     2

#Troubleshooting ofor the next step
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
GOwall.cluster3<- goseq(pwf.cluster3, GOref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
dim(GOwall.cluster3)
dim(GOwall.cluster3)
#[1] 14740     7

head(GOwall.cluster3)
#category over_represented_pvalue under_represented_pvalue
#2350 GO:0006259            5.962743e-20                        1
#2351 GO:0006260            1.003942e-15                        1
#7592 GO:0034061            3.834537e-14                        1
#1225 GO:0003964            1.737610e-13                        1
#2363 GO:0006278            2.796764e-13                        1
#4490 GO:0015074            1.567938e-12                        1
#numDEInCat numInCat                                   term
#2350        531     1042                  DNA metabolic process
#2351        305      580                        DNA replication
#7592        246      462                DNA polymerase activity
#1225        229      430   RNA-directed DNA polymerase activity
#2363        230      434 RNA-templated DNA biosynthetic process
#4490        148      259                        DNA integration

# Filter GO terms with p-values < 0.05
cluster3.GO.05 <- GOwall.cluster3$category[GOwall.cluster3$over_represented_pvalue < 0.05]

# Convert to dataframe
cluster3.GO.05 <- data.frame(category = cluster3.GO.05)

# Merge with the original GO results to get detailed information
cluster3.GO.05 <- merge(cluster3.GO.05, GOwall.cluster3, by = "category")

# Order by ontology, p-value, and number of DE genes
cluster3.GO.05 <- cluster3.GO.05[order(cluster3.GO.05$ontology, cluster3.GO.05$over_represented_pvalue, -cluster3.GO.05$numDEInCat),]

# Convert term column to factor
cluster3.GO.05$term <- as.factor(cluster3.GO.05$term)

# Display the dimensions of the filtered dataframe
dim(cluster3.GO.05)  # Number of significant GO terms
#516   7

# Save significant terms for cluster 6
write.csv(cluster3.GO.05, file = "GO_WGCNA_cluster3.csv", row.names = FALSE)


