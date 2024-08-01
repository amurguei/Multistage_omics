#GO Enrichment WGCNA data
#GO Enrichment WGCNA data
setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_GFF")
load("Environment_P_acuta.RData") #load all WCGNA data

#Prepare the dataframe
genes.GO <- as.data.frame(t(datExpr))
genes.GO <- cbind(gene_id = rownames(genes.GO), genes.GO)
rownames(genes.GO) <- NULL

head(genes.GO)
#gene_id SRR3051863 SRR3051864
#1 Pocillopora_acuta_HIv2___RNAseq.g23616.t1   7.500246   7.542322
#2  Pocillopora_acuta_HIv2___RNAseq.g4197.t1   7.216874   7.724165
#3 Pocillopora_acuta_HIv2___RNAseq.g12137.t1   7.980040   8.142403
#4  Pocillopora_acuta_HIv2___RNAseq.g4784.t1   7.527193   7.542322
#5      Pocillopora_acuta_HIv2___TS.g3788.t1   7.706131   7.788385
#6 Pocillopora_acuta_HIv2___RNAseq.g23764.t1   7.689690   7.712617
#SRR3051865 SRR3051866 SRR3051867 SRR3051868 SRR3051869 SRR3051870
#1   7.660746   7.216874   7.216874   7.468601   7.594604   7.584190
#2   7.522999   7.738497   7.216874   7.733712   7.406143   7.464883
#3   8.148667   7.519123   7.575558   7.776943   7.811606   7.810693
#4   7.624451   7.738497   7.645114   7.216874   7.216874   7.645390
#5   7.791389   7.451166   7.488312   7.342857   7.371449   7.216874
#6   7.539492   7.860950   7.773818   7.651782   7.578609   7.549287
#SRR3051871
#1   7.711538
#2   7.216874
#3   7.868900
#4   7.673264
#5   7.327999
#6   7.465113

### Select interesting Clusters (the ones with a pattern by lifestage)
#here I looked at the difference in mean module eigengene value (boxplot profiles)
#and significance in heatmap and decided to analyze cluster. I'll start with cluster 8

## Cluster 8 modules#check tomorrow with chatGPT

genes_clust8.GO <- genes.GO[,c(1:10)]
geneColor <- geneInfo %>% select(gene_id, moduleColor) #Make a dataframe containing just gene_id and moduleColor

module_cluster8 <- geneColor %>% filter(moduleColor == "salmon")

geneColor$gene_id <- as.factor(geneColor$gene_id) #Make factor for merge

genes_clust8.GO$gene_id <- as.factor(genes_clust8.GO$gene_id) #Make factor for merge
genes_clust8.GO <- merge(module_cluster8, genes_clust8.GO)


#I need a vector of gene lengths not present in the annotation file, can get it from the GFF 
## Generate a vector of gene lenghts 

#import file
gff <- rtracklayer::import("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_GFF/Pocillopora_acuta_HIv2.genes_fixed.gff3") #if this doesn't work, restart R and try again 

transcripts <- subset(gff, type == "transcript") #keep only transcripts 

transcripts_gr <- makeGRangesFromDataFrame(transcripts, keep.extra.columns=TRUE) #extract length information

transcript_lengths <- width(transcripts_gr) #isolate length of each transcript

seqnames<-transcripts_gr$ID #extract list of gene id 

lengths<-cbind(seqnames, transcript_lengths)

lengths<-as.data.frame(lengths) #convert to data frame

#Add in length to the annotation file.    

annotated_probes$length<-lengths$transcript_lengths[match(annotated_probes$gene, lengths$seqnames)]

which(is.na(annotated_probes$length)) #all genes have lengths 

#Back to GOseq, requires a vector of all genes and all differentially expressed genes. 

#Build a data frame that links the gene IDs, modules, and counts of expressed genes and the gene lengths.

GO.annot <- subset(annotated_probes, select= c(gene_id, GOs, length)) #Select only relevant information

GOref <- merge(genes.GO, GO.annot, by.x="gene_id")
head(GOref)

#GO.annot2 <- subset(annot, select= c(gene_id, GOs)) #Select only relevant information

#GOref2 <- merge(genes.GO, annot, by.x="gene_id")
#head(GOref2)

dim(GOref)
#[1] 20891    12
####### Cluster 8 ######
cluster8_gene.vector <- as.vector(genes_clust8.GO$gene_id)
cluster8_gene.vector=as.integer(GOref$gene_id %in% cluster8_gene.vector)
names(cluster8_gene.vector)=GOref$gene_id
head(cluster8_gene.vector)

cluster8_ID.vector <- as.character(GOref$gene_id) # Make ID vector
head(cluster8_ID.vector)
dim(cluster8_ID.vector)
#Null 

# Make length vector for cluster 8
cluster8_Length.vector <- GOref$length # Make length vector
cluster8_Length.vector <- as.numeric(cluster8_Length.vector)
head(cluster8_Length.vector)
# "9417"  "8248"  "18867" "12534" "13589" "6685" 

# Calculate Probability Weighting Function for cluster 8
pwf.cluster8 <- nullp(cluster8_gene.vector, cluster8_ID.vector, bias.data=cluster8_Length.vector) # Weight vector by length of gene

#Prepare GO term dataframe
GO.annot <- select(geneInfo, gene_id, Annotation.GO.ID)
splitted <- strsplit(as.character(GO.annot$Annotation.GO.ID), ";") #split into multiple GO ids
GO.terms <- data.frame(v1 = rep.int(GO.annot$gene_id, sapply(splitted, length)), v2 = unlist(splitted)) #list all genes with each of their GO terms in a single row
colnames(GO.terms) <- c("gene_id", "GO.ID")
head(GO.terms)
tail(GO.terms)

GO.terms$GO.ID<- as.character(GO.terms$GO.ID)
GO.terms$GO.ID <- replace_na(GO.terms$GO.ID, "unknown")
GO.terms$GO.ID <- as.factor(GO.terms$GO.ID)
GO.terms$gene_id <- as.factor(GO.terms$gene_id)
GO.terms$GO.ID <- gsub(" ", "", GO.terms$GO.ID)
dim(GO.terms)
#20891     2

#Troubleshooting once again for the next step

# Remove rows with invalid GO IDs
GO.terms <- GO.terms[!GO.terms$GO.ID %in% c("unknown", "-"), ]

# Split GO terms if they are in a single string separated by commas
GO.terms$GO.ID <- unlist(strsplit(as.character(GO.terms$GO.ID), ","))

# Reformat the dataframe
GO.terms <- data.frame(gene_id = rep(GO.terms$gene_id, sapply(strsplit(as.character(GO.terms$GO.ID), ","), length)),
                       GO.ID = unlist(strsplit(as.character(GO.terms$GO.ID), ",")))


##Find enriched GO terms, "selection-unbiased testing for category enrichment amongst significantly expressed genes for RNA-seq data"
GOwall.cluster8<- goseq(pwf.cluster8, GOref$gene_id, gene2cat=GO.terms, test.cats=c("GO:CC", "GO:BP", "GO:MF"), method="Wallenius", use_genes_without_cat=TRUE)
dim(GOwall.cluster5)
head(GOwall.cluster5)

# Filter GO terms with p-values < 0.05
cluster8.GO.05 <- GOwall.cluster8$category[GOwall.cluster8$over_represented_pvalue < 0.05]

# Convert to dataframe
cluster8.GO.05 <- data.frame(category = cluster8.GO.05)

# Merge with the original GO results to get detailed information
cluster8.GO.05 <- merge(cluster8.GO.05, GOwall.cluster8, by = "category")

# Order by ontology, p-value, and number of DE genes
cluster8.GO.05 <- cluster8.GO.05[order(cluster8.GO.05$ontology, cluster8.GO.05$over_represented_pvalue, -cluster8.GO.05$numDEInCat),]

# Convert term column to factor
cluster8.GO.05$term <- as.factor(cluster8.GO.05$term)

# Display the dimensions of the filtered dataframe
dim(cluster8.GO.05)  # Number of significant GO terms
#337   7

# Save significant terms for cluster 8
write.csv(cluster8.GO.05, file = "GO_WGCNA_cluster8.csv", row.names = FALSE)



