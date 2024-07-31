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
