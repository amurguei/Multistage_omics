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

