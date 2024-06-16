#Script for adding transcript and gene id into GFF file for alignment. Modified from: https://github.com/AHuffmyer/larval_symbiont_TPC/blob/main/scripts/bioinformatics/fix_gff_format.Rmd  

#Here, I'll be adding transcript_id= and gene_id= to 'gene' column because we needs that label to map our RNAseq data  

#Load libraries and data. 

#Load libraries
library(tidyverse)
#install.packages("R.utils")
library(R.utils)
```

#Load  gene gff file

gff <- read.csv("Multistage_omics/R scripts/M. capitata/New_genome/Montipora_capitata_HIv3.genes.gff3", header=FALSE, sep="\t")

#Rename columns

colnames(gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")

head(gff)

#Create transcript ID  

gff$transcript_id <- sub(";.*", "", gff$gene)
gff$transcript_id <- gsub("ID=", "", gff$transcript_id) #remove ID= 
gff$transcript_id <- gsub("Parent=", "", gff$transcript_id) #remove ID= 
head(gff)

#Create Parent ID 

gff$parent_id <- sub(".*Parent=", "", gff$gene)
gff$parent_id <- sub(";.*", "", gff$parent_id)
gff$parent_id <- gsub("ID=", "", gff$parent_id) #remove ID= 
head(gff)

#Now remove the transcript and parent ID separate columns.  

gff<-gff %>%
  select(!transcript_id)%>%
  select(!parent_id)

head(gff)

#Save file. Then upload this to Andromeda for use in bioinformatic steps.  

write.table(gff, file="Multistage_omics/R scripts/M. capitata/New_genome/Montipora_capitata_HIv3.genes_fixed.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

#gzip the file 
#gzip("Multistage_omics/R scripts/M. capitata/New_genome/Montipora_capitata_HIv3.genes_fixed.gff3")
