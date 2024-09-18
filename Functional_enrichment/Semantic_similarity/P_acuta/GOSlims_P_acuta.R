#Find GO Slim terms
library(goseq)
library(tidyverse)
library(GSEABase)
library(data.table)
library(ggplot2)
library(cowplot)
library(patchwork)
library(readr)
setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_GFF")
load("Environment_P_acuta") #load all WCGNA data

#Run GOslim to get broader categories
slim <- getOBOCollection("http://current.geneontology.org/ontology/subsets/goslim_generic.obo") #get GO database

#Load CSV files
larvae_GO <- read_csv("WGCNA_Pacu_larvaeGO.csv")
meta_GO <- read_csv("WGCNA_Pacu_metaGO.csv")
spat_GO <- read_csv("WGCNA_Pacu_spatGO.csv")

#Creating a new column to the code called "dir" prior to merging
larvae_GO$dir <- "larvae"
meta_GO$dir <- "meta" 
spat_GO$dir <- "spat" 
#Merging
all_GO <- bind_rows(larvae_GO, meta_GO, spat_GO)  #bind rows

## BP
BP_GO <- all_GO %>%
  filter(ontology=="BP")
BPGO_collection <- GOCollection(BP_GO$category) #Make library of query terms
slims_bp <- data.frame(goSlim(BPGO_collection, slim, "BP")) #Find common parent terms to slim down our list
slims_bp$category <- row.names(slims_bp) #save rownames as category

## MF
MF_GO <- all_GO %>%
  filter(ontology=="MF")
MFGO_collection <- GOCollection(MF_GO$category) #Make library of query terms
slims_mf <- data.frame(goSlim(MFGO_collection, slim, "MF")) #Find common parent terms to slim down our list
slims_mf$category <- row.names(slims_mf) #save rownames as category

## CC
CC_GO <- all_GO %>%
  filter(ontology=="CC")
CCGO_collection <- GOCollection(CC_GO$category) #Make library of query terms
slims_cc <- data.frame(goSlim(CCGO_collection, slim, "CC")) #Find common parent terms to slim down our list
slims_cc$category <- row.names(slims_cc) #save rownames as category

#Get mapped terms, using functions from Sam White's Biostars [post](https://support.bioconductor.org/p/128407/#128409).
#Write function mappedIds to get the query terms that mapped to the slim categories
mappedIds <-
  function(df, collection, OFFSPRING) #the command to run requires a dataframe of slim terms, like slims_MF above, your list of query terms, and the offspring from the GOCollection by goSlim
  {
    map <- as.list(OFFSPRING)[rownames(df)] # Subset GOcollection offspring by the rownames of your dataframe
    mapped <- lapply(map, intersect, ids(collection)) #Find the terms that intersect between the subset made above of your query terms and the GOids from the GO collection
    df[["go_terms"]] <- vapply(unname(mapped), paste, collapse = ";", character(1L)) #Add column "go_terms" with matching terms 
    df #show resulting dataframe
  }

#Run function for MF and BP terms
BPslim <- mappedIds(slims_bp, BPGO_collection, GOBPOFFSPRING)
MFslim <- mappedIds(slims_mf, MFGO_collection, GOMFOFFSPRING)
CCslim <- mappedIds(slims_cc, CCGO_collection, GOCCOFFSPRING)

#Remove duplicate matches, keeping the broader umbrella term
#BP
BPslim <- filter(BPslim, Count>0 & Term!="biological_process") #filter out empty slims and term "biological process"
BPsplitted <- strsplit(as.character(BPslim$go_terms), ";") #split into multiple GO ids
BPslimX <- data.frame(Term = rep.int(BPslim$Term, sapply(BPsplitted, length)), go_term = unlist(BPsplitted)) #list all
BPslimX <- merge(BPslimX, BPslim[,c(1,3:4)], by="Term") #Add back counts, term, and category info
BPslimX <- unique(setDT(BPslimX)[order(go_term, -Count)], by = "go_term") #remove duplicate offspring terms, keeping only those that appear in the larger umbrella term (larger Count number)
BPslim <- data.frame(slim_term=BPslimX$Term, slim_cat=BPslimX$category, category=BPslimX$go_term) #rename columns
head(BPslim)

#slim_term   slim_cat   category
#1              DNA recombination GO:0006310 GO:0000019
#2 carbohydrate metabolic process GO:0005975 GO:0000272
#3         mRNA metabolic process GO:0016071 GO:0000291
#4              DNA recombination GO:0006310 GO:0000335
#5              DNA recombination GO:0006310 GO:0000337
#6                     DNA repair GO:0006281 GO:0000731

#slim_term   slim_cat   category
#1                 DNA recombination GO:0006310 GO:0000018
#2         cytoskeleton organization GO:0007010 GO:0000022
#3  nitrogen cycle metabolic process GO:0071941 GO:0000050
#4                         signaling GO:0023052 GO:0000075
#5 sulfur compound metabolic process GO:0006790 GO:0000097
#6       DNA-templated transcription GO:0006351 GO:0000122

#MF
MFslim <- filter(MFslim, Count>0 & Term!="molecular_function") #filter out empty slims and term "molecular function"
MFsplitted <- strsplit(as.character(MFslim$go_terms), ";") #split into multiple GO ids
MFslimX <- data.frame(Term = rep.int(MFslim$Term, sapply(MFsplitted, length)), go_term = unlist(MFsplitted)) #list all
MFslimX <- merge(MFslimX, MFslim[,c(1,3:4)], by="Term")  #Add back counts, term, and category info
MFslimX <- unique(setDT(MFslimX)[order(go_term, -Count)], by = "go_term")  #remove duplicate offspring terms, keeping only
MFslim <- data.frame(slim_term=MFslimX$Term, slim_cat=MFslimX$category, category=MFslimX$go_term) #rename columns
head(MFslim)

#slim_term   slim_cat   category
#1 ATP-dependent activity GO:0140657 GO:0000146
#2   transporter activity GO:0005215 GO:0000295
#3            DNA binding GO:0003677 GO:0001147
#4     catalytic activity GO:0003824 GO:0001537
#5     catalytic activity GO:0003824 GO:0001681
#6     catalytic activity GO:0003824 GO:0002083

#CC
CCslim <- filter(CCslim, Count>0 & Term!="cellular component") #filter out empty slims and term "molecular function"
CCsplitted <- strsplit(as.character(CCslim$go_terms), ";") #split into multiple GO ids
CCslimX <- data.frame(Term = rep.int(CCslim$Term, sapply(CCsplitted, length)), go_term = unlist(CCsplitted)) #list all
CCslimX <- merge(CCslimX, CCslim[,c(1,3:4)], by="Term")  #Add back counts, term, and category info
CCslimX <- unique(setDT(CCslimX)[order(go_term, -Count)], by = "go_term")  #remove duplicate offspring terms, keeping only
CCslim <- data.frame(slim_term=CCslimX$Term, slim_cat=CCslimX$category, category=CCslimX$go_term) #rename columns
head(CCslim)
# head(CCslim)
#slim_term   slim_cat   category
#1       organelle GO:0043226 GO:0000123
#2       organelle GO:0043226 GO:0000124
#3       organelle GO:0043226 GO:0000139
#4       organelle GO:0043226 GO:0001725
#5 plasma membrane GO:0005886 GO:0001891
#6       organelle GO:0043226 GO:0002102

#Save slim info with GO enrichment info for heatmap dataframes.
GO.BP <- right_join(BPslim, filter(all_GO, ontology=="BP"), by="category") #add back GO enrichment info for each offspring term
GO.MF <- right_join(MFslim, filter(all_GO, ontology=="MF"), by="category") #add back GO enrichment info for each offspring term
GO.CC <- right_join(CCslim, filter(all_GO, ontology=="CC"), by="category")

GO.BP <- na.omit(GO.BP) 
GO.MF <- na.omit(GO.MF) 
GO.CC <- na.omit(GO.CC)

write.csv(GO.BP, "GOs per slim_list_BP.csv")
write.csv(GO.MF, "GOs per slim_list_MF.csv")
write.csv(GO.CC, "GOs per slim_list_CC.csv")
