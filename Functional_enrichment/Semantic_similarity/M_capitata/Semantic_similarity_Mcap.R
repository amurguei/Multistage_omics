rm(list=ls()) # removes all prior objects
library(tidyverse)
library(org.Hs.eg.db)
library(simplifyEnrichment)
#install.packages("magick")
library(magick)

setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/M. capitata/New_genome/fixed_gff")

#Testing semantic similarity by lifestage

larvaeMcapGO <- read_csv("WGCNA_Mcap_larvaeGO.csv", col_names = TRUE)
#
nrow(larvaeMcapGO) # total number of enriched terms. 335
nrow(filter(larvaeMcapGO, ontology == "BP")) # BP terms.182.
nrow(filter(larvaeMcapGO, ontology == "CC")) # CC terms.48.
nrow(filter(larvaeMcapGO, ontology == "MF")) # MF terms.87.

larvaeMcapGO_BP <- larvaeMcapGO %>% filter(ontology == "BP")
larvaeMcapGOBP <- larvaeMcapGO_BP$category

larvaeMcapGO_MF <- larvaeMcapGO %>% filter(ontology == "MF")
larvaeMcapGOMF <- larvaeMcapGO_MF$category

larvaeMcapGO_CC <- larvaeMcapGO %>% filter(ontology == "CC")
larvaeMcapGOCC <- larvaeMcapGO_CC$category

#Calculate a similarity matrix and save the output.
larvaeMcapGOdissBP = GO_similarity(larvaeMcapGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaeMcapGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaeMcapGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

larvaeMcapGOdissBP = GO_similarity(larvaeMcapGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaeMcapGOBPbigger_GOSEM.pdf", width = 10, height = 7)
simplifyGO(larvaeMcapGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=20)
dev.off()

larvaeMcapGOdissMF = GO_similarity(larvaeMcapGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaeMcapGOMF_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaeMcapGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

larvaeMcapGOdissCC = GO_similarity(larvaeMcapGOCC, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaeMcapGOCC_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaeMcapGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

#Metamorphosed

metaMcapGO <- read_csv("WGCNA_Mcap_metaGO.csv", col_names = TRUE)


metaMcapGO_BP <- metaMcapGO %>% filter(ontology == "BP")
metaMcapGOBP <- metaMcapGO_BP$category

metaMcapGO_MF <- metaMcapGO %>% filter(ontology == "MF")
metaMcapGOMF <- metaMcapGO_MF$category

metaMcapGO_CC <- metaMcapGO %>% filter(ontology == "CC")
metaMcapGOCC <- metaMcapGO_CC$category


#Calculate a similarity matrix and save the output.

nrow(metaMcapGO) # total number of enriched terms. 547
nrow(filter(metaMcapGO, ontology == "BP")) # BP terms.360.
nrow(filter(metaMcapGO, ontology == "CC")) # CC terms.72.
nrow(filter(metaMcapGO, ontology == "MF")) # MF terms.84.

metaMcapGOdissBP = GO_similarity(metaMcapGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaMcapGOBP_lesscrowde_GOSEM.pdf", width = 7, height = 5)
simplifyGO(metaMcapGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()


metaMcapGOdissBP = GO_similarity(metaMcapGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaMcapGOBP_lesscrowded_GOSEM.pdf", width = 7, height = 5)
simplifyGO(metaMcapGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=17)
dev.off()

metaMcapGOdissMF = GO_similarity(metaMcapGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaMcapGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(metaMcapGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

metaMcapGOdissCC = GO_similarity(metaMcapGO_CC, ont = "CC", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaMcapGOCC_GOSEM.pdf", width = 7, height = 5)
simplifyGO(metaMcapGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

#Spat


spatMcapGO <- read_csv("WGCNA_Mcap_spatGO.csv", col_names = TRUE)


spatMcapGO_BP <- spatMcapGO %>% filter(ontology == "BP")
spatMcapGOBP <- spatMcapGO_BP$category

spatMcapGO_MF <- spatMcapGO %>% filter(ontology == "MF")
spatMcapGOMF <- spatMcapGO_MF$category

spatMcapGO_CC <- spatMcapGO %>% filter(ontology == "CC")
spatMcapGOCC <- spatMcapGO_CC$category


#Calculate a similarity matrix and save the output.

nrow(spatMcapGO) # total number of enriched terms. 415
nrow(filter(spatMcapGO, ontology == "BP")) # BP terms.283.
nrow(filter(spatMcapGO, ontology == "CC")) # CC terms.47.
nrow(filter(spatMcapGO, ontology == "MF")) # MF terms.62.

spatMcapGOdissBP = GO_similarity(spatMcapGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatMcapGOBP_lesscrowded_GOSEM.pdf", width = 7, height = 5)
simplifyGO(spatMcapGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=17)
dev.off()

spatMcapGOdissBP = GO_similarity(spatMcapGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatMcapGOBP_evenlesscrowdedGOSEM.pdf", width = 10, height = 7)
simplifyGO(spatMcapGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=15)
dev.off()

spatMcapGOdissBP = GO_similarity(spatMcapGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "Semantic_similarity/spatMcapGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(spatMcapGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

spatMcapGOdissMF = GO_similarity(spatMcapGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatMcapGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(spatMcapGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

spatMcapGOdissCC = GO_similarity(spatMcapGOCC, ont = "CC", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatMcapGOCC_GOSEM.pdf", width = 7, height = 5)
simplifyGO(spatMcapGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()