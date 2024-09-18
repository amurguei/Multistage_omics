rm(list=ls()) # removes all prior objects
library(tidyverse)
library(org.Hs.eg.db)
library(simplifyEnrichment)
#install.packages("magick")
library(magick)

setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/P. acuta/New_genome_P_acuta/fixed_GFF")
load("Environment_P_acuta") #load all WCGNA data


#Testing semantic similarity by lifestage

larvaePacuGO <- read_csv("WGCNA_Pacu_larvaeGO.csv", col_names = TRUE)
#
nrow(larvaePacuGO) # total number of enriched terms. 428
nrow(filter(larvaePacuGO, ontology == "BP")) # BP terms.296.
nrow(filter(larvaePacuGO, ontology == "CC")) # CC terms.38.
nrow(filter(larvaePacuGO, ontology == "MF")) # MF terms.76.

larvaePacuGO_BP <- larvaePacuGO %>% filter(ontology == "BP")
larvaePacuGOBP <- larvaePacuGO_BP$category

larvaePacuGO_MF <- larvaePacuGO %>% filter(ontology == "MF")
larvaePacuGOMF <- larvaePacuGO_MF$category

larvaePacuGO_CC <- larvaePacuGO %>% filter(ontology == "CC")
larvaePacuGOCC <- larvaePacuGO_CC$category

#Calculate a similarity matrix and save the output.
larvaePacuGOdissBP = GO_similarity(larvaePacuGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaePacuGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaePacuGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

larvaePacuGOdissBP = GO_similarity(larvaePacuGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaePacuGOBP_GOSEM.pdf", width = 10, height = 7)
simplifyGO(larvaePacuGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=20)
dev.off()

larvaePacuGOdissMF = GO_similarity(larvaePacuGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaePacuGOMF_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaePacuGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()


larvaePacuGOdissCC = GO_similarity(larvaePacuGOCC, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaePacuGOCC_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaePacuGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

#Metamorphosed

metaPacuGO <- read_csv("WGCNA_Pacu_metaGO.csv", col_names = TRUE)


metaPacuGO_BP <- metaPacuGO %>% filter(ontology == "BP")
metaPacuGOBP <- metaPacuGO_BP$category

metaPacuGO_MF <- metaPacuGO %>% filter(ontology == "MF")
metaPacuGOMF <- metaPacuGO_MF$category

metaPacuGO_CC <- metaPacuGO %>% filter(ontology == "CC")
metaPacuGOCC <- metaPacuGO_CC$category


#Calculate a similarity matrix and save the output.

nrow(metaPacuGO) # total number of enriched terms. 547
nrow(filter(metaPacuGO, ontology == "BP")) # BP terms.360.
nrow(filter(metaPacuGO, ontology == "CC")) # CC terms.72.
nrow(filter(metaPacuGO, ontology == "MF")) # MF terms.84.

metaPacuGOdissBP = GO_similarity(metaPacuGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaPacuGOBP_lesscrowde_GOSEM.pdf", width = 7, height = 5)
simplifyGO(metaPacuGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()


metaPacuGOdissBP = GO_similarity(metaPacuGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaPacuGOBP_lesscrowded_GOSEM.pdf", width = 7, height = 5)
simplifyGO(metaPacuGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=17)
dev.off()

metaPacuGOdissMF = GO_similarity(metaPacuGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaPacuGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(metaPacuGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

metaPacuGOdissCC = GO_similarity(metaPacuGO_CC, ont = "CC", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaPacuGOCC_GOSEM.pdf", width = 7, height = 5)
simplifyGO(metaPacuGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

#Spat


spatPacuGO <- read_csv("WGCNA_Pacu_spatGO.csv", col_names = TRUE)


spatPacuGO_BP <- spatPacuGO %>% filter(ontology == "BP")
spatPacuGOBP <- spatPacuGO_BP$category

spatPacuGO_MF <- spatPacuGO %>% filter(ontology == "MF")
spatPacuGOMF <- spatPacuGO_MF$category

spatPacuGO_CC <- spatPacuGO %>% filter(ontology == "CC")
spatPacuGOCC <- spatPacuGO_CC$category


#Calculate a similarity matrix and save the output.

nrow(spatPacuGO) # total number of enriched terms. 415
nrow(filter(spatPacuGO, ontology == "BP")) # BP terms.283.
nrow(filter(spatPacuGO, ontology == "CC")) # CC terms.47.
nrow(filter(spatPacuGO, ontology == "MF")) # MF terms.62.

spatPacuGOdissBP = GO_similarity(spatPacuGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatPacuGOBP_lesscrowded_GOSEM.pdf", width = 7, height = 5)
simplifyGO(spatPacuGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=17)
dev.off()

spatPacuGOdissBP = GO_similarity(spatPacuGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatPacuGOBP_evenlesscrowdedGOSEM.pdf", width = 10, height = 7)
simplifyGO(spatPacuGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=15)
dev.off()

spatPacuGOdissBP = GO_similarity(spatPacuGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "Semantic_similarity/spatPacuGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(spatPacuGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

spatPacuGOdissMF = GO_similarity(spatPacuGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatPacuGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(spatPacuGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

spatPacuGOdissCC = GO_similarity(spatPacuGOCC, ont = "CC", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatPacuGOCC_GOSEM.pdf", width = 7, height = 5)
simplifyGO(spatPacuGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()