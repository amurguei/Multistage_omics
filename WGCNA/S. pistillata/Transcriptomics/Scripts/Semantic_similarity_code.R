rm(list=ls()) # removes all prior objects
library(tidyverse)
library(org.Hs.eg.db)
library(simplifyEnrichment)
#install.packages("magick")
library(magick)

setwd("E:/Users/amurgueitio/Documents/Multistage_omics/R scripts/S. pistillata")
#Testing semantic similarity by cluster

#WGCNA cluster 5

LarvaeSpisGO <- read_csv("GO_WGCNA_cluster5_larvaeup.csv", col_names = TRUE)
#
nrow(LarvaeSpisGO) # total number of enriched terms. 556.
nrow(filter(LarvaeSpisGO, ontology == "BP")) # BP terms.311.
nrow(filter(LarvaeSpisGO, ontology == "CC")) # CC terms.107.
nrow(filter(LarvaeSpisGO, ontology == "MF")) # MF terms.93.

LarvaeSpisGO_BP <- LarvaeSpisGO %>% filter(ontology == "BP")
LarvaeBP <- LarvaeSpisGO_BP$category

LarvaeSpisGO_MF <- LarvaeSpisGO %>% filter(ontology == "MF")
LarvaeMF <- LarvaeSpisGO_MF$category


#Calculate a similarity matrix and save the output.
LarvaedissBP = GO_similarity(LarvaeBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/LarvaeBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(LarvaedissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

pdf(file = "semantic_similarity/LarvaeBP_GOSEM.pdf", width = 10, height = 7)
simplifyGO(LarvaedissBP, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

LarvaedissMF = GO_similarity(LarvaeMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/LarvaeMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(LarvaedissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

#Meta up, cluster 8

Cluster8SpisGO <- read_csv("GO_WGCNA_cluster8_metaup.csv", col_names = TRUE)
#
nrow(Cluster8SpisGO) # total number of enriched terms. 556.
nrow(filter(Cluster8SpisGO, ontology == "BP")) # BP terms.311.
nrow(filter(Cluster8SpisGO, ontology == "CC")) # CC terms.107.
nrow(filter(Cluster8SpisGO, ontology == "MF")) # MF terms.93.

Cluster8SpisGO_BP <- Cluster8SpisGO %>% filter(ontology == "BP")
Cluster8SpisGOBP <- Cluster8SpisGO_BP$category

Cluster8SpisGO_MF <- LarvaeSpisGO %>% filter(ontology == "MF")
Cluster8SpisGOMF <- LarvaeSpisGO_MF$category


#Calculate a similarity matrix and save the output.
Cluster8SpisGOdissBP = GO_similarity(Cluster8SpisGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster8SpisGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(Cluster8SpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()


Cluster8SpisGOdissMF = GO_similarity(Cluster8SpisGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster8SpisGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(LarvaedissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()


#Meta up, cluster 8

Cluster8SpisGO <- read_csv("GO_WGCNA_cluster8_metaup.csv", col_names = TRUE)
#
nrow(Cluster8SpisGO) # total number of enriched terms. 1011
nrow(filter(Cluster8SpisGO, ontology == "BP")) # BP terms.650.
nrow(filter(Cluster8SpisGO, ontology == "CC")) # CC terms.129.
nrow(filter(Cluster8SpisGO, ontology == "MF")) # MF terms.138.

Cluster8SpisGO_BP <- Cluster8SpisGO %>% filter(ontology == "BP")
Cluster8SpisGOBP <- Cluster8SpisGO_BP$category

Cluster8SpisGO_MF <- Cluster8SpisGO %>% filter(ontology == "MF")
Cluster8SpisGOMF <- Cluster8SpisGO_MF$category

#Calculate a similarity matrix and save the output.
Cluster8SpisGOdissBP = GO_similarity(Cluster8SpisGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster8SpisGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(Cluster8SpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()
#Calculate a similarity matrix and save the output.
Cluster8SpisGOdissBP = GO_similarity(Cluster8SpisGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster8SpisGOBP_GOSEM_increasedwidth.pdf", width = 10, height = 7)
simplifyGO(Cluster8SpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()
#Max words increased layout
simplifyGO(Cluster8SpisGOdissBP, word_cloud_layout = "circle", word_cloud_grob_param = list(max_width = 50), max_words = 20)
pdf(file = "semantic_similarity/Cluster8SpisGOBP_GOSEM_max.pdf", width = 10, height = 7)
simplifyGO(Cluster8SpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=25)
dev.off()

Cluster8SpisGOdissMF = GO_similarity(Cluster8SpisGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster8SpisGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(Cluster8SpisGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

#Spat up, cluster 3

Cluster3SpisGO <- read_csv("GO_WGCNA_cluster3.csv", col_names = TRUE)
#
nrow(Cluster3SpisGO) # total number of enriched terms. 516
nrow(filter(Cluster3SpisGO, ontology == "BP")) # BP terms.321.
nrow(filter(Cluster3SpisGO, ontology == "CC")) # CC terms.30.
nrow(filter(Cluster3SpisGO, ontology == "MF")) # MF terms.117.

Cluster3SpisGO_BP <- Cluster3SpisGO %>% filter(ontology == "BP")
Cluster3SpisGOBP <- Cluster3SpisGO_BP$category

Cluster3SpisGO_MF <- Cluster3SpisGO %>% filter(ontology == "MF")
Cluster3SpisGOMF <- Cluster3SpisGO_MF$category

#Calculate a similarity matrix and save the output.
Cluster3SpisGOdissBP = GO_similarity(Cluster3SpisGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster3SpisGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(Cluster3SpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

#Checking if this works
# Calculate GO similarity matrix
Cluster3SpisGOdissBP = GO_similarity(Cluster3SpisGOBP, ont = "BP", db = "org.Hs.eg.db")

# Set up PDF for saving the word cloud with increased dimensions
pdf(file = "semantic_similarity/Cluster3SpisGOBP_modparam_GOSEM.pdf", width = 10, height = 7)

# Create the word cloud with customized parameters to avoid overlapping
simplifyGO(
  Cluster3SpisGOdissBP, 
  word_cloud_grob_param = list(max_width = 60), 
  max_words = 17 # Reduce the number of words for better readability
)

# Close the PDF device
dev.off()

#Calculate a similarity matrix and save the output.
Cluster3SpisGOdissBP = GO_similarity(Cluster3SpisGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster3SpisGOBP_GOSEM_increasedwidth.pdf", width = 10, height = 7)
simplifyGO(Cluster3SpisGOdissBP, word_cloud_grob_param = list(max_width = 55), max_words=20)
dev.off()

Cluster3SpisGOdissMF = GO_similarity(Cluster3SpisGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster3SpisGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(Cluster3SpisGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

#Testing semantic similarity by lifestage

larvaeSpisGO <- read_csv("WGCNA_Spis_larvaeGO.csv", col_names = TRUE)
#
nrow(larvaeSpisGO) # total number of enriched terms. 914
nrow(filter(larvaeSpisGO, ontology == "BP")) # BP terms.569.
nrow(filter(larvaeSpisGO, ontology == "CC")) # CC terms.137.
nrow(filter(larvaeSpisGO, ontology == "MF")) # MF terms.140.

Cluster8SpisGO_BP <- Cluster8SpisGO %>% filter(ontology == "BP")
Cluster8SpisGOBP <- Cluster8SpisGO_BP$category

Cluster8SpisGO_MF <- LarvaeSpisGO %>% filter(ontology == "MF")
Cluster8SpisGOMF <- LarvaeSpisGO_MF$category


#Calculate a similarity matrix and save the output.
Cluster8SpisGOdissBP = GO_similarity(Cluster8SpisGOBP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster8SpisGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(Cluster8SpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()


Cluster8SpisGOdissMF = GO_similarity(Cluster8SpisGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/Cluster8SpisGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(LarvaedissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()




