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

larvaeSpisGO_BP <- larvaeSpisGO %>% filter(ontology == "BP")
larvaeSpisGO_BP <- larvaeSpisGO_BP$category

larvaeSpisGO_MF <- larvaeSpisGO %>% filter(ontology == "MF")
larvaeSpisGOMF <- larvaeSpisGO_MF$category

larvaeSpisGO_CC <- larvaeSpisGO %>% filter(ontology == "CC")
larvaeSpisGOCC <- larvaeSpisGO_CC$category


#Calculate a similarity matrix and save the output.
larvaeSpisGOdissBP = GO_similarity(larvaeSpisGO_BP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaeSpisGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaeSpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()


#Calculate a similarity matrix and save the output.
larvaeSpisGOdissBP = GO_similarity(larvaeSpisGO_BP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaeSpisGOBP_lesscrowded_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaeSpisGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=17)
dev.off()

larvaeSpisGOdissMF = GO_similarity(larvaeSpisGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaeSpisGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(larvaeSpisGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

larvaeSpisGOdissCC = GO_similarity(larvaeSpisGO_CC, ont = "CC", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/larvaeSpisGOCC_GOSEM.pdf", width = 7, height = 5)
simplifyGO(larvaeSpisGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

metaSpisGO <- read_csv("WGCNA_Spis_metaGO.csv", col_names = TRUE)
#
nrow(metaSpisGO) # total number of enriched terms. 895
nrow(filter(metaSpisGO, ontology == "BP")) # BP terms.579.
nrow(filter(metaSpisGO, ontology == "CC")) # CC terms.79.
nrow(filter(metaSpisGO, ontology == "MF")) # MF terms.157.

metaSpisGO_BP <- metaSpisGO %>% filter(ontology == "BP")
metaSpisGOBP <- metaSpisGO_BP$category

metaSpisGO_MF <- metaSpisGO %>% filter(ontology == "MF")
metaSpisGOMF <- metaSpisGO_MF$category

metaSpisGO_CC <- metaSpisGO %>% filter(ontology == "CC")
metaSpisGOCC <- metaSpisGO_CC$category

#Calculate a similarity matrix and save the output.
metaSpisGOdissBP = GO_similarity(metaSpisGO_BP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaSpisGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(metaSpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

#Calculate a similarity matrix and save the output.
metaSpisGOdissBP = GO_similarity(metaSpisGO_BP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaSpisGOBP_lesscrowded_GOSEM.pdf", width = 7, height = 5)
simplifyGO(metaSpisGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=17)
dev.off()


#Calculate a similarity matrix and save the output.
metaSpisGOdissBP = GO_similarity(metaSpisGO_BP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaSpisGOBP_evenlesscrowded_GOSEM.pdf", width = 10, height = 7)
simplifyGO(metaSpisGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=13)
dev.off()

metaSpisGOdissMF = GO_similarity(metaSpisGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaSpisGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(metaSpisGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

metaSpisGOdissCC = GO_similarity(metaSpisGOCC, ont = "CC", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/metaSpisGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(metaSpisGOdissCC, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

spatSpisGO <- read_csv("WGCNA_Spis_spatGO.csv", col_names = TRUE)
#
nrow(spatSpisGO) # total number of enriched terms. 281
nrow(filter(spatSpisGO, ontology == "BP")) # BP terms.178.
nrow(filter(spatSpisGO, ontology == "CC")) # CC terms.29.
nrow(filter(spatSpisGO, ontology == "MF")) # MF terms.55.

spatSpisGO_BP <- spatSpisGO %>% filter(ontology == "BP")
spatSpisGOBP <- spatSpisGO_BP$category

spatSpisGO_MF <- spatSpisGO %>% filter(ontology == "MF")
spatSpisGOMF <- spatSpisGO_MF$category

spatSpisGO_CC <- spatSpisGO %>% filter(ontology == "CC")
spatSpisGOCC <- spatSpisGO_CC$category

#Calculate a similarity matrix and save the output.
spatSpisGOdissBP = GO_similarity(spatSpisGO_BP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatSpisGOBP_GOSEM.pdf", width = 7, height = 5)
simplifyGO(spatSpisGOdissBP, word_cloud_grob_param = list(max_width = 50), max_words=20)
dev.off()

#Calculate a similarity matrix and save the output.
spatSpisGOdissBP = GO_similarity(spatSpisGO_BP, ont = "BP", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatSpisGOBP_lesscrowded_GOSEM.pdf", width = 7, height = 5)
simplifyGO(spatSpisGOdissBP, word_cloud_grob_param = list(max_width = 60), max_words=17)
dev.off()

spatSpisGOdissMF = GO_similarity(spatSpisGOMF, ont = "MF", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatSpisGOMF_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(spatSpisGOdissMF, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

spatSpisGOdissCC = GO_similarity(spatSpisGOCC, ont = "CC", db = "org.Hs.eg.db")
pdf(file = "semantic_similarity/spatSpisGOCC_GOSEM_Increasedwidth.pdf", width = 10, height = 7)
simplifyGO(spatSpisGOdissCC, word_cloud_grob_param = list(max_width = 50), max_words = 20)
dev.off()

#Barplots per lifestage

#Larvae
LarvaeSpisGO_BP %>% 
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  ggplot(aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalueover_represented_pvalue") +
  ggtitle(expression("Larvae"~italic("S. pistillata"))) + # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank(),#Set the plot background
        legend.position="none")

#Only 10 top 
LarvaeSpisGO_BP %>%
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  top_n(-10, over_represented_pvalue) %>%  # Select top 10 terms with smallest p-values
  ggplot(aes(x=term, y=over_represented_pvalue)) +
  geom_segment(aes(x=term, xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalue") +
  ggtitle(expression("Larvae"~italic("S. pistillata"))) +  # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() +  # Set background color
  theme(panel.border = element_blank(),  # Set border
        panel.grid.major = element_blank(),  # Set major gridlines
        panel.grid.minor = element_blank(),  # Set minor gridlines
        axis.line = element_line(colour = "black"),  # Set axes color
        plot.background=element_blank(),  # Set the plot background
        legend.position="none")
#Top20
LarvaeSpisGO_BP %>%
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  top_n(-20, over_represented_pvalue) %>%  # Select top 10 terms with smallest p-values
  ggplot(aes(x=term, y=over_represented_pvalue)) +
  geom_segment(aes(x=term, xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalue") +
  ggtitle(expression("Larvae"~italic("S. pistillata"))) +  # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() +  # Set background color
  theme(panel.border = element_blank(),  # Set border
        panel.grid.major = element_blank(),  # Set major gridlines
        panel.grid.minor = element_blank(),  # Set minor gridlines
        axis.line = element_line(colour = "black"),  # Set axes color
        plot.background=element_blank(),  # Set the plot background
        legend.position="none")

#Meta

metaSpisGO_BP %>% 
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  ggplot(aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalueover_represented_pvalue") +
  ggtitle(expression("Metamorphosed"~italic("S. pistillata"))) + # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank(),#Set the plot background
        legend.position="none")


metaSpisGO_BP %>%
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  top_n(-10, over_represented_pvalue) %>%  # Select top 10 terms with smallest p-values
  ggplot(aes(x=term, y=over_represented_pvalue)) +
  geom_segment(aes(x=term, xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalue") +
  ggtitle(expression("meta"~italic("S. pistillata"))) +  # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() +  # Set background color
  theme(panel.border = element_blank(),  # Set border
        panel.grid.major = element_blank(),  # Set major gridlines
        panel.grid.minor = element_blank(),  # Set minor gridlines
        axis.line = element_line(colour = "black"),  # Set axes color
        plot.background=element_blank(),  # Set the plot background
        legend.position="none")


metaSpisGO_BP %>%
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  top_n(-20, over_represented_pvalue) %>%  # Select top 10 terms with smallest p-values
  ggplot(aes(x=term, y=over_represented_pvalue)) +
  geom_segment(aes(x=term, xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalue") +
  ggtitle(expression("meta"~italic("S. pistillata"))) +  # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() +  # Set background color
  theme(panel.border = element_blank(),  # Set border
        panel.grid.major = element_blank(),  # Set major gridlines
        panel.grid.minor = element_blank(),  # Set minor gridlines
        axis.line = element_line(colour = "black"),  # Set axes color
        plot.background=element_blank(),  # Set the plot background
        legend.position="none")

#Spat

spatSpisGO_BP %>% 
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  ggplot(aes(x=term, y=over_represented_pvalue) ) +
  geom_segment( aes(x=term ,xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalueover_represented_pvalue") +
  ggtitle(expression("Spat"~italic("S. pistillata"))) + # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank(),#Set the plot background
        legend.position="none")


spatSpisGO_BP %>%
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  top_n(-10, over_represented_pvalue) %>%  # Select top 10 terms with smallest p-values
  ggplot(aes(x=term, y=over_represented_pvalue)) +
  geom_segment(aes(x=term, xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalue") +
  ggtitle(expression("Spat"~italic("S. pistillata"))) +  # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() +  # Set background color
  theme(panel.border = element_blank(),  # Set border
        panel.grid.major = element_blank(),  # Set major gridlines
        panel.grid.minor = element_blank(),  # Set minor gridlines
        axis.line = element_line(colour = "black"),  # Set axes color
        plot.background=element_blank(),  # Set the plot background
        legend.position="none")


spatSpisGO_BP %>%
  mutate(term = fct_reorder(term, over_represented_pvalue, .desc = TRUE)) %>%
  top_n(-20, over_represented_pvalue) %>%  # Select top 10 terms with smallest p-values
  ggplot(aes(x=term, y=over_represented_pvalue)) +
  geom_segment(aes(x=term, xend=term, y=0, yend=over_represented_pvalue), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("over_represented_pvalue") +
  ggtitle(expression("Spat"~italic("S. pistillata"))) +  # Add a main title with italicized species name
  theme(plot.title = element_text(face = 'bold',
                                  size = 12,
                                  hjust = 0)) +
  theme_bw() +  # Set background color
  theme(panel.border = element_blank(),  # Set border
        panel.grid.major = element_blank(),  # Set major gridlines
        panel.grid.minor = element_blank(),  # Set minor gridlines
        axis.line = element_line(colour = "black"),  # Set axes color
        plot.background=element_blank(),  # Set the plot background
        legend.position="none")



