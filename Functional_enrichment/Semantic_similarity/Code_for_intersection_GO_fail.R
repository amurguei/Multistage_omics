

#I'll try intersecting GO terms in these two, see what happens

library(simplifyEnrichment)


WGCNA_Spis_larvaeGO <- read_csv("GitHub/Multistage_omics/Functional_enrichment/Semantic_similarity/S_pistillata/Outputs/WGCNA_Spis_larvaeGO.csv")
WGCNA_Pacu_larvaeGO <- read_csv("GitHub/Multistage_omics/Functional_enrichment/Semantic_similarity/P_acuta/WGCNA_Pacu_larvaeGO.csv")
WGCNA_Mcap_larvaeGO <- read_csv("GitHub/Multistage_omics/Functional_enrichment/Semantic_similarity/M_capitata/WGCNA_Mcap_larvaeGO.csv")

# Find intersection of GO terms
common_GO_terms <- intersect(WGCNA_Spis_larvaeGO$category, WGCNA_Pacu_larvaeGO$category)

core_stage1 <- Reduce(intersect, list(WGCNA_Spis_larvaeGO$category, 
                                      WGCNA_Pacu_larvaeGO$category, 
                                      WGCNA_Mcap_larvaeGO$category))
#Nothing here

# Print the common GO terms
print(common_GO_terms)

# Install if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GOSemSim")

# Load the package
library(GOSemSim)

# Load GO database (change species if needed)
hsGO <- godata('org.Hs.eg.db', ont="BP")  # BP = Biological Process, change to "MF" or "CC" if needed

# Compute similarity between all GO terms from both datasets
# Compute semantic similarity between the two sets of GO terms
sim_matrix <- mgoSim(WGCNA_Spis_larvaeGO$category, 
                     WGCNA_Pacu_larvaeGO$category, 
                     semData=hsGO, 
                     measure="Wang", 
                     combine=NULL)  # "NULL" keeps pairwise similarity scores

# Check output
print(sim_matrix)

# Optionally, find the highest similarity values
highly_similar <- which(sim_matrix > 0.7, arr.ind = TRUE)  # Adjust threshold if needed
print(highly_similar)

#Now testing with the 3 species: 

# Load the necessary package
library(GOSemSim)
library(readr)

# Load GO enrichment results for three species
WGCNA_Spis_larvaeGO <- read_csv("GitHub/Multistage_omics/Functional_enrichment/Semantic_similarity/S_pistillata/Outputs/WGCNA_Spis_larvaeGO.csv")
WGCNA_Pacu_larvaeGO <- read_csv("GitHub/Multistage_omics/Functional_enrichment/Semantic_similarity/P_acuta/WGCNA_Pacu_larvaeGO.csv")
WGCNA_Mcap_larvaeGO <- read_csv("GitHub/Multistage_omics/Functional_enrichment/Semantic_similarity/M_capitata/WGCNA_Mcap_larvaeGO.csv")

# Load GO database (adjust species if needed)
hsGO <- godata('org.Hs.eg.db', ont="BP")  # "BP" = Biological Process

# Compute pairwise semantic similarity matrices
sim_Spis_Pacu <- mgoSim(WGCNA_Spis_larvaeGO$category, 
                        WGCNA_Pacu_larvaeGO$category, 
                        semData=hsGO, 
                        measure="Wang", 
                        combine=NULL)

sim_Spis_Mcap <- mgoSim(WGCNA_Spis_larvaeGO$category, 
                        WGCNA_Mcap_larvaeGO$category, 
                        semData=hsGO, 
                        measure="Wang", 
                        combine=NULL)

sim_Pacu_Mcap <- mgoSim(WGCNA_Pacu_larvaeGO$category, 
                        WGCNA_Mcap_larvaeGO$category, 
                        semData=hsGO, 
                        measure="Wang", 
                        combine=NULL)

# Find highly similar GO terms (threshold can be adjusted)
threshold <- 0.7

similar_Spis_Pacu <- which(sim_Spis_Pacu > threshold, arr.ind = TRUE)
similar_Spis_Mcap <- which(sim_Spis_Mcap > threshold, arr.ind = TRUE)
similar_Pacu_Mcap <- which(sim_Pacu_Mcap > threshold, arr.ind = TRUE)

print(similar_Spis_Pacu)  # For S. pistillata and P. acuta
print(similar_Spis_Mcap)  # For S. pistillata and M. capitata
print(similar_Pacu_Mcap)  # For P. acuta and M. capitata


# Extract the actual GO terms that meet the threshold
core_GO_terms <- intersect(
  intersect(WGCNA_Spis_larvaeGO$category[similar_Spis_Pacu[, 1]], 
            WGCNA_Spis_larvaeGO$category[similar_Spis_Mcap[, 1]]),
  WGCNA_Pacu_larvaeGO$category[similar_Pacu_Mcap[, 1]]
)

# Print core GO terms shared across all three species
print(core_GO_terms)




# Load GO database (adjust species)
hsGO <- godata('org.Hs.eg.db', ont="BP") 

# Combine all GO terms from the three species
all_GO_terms <- unique(c(WGCNA_Spis_larvaeGO$category, 
                         WGCNA_Pacu_larvaeGO$category, 
                         WGCNA_Mcap_larvaeGO$category))

all_GO_terms <- as.character(all_GO_terms)  # Ensure it's a character vector


# Compute similarity matrix
go_sim_matrix <- GO_similarity(all_GO_terms, measure = "Wang", sem_data = hsGO)

# Perform clustering & visualize
simplifyGO(go_sim_matrix)







library(GOSemSim)

# Load GO database (adjust species)
hsGO <- godata('org.Hs.eg.db', ont = "BP") 

# Combine all GO terms from the three species
all_GO_terms <- unique(c(WGCNA_Spis_larvaeGO$category, 
                         WGCNA_Pacu_larvaeGO$category, 
                         WGCNA_Mcap_larvaeGO$category))

all_GO_terms <- as.character(all_GO_terms)  # Ensure it's a character vector

# Compute pairwise similarity for all GO terms
go_sim_matrix <- mgoSim(all_GO_terms, all_GO_terms, semData = hsGO, measure = "Wang", combine = NULL)

# Convert to a distance matrix (for clustering)
go_dist <- as.dist(1 - go_sim_matrix)  

# Perform hierarchical clustering
hclust_obj <- hclust(go_dist, method = "average")

# Cut the tree to define clusters (adjust k as needed)
clusters <- cutree(hclust_obj, k = 10)  # Adjust k based on visualization

# Assign clusters to GO terms
clustered_GO <- data.frame(GO_term = all_GO_terms, cluster = clusters)

# Check clusters that contain GO terms from all species
library(dplyr)
shared_clusters <- clustered_GO %>%
  group_by(cluster) %>%
  filter(any(GO_term %in% WGCNA_Spis_larvaeGO$category) & 
           any(GO_term %in% WGCNA_Pacu_larvaeGO$category) & 
           any(GO_term %in% WGCNA_Mcap_larvaeGO$category))

# Extract representative GO terms for shared clusters
core_GO_terms <- unique(shared_clusters$GO_term)

print(core_GO_terms)  # Final core GO terms shared across species

# Save GO terms for REVIGO
write.table(core_GO_terms, file = "GO_terms_for_REVIGO.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

getwd()


library(simplifyEnrichment)

# Compute similarity again didn´t like this one...
go_sim_matrix <- mgoSim(all_GO_terms, all_GO_terms, semData = hsGO, measure = "Wang", combine = NULL)

# Run `simplifyEnrichment`
simplifyEnrichment(go_sim_matrix, method = "binary_cut")



#Trying with Bingo
# Install Bingo from Bioconductor (if you haven't already)
if (!requireNamespace("Bingo", quietly = TRUE)) {
  BiocManager::install("Bingo")
}
# Load Bingo
library(Bingo)
