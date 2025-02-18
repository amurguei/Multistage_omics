# 18 November 2024
#Based on the work of Matthew Haas https://github.com/MatthewHaas/genome_assembly/blob/master/venn_diagram_orthogroups_with_rice_relatives.R#L1
# Purpose of this code is to make a venn diagram for overlapping Orthogroups for the NWR genome paper
# This script uses the updated OrthoFinder results that includes Z. latifolia and other Oryza species.
# Done in RStudio

# Set working directory
setwd("C:/Users/amurgueitio/Documents/Multistage_omics/Orthogroups4sp")

# Load required packages
library(data.table)
library(VennDiagram)
# Since we wanted to add commas as the separator for thousands, I modified the venn.diagram function.
# We need to use source() to load the modified function.
source("C:/Users/amurgueitio/Documents/GitHub/Multistage_omics/Orthogroups/Scripts/modified_venn_diagram.R")

# Read in data
x <- fread("C:/Users/amurgueitio/Documents/Multistage_Omics/Orthogroups4sp/Orthogroups.tsv")

# Find all rows for which each species has at least one gene present in the orthogroup
Montipora_capitata_HIv3.genes.pep <- x[Montipora_capitata_HIv3.genes.pep > 0]
Pocillopora_acuta_HIv2.genes.pep <- x[Pocillopora_acuta_HIv2.genes.pep > 0]
Spis.genome.annotation.pep.longest <- x[Spis.genome.annotation.pep.longest > 0]
Acropora_tenuis_0.11.maker_post_001.proteins <- x[Acropora_tenuis_0.11.maker_post_001.proteins >0]

# Get the orthogroup IDS for each species.
Montipora_capitata_HIv3.genes.pep_orthogroups <- Montipora_capitata_HIv3.genes.pep$Orthogroup
Pocillopora_acuta_HIv2.genes.pep_orthogroups <- Pocillopora_acuta_HIv2.genes.pep$Orthogroup
Spis.genome.annotation.pep.longest_orthogroups <- Spis.genome.annotation.pep.longest$Orthogroup
Acropora_tenuis_0.11.maker_post_001.proteins_orthogroups <- Acropora_tenuis_0.11.maker_post_001.proteins$Orthogroup

# Save data
save(Montipora_capitata_HIv3.genes.pep_orthogroups, Pocillopora_acuta_HIv2.genes.pep_orthogroups, Spis.genome.annotation.pep.longest_orthogroups, Acropora_tenuis_0.11.maker_post_001.proteins_orthogroups, file="orthogroup_names_for_venn_diagram_update_trial.Rdata")

colors= c("hotpink1", "darkmagenta", "cyan4","steelblue3")

# Load the necessary library
library(VennDiagram)

# Define the orthogroup data for the three species
venn.diagram(
  x = list(
    Montipora_capitata = Montipora_capitata_HIv3.genes.pep_orthogroups,
    Pocillopora_acuta = Pocillopora_acuta_HIv2.genes.pep_orthogroups,
    Stylophora_pistillata = Spis.genome.annotation.pep.longest_orthogroups,
    Acropora_tenuis = Acropora_tenuis_0.11.maker_post_001.proteins_orthogroups
  ),
  category.names = c(
    expression(atop(italic("Montipora capitata"), plain(16410))),
    expression(atop(italic("Pocillopora acuta"), plain(15456))),
    expression(atop(italic("Stylophora pistillata"), plain(13666))),
    expression(atop(italic("Acropora tenuis"), plain(13324)))
  ),
  filename = "venn_diagram_coral_species.png",
  output = TRUE,
  imagetype = "png",
  scaled = FALSE,
  col = "black",
  fill = colors,
  cat.col = colors,
  cat.cex = 0.95,
  cat.dist = 0.35,
  margin = 0.30,
  euler.d = FALSE
)


#Another version with numbers and bigger circles
colors = c("hotpink1", "darkmagenta", "cyan4", "steelblue3")

# Ensure the input data for orthogroups are treated as lists
Montipora_capitata_orthogroups <- unique(Montipora_capitata_HIv3.genes.pep_orthogroups)
Pocillopora_acuta_orthogroups <- unique(Pocillopora_acuta_HIv2.genes.pep_orthogroups)
Stylophora_pistillata_orthogroups <- unique(Spis.genome.annotation.pep.longest_orthogroups)
Acropora_tenuis_0.11.maker_post_001.proteins_orthogroups <- unique(Acropora_tenuis_0.11.maker_post_001.proteins_orthogroups)

# Check lengths
length_Montipora <- length(Montipora_capitata_orthogroups)
length_Pocillopora <- length(Pocillopora_acuta_orthogroups)
length_Stylophora <- length(Stylophora_pistillata_orthogroups)
length_Acropora <- length(Acropora_tenuis_0.11.maker_post_001.proteins_orthogroups)

# Generate the Venn diagram
venn.diagram(
  x = list(
    Montipora_capitata = Montipora_capitata_orthogroups,
    Pocillopora_acuta = Pocillopora_acuta_orthogroups,
    Stylophora_pistillata = Stylophora_pistillata_orthogroups
  ),
  category.names = c(
    expression(atop(italic("Montipora capitata"), plain(format(length_Montipora, big.mark = ",")))),
    expression(atop(italic("Pocillopora acuta"), plain(format(length_Pocillopora, big.mark = ",")))),
    expression(atop(italic("Stylophora pistillata"), plain(format(length_Stylophora, big.mark = ",")))),
    expression(atop(italic("Acropora tenuis"), plain(format(length_Acropora, big.mark = ","))))
  ),
  filename = "venn_diagram_mod.png",
  output = TRUE,
  imagetype = "png",
  scaled = FALSE,
  col = "black",
  fill = colors,
  cat.col = colors,
  cat.cex = 0.95,
  cat.dist = 0.35,
  margin = 0.30,
  euler.d = FALSE
)


str(Montipora_capitata_HIv3.genes.pep_orthogroups)
str(Pocillopora_acuta_HIv2.genes.pep_orthogroups)
str(Spis.genome.annotation.pep.longest_orthogroups)
