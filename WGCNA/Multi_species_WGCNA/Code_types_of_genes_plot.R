# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)

# Set working directory
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
install.packages("ggrepel")  # Only run if not installed
library(ggrepel)


# Load your data
df <- read_excel("Single_copy_genes_classified.xlsx")

# Rename "No" to "All_other_genes"

df <- df %>%
  mutate(SOMP_Status = case_when(
    SOMP_Status == "No" ~ "All_other_genes",
    SOMP_Status == "SOMP" ~ "Putative_SOMPs",
    TRUE ~ SOMP_Status
  ))
# Count the number of genes per SOMP_status
status_counts <- df %>%
  count(SOMP_Status)

# Create the pie chart
# Set working directory
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")

library(ggplot2)
library(dplyr)
library(readxl)
library(scales)
library(ggrepel)

# Load data
somp_data <- read_excel("Single_copy_genes_classified.xlsx")

# Define colors
color_palette <- c(
  "All_other_genes" = "#b8e186",
  "Putative_SOMPs" = "#d01c8b",
  "Other_Biomin_toolkit" = "#f1b6da",
  "Unsure" = "#4dac26"
)

# Define readable labels for the legend
label_names <- c(
  "All_other_genes" = "All other genes",
  "Putative_SOMPs" = "Putative SOMPs",
  "Other_Biomin_toolkit" = "Other biomineralization genes",
  "Unsure" = "Unsure"
)

# Prepare the data
somp_counts <- somp_data %>%
  mutate(SOMP_Status = case_when(
    SOMP_Status == "No" ~ "All_other_genes",
    SOMP_Status == "SOMP" ~ "Putative_SOMPs",
    TRUE ~ SOMP_Status
  )) %>%
  count(SOMP_Status) %>%
  mutate(perc = n / sum(n),
         label = paste0(scales::percent(perc, accuracy = 0.1)),
         Category = label_names[SOMP_Status])

# Calculate positions
somp_counts <- somp_counts %>%
  arrange(desc(SOMP_Status)) %>%
  mutate(ypos = cumsum(perc) - 0.5 * perc)

# Plot with only percentage labels (fixed alignment)
ggplot(somp_counts, aes(x = "", y = perc, fill = SOMP_Status)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(y = ypos, label = label),
                   size = 4,
                   nudge_x = 1,
                   show.legend = FALSE,
                   direction = "y",
                   segment.size = 0.5,
                   box.padding = 0.5,
                   xlim = c(1.2, NA)) +
  scale_fill_manual(values = color_palette, labels = label_names, name = "Gene Category") +
  theme_void() +
  ggtitle("Gene Proportions by SOMP Category")







#SOMP + Other Biomin genes

library(readxl)
library(dplyr)
library(ggplot2)

# Load data
data <- read_excel("Single_copy_genes_classified.xlsx")

# Clean Type_S_B
data <- data %>%
  mutate(Type_S_B = trimws(Type_S_B)) %>% 
  filter(!is.na(Type_S_B))

# Standardize SOMP labels
data <- data %>%
  mutate(SOMP_Status = recode(SOMP_Status,
                              "SOMP" = "Putative_SOMPs",
                              "No" = "All_other_genes"))

# 1. Pie chart for Putative_SOMPs
somp_data <- data %>%
  filter(SOMP_Status == "Putative_SOMPs") %>%
  count(Type_S_B) %>%
  mutate(perc = n / sum(n),
         label = paste0(scales::percent(perc)))

ggplot(somp_data, aes(x = "", y = perc, fill = Type_S_B)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 4) +
  theme_void() +
  ggtitle("Putative SOMPs – Type Distribution")

# 2. Pie chart for Other_Biomin_toolkit
biomin_data <- data %>%
  filter(SOMP_Status == "Other_Biomin_toolkit") %>%
  count(Type_S_B) %>%
  mutate(perc = n / sum(n),
         label = paste0(scales::percent(perc)))

ggplot(biomin_data, aes(x = "", y = perc, fill = Type_S_B)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 4) +
  theme_void() +
  ggtitle("Other Biomineralization Toolkit – Type Distribution")


library(tidyverse)
library(ggrepel)
library(readxl)

# Set your working directory
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")

# Load data
data <- read_excel("Single_copy_genes_classified.xlsx")

# Filter to only Putative SOMPs
somp_data <- data %>%
  filter(SOMP_Status == "Putative_SOMPs")

# Clean and merge Type_S_B categories
somp_data <- somp_data %>%
  mutate(Type_S_B = case_when(
    Type_S_B %in% c("CA/ Mucin", "CA_Mucin") ~ "Carbonic anhydrase / Mucin",
    Type_S_B %in% c("Sushi_Von Willebrand_EGF", "Sushi, Von Willebrand, EGF") ~ "Sushi, Von Willebrand, EGF",
    TRUE ~ Type_S_B
  ))

# Count and calculate percentages
type_counts <- somp_data %>%
  count(Type_S_B) %>%
  mutate(perc = n / sum(n),
         label = paste0(scales::percent(perc, accuracy = 0.1)),
         ypos = cumsum(perc) - 0.5 * perc)

# Plot
ggplot(type_counts, aes(x = "", y = perc, fill = Type_S_B)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_label_repel(aes(y = ypos, label = label), 
                   size = 4, show.legend = FALSE, nudge_x = 1.5, segment.size = 0.5) +
  labs(title = "Putative SOMPs – Type Distribution", fill = "Gene Type") +
  theme_void()


library(tidyverse)
library(ggrepel)
library(readxl)

# Load your data
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")
data <- read_excel("Single_copy_genes_classified.xlsx")

# Filter for Putative SOMPs, remove NA and standardize category names
somp_data <- data %>%
  filter(SOMP_Status == "SOMP", !is.na(Type_S_B)) %>%
  mutate(Type_S_B = case_when(
    Type_S_B %in% c("CA/ Mucin", "CA_Mucin") ~ "Carbonic anhydrase / Mucin",
    Type_S_B %in% c("Sushi_Von Willebrand_EGF", "Sushi, Von Willebrand, EGF") ~ "Sushi, Von Willebrand, EGF",
    TRUE ~ Type_S_B
  ))

# Summarize counts and calculate positions
type_counts <- somp_data %>%
  count(Type_S_B) %>%
  mutate(perc = n / sum(n),
         label = paste0(scales::percent(perc, accuracy = 0.1))) %>%
  arrange(desc(perc)) %>%
  mutate(ypos = cumsum(perc) - 0.5 * perc)

# Plot
ggplot(type_counts, aes(x = "", y = perc, fill = Type_S_B)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 4) +
  labs(title = "Putative SOMPs – Type Distribution", fill = "Gene Type") +
  theme_void()


library(readxl)
library(dplyr)
library(ggplot2)

# Distinct and color-blind-friendly palette
# Color-blind-friendly palette
color_palette <- c(
  "BMP" = "#E69F00",         # Orange
  "Tolloid_like" = "#56B4E9", # Sky Blue
  "Wnt_Protein" = "#009E73"   # Bluish Green
)



# Load data
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")
data <- read_excel("Single_copy_genes_classified.xlsx")

# Filter for Other_Biomin_toolkit genes and clean categories
biomin_data <- data %>%
  filter(SOMP_Status == "Other_Biomin_toolkit", !is.na(Type_S_B)) %>%
  mutate(Type_S_B = case_when(
    Type_S_B %in% c("CA/ Mucin", "CA_Mucin") ~ "Carbonic anhydrase / Mucin",
    Type_S_B %in% c("Sushi_Von Willebrand_EGF", "Sushi, Von Willebrand, EGF") ~ "Sushi, Von Willebrand, EGF",
    TRUE ~ Type_S_B
  ))

# Summarize counts and compute percentages
type_counts <- biomin_data %>%
  count(Type_S_B) %>%
  mutate(perc = n / sum(n),
         label = paste0(scales::percent(perc, accuracy = 0.1)))

# Plot
ggplot(type_counts, aes(x = "", y = perc, fill = Type_S_B)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 4) +
  scale_fill_manual(values = color_palette) +
  labs(title = "Other Biomineralization Genes – Type Distribution", fill = "Gene Type") +
  theme_void()





library(tidyverse)

# Load data
gcount_annot <- read_csv("gcount_raw_with_annotations.csv")

# Define species keywords to remove from descriptions
species_keywords <- c("homo sapiens", "mus musculus", "gallus gallus", 
                      "xenopus", "danio rerio", "rattus norvegicus", 
                      "nematostella vectensis", "strongylocentrotus purpuratus",
                      "branchiostoma floridae", "amphimedon queenslandica", 
                      "ciona intestinalis")

# Clean hit descriptions directly in the data frame
gcount_annot <- gcount_annot %>%
  mutate(
    description_clean = str_to_lower(`Hit description`)
  )

# Remove species names
for (species in species_keywords) {
  gcount_annot$description_clean <- str_replace_all(gcount_annot$description_clean, fixed(species), "")
}

# Define keyword-to-category mapping as a named vector
category_mapping <- setNames(
  c(
    "Signaling (Kinases)",
    "Signaling (Receptors)",
    "Binding Proteins",
    "Transcription Regulation",
    "Metabolism (Mitochondrial)",
    "Protein Degradation",
    "Protein Degradation",
    "Protein Complex Components",
    "Protein Complex Components",
    "Structural Proteins",
    "Structural Proteins",
    "Structural Proteins",
    "Transcription Regulation",
    "Homologous Proteins",
    "Homologous Proteins",
    "General/Associated Functions",
    "Unknown/Predicted",
    "Unknown/Predicted",
    "General Proteins"
  ),
  c(
    "kinase",
    "receptor",
    "binding",
    "transcription",
    "mitochondrial",
    "ubiquitin",
    "ligase",
    "subunit",
    "complex",
    "domain",
    "containing",
    "repeat",
    "factor",
    "homolog",
    "like",
    "associated",
    "predicted",
    "uncharacterized",
    "protein"
  )
)

# Assign category based on keyword matching
gcount_annot$category <- sapply(gcount_annot$description_clean, function(d) {
  keyword <- names(category_mapping)[which(str_detect(d, names(category_mapping)))]
  if (length(keyword) > 0) category_mapping[keyword[1]] else "Other"
})

# Extract species from GeneID_Spis
gcount_annot <- gcount_annot %>%
  mutate(Species = case_when(
    str_detect(GeneID_Spis, "Spis") ~ "S. pistillata",
    str_detect(GeneID_Spis, "Pacu") ~ "P. acuta",
    str_detect(GeneID_Spis, "Mcap") ~ "M. capitata",
    str_detect(GeneID_Spis, "Aten") ~ "A. tenuis",
    TRUE ~ "Unknown"
  ))

# Summarize: count of categories per species
species_category_counts <- gcount_annot %>%
  group_by(Species, category) %>%
  summarize(count = n(), .groups = "drop")

# Plot
ggplot(species_category_counts, aes(x = category, y = count, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    y = "Gene Count",
    x = "Functional Category",
    title = "Gene Functional Profiles by Category and Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
library(tidyverse)

# Load data
gcount_annot <- read_csv("gcount_raw_with_annotations.csv")

# Define species keywords to remove from descriptions
species_keywords <- c("homo sapiens", "mus musculus", "gallus gallus", 
                      "xenopus", "danio rerio", "rattus norvegicus", 
                      "nematostella vectensis", "strongylocentrotus purpuratus",
                      "branchiostoma floridae", "amphimedon queenslandica", 
                      "ciona intestinalis")

# Clean descriptions and remove species names
gcount_annot <- gcount_annot %>%
  mutate(
    description_clean = str_to_lower(`Hit description`)
  )

for (species in species_keywords) {
  gcount_annot$description_clean <- str_replace_all(
    gcount_annot$description_clean,
    fixed(species),
    ""
  )
}

# Define keyword-to-category mapping (named vector)
category_mapping <- setNames(
  c(
    "Signaling (Kinases)", "Signaling (Receptors)", "Binding Proteins",
    "Transcription Regulation", "Metabolism (Mitochondrial)",
    "Protein Degradation", "Protein Degradation",
    "Protein Complex Components", "Protein Complex Components",
    "Structural Proteins", "Structural Proteins", "Structural Proteins",
    "Transcription Regulation", "Homologous Proteins", "Homologous Proteins",
    "General/Associated Functions", "Unknown/Predicted", "Unknown/Predicted",
    "General Proteins"
  ),
  c(
    "kinase", "receptor", "binding",
    "transcription", "mitochondrial",
    "ubiquitin", "ligase",
    "subunit", "complex",
    "domain", "containing", "repeat",
    "factor", "homolog", "like",
    "associated", "predicted", "uncharacterized",
    "protein"
  )
)

# Order categories by total count (for cleaner plot)
species_category_counts <- species_category_counts %>%
  mutate(category = fct_reorder(category, count, .fun = sum, .desc = TRUE))

# Plot: horizontal side-by-side bar plot
ggplot(species_category_counts, aes(x = category, y = count, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Functional Category",
    y = "Gene Count",
    title = "Gene Functional Profiles by Category and Species"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))





library(tidyverse)

# Load data
gcount_annot <- read_csv("gcount_raw_with_annotations.csv")

# Species names to remove from descriptions
species_keywords <- c("homo sapiens", "mus musculus", "gallus gallus", 
                      "xenopus", "danio rerio", "rattus norvegicus", 
                      "nematostella vectensis", "strongylocentrotus purpuratus",
                      "branchiostoma floridae", "amphimedon queenslandica", 
                      "ciona intestinalis")

# Clean descriptions
gcount_annot <- gcount_annot %>%
  mutate(description_clean = str_to_lower(`Hit description`))

for (species in species_keywords) {
  gcount_annot$description_clean <- str_replace_all(
    gcount_annot$description_clean,
    fixed(species),
    ""
  )
}

# Keyword-to-category mapping
category_mapping <- setNames(
  c(
    "Signaling (Kinases)", "Signaling (Receptors)", "Binding Proteins",
    "Transcription Regulation", "Metabolism (Mitochondrial)",
    "Protein Degradation", "Protein Degradation",
    "Protein Complex Components", "Protein Complex Components",
    "Structural Proteins", "Structural Proteins", "Structural Proteins",
    "Transcription Regulation", "Homologous Proteins", "Homologous Proteins",
    "General/Associated Functions", "Unknown/Predicted", "Unknown/Predicted",
    "General Proteins"
  ),
  c(
    "kinase", "receptor", "binding",
    "transcription", "mitochondrial",
    "ubiquitin", "ligase",
    "subunit", "complex",
    "domain", "containing", "repeat",
    "factor", "homolog", "like",
    "associated", "predicted", "uncharacterized",
    "protein"
  )
)

# Assign one category per description
gcount_annot$category <- sapply(gcount_annot$description_clean, function(d) {
  keyword <- names(category_mapping)[which(str_detect(d, names(category_mapping)))]
  if (length(keyword) > 0) category_mapping[keyword[1]] else "Other"
})


  mutate(category = fct_inorder(category))  # preserves the arranged order

# Count and reorder by count descending
category_counts <- gcount_annot %>%
  count(category) %>%
  mutate(category = fct_reorder(category, n, .desc = TRUE))  # FIXED ORDERING

# Plot
ggplot(category_counts, aes(x = category, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Functional Category",
    y = "Number of Genes",
    title = "Spis-Annotated Orthogroup Functional Categories"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))


library(tidyverse)

# Load data
gcount_annot <- read_csv("gcount_raw_with_annotations.csv")

# Species names to remove from descriptions
species_keywords <- c("homo sapiens", "mus musculus", "gallus gallus", 
                      "xenopus", "danio rerio", "rattus norvegicus", 
                      "nematostella vectensis", "strongylocentrotus purpuratus",
                      "branchiostoma floridae", "amphimedon queenslandica", 
                      "ciona intestinalis")

# Clean descriptions
gcount_annot <- gcount_annot %>%
  mutate(description_clean = str_to_lower(`Hit description`))

for (species in species_keywords) {
  gcount_annot$description_clean <- str_replace_all(
    gcount_annot$description_clean,
    fixed(species),
    ""
  )
}

# Keyword-to-category mapping
category_mapping <- setNames(
  c(
    "Signaling (Kinases)", "Signaling (Receptors)", "Binding Proteins",
    "Transcription Regulation", "Metabolism (Mitochondrial)",
    "Protein Degradation", "Protein Degradation",
    "Protein Complex Components", "Protein Complex Components",
    "Structural Proteins", "Structural Proteins", "Structural Proteins",
    "Transcription Regulation", "Homologous Proteins", "Homologous Proteins",
    "General/Associated Functions", "Unknown/Predicted", "Unknown/Predicted",
    "General Proteins"
  ),
  c(
    "kinase", "receptor", "binding",
    "transcription", "mitochondrial",
    "ubiquitin", "ligase",
    "subunit", "complex",
    "domain", "containing", "repeat",
    "factor", "homolog", "like",
    "associated", "predicted", "uncharacterized",
    "protein"
  )
)

# Assign one category per description
gcount_annot$category <- sapply(gcount_annot$description_clean, function(d) {
  keyword <- names(category_mapping)[which(str_detect(d, names(category_mapping)))]
  if (length(keyword) > 0) category_mapping[keyword[1]] else "Other"
})

# Count and manually set factor levels in descending order
category_counts <- gcount_annot %>%
  count(category) %>%
  arrange(desc(n)) %>%
  mutate(category = factor(category, levels = unique(category)))

# Plot
ggplot(category_counts, aes(x = category, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Functional Category",
    y = "Number of Genes",
    title = "Spis-Annotated Orthogroup Functional Categories"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))



library(tidyverse)

# Load data
gcount_annot <- read_csv("gcount_raw_with_annotations.csv")

# Species names to remove from descriptions
species_keywords <- c("homo sapiens", "mus musculus", "gallus gallus", 
                      "xenopus", "danio rerio", "rattus norvegicus", 
                      "nematostella vectensis", "strongylocentrotus purpuratus",
                      "branchiostoma floridae", "amphimedon queenslandica", 
                      "ciona intestinalis")

# Clean descriptions
gcount_annot <- gcount_annot %>%
  mutate(description_clean = str_to_lower(`Hit description`))

for (species in species_keywords) {
  gcount_annot$description_clean <- str_replace_all(
    gcount_annot$description_clean,
    fixed(species),
    ""
  )
}

# Keyword-to-category mapping
category_mapping <- setNames(
  c(
    "Signaling (Kinases)", "Signaling (Receptors)", "Binding Proteins",
    "Transcription Regulation", "Metabolism (Mitochondrial)",
    "Protein Degradation", "Protein Degradation",
    "Protein Complex Components", "Protein Complex Components",
    "Structural Proteins", "Structural Proteins", "Structural Proteins",
    "Transcription Regulation", "Homologous Proteins", "Homologous Proteins",
    "General/Associated Functions", "Unknown/Predicted", "Unknown/Predicted",
    "General Proteins"
  ),
  c(
    "kinase", "receptor", "binding",
    "transcription", "mitochondrial",
    "ubiquitin", "ligase",
    "subunit", "complex",
    "domain", "containing", "repeat",
    "factor", "homolog", "like",
    "associated", "predicted", "uncharacterized",
    "protein"
  )
)

# Assign one category per description
gcount_annot <- gcount_annot %>%
  rowwise() %>%
  mutate(
    category = {
      keyword <- names(category_mapping)[str_detect(description_clean, names(category_mapping))]
      if (length(keyword) > 0) category_mapping[[keyword[1]]] else "Other"
    }
  ) %>%
  ungroup()

# Count categories
category_counts <- gcount_annot %>%
  count(category) %>%
  arrange(desc(n))

# Force factor levels based on sorted counts
category_counts <- category_counts %>%
  mutate(category = factor(category, levels = category))

# Plot
ggplot(category_counts, aes(x = category, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Functional Category",
    y = "Number of Genes",
    title = "Single-Copy Orthologues by Categories"
  ) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1))

