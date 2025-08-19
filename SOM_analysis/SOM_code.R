library(readr)
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/SOM_analysis")

file.exists("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/SOM_analysis")

SOM_R <- read_csv("SOM_R_dataset.csv")


# Calculate Z-scores and store them in a new column 'Z_scores'
library(dplyr)
library(ggplot2)
SOM_R <- SOM_R %>%
  group_by(Gene) %>%  
  mutate(Z_scores = (Expression - mean(Expression)) / sd(Expression))

# View the updated dataset
head(SOM_R)

# Calculate mean Z-scores for each species at each life stage
SOM_R_summary <- SOM_R %>%
  group_by(Species, Life_stage) %>%
  summarize(mean_Z_score = mean(Z_scores), .groups = "drop")

# Create the plot
ggplot(SOM_R_summary, aes(x = Life_stage, y = mean_Z_score, color = Species, group = Species)) +
  geom_line(size = 1) +                       # Line plot
  geom_point(size = 3) +                      # Add points for clarity
  theme_minimal() +                           # Clean theme
  labs(x = "Life Stage", y = "Mean Z-score", 
       title = "Normalized Gene Expression Across Life Stages in Coral Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Load required libraries
library(ggplot2)
library(dplyr)

# Summarize the data by species and life stage
plot_data <- SOM_R %>%
  group_by(Species, Life_stage) %>%
  summarize(
    Mean_Z_score = mean(Z_scores, na.rm = TRUE),
    SE_Z_score = sd(Z_scores, na.rm = TRUE) / sqrt(n())
  )

# Create the plot
ggplot(plot_data, aes(x = Life_stage, y = Mean_Z_score, color = Species, shape = Species, group = Species)) +
  geom_point(size = 3) +  # Add points
  geom_line(linewidth = 1) +  # Add lines connecting points
  geom_errorbar(aes(ymin = Mean_Z_score - SE_Z_score, ymax = Mean_Z_score + SE_Z_score), width = 0.2) +  # Add error bars
  geom_smooth(se = FALSE, method = "loess", linetype = "dashed", linewidth = 0.8) +  # Add smooth trend lines
  scale_color_brewer(palette = "Set1") +  # Use a colorblind-friendly palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Normalized Gene Expression Across Life Stages in Coral Species",
    x = "Life Stage",
    y = "Mean Z-Score"
  )

ggplot(plot_data, aes(x = Life_stage, y = Mean_Z_score, color = Species, shape = Species, group = Species)) +
  geom_point(size = 3) +  
  geom_line(linewidth = 1) +  
  geom_errorbar(aes(ymin = Mean_Z_score - SE_Z_score, ymax = Mean_Z_score + SE_Z_score), width = 0.2) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Normalized Gene Expression Across Life Stages in Coral Species",
    x = "Life Stage",
    y = "Mean Z-Score"
  )


ggplot(plot_data, aes(x = Life_stage, y = Mean_Z_score, color = Species, shape = Species, group = Species)) +
  geom_point(size = 3) +  
  geom_line(linewidth = 1) +  
  geom_smooth(se = FALSE, method = "loess", linetype = "dashed", linewidth = 0.8) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "CARP1",
    x = "Life Stage",
    y = "Mean Z-Score"
  )

#here starts the code for each protein: 
library(readr)
library(dplyr)
library(ggplot2)


# Bring in dataset
SOM_R <- read_csv("GitHub/Multistage_omics/SOM_analysis/SOM_R_dataset.csv")


# Remove empty columns and rows
SOM_R <- SOM_R %>% select(where(~ !all(is.na(.))))
SOM_R <- SOM_R %>% filter(if_any(everything(), ~ !is.na(.)))

#Calculate Z scores for every gene 
SOM_R <- SOM_R %>%
  group_by(Protein, Gene) %>%  
  mutate(Z_scores = (Expression - mean(Expression, na.rm = TRUE)) / sd(Expression, na.rm = TRUE))

#Calculate mean expression for every life stage of every species per protein
SOM_R_summary <- SOM_R %>%
  group_by(Protein, Species, Life_stage) %>%
  summarize(mean_expression = mean(Expression, na.rm = TRUE),
            mean_Z_score = mean(Z_scores, na.rm = TRUE),
            .groups = "drop")

ggplot(SOM_R_summary, aes(x = Life_stage, y = mean_Z_score, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Protein, scales = "free") +  # Facet by Protein
  theme_minimal() +
  labs(title = "Mean Z-scores per Life Stage and Species",
       x = "Life Stage", 
       y = "Mean Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggplot2)

ggplot(SOM_R_summary, aes(x = Life_stage, y = mean_Z_score, color = Species, group = Species)) +
  geom_line(size = 1) +  # Line for each species
  geom_point(size = 3) +  # Points at each life stage
  facet_wrap(~ Protein, scales = "free_y") +  # Separate plots for each protein
  theme_minimal() +
  labs(title = "Mean Z-scores Across Life Stages by Species",
       x = "Life Stage", 
       y = "Mean Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Family

library(dplyr)
library(ggplot2)

# Remove empty columns and rows
SOM_R <- SOM_R %>% select(where(~ !all(is.na(.))))
SOM_R <- SOM_R %>% filter(if_any(everything(), ~ !is.na(.)))

# Add Family column based on Species
SOM_R <- SOM_R %>%
  mutate(Family = case_when(
    Species %in% c("Acropora_tenuis", "Montipora_capitata") ~ "Acroporidae",
    Species %in% c("Stylophora_pistillata", "Pocillopora_acuta") ~ "Pocilloporidae",
    TRUE ~ "Other"  # Optional: for species not listed
  ))

# Calculate Z scores for every gene
SOM_R <- SOM_R %>%
  group_by(Protein, Gene) %>%  
  mutate(Z_scores = (Expression - mean(Expression, na.rm = TRUE)) / sd(Expression, na.rm = TRUE))


# Calculate mean expression for every life stage of every species per protein
SOM_R_summary <- SOM_R %>%
  group_by(Protein, Species, Family, Life_stage, Gene) %>%
  summarize(mean_expression = mean(Expression, na.rm = TRUE),
            mean_Z_score = mean(Z_scores, na.rm = TRUE),
            .groups = "drop")

# Define colors manually so species from the same family share a color
family_colors <- c("Acroporidae" = "blue", "Pocilloporidae" = "red")

# Bar plot with species information but colored by family
ggplot(SOM_R_summary, aes(x = Life_stage, y = mean_Z_score, fill = Family)) +
  geom_bar(stat = "identity", position = "dodge", aes(pattern = Species)) +  # Retain species info via pattern
  facet_wrap(~ Protein, scales = "free") +  
  scale_fill_manual(values = family_colors) +  # Assign family-based colors
  theme_minimal() +
  labs(title = "Mean Z-scores per Life Stage and Species (Colored by Family)",
       x = "Life Stage", 
       y = "Mean Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Line plot with species as separate lines but colored by family
ggplot(SOM_R_summary, aes(x = Life_stage, y = mean_Z_score, color = Family, group = Species, linetype = Species)) +
  geom_line(size = 1) +  
  geom_point(size = 3) +  
  facet_wrap(~ Protein, scales = "free_y") +  
  scale_color_manual(values = family_colors) +  # Assign colors by family
  theme_minimal() +
  labs(title = "Mean Z-scores Across Life Stages by Species (Colored by Family)",
       x = "Life Stage", 
       y = "Mean Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


species_palette <- c(
  "Montipora_capitata" = "hotpink1",
  "Pocillopora_acuta" = "darkmagenta",
  "Acropora_tenuis" = "cyan4",
  "Stylophora_pistillata" = "chartreuse3"
)

ggplot(SOM_R_summary, aes(x = Life_stage, y = mean_Z_score, 
                          color = Species, group = interaction(Species, Gene))) +
  geom_line(size = 1) + 
  geom_point(size = 2) +  
  facet_wrap(~ Protein, scales = "free_y") +  
  theme_minimal() +
  labs(title = "Mean Z-scores Across Life Stages by SOM Protein-Coding Gene",
       x = "Life Stage", 
       y = "Mean Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = species_palette)  # Apply custom colors


#Now with same colour and different lines

family_palette <- c(
  "Acroporidae" = "#4dac26",
  "Pocilloporidae" = "#d01c8b"
)


species_linetypes <- c(
  "Montipora_capitata" = "solid",
  "Acropora_tenuis" = "dashed",
  "Pocillopora_acuta" = "dotted",
  "Stylophora_pistillata" = "dotdash"
)

ggplot(SOM_R_summary, aes(x = Life_stage, y = mean_Z_score, 
                          color = Family, group = interaction(Species, Gene))) +
  geom_line(size = 1, aes(linetype = Species)) +  # Different linetypes for species
  geom_point(size = 2) +  
  facet_wrap(~ Protein, scales = "free_y") +  
  theme_minimal() +
  labs(title = "Mean Z-scores Across Life Stages by SOM Protein-Coding Gene",
       x = "Life Stage", 
       y = "Mean Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = family_palette) +  # Apply family colors
  scale_linetype_manual(values = species_linetypes)  # Apply species linetypes



# Subset for selected proteins
subset_proteins <- c("MAM and LDL-receptor class A domain-containing 1", 
                     "CARP2", "Cadherin", "Vitellogenin")

SOM_R_subset <- SOM_R_summary %>% 
  dplyr::filter(Protein %in% subset_proteins)

# Plot with vertical facets
ggplot(SOM_R_subset, aes(x = Life_stage, y = mean_Z_score, 
                         color = Family, group = interaction(Species, Gene))) +
  geom_line(size = 1, aes(linetype = Species)) +
  geom_point(size = 2) +  
  facet_wrap(~ Protein, scales = "free_y", ncol = 1) +  # Vertical facets
  theme_minimal() +
  labs(title = "Mean Z-scores Across Life Stages by Selected SOM Proteins",
       x = "Life Stage", 
       y = "Mean Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = family_palette) +
  scale_linetype_manual(values = species_linetypes)




library(dplyr)
library(ggplot2)

# Define visual styles
family_palette <- c(
  "Acroporidae" = "#4dac26",
  "Pocilloporidae" = "#d01c8b"
)

species_linetypes <- c(
  "Montipora_capitata" = "solid",
  "Acropora_tenuis" = "dashed",
  "Pocillopora_acuta" = "dotted",
  "Stylophora_pistillata" = "dotdash"
)

# Define target proteins
subset_proteins <- c(
  "MAM and LDL-receptor class A domain-containing 1", 
  "CARP2", 
  "Vitellogenin",
  "Cadherin"
)

#Creating subset
SOM_R_subset <- SOM_R_summary %>%
  filter(Protein %in% subset_proteins) %>%
  mutate(
    Protein = factor(Protein, levels = subset_proteins),  # set facet order
    Life_stage = recode(
      Life_stage,
      "Larvae" = "Stage I",
      "Metamorphosed" = "Stage II",
      "Spat" = "Stage III"
    ),
    Life_stage = factor(Life_stage, levels = c("Stage I", "Stage II", "Stage III"))
  )


# Plot
plot <- ggplot(SOM_R_subset, aes(x = Life_stage, y = mean_Z_score, 
                                 color = Family, group = interaction(Species, Gene))) +
  geom_line(size = 1, aes(linetype = Species)) +
  geom_point(size = 2) +  
  facet_wrap(~ Protein, scales = "free_y", ncol = 1) +  # Vertical facets
  theme_minimal() +
  labs(title = "Mean Z-scores Across Life Stages for Selected SOM Proteins",
       x = "Life Stage", 
       y = "Mean Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = family_palette) +
  scale_linetype_manual(values = species_linetypes)

# Display
print(plot)

# Save to PNG
ggsave("SOM_subset_staged.png", plot, width = 6, height = 10, dpi = 300)

