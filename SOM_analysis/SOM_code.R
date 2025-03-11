SOM_R <- read_csv("C:/Users/amurg/Downloads/SOM_R_dataset - Sheet1.csv")

# Calculate Z-scores and store them in a new column 'Z_scores'
library(dplyr)

SOM_R <- SOM_R %>%
  group_by(Gene) %>%  
  mutate(Z_scores = (Expression - mean(Expression)) / sd(Expression))

# View the updated dataset
head(SOM_R)
