# 1. Load required libraries
library(edgeR)
library(limma)
library(variancePartition)
library(tidyverse)
library(genefilter)

# 2. Load expression data
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs")
gcount <- read.csv("Mcap_transcript_count_matrix.csv", row.names = "gene_id")
treatmentinfo <- read_csv("5-Mcap-SampleInfo.csv")

# 3. Filter low-expression genes
filt <- filterfun(pOverA(0.33, 10))  # At least 10 counts in 3/9 samples
gfilt <- genefilter(gcount, filt)
gcount_filt <- gcount[gfilt, ]

# 4. Add additional metadata

## 4a. Alignment rate
alignment_table <- tribble(
  ~SampleID, ~SRR, ~AlignRate,
  "AH1", "SRR14864072", 73.0,
  "AH2", "SRR14864071", 80.8,
  "AH3", "SRR14864070", 82.1,
  "AH4", "SRR14864069", 79.9,
  "AH5", "SRR14864068", 71.9,
  "AH6", "SRR14864067", 73.0,
  "AH7", "SRR14864066", 47.3,
  "AH8", "SRR14864065", 65.0,
  "AH9", "SRR14864064", 63.4
)
treatmentinfo <- left_join(treatmentinfo, alignment_table, by = c("sampleID" = "SampleID"))

## 4b. TIN scores
tin_data <- read.table("multiqc_tin.txt", header = TRUE)
colnames(tin_data) <- c("SRR", "TIN")
tin_data <- left_join(tin_data, alignment_table, by = "SRR")
treatmentinfo <- left_join(treatmentinfo, tin_data %>% select(SampleID, TIN),
                           by = c("sampleID" = "SampleID"))

## 4c. RNA metrics
rna_quality <- tribble(
  ~SampleID, ~TotalRNA, ~RIN, ~RNA_260_230, ~RNA_260_280,
  "AH1", 28.7, 8.4, 0.3, 1.88,
  "AH2", 45.4, 9.2, 0.37, 2.06,
  "AH3", 27.8, 8.5, 0.9, 1.75,
  "AH4", 32.3, 8.6, 1.16, 1.85,
  "AH5", 32.3, 8.9, 0.65, 1.90,
  "AH6", 30.8, 8.7, 1.18, 1.87,
  "AH7", 26.4, 6.1, 1.05, 1.87,
  "AH8", 56.3, 7.4, 1.28, 1.96,
  "AH9", 41.8, 6.8, 0.22, 1.99
)
treatmentinfo <- left_join(treatmentinfo, rna_quality, by = c("sampleID" = "SampleID"))

# 5. Scale technical covariates
treatmentinfo <- treatmentinfo %>%
  mutate(
    TIN_scaled = scale(TIN),
    RIN_scaled = scale(RIN),
    RNA_260_230_scaled = scale(RNA_260_230),
    RNA_260_280_scaled = scale(RNA_260_280),
    TotalRNA_scaled = scale(TotalRNA),
    AlignRate_scaled = scale(AlignRate)
  )

# 6. Set rownames of metadata to match columns of expression data
stopifnot(all(colnames(gcount_filt) == treatmentinfo$sampleID))
treatmentinfo <- as.data.frame(treatmentinfo)
rownames(treatmentinfo) <- treatmentinfo$sampleID

# 7. Ensure timepoint is a factor with >1 level
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I", "II", "III"))
stopifnot(length(levels(treatmentinfo$timepoint)) > 1)

# 8. Create DGEList object and normalize
dge <- DGEList(counts = gcount_filt)
dge <- calcNormFactors(dge)
dge$counts
dege as dataframe

# 9. Filter genes based on expression
fixed_formula <- ~ timepoint + TIN_scaled + RNA_260_230_scaled + RNA_260_280_scaled + TotalRNA_scaled
keep <- filterByExpr(dge, design = model.matrix(fixed_formula, data = treatmentinfo))
dge <- dge[keep, ]

# 10. Define formula with random effect
form <- ~ (1|timepoint) + TIN_scaled + RNA_260_230_scaled + RNA_260_280_scaled + TotalRNA_scaled

# 11. Run voom with dream weights
vobj <- voomWithDreamWeights(dge, formula = form, data = treatmentinfo)

# 12. Fit dream model (use fully qualified function to avoid ambiguity)
fit <- variancePartition::dream(vobj, formula = form, data = treatmentinfo)

# 13. Check class to confirm it's a Dream object
print(class(fit))  # should include "Dream"

# 14. Extract and plot variance fractions
vp <- fitExtractVarPartModel(fit)
plotVarPart(vp)


---------------------------------------------------------------------------------



# 1. Load required libraries
library(edgeR)
library(limma)
library(variancePartition)
library(tidyverse)
library(genefilter)


setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs")

# 2. Load expression data
gcount <- read.csv("Mcap_transcript_count_matrix.csv", row.names = "gene_id")
treatmentinfo <- read_csv("5-Mcap-SampleInfo.csv")

# 3. Filter low-expression genes
filt <- filterfun(pOverA(0.33, 10))  # At least 10 counts in 3/9 samples
gfilt <- genefilter(gcount, filt)
gcount_filt <- gcount[gfilt, ]

# 4. Add additional metadata: alignment rate, TIN, RNA quality metrics

# 4a. Add alignment rate
alignment_table <- tribble(
  ~SampleID, ~SRR, ~AlignRate,
  "AH1", "SRR14864072", 73.0,
  "AH2", "SRR14864071", 80.8,
  "AH3", "SRR14864070", 82.1,
  "AH4", "SRR14864069", 79.9,
  "AH5", "SRR14864068", 71.9,
  "AH6", "SRR14864067", 73.0,
  "AH7", "SRR14864066", 47.3,
  "AH8", "SRR14864065", 65.0,
  "AH9", "SRR14864064", 63.4
)

treatmentinfo <- left_join(treatmentinfo, alignment_table, by = c("sampleID" = "SampleID"))

# 4b. Add TIN
tin_data <- read.table("multiqc_tin.txt", header = TRUE)
colnames(tin_data) <- c("SRR", "TIN")

tin_data <- left_join(tin_data, alignment_table, by = "SRR")
treatmentinfo <- left_join(treatmentinfo, tin_data %>% select(SampleID, TIN),
                           by = c("sampleID" = "SampleID"))

# 4c. Add RNA metrics
rna_quality <- tribble(
  ~SampleID, ~TotalRNA, ~RIN, ~RNA_260_230, ~RNA_260_280,
  "AH1", 28.7, 8.4, 0.3, 1.88,
  "AH2", 45.4, 9.2, 0.37, 2.06,
  "AH3", 27.8, 8.5, 0.9, 1.75,
  "AH4", 32.3, 8.6, 1.16, 1.85,
  "AH5", 32.3, 8.9, 0.65, 1.90,
  "AH6", 30.8, 8.7, 1.18, 1.87,
  "AH7", 26.4, 6.1, 1.05, 1.87,
  "AH8", 56.3, 7.4, 1.28, 1.96,
  "AH9", 41.8, 6.8, 0.22, 1.99
)

treatmentinfo <- left_join(treatmentinfo, rna_quality, by = c("sampleID" = "SampleID"))

# 5. Scale technical covariates
treatmentinfo <- treatmentinfo %>%
  mutate(
    TIN_scaled = scale(TIN),
    RIN_scaled = scale(RIN),
    RNA_260_230_scaled = scale(RNA_260_230),
    RNA_260_280_scaled = scale(RNA_260_280),
    TotalRNA_scaled = scale(TotalRNA),
    AlignRate_scaled = scale(AlignRate)
  )

# 6. Set rownames of metadata to match columns of expression data
stopifnot(all(colnames(gcount_filt) == treatmentinfo$sampleID))  # safety check
treatmentinfo <- as.data.frame(treatmentinfo)
rownames(treatmentinfo) <- treatmentinfo$sampleID

# Ensure timepoint is a factor
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I", "II", "III"))
rownames(treatmentinfo) <- treatmentinfo$sampleID  # required by dream

# Filter formula (only fixed effects)
filter_formula <- ~ timepoint + TIN_scaled + RNA_260_230_scaled + RNA_260_280_scaled + TotalRNA_scaled
keep <- filterByExpr(dge, design = model.matrix(filter_formula, data = treatmentinfo))
dge <- dge[keep, ]

# Define variancePartition formula with random effects
form <- ~ (1|timepoint) + TIN_scaled + RNA_260_230_scaled + RNA_260_280_scaled + TotalRNA_scaled

# Voom with precision weights
vobj <- voomWithDreamWeights(dge, formula = form, data = treatmentinfo)

# Fit the dream model
fit <- dream(vobj, formula = form, data = treatmentinfo)

# Confirm fit class
class(fit)  # should be "Dream"

# Extract variance fractions
vp <- fitExtractVarPartModel(fit)

# Plot variance fractions
plotVarPart(vp)



# 7. Create DGEList object and normalize
dge <- DGEList(counts = gcount_filt)
dge <- calcNormFactors(dge)

treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I", "II", "III"))

table(treatmentinfo$timepoint)
str(treatmentinfo$timepoint)


# 8. Filter genes for sufficient expression
form <- ~ (1|timepoint) + TIN_scaled + RNA_260_230_scaled + RNA_260_280_scaled + TotalRNA_scaled
keep <- filterByExpr(dge, design = model.matrix(~ timepoint + TIN_scaled + RNA_260_230_scaled + RNA_260_280_scaled + TotalRNA_scaled, data = treatmentinfo))
dge <- dge[keep, ]

# 9. Run voom with dream weights
vobj <- voomWithDreamWeights(dge, formula = form, data = treatmentinfo)

# 10. Fit the dream model
fit <- dream(vobj, formula = form, data = treatmentinfo)

# 11. Extract and plot variance fractions
vp <- fitExtractVarPartModel(fit)
plotVarPart(vp)
