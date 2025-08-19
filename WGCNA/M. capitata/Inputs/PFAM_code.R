
# Load libraries
library(tidyverse)
library(stringr)
library(readr)


setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs")

# Load eggNOG-mapper output (tab-separated, skip hash comments if any)
eggnog<- read_delim("Montipora_capitata_HIv3.genes.EggNog_results.txt", 
                                                           delim = "\t", escape_double = FALSE, 
                                                           trim_ws = TRUE)
somp_domains <- c(
  "PF00194",  # Carbonic_anhydrase
  "PF00431",  # CUB
  "PF00092",  # VWA (von Willebrand factor type A)
  "PF00093",  # VWFC (von Willebrand factor type C)
  "PF00090",  # Thrombospondin type 1
  "PF14649",  # Galaxin
  "PF00064",  # Cadherin
  "PF00059",  # EGF-like domain
  "PF00089",  # Lectin_C
  "PF05762",  # Mucin-like
  "PF00136",  # Kazal-type serine protease inhibitor
  "PF01456",  # C1q domain
  "PF07679",  # Low-complexity acidic domain
  "PF01347",  # Vitellogenin_N
  "PF10911",  # Vitellogenin_open_beta
  "PF00655",  # von Willebrand factor type D (in Vitellogenin)
  "PF07496"   # DUF1943 (domain of unknown function, sometimes in Vtg)
)

#this didn't work..
#somp_candidates <- eggnog %>%
#  filter(str_detect(PFAMs, str_c(somp_domains, collapse = "|")))
#colnames(eggnog)
#colnames(eggnog)
#[1] "#query"         "seed_ortholog"  "evalue"         "score"         
#[5] "eggNOG_OGs"     "max_annot_lvl"  "COG_category"   "Description"   
#[9] "Preferred_name" "GOs"            "EC"             "KEGG_ko"       
#[13] "KEGG_Pathway"   "KEGG_Module"    "KEGG_Reaction"  "KEGG_rclass"   
#[17] "BRITE"          "KEGG_TC"        "CAZy"           "BiGG_Reaction" 
#[21] "PFAMs"

eggnog$PFAMs[1:10]

#[1] "7tm_1"                  "L_HMGIC_fpl"            "SPRY"                  
#[4] "ZZ"                     "Ldl_recept_a"           "-"                     
#[7] "-"                      "Peptidase_C48"          "DUF1298,WES_acyltransf"
#[10] "RVT_1,rve,zf-H2C2" 

# Get a list of all unique domain names in PFAMs
eggnog %>%
  filter(!is.na(PFAMs)) %>%
  separate_rows(PFAMs, sep = ";") %>%
  count(PFAMs, sort = TRUE)

#PFAMs                         n
#<chr>                     <int>
#  1 -                          3451
#2 RVT_1,rve                   636
#3 RVT_1                       612
#4 Exo_endo_phos_2,RVT_1       419
#5 rve                         354
#6 7tm_1                       336
#7 DDE_Tnp_4                   266
#8 DUF1759,Peptidase_A17,rve   232
#9 RVT_1,rve,zf-H2C2           206
#10 DDE_Tnp_4,HTH_Tnp_4,THAP    180
# i 5,343 more rows


# Define known SOMP-related domain names
somp_domains <- c(
  "Carbonic_anhydrase", "CUB", "VWA", "VWFC", "Thrombospondin", "Galaxin",
  "Cadherin", "EGF", "Lectin_C", "Mucin", "Kazal", "C1q", "Low_complexity",
  "Vitellogenin", "Vitellogenin_N", "Vitellogenin_open_beta",
  "VWD", "DUF1943"
)

#Clean the data — remove blanks and NAs
eggnog_clean <- eggnog %>%
  filter(!is.na(PFAMs) & PFAMs != "-")

# Split comma-separated PFAMs into separate rows
eggnog_split <- eggnog_clean %>%
  separate_rows(PFAMs, sep = ",")

# Trim whitespace if needed
eggnog_split$PFAMs <- str_trim(eggnog_split$PFAMs)

# Filter to SOMP-related domain hits
somp_candidates <- eggnog_split %>%
  filter(PFAMs %in% somp_domains)

table(somp_candidates$PFAMs)

#C1q Cadherin      CUB      EGF Lectin_C      VWA      VWD 
#12       18       90      110      194      108       20 



# Define PFAMs of interest
somp_pfams <- c("VWA", "VWD", "CUB", "EGF", "Lectin_C", "C1q", "Kazal_1", "ZP")

# Filter rows that contain any of those PFAMs
library(dplyr)
somp_candidates <- eggnog %>%
  filter(str_detect(PFAMs, paste(somp_pfams, collapse = "|")))

# Optional: check summary
table(somp_candidates$PFAMs)

library(tidyverse)

somp_candidates %>%
  separate_rows(PFAMs, sep = ",") %>%
  count(PFAMs, sort = TRUE)


library(dplyr)
library(stringr)

dev_go_hits <- eggnog %>%
  filter(str_detect(GOs, regex("develop|morphogenesis|embryo|patterning", ignore_case = TRUE)))

dev_desc_hits <- eggnog %>%
  filter(str_detect(Description, regex("develop|morphogenesis|patterning|embryo", ignore_case = TRUE)))

dev_pfams <- c("Homeobox", "Homeodomain", "Forkhead", "T-box", "Wnt", "Frizzled", "Notch", "bHLH", "Hedgehog", "DSL", "ANK", "PAX", "Zinc_finger")

dev_pfam_hits <- eggnog %>%
  filter(str_detect(PFAMs, paste(dev_pfams, collapse = "|")))

dev_genes <- bind_rows(dev_go_hits, dev_desc_hits, dev_pfam_hits) %>%
  distinct(`#query`, .keep_all = TRUE)

