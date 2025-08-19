# ============================
# Montipora capitata — DEGs per life stage (vs Others + Pairwise) & GO enrichment via EggNOG
# ============================

# ---------- Packages ----------
# install.packages(c("tidyverse","data.table"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2","clusterProfiler","GO.db"))

library(tidyverse)
library(data.table)
library(DESeq2)
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)
setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------- User settings ----------
sample_file <- "5-Mcap-SampleInfo.csv"
counts_file <- "Mcap_transcript_count_matrix.csv"
eggnog_file <- "Montipora_capitata_HIv3.genes.EggNog_results.txt"

life_stage_col <- "timepoint"   # <-- your sample-info column name

# DEG thresholds
padj_cutoff <- 0.05
lfc_cutoff  <- 1

# If EggNOG headers don't autodetect, set these (otherwise leave NULL)
id_col_override <- NULL      # e.g., "query"
go_col_override <- NULL      # e.g., "GOs"

# ---------- Output dirs ----------
out_dir <- "results_mcap_DEG_enrich"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "pairwise"), showWarnings = FALSE)
dir.create(file.path(out_dir, "vs_others"), showWarnings = FALSE)
dir.create(file.path(out_dir, "enrichment"), showWarnings = FALSE)

# ---------- Load data ----------
message("Loading sample info and count matrix ...")
sample_info <- read.csv(sample_file, row.names = 1, check.names = FALSE)
counts      <- read.csv(counts_file, row.names = 1, check.names = FALSE)

# sanity checks
stopifnot(all(colnames(counts) %in% rownames(sample_info)))
sample_info <- sample_info[colnames(counts), , drop = FALSE]
stopifnot(life_stage_col %in% colnames(sample_info))

# factor with exactly 3 levels
sample_info[[life_stage_col]] <- factor(sample_info[[life_stage_col]])
stages <- levels(sample_info[[life_stage_col]])
if (length(stages) != 3) {
  stop("Expected exactly 3 levels in '", life_stage_col, "'. Found: ", paste(stages, collapse = ", "))
}

# ---------- Fit DESeq once with unified design ----------
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = sample_info,
  design    = as.formula(paste0("~ ", life_stage_col))
)

# optional prefilter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds, quiet = TRUE)  # size factors + dispersions + Wald test
stages <- levels(colData(dds)[[life_stage_col]])

# ---------- Pairwise contrasts (A vs B, etc.) ----------
message("Running pairwise DESeq2 contrasts ...")
pairwise_res <- list()
pairwise_deg <- list()

pairs <- t(combn(stages, 2))
for (i in seq_len(nrow(pairs))) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  
  # a vs b
  rab <- results(dds, contrast = c(life_stage_col, a, b)) %>%
    as.data.frame() %>% rownames_to_column("feature_id") %>% arrange(padj)
  dab <- dplyr::filter(rab, !is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  key_ab <- paste(a, "vs", b, sep = "_")
  pairwise_res[[key_ab]] <- rab
  pairwise_deg[[key_ab]] <- dab
  write.csv(rab, file.path(out_dir, "pairwise", paste0("DESeq2_results_", key_ab, ".csv")), row.names = FALSE)
  write.csv(dab, file.path(out_dir, "pairwise", paste0("DEGs_", key_ab, "_padj", padj_cutoff, "_lfc", lfc_cutoff, ".csv")), row.names = FALSE)
  message(key_ab, ": ", nrow(dab), " DEGs")
  # b vs a
  rba <- results(dds, contrast = c(life_stage_col, b, a)) %>%
    as.data.frame() %>% rownames_to_column("feature_id") %>% arrange(padj)
  dba <- dplyr::filter(rba, !is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  key_ba <- paste(b, "vs", a, sep = "_")
  pairwise_res[[key_ba]] <- rba
  pairwise_deg[[key_ba]] <- dba
  write.csv(rba, file.path(out_dir, "pairwise", paste0("DESeq2_results_", key_ba, ".csv")), row.names = FALSE)
  write.csv(dba, file.path(out_dir, "pairwise", paste0("DEGs_", key_ba, "_padj", padj_cutoff, "_lfc", lfc_cutoff, ".csv")), row.names = FALSE)
  message(key_ba, ": ", nrow(dba), " DEGs")
}

#I_vs_II: 249 DEGs
#II_vs_I: 249 DEGs
#I_vs_III: 12165 DEGs
#III_vs_I: 12165 DEGs
#II_vs_III: 12512 DEGs
#III_vs_II: 12512 DEGs

# ---------- "Stage vs Others" using list-contrasts ----------
# A vs geometric mean(B,C) on log2 scale: contrast = [coef(A) - 0.5*coef(B) - 0.5*coef(C)]
# Implementation detail: we rotate the reference to one of the "others" so coefficient names exist, then combine coefs.
message("Running 'stage vs others' list-contrasts ...")
vs_others_res <- list()
vs_others_deg <- list()

for (a in stages) {
  others <- setdiff(stages, a)
  ref <- others[1]
  
  # relevel reference and recompute Wald test coefs without refitting size factors/dispersion
  dds_ref <- dds
  colData(dds_ref)[[life_stage_col]] <- relevel(colData(dds_ref)[[life_stage_col]], ref = ref)
  dds_ref <- nbinomWaldTest(dds_ref)  # reuse size factors/dispersion
  
  coef_a_ref <- paste0(life_stage_col, "_", a, "_vs_", ref)
  coef_o_ref <- paste0(life_stage_col, "_", others[2], "_vs_", ref)
  
  if (!all(c(coef_a_ref, coef_o_ref) %in% resultsNames(dds_ref))) {
    stop("Could not find coefficients: ", coef_a_ref, " and/or ", coef_o_ref,
         ". Available: ", paste(resultsNames(dds_ref), collapse = ", "))
  }   
  
  r <- results(
    dds_ref,
    contrast   = list(c(coef_a_ref, coef_o_ref)),
    listValues = c(1, -0.5)
  ) %>% as.data.frame() %>% rownames_to_column("feature_id") %>% arrange(padj)
  
  d <- dplyr::filter(r, !is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  
  key <- paste0(a, "_vs_others")
  vs_others_res[[key]] <- r
  vs_others_deg[[key]] <- d
  
  write.csv(r, file.path(out_dir, "vs_others", paste0("DESeq2_results_", key, ".csv")), row.names = FALSE)
  write.csv(d, file.path(out_dir, "vs_others", paste0("DEGs_", key, "_padj", padj_cutoff, "_lfc", lfc_cutoff, ".csv")), row.names = FALSE)
  message(key, ": ", nrow(d), " DEGs")
}

# ---------- Parse EggNOG for GO mappings ----------
message("Parsing EggNOG for GO mappings ...")

# Read file
egg <- data.table::fread(eggnog_file, sep = "\t", header = TRUE, data.table = FALSE, quote = "")

# Force the column names we want
id_col <- "#query"
go_col <- "GOs"

if (!id_col %in% colnames(egg)) {
  stop("ID column '#query' not found. Found: ", paste(colnames(egg), collapse = ", "))
}
if (!go_col %in% colnames(egg)) {
  stop("GO column 'GOs' not found. Found: ", paste(colnames(egg), collapse = ", "))
}

# Select and rename
egg2 <- egg %>%
  dplyr::select(all_of(c(id_col, go_col))) %>%
  dplyr::rename(feature_id = all_of(id_col),
                GO_raw     = all_of(go_col)) %>%
  dplyr::mutate(GO_raw = ifelse(is.na(GO_raw), "", GO_raw))

# Split GOs by comma/semicolon/space
split_go <- function(x) {
  parts <- unlist(strsplit(x, "[,; ]+"))
  unique(parts[nzchar(parts)])
}

# Build TERM2GENE
TERM2GENE <- egg2 %>%
  dplyr::mutate(GO = purrr::map(GO_raw, split_go)) %>%
  dplyr::select(feature_id, GO) %>%
  tidyr::unnest(GO) %>%
  dplyr::filter(grepl("^GO:\\d{7}$", GO)) %>%
  dplyr::distinct(GO, feature_id)

# Build TERM2NAME for descriptions
all_go <- unique(TERM2GENE$GO)
go_terms <- AnnotationDbi::select(GO.db, keys = all_go, columns = c("TERM"), keytype = "GOID")
TERM2NAME <- as_tibble(go_terms) %>%
  dplyr::rename(GO = GOID, Description = TERM) %>%
  dplyr::distinct()

# Background universe = all tested features with GO mapping
universe <- intersect(rownames(dds), TERM2GENE$feature_id)

# ---------- Function for enrichment ----------
do_enrich <- function(gene_ids, label) {
  gene_set <- intersect(gene_ids, TERM2GENE$feature_id)
  if (length(gene_set) < 3) {
    message(label, ": Not enough mapped DEGs with GO terms (n=", length(gene_set), "). Skipping.")
    return(NULL)
  }
  eg <- enricher(
    gene         = gene_set,
    TERM2GENE    = TERM2GENE,
    TERM2NAME    = TERM2NAME,
    universe     = universe,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 5
  )
  if (!is.null(eg) && nrow(as.data.frame(eg)) > 0) {
    write.csv(as.data.frame(eg),
              file = file.path(out_dir, "enrichment",
                               paste0("GO_", gsub("[^A-Za-z0-9_]+","_", label), ".csv")),
              row.names = FALSE)
    # Optional: dotplot
    try({
      p <- dotplot(eg, showCategory = 20)
      ggsave(file.path(out_dir, "enrichment",
                       paste0("GO_dotplot_", gsub("[^A-Za-z0-9_]+","_", label), ".pdf")),
             p, width = 8, height = 6)
      ggsave(file.path(out_dir, "enrichment",
                       paste0("GO_dotplot_", gsub("[^A-Za-z0-9_]+","_", label), ".png")),
             p, width = 8, height = 6, dpi = 300)
    }, silent = TRUE)
  } else {
    message(label, ": No significant GO terms.")
  }
  eg
}

# ---------- Enrichment for pairwise ----------
message("GO enrichment for pairwise DEG sets ...")
enrich_pairwise <- list()
for (nm in names(pairwise_deg)) {
  genes <- pairwise_deg[[nm]]$feature_id
  enrich_pairwise[[nm]] <- do_enrich(genes, paste0("Pairwise_", nm, "_DEGs"))
}

# ---------- Enrichment for vs-others ----------
message("GO enrichment for 'vs others' DEG sets ...")
enrich_vs_others <- list()
for (nm in names(vs_others_deg)) {
  genes <- vs_others_deg[[nm]]$feature_id
  enrich_vs_others[[nm]] <- do_enrich(genes, paste0("VsOthers_", nm, "_DEGs"))
}

# ---------- Summary table ----------
summary_df <- tibble(
  Contrast   = c(names(pairwise_deg), names(vs_others_deg)),
  n_DEGs     = c(map_int(pairwise_deg, nrow), map_int(vs_others_deg, nrow)),
  n_GO_terms = c(
    map_int(names(pairwise_deg), ~ { eg <- enrich_pairwise[[.x]]; if (is.null(eg)) 0 else nrow(as.data.frame(eg)) }),
    map_int(names(vs_others_deg), ~ { eg <- enrich_vs_others[[.x]]; if (is.null(eg)) 0 else nrow(as.data.frame(eg)) })
  )
)
write.csv(summary_df, file.path(out_dir, "summary_counts.csv"), row.names = FALSE)

message("GO enrichment complete. Results are in: ", normalizePath(file.path(out_dir, "enrichment")))
