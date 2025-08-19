suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(AnnotationDbi)
})

options(stringsAsFactors = FALSE)

# =========================
# Helpers
# =========================

# Build a length vector from a GFF/GTF via GenomicFeatures.
# level = "gene" or "tx" (transcript). If NULL, we will try both and pick the best match.
build_length_vector <- function(gff_file, level = c("gene","tx")) {
  if (is.null(level)) stop("Internal: 'level' must be supplied to build_length_vector().")
  level <- match.arg(level)
  txdb  <- GenomicFeatures::makeTxDbFromGFF(gff_file, format = tools::file_ext(gff_file))
  
  if (level == "gene") {
    ex_by_gene   <- exonsBy(txdb, by = "gene")
    gene_lengths <- sum(width(reduce(ex_by_gene)))
    len_vec      <- as.numeric(gene_lengths)
    names(len_vec) <- names(gene_lengths)          # TxDb gene IDs (from GFF)
  } else {
    ex_by_tx   <- exonsBy(txdb, by = "tx")
    tx_lengths <- sum(width(reduce(ex_by_tx)))
    # Try to map TXID -> TXNAME (public transcript IDs). If absent, keep TXID.
    tx_map <- AnnotationDbi::select(txdb, keys = names(tx_lengths),
                                    keytype = "TXID", columns = "TXNAME")
    txname_by_txid <- setNames(tx_map$TXNAME, tx_map$TXID)
    nm <- unname(txname_by_txid[names(tx_lengths)])
    names(tx_lengths) <- ifelse(is.na(nm) | nm == "", names(tx_lengths), nm)
    len_vec <- as.numeric(tx_lengths)
    names(len_vec) <- names(tx_lengths)
  }
  
  keep <- !is.na(names(len_vec)) & names(len_vec) != ""
  len_vec[keep]
}

# Read a count matrix and return a numeric matrix with rownames from an ID column
read_counts <- function(counts_file,
                        id_col_candidates = c("gene_id","Geneid","transcript_id","ID")) {
  df <- readr::read_csv(counts_file, show_col_types = FALSE)
  id_col <- intersect(id_col_candidates, names(df))
  if (length(id_col) == 0) {
    stop("Could not find an ID column in: ", counts_file,
         " (looked for: ", paste(id_col_candidates, collapse = ", "), ")")
  }
  id_col <- id_col[1]
  df <- as.data.frame(df)
  rownames(df) <- as.character(df[[id_col]])
  df[[id_col]] <- NULL
  
  # Coerce remaining columns to integers
  df[] <- lapply(df, function(x) {
    if (is.logical(x)) as.integer(x) else suppressWarnings(as.integer(round(x)))
  })
  
  # Drop rows that became NA entirely
  df <- df[apply(!is.na(df), 1, any), , drop = FALSE]
  df[is.na(df)] <- 0L
  as.matrix(df)
}

# Compute TPM. Removes 0-length rows (with a warning).
counts_to_tpm <- function(counts_mat, lengths_vec) {
  stopifnot(is.numeric(lengths_vec))
  L <- lengths_vec[rownames(counts_mat)]
  if (any(is.na(L))) {
    stop("NA in matched lengths. Check ID overlap.")
  }
  if (any(L == 0)) {
    warning("Zero lengths detected. Dropping those rows before TPM.")
    keep <- L > 0
    counts_mat <- counts_mat[keep, , drop = FALSE]
    L <- L[keep]
  }
  rate <- sweep(counts_mat, 1, L, "/")
  lib <- colSums(rate)
  if (any(lib == 0)) {
    zero_cols <- colnames(counts_mat)[lib == 0]
    warning("Some samples have zero library size after length-scaling: ",
            paste(zero_cols, collapse = ", "))
    lib[lib == 0] <- NA_real_
  }
  t(t(rate) / lib) * 1e6
}

# Try both gene- and transcript-length vectors and choose the one with better overlap
choose_best_lengths <- function(count_rows, gff_file) {
  message("  • Building gene-lengths…")
  gene_len <- build_length_vector(gff_file, "gene")
  overlap_gene <- sum(count_rows %in% names(gene_len))
  
  message("  • Building transcript-lengths…")
  tx_len   <- build_length_vector(gff_file, "tx")
  overlap_tx <- sum(count_rows %in% names(tx_len))
  
  message(sprintf("    Overlap (gene): %d   Overlap (tx): %d", overlap_gene, overlap_tx))
  
  if (overlap_tx >= overlap_gene && overlap_tx > 0) {
    list(level = "tx", lengths = tx_len, overlap = overlap_tx)
  } else if (overlap_gene > 0) {
    list(level = "gene", lengths = gene_len, overlap = overlap_gene)
  } else {
    stop("No overlap with either gene or transcript IDs. ",
         "Check whether row IDs match the GFF (names/versions).")
  }
}

# Main driver that:
#  - builds/loads length vector (from gff_file OR lengths_file),
#  - picks level automatically if needed,
#  - computes TPM, and writes <basename>__TPM.csv
process_dataset <- function(counts_file,
                            gff_file = NULL,
                            lengths_file = NULL,     # CSV with columns: gene_id, average_length (or length)
                            level = NULL,            # "gene" or "tx" or NULL to auto-choose
                            out_suffix = "__TPM.csv",
                            id_col_candidates = c("gene_id","Geneid","transcript_id","ID")) {
  
  message("\n==========")
  message("Processing: ", counts_file)
  cm <- read_counts(counts_file, id_col_candidates = id_col_candidates)
  feature_ids <- rownames(cm)
  
  # Pick/prepare lengths
  if (!is.null(lengths_file)) {
    message("  • Using provided lengths file: ", lengths_file)
    len_df <- readr::read_csv(lengths_file, show_col_types = FALSE)
    names(len_df) <- tolower(names(len_df))
    # Accept a few common column names
    id_col  <- if ("gene_id" %in% names(len_df)) "gene_id" else
      if ("transcript_id" %in% names(len_df)) "transcript_id" else
        stop("Lengths file needs a 'gene_id' or 'transcript_id' column.")
    len_col <- if ("average_length" %in% names(len_df)) "average_length" else
      if ("length" %in% names(len_df)) "length" else
        stop("Lengths file needs a 'average_length' or 'length' column.")
    length_vec <- len_df[[len_col]]
    names(length_vec) <- as.character(len_df[[id_col]])
    length_vec <- as.numeric(length_vec)
    names(length_vec) <- as.character(names(length_vec))
    level_used <- if (id_col == "gene_id") "gene" else "tx"
  } else if (!is.null(gff_file)) {
    if (is.null(level)) {
      message("  • Detecting whether counts are gene- or transcript-level…")
      pick <- choose_best_lengths(feature_ids, gff_file)
      level_used  <- pick$level
      length_vec  <- pick$lengths
      message("    → Selected: ", level_used)
    } else {
      message("  • Building ", level, "-lengths from GFF: ", gff_file)
      length_vec <- build_length_vector(gff_file, level)
      level_used <- level
    }
  } else {
    stop("Provide either gff_file or lengths_file.")
  }
  
  # Match & report overlap
  shared <- intersect(feature_ids, names(length_vec))
  pct <- 100 * length(shared) / length(feature_ids)
  message(sprintf("  • Matched %d / %d rows (%.1f%%) to %s lengths.",
                  length(shared), length(feature_ids), pct, level_used))
  if (length(shared) == 0) stop("No features matched between counts and lengths.")
  
  cm_match <- cm[shared, , drop = FALSE]
  len_match <- length_vec[shared]
  
  # TPM
  tpm <- counts_to_tpm(cm_match, len_match)
  
  # Write
  out_file <- file.path(dirname(counts_file),
                        paste0(tools::file_path_sans_ext(basename(counts_file)), out_suffix))
  out_df <- data.frame(feature_id = rownames(tpm), tpm, check.names = FALSE)
  readr::write_csv(out_df, out_file)
  message("  • Wrote: ", out_file)
  
  invisible(list(out_file = out_file,
                 matched_rows = length(shared),
                 level = level_used))
}

# =========================
# CONFIG: your four datasets
# =========================
# 👉 EDIT THESE PATHS to match your folders

cfg <- list(
  # 1) Montipora capitata — transcript-level counts
  Mcap = list(
    counts_file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs/Mcap_transcript_count_matrix.csv",
    gff_file    = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs/Montipora_capitata_HIv3.genes.gff3",
    level       = "tx"      # known tx-level
  ),
  
  # 2) Pocillopora acuta — your file name says "gene", but earlier code matched transcripts.
  #    We'll let the script auto-detect; if you know for sure, set level="gene" or "tx".
  Pacu = list(
    counts_file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/P. acuta/Input/Pacu_gene_count_matrix_newGFF.csv",
    gff_file    = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/P. acuta/Input/Pocillopora_acuta_HIv2.genes_fixed.gff3",
    level       = NULL      # auto-choose best overlap
  ),
  
  # 3) Stylophora pistillata — gene-level counts with precomputed gene lengths CSV
  Spis = list(
    counts_file  = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Transcriptomics/Inputs/4-Spis-GeneCountMatrix.csv",
    lengths_file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Transcriptomics/Inputs/averaged_gene_lengths.csv"
    # level is inferred from the lengths_file (gene_id)
  ),
  
  # 4) Acropora tenuis — gene-level counts + GFF
  Aten = list(
    counts_file = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis/gene_count_matrix_Acropora.csv",
    gff_file    = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis/aten_0.11.maker_post_001.genes.gff",
    level       = "gene"    # known gene-level
  )
)

# =========================
# RUN
# =========================
results <- lapply(cfg, function(x) {
  do.call(process_dataset, x)
})

# Quick summary
cat("\nSummary:\n")
for (nm in names(results)) {
  r <- results[[nm]]
  cat(sprintf("  %-4s → %s (level: %s, matched rows: %d)\n",
              nm, r$out_file, r$level, r$matched_rows))
}
cat("Done.\n")



library(readr)
library(dplyr)

# Paths (edit if different)
counts_file  <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Transcriptomics/Inputs/4-Spis-GeneCountMatrix.csv"
lengths_file <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Transcriptomics/Inputs/averaged_gene_lengths.csv"

# Read counts (find an ID column and make it rownames)
read_counts_simple <- function(f, id_cols = c("gene_id","Geneid","transcript_id","ID","geneid","gene")) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  id_col <- intersect(id_cols, names(df))
  if (!length(id_col)) stop("No ID column found in counts file.")
  id_col <- id_col[1]
  df <- as.data.frame(df)
  rownames(df) <- as.character(df[[id_col]])
  df[[id_col]] <- NULL
  df[] <- lapply(df, function(x) suppressWarnings(as.integer(round(x))))
  df[is.na(df)] <- 0L
  as.data.frame(df)
}

cm <- read_counts_simple(counts_file)

# Peek at row IDs in counts
head(rownames(cm), 20)
length(rownames(cm))

# Read lengths; try to auto-detect id/length columns
len_df <- readr::read_csv(lengths_file, show_col_types = FALSE)
names(len_df) <- tolower(names(len_df))
candidate_id_cols  <- names(len_df)[grepl("gene|transcript|id", names(len_df))]
candidate_len_cols <- names(len_df)[grepl("length", names(len_df))]

candidate_id_cols
candidate_len_cols

# Print a few IDs from lengths file
head(len_df[[candidate_id_cols[1]]], 20)
nrow(len_df)

cm <- read_counts_simple(counts_file)

# Build a *named* length vector without dropping names
len_df <- readr::read_csv(lengths_file, show_col_types = FALSE)
names(len_df) <- tolower(names(len_df))
length_vec <- setNames(as.numeric(len_df$average_length), as.character(len_df$gene_id))

# Check overlap (should be 25,769)
shared <- intersect(rownames(cm), names(length_vec))
length(shared)  # expect 25769

counts_to_tpm <- function(counts_mat, lengths_vec) {
  L <- lengths_vec[rownames(counts_mat)]
  keep <- !is.na(L) & L > 0
  counts_mat <- as.matrix(counts_mat[keep, , drop = FALSE])
  L <- L[keep]
  rate <- sweep(counts_mat, 1, L, "/")
  t(t(rate) / colSums(rate)) * 1e6
}

tpm <- counts_to_tpm(cm[shared, , drop = FALSE], length_vec[shared])

out_file <- sub("\\.csv$", "__TPM.csv", counts_file)
readr::write_csv(data.frame(feature_id = rownames(tpm), tpm, check.names = FALSE), out_file)
message("Wrote: ", out_file)


suppressPackageStartupMessages({
  library(readr); library(dplyr)
  library(GenomicFeatures); library(GenomicRanges)
})

# ---- paths (edit if needed)
counts_file <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis/gene_count_matrix_Acropora.csv"
gff_file    <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis/aten_0.11.maker_post_001.genes.gff"

# ---- helpers
read_counts_simple <- function(f, id_cols = c("gene_id","Geneid","transcript_id","ID","geneid","gene")) {
  df <- readr::read_csv(f, show_col_types = FALSE)
  id_col <- intersect(id_cols, names(df))[1]
  df <- as.data.frame(df)
  rownames(df) <- as.character(df[[id_col]])
  df[[id_col]] <- NULL
  df[] <- lapply(df, function(x) suppressWarnings(as.integer(round(x))))
  df[is.na(df)] <- 0L
  as.matrix(df)
}

counts_to_tpm <- function(counts_mat, lengths_vec) {
  L <- lengths_vec[rownames(counts_mat)]
  keep <- !is.na(L) & L > 0
  counts_mat <- counts_mat[keep, , drop = FALSE]
  L <- L[keep]
  rate <- sweep(counts_mat, 1, L, "/")
  t(t(rate) / colSums(rate)) * 1e6
}

# ---- build gene lengths from GFF (Aten is gene-level)
txdb <- GenomicFeatures::makeTxDbFromGFF(gff_file, format = tools::file_ext(gff_file))
ex_by_gene <- exonsBy(txdb, by = "gene")
gene_lengths <- sum(width(reduce(ex_by_gene)))
length_vec <- setNames(as.numeric(gene_lengths), names(gene_lengths))

# ---- read counts and match IDs
cm <- read_counts_simple(counts_file)
shared <- intersect(rownames(cm), names(length_vec))
stopifnot(length(shared) > 0)

# ---- TPM + LogTPM
tpm <- counts_to_tpm(cm[shared, , drop = FALSE], length_vec[shared])
logtpm <- log2(tpm + 1)

# ---- write both
out_tpm    <- sub("\\.csv$", "__TPM.csv", counts_file)
out_logtpm <- sub("\\.csv$", "__LogTPM2.csv", counts_file)  # <- your requested name

write_csv(data.frame(feature_id = rownames(tpm), tpm, check.names = FALSE), out_tpm)
write_csv(data.frame(feature_id = rownames(logtpm), logtpm, check.names = FALSE), out_logtpm)

message("Wrote: ", out_tpm)
message("Wrote: ", out_logtpm)


suppressPackageStartupMessages({ library(readr) })

# Convert a single TPM file → LogTPM2
tpm_to_logtpm2 <- function(tpm_file, id_cols = c("feature_id","gene_id","transcript_id","ID")) {
  message("Reading: ", tpm_file)
  df <- readr::read_csv(tpm_file, show_col_types = FALSE)
  
  # Find the ID column (first match wins)
  id_col <- intersect(id_cols, names(df))
  if (!length(id_col)) stop("No ID column found in TPM file: ", tpm_file)
  id_col <- id_col[1]
  
  ids <- as.character(df[[id_col]])
  mat <- as.matrix(df[, setdiff(names(df), id_col), drop = FALSE])
  
  # Ensure numeric
  storage.mode(mat) <- "numeric"
  
  # Compute log2(TPM + 1)
  logtpm <- log2(mat + 1)
  rownames(logtpm) <- ids
  
  out_file <- sub("__TPM\\.csv$", "__LogTPM2.csv", tpm_file)
  readr::write_csv(
    data.frame(feature_id = rownames(logtpm), logtpm, check.names = FALSE),
    out_file
  )
  message("Wrote: ", out_file)
  invisible(out_file)
}

#
tpm_files <- c(
  "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs/Mcap_transcript_count_matrix__TPM.csv",
  "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/P. acuta/Input/Pacu_gene_count_matrix_newGFF__TPM.csv",
  "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Transcriptomics/Inputs/4-Spis-GeneCountMatrix__TPM.csv"
)

invisible(lapply(tpm_files, tpm_to_logtpm2))






############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(GenomicFeatures)
  library(GenomicRanges)
})

# ── Paths (EDIT if yours differ) ───────────────────────────
counts_path <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis/gene_count_matrix_Acropora.csv"
gff_path    <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis/aten_0.11.maker_post_001.genes.gff"
out_dir     <- "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis"

# samples to remove (drop ASAP)
remove_samples <- c("DRR318292", "DRR318296", "DRR318290")

# ── 1) Load counts, drop genomic cols, set rownames, remove bad samples ─
df <- read_csv(counts_path, show_col_types = FALSE)

# normalize ID column name
if ("Geneid" %in% names(df)) names(df)[names(df) == "Geneid"] <- "gene_id"
stopifnot("gene_id" %in% names(df))

# drop genomic columns if present
drop_cols <- intersect(c("Chr","Start","End","Strand","Length"), names(df))
df <- dplyr::select(df, -all_of(drop_cols))


# make matrix: rows=genes, cols=samples
ids <- as.character(df$gene_id)
mat <- as.matrix(df[, setdiff(names(df), "gene_id"), drop = FALSE])
storage.mode(mat) <- "integer"
rownames(mat) <- ids
mat[is.na(mat)] <- 0L

# remove unwanted samples *before* normalization
keep_cols <- setdiff(colnames(mat), remove_samples)
mat <- mat[, keep_cols, drop = FALSE]

cat("Counts matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")
cat("First samples:", paste(head(colnames(mat), 6), collapse=", "), "\n")
cat("First genes:", paste(head(rownames(mat), 6), collapse=", "), "\n")

# ── 2) Build gene lengths from GFF and check overlap ─────────
txdb <- makeTxDbFromGFF(gff_path, format = tools::file_ext(gff_path))
ex_by_gene   <- exonsBy(txdb, by = "gene")
gene_lengths <- sum(width(reduce(ex_by_gene)))
len_vec      <- setNames(as.numeric(gene_lengths), names(gene_lengths))

shared <- intersect(rownames(mat), names(len_vec))
cat("Length overlap:", length(shared), "/", nrow(mat),
    sprintf("(%.1f%%)\n", 100*length(shared)/nrow(mat)))

if (!length(shared)) stop("No overlap between counts and gene lengths.")

# ── 3) TPM (only shared rows) ───────────────────────────────
counts_to_tpm <- function(counts_mat, lengths_vec) {
  L <- lengths_vec[rownames(counts_mat)]
  keep <- !is.na(L) & L > 0
  counts_mat <- counts_mat[keep, , drop = FALSE]
  L <- L[keep]
  rate <- sweep(counts_mat, 1, L, "/")
  t(t(rate) / colSums(rate)) * 1e6
}

tpm <- counts_to_tpm(mat[shared, , drop = FALSE], len_vec[shared])

cat("TPM dims:", paste(dim(tpm), collapse=" x "), "\n")
cat("TPM quick summary (first 100x5):\n")
print(summary(as.numeric(tpm[seq_len(min(100, nrow(tpm))), seq_len(min(5, ncol(tpm)))])))
stopifnot(!anyNA(tpm), !any(colSums(tpm) == 0))

# ── 4) Log2(TPM+1) ─────────────────────────────────────────
logtpm <- log2(tpm + 1)
stopifnot(all(colnames(logtpm) == colnames(tpm)))

# ── 5) Write clean outputs with a proper feature_id column ─
tpm_out <- file.path(out_dir, "gene_count_matrix_Acropora__TPM_CLEAN_filtered.csv")
log_out <- file.path(out_dir, "gene_count_matrix_Acropora__LogTPM2_CLEAN_filtered.csv")

write_csv(data.frame(feature_id = rownames(tpm),   tpm,    check.names = FALSE), tpm_out)
write_csv(data.frame(feature_id = rownames(logtpm), logtpm, check.names = FALSE), log_out)

# round-trip peek
peek <- read_csv(log_out, n_max = 3, show_col_types = FALSE)
cat("Wrote:\n  ", tpm_out, "\n  ", log_out, "\n")
cat("LogTPM2 columns peek:\n", paste(names(peek), collapse=" | "), "\n")

#Making a multi-species matrix
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ─────────────────────────────────────────────────────────────
# 0) Paths
# ─────────────────────────────────────────────────────────────

# Log2(TPM+1) matrices you just wrote
logtpm_files <- list(
  Aten = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/A_tenuis/gene_count_matrix_Acropora__LogTPM2.csv",
  Mcap = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/M. capitata/Inputs/Mcap_transcript_count_matrix__LogTPM2.csv",
  Pacu = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/P. acuta/Input/Pacu_gene_count_matrix_newGFF__LogTPM2.csv",
  Spis = "C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/S. pistillata/Transcriptomics/Inputs/4-Spis-GeneCountMatrix__LogTPM2.csv"
)

# Column names in orthogroups file that hold gene IDs per species
ortho_cols <- c(
  Aten = "Acropora_tenuis_0.11.maker_post_001.proteins",
  Mcap = "Montipora_capitata_HIv3.genes.pep",
  Pacu = "Pocillopora_acuta_HIv2.genes.pep",
  Spis = "Spis.genome.annotation.pep.longest"
)

# If your stylophora column with no suffix exists, prefer it:
# (we’ll auto-detect; leave as-is)
# ortho_cols["Spis"] <- "Spis.genome.annotation.pep.longest_no_suffixes"

# ─────────────────────────────────────────────────────────────
# 1) Helpers
# ─────────────────────────────────────────────────────────────
read_logtpm <- function(f, id_cols = c("feature_id","gene_id","transcript_id","ID")) {
  stopifnot(file.exists(f))
  df <- readr::read_csv(f, show_col_types = FALSE)
  idcol <- intersect(id_cols, names(df))
  if (!length(idcol)) stop("No ID column found in: ", f)
  idcol <- idcol[1]
  ids <- as.character(df[[idcol]])
  mat <- as.matrix(df[, setdiff(names(df), idcol), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  mat
}

# Try to resolve known naming quirks per species
resolve_id <- function(id, species, rowname_set) {
  cand <- character(0)
  if (species == "Spis") {
    # orthogroups often store "Spis1234" while matrix might have "SpisGene1234" or transcript suffixes
    base <- id
    # If the orthogroup has transcript suffix (.t1a etc.), also try without
    base <- sub("\\.t[0-9]+[a-z]?$", "", base)
    cand <- c(
      base,
      sub("^Spis", "SpisGene", base),         # add SpisGene prefix
      paste0(base, ".t1"),                    # common default isoform
      paste0(sub("^Spis", "SpisGene", base), ".t1")
    )
  } else if (species == "Aten") {
    # your counts often need ".m1" appended
    cand <- c(id, paste0(id, ".m1"))
  } else {
    # Mcap & Pacu: try as-is and minor cleanups
    base <- id
    cand <- c(
      base,
      sub("\\s+", "", base),                  # strip spaces
      sub("\\|.*$", "", base)                 # drop pipe annotations if any
    )
  }
  # return the first candidate that exists, else NA
  hit <- cand[cand %in% rowname_set]
  if (length(hit)) hit[1] else NA_character_
}

# Build per-species OG-indexed expression matrix
build_og_matrix <- function(orthos, species, logtpm_mat, ortho_col_name) {
  if (!ortho_col_name %in% names(orthos)) {
    stop("Column '", ortho_col_name, "' not in orthogroups table.")
  }
  og_ids <- orthos$Orthogroup
  g_ids  <- as.character(orthos[[ortho_col_name]])
  
  # If the orthogroup column has comma-separated multiple IDs, keep single-copy only
  # (split and drop rows with != 1 gene)
  split_ids <- strsplit(ifelse(is.na(g_ids), "", g_ids), ",")
  split_ids <- lapply(split_ids, function(x) unique(trimws(x[x != ""])))
  keep <- vapply(split_ids, length, integer(1)) == 1
  og_ids <- og_ids[keep]
  g_ids  <- vapply(split_ids[keep], `[`, character(1), 1)
  
  # Resolve to matrix rownames
  rnames <- rownames(logtpm_mat)
  resolved <- vapply(g_ids, resolve_id, character(1), species = species, rowname_set = rnames)
  
  # Report mapping stats
  message(sprintf("[%s] %d OGs, %d resolved (%.1f%%)",
                  species, length(og_ids), sum(!is.na(resolved)),
                  100*mean(!is.na(resolved))))
  
  # Build matrix with OG rows; missing genes become NA
  res_idx <- match(resolved, rnames)
  mat <- matrix(NA_real_, nrow = length(og_ids), ncol = ncol(logtpm_mat),
                dimnames = list(og_ids, colnames(logtpm_mat)))
  ok <- !is.na(res_idx)
  if (any(ok)) mat[ok, ] <- logtpm_mat[res_idx[ok], , drop = FALSE]
  attr(mat, "GeneID") <- setNames(resolved, og_ids)
  mat
}

# ─────────────────────────────────────────────────────────────
# 2) Load data
# ─────────────────────────────────────────────────────────────
orthos <- readr::read_csv(orthos_file, show_col_types = FALSE)

# Prefer a no-suffix Spis column if present
if ("Spis.genome.annotation.pep.longest_no_suffixes" %in% names(orthos)) {
  ortho_cols["Spis"] <- "Spis.genome.annotation.pep.longest_no_suffixes"
}
#.m1 in Aten
expr <- lapply(logtpm_files, read_logtpm)

resolve_id <- function(id, species, rowname_set){
  cand <- character(0)
  if (species == "Spis"){
    base <- sub("\\.t[0-9]+[a-z]?$", "", id)
    cand <- c(base, sub("^Spis","SpisGene", base), paste0(base, ".t1"), paste0(sub("^Spis","SpisGene", base), ".t1"))
  } else if (species == "Aten"){
    # normalize
    base <- tolower(gsub("\\s+", "", id))
    # main fix: strip a trailing ".m1" if present
    c1 <- sub("\\.m1$", "", base)
    cand <- c(c1, base)
  } else {
    base <- tolower(gsub("\\s+", "", id))
    cand <- c(base, sub("\\|.*$", "", base))
  }
  # case-insensitive match against matrix rownames
  rset <- tolower(rowname_set)
  hit  <- cand[ match(cand, rset, nomatch = 0) > 0 ]
  if (length(hit)) {
    return(rowname_set[ match(hit[1], rset) ]) # return original-case rowname
  } else NA_character_
}

# ─────────────────────────────────────────────────────────────
# 3) Build per-species OG matrices (LogTPM2)
# ─────────────────────────────────────────────────────────────
og_mats <- list(
  Aten = build_og_matrix(orthos, "Aten", expr$Aten, ortho_cols["Aten"]),
  Mcap = build_og_matrix(orthos, "Mcap", expr$Mcap, ortho_cols["Mcap"]),
  Pacu = build_og_matrix(orthos, "Pacu", expr$Pacu, ortho_cols["Pacu"]),
  Spis = build_og_matrix(orthos, "Spis", expr$Spis, ortho_cols["Spis"])
)

# Save the per-species wide matrices (rows = Orthogroup, cols = samples)
out_dir <- dirname(orthos_file)
for (sp in names(og_mats)) {
  out_path <- file.path(out_dir, paste0("OG_LogTPM2_", sp, ".csv"))
  gid <- attr(og_mats[[sp]], "GeneID")
  out_df <- cbind(Orthogroup = rownames(og_mats[[sp]]),
                  GeneID = unname(gid[rownames(og_mats[[sp]])]),
                  as.data.frame(og_mats[[sp]], check.names = FALSE))
  readr::write_csv(out_df, out_path)
  message("Wrote: ", out_path)
}

# ─────────────────────────────────────────────────────────────
# 4) Intersect to OGs present in ALL 4 species (optional, handy for downstream)
# ─────────────────────────────────────────────────────────────
common_ogs <- Reduce(intersect, lapply(og_mats, rownames))
message("Orthogroups present in ALL 4 species: ", length(common_ogs))

og_mats_all4 <- lapply(og_mats, function(m) m[common_ogs, , drop = FALSE])

# Save intersected per-species matrices
for (sp in names(og_mats_all4)) {
  out_path <- file.path(out_dir, paste0("OG_LogTPM2_", sp, "_ALL4.csv"))
  gid <- attr(og_mats[[sp]], "GeneID")
  out_df <- cbind(Orthogroup = rownames(og_mats_all4[[sp]]),
                  GeneID = unname(gid[rownames(og_mats_all4[[sp]])]),
                  as.data.frame(og_mats_all4[[sp]], check.names = FALSE))
  readr::write_csv(out_df, out_path)
  message("Wrote: ", out_path)
}

# ─────────────────────────────────────────────────────────────
# 5) Long / tidy table (all species together)
# ─────────────────────────────────────────────────────────────
long_list <- list()
for (sp in names(og_mats)) {
  m <- og_mats[[sp]]
  gid <- attr(m, "GeneID")
  df <- as.data.frame(m)
  df$Orthogroup <- rownames(m)
  df$GeneID <- unname(gid[rownames(m)])
  df <- df |> relocate(Orthogroup, GeneID)
  long <- df |>
    pivot_longer(cols = -c(Orthogroup, GeneID),
                 names_to = "Sample", values_to = "LogTPM2") |>
    mutate(Species = sp, .before = 1)
  long_list[[sp]] <- long
}
long_all <- bind_rows(long_list)

# Save the long format
out_long <- file.path(out_dir, "OG_LogTPM2_long_all_species.csv")
readr::write_csv(long_all, out_long)
message("Wrote: ", out_long)


# 1) Pick the correct Spis column from your orthogroups table
spis_col <- if ("Spis.genome.annotation.pep.longest_no_suffixes" %in% names(orthos))
  "Spis.genome.annotation.pep.longest_no_suffixes" else
    "Spis.genome.annotation.pep.longest"

# 2) Take single-copy only (exactly one ID per OG)
spis_ids_raw <- orthos[[spis_col]]
split_ids <- strsplit(ifelse(is.na(spis_ids_raw), "", spis_ids_raw), ",")
split_ids <- lapply(split_ids, function(x) unique(trimws(x[x != ""])))
keep_one  <- vapply(split_ids, length, integer(1)) == 1
og_ids_sc <- orthos$Orthogroup[keep_one]
spis_ids  <- vapply(split_ids[keep_one], `[`, character(1), 1)

# 3) Peek a few and compare to the expression matrix rownames
cat("Spis OG examples:\n"); print(head(spis_ids, 10))
cat("Spis matrix rownames:\n"); print(head(rownames(expr$Spis), 10))

# Replace your resolve_id() with this version (Aten/Mcap/Pacu unchanged, Spis fixed)
resolve_id <- function(id, species, rowname_set){
  norm <- function(x){
    x <- tolower(trimws(x))
    gsub("\\s+", "", x)
  }
  rset_lc <- norm(rowname_set)
  
  if (species == "Spis") {
    base <- norm(id)                                   # e.g., "spis10345"
    base_no_t <- sub("\\.t[0-9]+[a-z]?$", "", base)    # drop any ".t1", ".t2a" if present
    
    # Build candidates in order of likelihood:
    cand <- unique(c(
      base,                                            # "spis10345"
      sub("^spis", "spisgene", base),                  # "spisgene10345"
      base_no_t,                                       # "spis10345" (if OG had .t1)
      sub("^spis", "spisgene", base_no_t),             # "spisgene10345"
      paste0(sub("^spis", "spisgene", base_no_t), ".t1"),   # "spisgene10345.t1"
      paste0(base_no_t, ".t1")                              # "spis10345.t1" (just in case)
    ))
    
    m <- match(cand, rset_lc, nomatch = 0)
    if (any(m > 0)) return(rowname_set[m[which(m > 0)[1]]])
    return(NA_character_)
  }
  
  if (species == "Aten") {
    base <- norm(id)
    base_fix <- sub("\\.m1$", "", base)   # strip the *extra* trailing ".m1" from OG side
    cand <- unique(c(base_fix, base))
    m <- match(cand, rset_lc, nomatch = 0)
    if (any(m > 0)) return(rowname_set[m[which(m > 0)[1]]])
    return(NA_character_)
  }
  
  # Mcap & Pacu
  base <- norm(id)
  cand <- c(base, sub("\\|.*$", "", base))
  m <- match(cand, rset_lc, nomatch = 0)
  if (any(m > 0)) return(rowname_set[m[which(m > 0)[1]]])
  NA_character_
}

# if you created a no-suffix column earlier, prefer it; else use the default
spis_col <- if ("Spis.genome.annotation.pep.longest_no_suffixes" %in% names(orthos))
  "Spis.genome.annotation.pep.longest_no_suffixes" else
    "Spis.genome.annotation.pep.longest"

og_mats <- list(
  Aten = build_og_matrix(orthos, "Aten", expr$Aten, ortho_cols[["Aten"]]),
  Mcap = build_og_matrix(orthos, "Mcap", expr$Mcap, ortho_cols[["Mcap"]]),
  Pacu = build_og_matrix(orthos, "Pacu", expr$Pacu, ortho_cols[["Pacu"]]),
  Spis = build_og_matrix(orthos, "Spis", expr$Spis, spis_col)
)

# show first 10 OGs with their GeneID for each species
for (sp in names(og_mats)) {
  cat("\n==", sp, "==\n")
  gid <- attr(og_mats[[sp]], "GeneID")
  print(head(data.frame(Orthogroup = names(gid), GeneID = unname(gid)), 10))
}

out_dir <- dirname(orthos_file)
for (sp in names(og_mats)) {
  out_path <- file.path(out_dir, paste0("OG_LogTPM2_", sp, ".csv"))
  gid <- attr(og_mats[[sp]], "GeneID")
  out_df <- cbind(
    Orthogroup = rownames(og_mats[[sp]]),
    GeneID     = unname(gid[rownames(og_mats[[sp]])]),
    as.data.frame(og_mats[[sp]], check.names = FALSE)
  )
  readr::write_csv(out_df, out_path)
  message("Wrote: ", out_path)
}



suppressPackageStartupMessages({ library(readr); library(dplyr) })

# --- 0) helper to prefix sample names to avoid collisions ---
prefix_cols <- function(mat, prefix){
  colnames(mat) <- paste0(prefix, "__", colnames(mat))
  mat
}

# --- 1) prefix columns per species (prevents same sample names clashing) ---
og_mats$Aten <- prefix_cols(og_mats$Aten, "Aten")
og_mats$Mcap <- prefix_cols(og_mats$Mcap, "Mcap")
og_mats$Pacu <- prefix_cols(og_mats$Pacu, "Pacu")
og_mats$Spis <- prefix_cols(og_mats$Spis, "Spis")

# --- 2) ALL4 (intersection) ---
common_ogs <- Reduce(intersect, lapply(og_mats, rownames))
message("Common OGs across all four species: ", length(common_ogs))

# subset and order identically
A <- og_mats$Aten[common_ogs, , drop = FALSE]
M <- og_mats$Mcap[common_ogs, , drop = FALSE]
P <- og_mats$Pacu[common_ogs, , drop = FALSE]
S <- og_mats$Spis[common_ogs, , drop = FALSE]

# sanity
stopifnot(identical(rownames(A), rownames(M)),
          identical(rownames(A), rownames(P)),
          identical(rownames(A), rownames(S)))

# bind columns
merged_all4 <- cbind(A, M, P, S)
merged_all4_df <- data.frame(Orthogroup = rownames(merged_all4),
                             merged_all4,
                             check.names = FALSE)

# write
out_dir <- dirname(orthos_file)
out_all4 <- file.path(out_dir, "OG_LogTPM2_ALL4_merged_wide.csv")
readr::write_csv(merged_all4_df, out_all4)
message("Wrote: ", out_all4)

# --- 3) UNION (all OGs seen in any species; NAs where missing) ---
# turn each into a data.frame with Orthogroup column for merges
dfA <- data.frame(Orthogroup = rownames(og_mats$Aten), og_mats$Aten, check.names = FALSE)
dfM <- data.frame(Orthogroup = rownames(og_mats$Mcap), og_mats$Mcap, check.names = FALSE)
dfP <- data.frame(Orthogroup = rownames(og_mats$Pacu), og_mats$Pacu, check.names = FALSE)
dfS <- data.frame(Orthogroup = rownames(og_mats$Spis), og_mats$Spis, check.names = FALSE)

merged_union <- Reduce(function(x, y) merge(x, y, by = "Orthogroup", all = TRUE),
                       list(dfA, dfM, dfP, dfS))

out_union <- file.path(out_dir, "OG_LogTPM2_UNION_merged_wide.csv")
readr::write_csv(merged_union, out_union)
message("Wrote: ", out_union)

# --- 4) PCA example on ALL4 (samples as rows, features as columns) ---
# Use the numeric block only (drop Orthogroup)
mat_all4 <- as.matrix(merged_all4_df[ , -1, drop = FALSE])
# transpose so rows = samples
pca_all4 <- prcomp(t(mat_all4), scale. = TRUE)

# quick check of explained variance in first 2 PCs
pv <- round(100 * summary(pca_all4)$importance[2, 1:2], 1)
message(sprintf("PC1=%.1f%%, PC2=%.1f%%", pv[1], pv[2]))





# ---------------------------------------
# PCA on TPM (no batch correction)
# ---------------------------------------

library(DESeq2)
library(ggplot2)
library(dplyr)
library(genefilter)

setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")

treatmentinfo <- read_csv("treatmentinfo_merged.csv")

# build DESeq2 object
treatmentinfo$timepoint <- factor(treatmentinfo$timepoint, levels = c("I","II","III"))
treatmentinfo$Species   <- factor(treatmentinfo$Species, 
                                  levels = c("Acropora_tenuis", "Montipora_capitata", 
                                             "Pocillopora_acuta","Stylophora_pistillata"))


# Set filter values for PoverA
filt <- filterfun(pOverA(0.083, 10))

# Apply filter to count data only
gfilt <- genefilter(gcount, filt)


gdds_nobatch <- DESeqDataSetFromMatrix(countData = gcount_filt,
                                       colData   = treatmentinfo,
                                       design    = ~ Species * timepoint)

# VST transform (normalizes library size, stabilizes variance)
gvst_nobatch <- vst(gdds_nobatch, blind = FALSE)

# PCA
pca_nobatch <- prcomp(t(assay(gvst_nobatch)))

# Make dataframe for plotting
pca_df_nobatch <- as.data.frame(pca_nobatch$x)
pca_df_nobatch$timepoint <- treatmentinfo$timepoint
pca_df_nobatch$Species   <- treatmentinfo$Species

# Color-blind-friendly palette
timepoint_palette <- c(
  "I"   = "#e9a3c9",         # pink
  "II"  = "lightgoldenrod2", # yellow
  "III" = "#a1d76a"          # green
)

plot_pca_cb <- function(pca_df, pca_obj, title = "Global Gene Expression by Life Stage and Species") {
  percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))[1:2]
  
  ggplot(pca_df, aes(x = PC1, y = PC2)) +
    # ellipses by timepoint
    stat_ellipse(aes(group = timepoint, fill = timepoint),
                 geom = "polygon", alpha = 0.3, color = NA, level = 0.95) +
    # points
    geom_point(aes(fill = timepoint, shape = Species),
               size = 5, color = "black", stroke = 0.8) +
    scale_fill_manual(values = timepoint_palette, name = "Timepoint") +
    scale_shape_manual(
      values = c(
        "Acropora_tenuis"      = 21,
        "Montipora_capitata"   = 22,
        "Pocillopora_acuta"    = 24,
        "Stylophora_pistillata"= 23
      ),
      labels = c(
        "Acropora_tenuis"      = expression(italic("A. tenuis")),
        "Montipora_capitata"   = expression(italic("M. capitata")),
        "Pocillopora_acuta"    = expression(italic("P. acuta")),
        "Stylophora_pistillata"= expression(italic("S. pistillata"))
      ),
      name = "Species"
    ) +
    guides(
      fill  = guide_legend(override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(override.aes = list(fill = "grey80", color = "black"))
    ) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(title) +
    coord_fixed() +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid   = element_blank(),
      plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11)
    )
}

# Plot
p_nobatch <- plot_pca_cb(pca_df_nobatch, pca_nobatch)
print(p_nobatch)

# Save
ggsave("PCA_Global_Expression_noBatch.png", plot = p_nobatch,
       width = 7, height = 5.5, dpi = 600)


# ─────────────────────────────────────────────────────────────
# PCA on merged TPM matrix (no VST, no batch correction)
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
})

setwd("C:/Users/amurg/OneDrive/Documentos/GitHub/Multistage_omics/WGCNA/Multi_species_WGCNA")

# 1) Files
tpm_file <- "OG_LogTPM2_UNION_merged_wide.csv"        # <-- SET THIS to your merged TPM *ALL4* wide file
meta_file <- "treatmentinfo_merged.csv"

# --- Clean sample names in expression matrix ---
samples <- colnames(og_mat)
samples_clean <- sub("^[A-Za-z]+__", "", samples)  # strip prefixes like Aten__, Mcap__, Pacu__, Spis__


# --- Match metadata ---
if (!"sampleID" %in% colnames(tinfo)) {
  stop("Metadata file missing 'sampleID' column. Found: ", paste(colnames(tinfo), collapse = ", "))
}

idx <- match(samples_clean, tinfo$sampleID)

# Check for unmatched samples
if (any(is.na(idx))) {
  stop("Unmatched samples: ", paste(samples_clean[is.na(idx)], collapse=", "))
}

# Reorder metadata to match expression matrix
tinfo <- tinfo[idx, ]

# Confirm alignment
stopifnot(all(samples_clean == tinfo$sampleID))



## ----------------------------
## Continue with in-memory objs
## ----------------------------

# 0) Pick your expression matrix object that's currently in memory
#    (use the first one that exists)
if (exists("og_mat_f")) {
  expr <- og_mat_f
} else if (exists("og_mat")) {
  expr <- og_mat
} else if (exists("og_wide")) {
  expr <- og_wide
} else {
  stop("No expression matrix found in memory (expected og_mat_f / og_mat / og_wide).")
}

# 1) Clean sample names (strip 'Aten__', 'Mcap__', 'Pacu__', 'Spis__')
samples_full  <- colnames(expr)
samples_clean <- sub("^.*__", "", samples_full)
colnames(expr) <- samples_clean

# 2) Get metadata already in memory
if (exists("tinfo")) {
  meta <- tinfo
} else if (exists("treatmentinfo")) {
  meta <- treatmentinfo
} else {
  stop("No metadata object found in memory (expected tinfo or treatmentinfo).")
}

# 3) Align metadata to expression columns by `sampleID`
if (!"sampleID" %in% colnames(meta)) {
  stop("Your metadata object must contain a 'sampleID' column. Found: ",
       paste(colnames(meta), collapse=", "))
}

idx <- match(samples_clean, meta$sampleID)
if (any(is.na(idx))) {
  stop("Unmatched samples in metadata: ", paste(samples_clean[is.na(idx)], collapse=", "))
}
ti <- meta[idx, , drop = FALSE]

# 4) Ensure factors for plotting
if (!all(c("timepoint","Species") %in% names(ti))) {
  stop("Metadata must include 'timepoint' and 'Species' columns.")
}
ti$timepoint <- factor(ti$timepoint, levels = c("I","II","III"))
ti$Species   <- factor(ti$Species,
                       levels = c("Acropora_tenuis","Montipora_capitata",
                                  "Pocillopora_acuta","Stylophora_pistillata"))

# 5) Confirm matrix is LogTPM2; if it looks like raw TPM, log it
x <- as.matrix(expr)
if (max(x, na.rm = TRUE) > 50) {  # heuristic: likely raw TPM
  message("Detected large values; applying log2(TPM+1).")
  x <- log2(x + 1)
}

# final guard
if (anyNA(x)) {
  bad <- rownames(x)[rowSums(is.na(x)) > 0][1]
  stop("Found NA values (e.g., row ", bad, "). Fix before PCA.")
}

# 6) PCA (samples as rows)
pca_obj <- prcomp(t(x), scale. = TRUE)

pca_df <- as.data.frame(pca_obj$x)
pca_df$timepoint <- ti$timepoint
pca_df$Species   <- ti$Species

# 7) Plot PC1–PC2
suppressPackageStartupMessages(library(ggplot2))
timepoint_palette <- c("I"="#e9a3c9", "II"="lightgoldenrod2", "III"="#a1d76a")

plot_pca_cb <- function(pca_df, pca_obj, title = "Global Gene Expression (LogTPM2)") {
  percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))[1:2]
  ggplot(pca_df, aes(PC1, PC2)) +
    stat_ellipse(aes(group = timepoint, fill = timepoint),
                 geom = "polygon", alpha = 0.30, color = NA, level = 0.95) +
    geom_point(aes(fill = timepoint, shape = Species),
               size = 5, color = "black", stroke = 0.8) +
    scale_fill_manual(values = timepoint_palette, name = "Timepoint") +
    scale_shape_manual(
      values = c("Acropora_tenuis"=21,"Montipora_capitata"=22,"Pocillopora_acuta"=24,"Stylophora_pistillata"=23),
      labels = c("Acropora_tenuis"=expression(italic("A. tenuis")),
                 "Montipora_capitata"=expression(italic("M. capitata")),
                 "Pocillopora_acuta"=expression(italic("P. acuta")),
                 "Stylophora_pistillata"=expression(italic("S. pistillata"))),
      name = "Species"
    ) +
    guides(fill  = guide_legend(override.aes = list(shape = 21, color = "black")),
           shape = guide_legend(override.aes = list(fill = "grey80", color = "black"))) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(title) +
    coord_fixed() +
    theme_minimal(base_size = 14) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.title = element_text(size = 12),
          legend.text  = element_text(size = 11))
}

p <- plot_pca_cb(pca_df, pca_obj)
print(p)

ggsave("PCA_LogTPM2_PC1_PC2.png", plot = p, width = 7, height = 5.5, dpi = 600)

# QC
pv <- round(100 * summary(pca_obj)$importance[2, 1:4], 1)
cat("Explained variance (%): ", paste(sprintf("PC%d=%0.1f", 1:4, pv), collapse=", "), "\n")
cat("Matrix dims used (OGs x samples): ", nrow(x), " x ", ncol(x), "\n")

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(ggplot2)
  library(gridExtra)  # for side-by-side plots
})

# log_mat = log2(TPM+1) matrix (genes x samples)
# tinfo   = metadata aligned to columns

## ------------------------
## 1) PCA on raw log_mat
## ------------------------
pca_raw <- prcomp(t(log_mat), scale. = TRUE)

plot_pca <- function(pca_obj, meta, title) {
  percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))[1:2]
  df <- as.data.frame(pca_obj$x)
  df$timepoint <- meta$timepoint
  df$Species   <- meta$Species
  
  ggplot(df, aes(PC1, PC2)) +
    stat_ellipse(aes(group = timepoint, fill = timepoint),
                 geom = "polygon", alpha = 0.3, color = NA) +
    geom_point(aes(fill = timepoint, shape = Species),
               size = 4, color = "black", stroke = 0.8) +
    scale_fill_manual(
      values = c("I"="#e9a3c9","II"="lightgoldenrod2","III"="#a1d76a"),
      name = "Timepoint"
    ) +
    scale_shape_manual(
      values = c("Acropora_tenuis"=21,"Montipora_capitata"=22,
                 "Pocillopora_acuta"=24,"Stylophora_pistillata"=23),
      name = "Species"
    ) +
    labs(x = paste0("PC1 (", percentVar[1], "%)"),
         y = paste0("PC2 (", percentVar[2], "%)"),
         title = title) +
    coord_fixed() +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

p1 <- plot_pca(pca_raw, tinfo, "log2(TPM+1)")

## ------------------------
## 2) PCA after limma::removeBatchEffect
## ------------------------
x_limma <- removeBatchEffect(log_mat, batch = tinfo$Species)
pca_limma <- prcomp(t(x_limma), scale. = TRUE)
p2 <- plot_pca(pca_limma, tinfo, "Log2(TPM+1) Limma")

## ------------------------
## 3) PCA after ComBat
## ------------------------
batch <- tinfo$Species
mod   <- model.matrix(~ timepoint, data = tinfo)  # preserve timepoint effect
x_combat <- ComBat(dat = log_mat, batch = batch, mod = mod,
                   par.prior = TRUE, prior.plots = FALSE)
pca_combat <- prcomp(t(x_combat), scale. = TRUE)
p3 <- plot_pca(pca_combat, tinfo, "Log2+1 PCA Combat")

## ------------------------
## Arrange side by side
## ------------------------
grid.arrange(p1, p2, p3, ncol = 3)

plot(p1)
plot(p2)
plot(p3)


suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(gridExtra)
})

# Make sure factors are set
tinfo$timepoint <- factor(tinfo$timepoint, levels = c("I","II","III"))
tinfo$Species   <- factor(tinfo$Species,
                          levels = c("Acropora_tenuis","Montipora_capitata",
                                     "Pocillopora_acuta","Stylophora_pistillata"))

# 1) RAW (no correction)
pca_raw <- prcomp(t(log_mat), scale. = TRUE)
p_raw   <- plot_pca(pca_raw, tinfo, "Raw log2(TPM+1)")

# 2) limma::removeBatchEffect (remove Species)
x_limma   <- removeBatchEffect(log_mat, batch = tinfo$Species)
pca_limma <- prcomp(t(x_limma), scale. = TRUE)
p_limma   <- plot_pca(pca_limma, tinfo, "After removeBatchEffect (Species)")

# 3) ComBat (remove Species, preserve timepoint)
batch      <- tinfo$Species
mod        <- model.matrix(~ timepoint, data = tinfo)
x_combat   <- ComBat(dat = log_mat, batch = batch, mod = mod,
                     par.prior = TRUE, prior.plots = FALSE)
pca_combat <- prcomp(t(x_combat), scale. = TRUE)
p_combat   <- plot_pca(pca_combat, tinfo, "After ComBat (Species; keep timepoint)")

# Show side by side
gridExtra::grid.arrange(p_raw, p_limma, p_combat, ncol = 3)

# Save individually
ggplot2::ggsave("PCA_LogTPM2_raw.png",     p_raw,    width = 7, height = 5.5, dpi = 600)
ggplot2::ggsave("PCA_LogTPM2_limma.png",   p_limma,  width = 7, height = 5.5, dpi = 600)
ggplot2::ggsave("PCA_LogTPM2_ComBat.png",  p_combat, width = 7, height = 5.5, dpi = 600)

# (Optional) export adjusted matrices for downstream use
# write.csv(x_limma,  "LogTPM2_removeBatchEffect_Species.csv")
# write.csv(x_combat, "LogTPM2_ComBat_Species_keepTimepoint.csv")

tinfo$timepoint <- factor(trimws(tinfo$timepoint), levels = c("I","II","III"))

plot_pca <- function(pca_obj, meta, title) {
  percentVar <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)))[1:2]
  df <- as.data.frame(pca_obj$x)
  df$timepoint <- factor(trimws(meta$timepoint), levels = c("I","II","III"))
  df$Species   <- meta$Species
  
  ggplot(df, aes(PC1, PC2)) +
    stat_ellipse(aes(group = timepoint, fill = timepoint),
                 geom = "polygon", alpha = 0.30, color = NA) +
    geom_point(aes(fill = timepoint, shape = Species),
               size = 4, color = "black", stroke = 0.8) +
    scale_fill_manual(values = c("I"="#e9a3c9","II"="lightgoldenrod2","III"="#a1d76a"),
                      name = "Timepoint", drop = FALSE) +
    scale_shape_manual(values = c("Acropora_tenuis"=21,"Montipora_capitata"=22,
                                  "Pocillopora_acuta"=24,"Stylophora_pistillata"=23),
                       name = "Species") +
    labs(x = paste0("PC1 (", percentVar[1], "%)"),
         y = paste0("PC2 (", percentVar[2], "%)"),
         title = title) +
    coord_fixed() +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}


# --- PCA on raw ---
pca_raw <- prcomp(t(log_mat), scale. = TRUE)
p1 <- plot_pca(pca_raw, tinfo, "Raw logTPM")

plot (p1)
# --- PCA on ComBat-corrected ---
pca_combat <- prcomp(t(combat_mat), scale. = TRUE)
p2 <- plot_pca(pca_combat, tinfo, "ComBat-corrected")
plot(p2)
# --- PCA on limma removeBatchEffect corrected ---
pca_limma <- prcomp(t(limma_mat), scale. = TRUE)
p3 <- plot_pca(pca_limma, tinfo, "limma::removeBatchEffect-corrected")
plot(p3)
# View them individually
print(p1)
print(p2)
print(p3)

# Or arrange side-by-side
library(patchwork)
(p1 | p2 | p3)
