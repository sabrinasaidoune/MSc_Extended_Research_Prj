# =============================================================================
# 08_multiomics_integration.R
# Mouse Ovary Multi-Omics Integration Pipeline
# Integrates DESeq2 transcriptomics (young vs. old, n=4/group) with
# Spectronaut proteomics (6 timepoints, n=5/timepoint)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# =============================================================================
# SECTION 1: Configuration and File Checks
# =============================================================================

# Input file paths
SOLUBLE_PROT_PATH  <- "proteomics/Singh-Soluble-20747-7/tables/Table-17.txt"
CELLULAR_PROT_PATH <- "proteomics/Singh-Cellular-20713-4-V2/tables/Table-19.txt"
RNA_PATH           <- "DESeq2_results_expr_filtered.csv"

# Output paths
OUTPUT_DIR <- "results"
TABLES_DIR <- file.path(OUTPUT_DIR, "tables")
PLOTS_DIR  <- file.path(OUTPUT_DIR, "plots")

# Expected timepoints in proteomics data
TIMEPOINTS <- c("D14", "D21", "M03", "M06", "M12", "M18")

# Create output directories if they do not exist
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR,  recursive = TRUE, showWarnings = FALSE)

# Validate required input files
required_files <- c(SOLUBLE_PROT_PATH, CELLULAR_PROT_PATH, RNA_PATH)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Missing required input files:\n  ", paste(missing, collapse = "\n  "))
}
message("✔ All input files found")

# =============================================================================
# SECTION 2: Proteomics Cleaning Function
# Handles both soluble and cellular fractions identically.
# =============================================================================

#' Clean and QC-filter a raw Spectronaut proteomics table.
#'
#' @param prot_path  Path to the tab-delimited Spectronaut output.
#' @param fraction   Character label for this fraction ("soluble" or "cellular").
#' @return           A tibble with cleaned, QC-filtered proteomics data.
clean_proteomics <- function(prot_path, fraction) {

  message(sprintf("\n--- Cleaning %s proteomics ---", fraction))

  prot_raw <- read_tsv(prot_path, show_col_types = FALSE)
  message(sprintf("  Loaded %d proteins (raw)", nrow(prot_raw)))

  # --- Identify column groups ------------------------------------------------

  # Metadata / annotation columns that must be present
  meta_cols <- c(
    "PG.ProteinGroups", "PG.Genes",
    "PG.ProteinDescriptions", "PG.UniProtIds",
    "PG.MolecularWeight",
    "PG.NrOfPrecursorsIdentified (Experiment-wide)",
    "GOMF name", "GOBP slim name", "KEGG name", "GOCC slim name"
  )

  # Validate required metadata columns exist
  missing_meta <- setdiff(meta_cols, names(prot_raw))
  if (length(missing_meta) > 0) {
    warning(sprintf("  [%s] Missing expected metadata columns: %s",
                    fraction, paste(missing_meta, collapse = ", ")))
    meta_cols <- intersect(meta_cols, names(prot_raw))
  }

  # Per-sample log2 abundance columns (e.g. "D14.R1.Quantity.log2")
  log2_cols <- grep("\\.Quantity\\.log2$", names(prot_raw), value = TRUE)
  if (length(log2_cols) == 0) stop(sprintf("[%s] No log2 abundance columns found.", fraction))

  # Upstream statistical result columns (comparisons pre-computed by Spectronaut)
  stat_cols <- grep(
    "^(p\\.value|adj\\.p\\.value|log2\\.fold\\.change|is\\.significant|t\\.statistic)\\.",
    names(prot_raw), value = TRUE
  )

  # Per-sample precursor counts used for quantification QC
  precursor_quant_cols <- grep("\\.NrOfPrecursorsUsedForQuantification$",
                               names(prot_raw), value = TRUE)

  # --- Select, extract, and filter -------------------------------------------

  prot_clean <- prot_raw %>%
    select(all_of(c(meta_cols, log2_cols, stat_cols))) %>%

    # Extract first (primary) entry from semicolon-delimited fields
    # str_extract("^[^;]+") is more efficient than str_split_fixed when only
    # the first token is needed — no full matrix allocation.
    mutate(
      gene_symbol         = na_if(str_trim(str_extract(PG.Genes, "^[^;]+")), ""),
      uniprot_id          = str_extract(PG.UniProtIds, "^[^;]+"),
      protein_description = str_extract(PG.ProteinDescriptions, "^[^;]+"),
      mol_weight_da       = suppressWarnings(
                              as.numeric(str_extract(PG.MolecularWeight, "^[^;]+"))
                            )
    ) %>%

    # Replace NaN with NA in abundance columns (NaN arises from log2(0))
    mutate(across(all_of(log2_cols), ~ na_if(., NaN))) %>%

    # --- QC metrics ---
    mutate(
      n_samples_detected   = rowSums(!is.na(across(all_of(log2_cols)))),
      pct_samples_detected = round(100 * n_samples_detected / length(log2_cols), 1),

      # Median precursors across samples; NA-safe; treats 0 as missing.
      # Computed outside mutate() for clarity, then attached as a column.
      median_precursors_used = if (length(precursor_quant_cols) > 0)
        apply(prot_raw[, precursor_quant_cols, drop = FALSE], 1,
              function(x) median(x[x > 0], na.rm = TRUE))
      else NA_real_
    )

  # Log proteins removed by precursor filter
  n_before <- nrow(prot_clean)
  prot_clean <- prot_clean %>%
    filter(is.na(median_precursors_used) | median_precursors_used >= 2)
  message(sprintf("  Precursor filter (>=2): %d -> %d proteins",
                  n_before, nrow(prot_clean)))

  # --- Per-timepoint mean log2 abundance ------------------------------------
  for (tp in TIMEPOINTS) {
    cols_tp <- grep(paste0("^", tp, "\\..*\\.Quantity\\.log2$"),
                    names(prot_clean), value = TRUE)
    if (length(cols_tp) > 0) {
      prot_clean[[paste0("mean_log2_", tp)]] <-
        rowMeans(prot_clean[, cols_tp, drop = FALSE], na.rm = TRUE)
    } else {
      warning(sprintf("  [%s] No log2 columns found for timepoint %s", fraction, tp))
      prot_clean[[paste0("mean_log2_", tp)]] <- NA_real_
    }
  }

  # --- Final column selection and renaming ----------------------------------
  prot_clean <- prot_clean %>%
    rename(
      protein_groups        = PG.ProteinGroups,
      n_precursors_exp      = `PG.NrOfPrecursorsIdentified (Experiment-wide)`,
      go_molecular_function = `GOMF name`,
      go_biological_process = `GOBP slim name`,
      go_cellular_component = `GOCC slim name`,
      kegg_pathway          = `KEGG name`
    ) %>%
    # Add fraction label for provenance tracking
    mutate(fraction = fraction) %>%
    select(
      fraction,
      gene_symbol, uniprot_id, protein_groups, protein_description,
      mol_weight_da, n_precursors_exp,
      n_samples_detected, pct_samples_detected, median_precursors_used,
      starts_with("mean_log2_"),
      all_of(log2_cols),
      all_of(stat_cols),
      go_molecular_function, go_biological_process,
      go_cellular_component, kegg_pathway
    )

  message(sprintf("  ✔ Final %s proteomics: %d proteins", fraction, nrow(prot_clean)))
  prot_clean
}

# =============================================================================
# SECTION 3: Load and Clean Both Proteomics Fractions
# =============================================================================

prot_soluble  <- clean_proteomics(SOLUBLE_PROT_PATH,  "soluble")
prot_cellular <- clean_proteomics(CELLULAR_PROT_PATH, "cellular")

write_csv(prot_soluble,  file.path(TABLES_DIR, "soluble_proteomics_clean.csv"),  na = "")
write_csv(prot_cellular, file.path(TABLES_DIR, "cellular_proteomics_clean.csv"), na = "")
message("\n✔ Wrote soluble_proteomics_clean.csv and cellular_proteomics_clean.csv")

# =============================================================================
# SECTION 4: Load and Prepare RNA-seq (DESeq2) Data
# =============================================================================

message("\n--- Loading RNA-seq data ---")

rna <- read_csv(RNA_PATH, show_col_types = FALSE)

# Normalise column names: accept either 'mgi_symbol' or 'gene_symbol'
if ("mgi_symbol" %in% names(rna) && !"gene_symbol" %in% names(rna)) {
  rna <- rename(rna, gene_symbol = mgi_symbol)
}

# Remove row-index column if present (artifact of write.csv with row.names)
rna <- select(rna, -any_of("...1"))

# Validate required columns
rna_required <- c("gene_symbol", "baseMean", "log2FoldChange", "padj")
missing_rna <- setdiff(rna_required, names(rna))
if (length(missing_rna) > 0) {
  stop("RNA data missing required columns: ", paste(missing_rna, collapse = ", "))
}

# Standardise RNA column names for clarity in merged output
rna <- rna %>%
  rename(
    rna_baseMean = baseMean,
    rna_log2FC   = log2FoldChange,
    rna_padj     = padj
  )
# Conditionally rename optional columns that may or may not be present
if ("lfcSE"          %in% names(rna)) rna <- rename(rna, rna_lfcSE       = lfcSE)
if ("pvalue"         %in% names(rna)) rna <- rename(rna, rna_pvalue      = pvalue)
if ("ensembl_gene_id" %in% names(rna)) rna <- rename(rna, rna_ensembl_id = ensembl_gene_id)

message(sprintf("  ✔ RNA data: %d genes loaded", nrow(rna)))

# =============================================================================
# SECTION 5: Integration Helper and Significance Functions
# =============================================================================

#' Determine which proteins are significant in at least one upstream comparison.
#'
#' Spectronaut exports boolean `is.significant.*` columns for each pairwise
#' comparison. This function returns TRUE if any of those columns is TRUE,
#' avoiding hard-coded comparison names.
#'
#' @param df  A proteomics data frame with is.significant.* columns.
#' @return    Logical vector of length nrow(df).
any_significant_protein <- function(df) {
  sig_cols <- grep("^is\\.significant\\.", names(df), value = TRUE)
  if (length(sig_cols) == 0) {
    message("  No is.significant.* columns found; treating all proteins as significant")
    return(rep(TRUE, nrow(df)))
  }
  # Convert to logical matrix, coercing character "TRUE"/"FALSE" if needed
  sig_mat <- df[, sig_cols, drop = FALSE] %>%
    mutate(across(everything(), ~ as.logical(.))) %>%
    as.matrix()
  rowSums(sig_mat, na.rm = TRUE) >= 1
}

#' Integrate a proteomics fraction with RNA-seq results.
#'
#' Steps:
#'   1. Drop proteins without gene symbols.
#'   2. Inner-join on gene_symbol (keeps only genes present in both datasets).
#'   3. Deduplicate: if multiple RNA rows match a gene, keep the one with
#'      the lowest adjusted p-value.
#'   4. Apply significance filters:
#'        - Protein: significant in >= 1 Spectronaut comparison.
#'        - RNA: rna_padj < 0.05.
#'
#' @param prot_df    Cleaned proteomics data frame.
#' @param rna_df     Cleaned RNA-seq data frame.
#' @param fraction   Label used in log messages.
#' @return           Integrated, significance-filtered data frame.
integrate_omics <- function(prot_df, rna_df, fraction) {

  message(sprintf("\n--- Integrating %s fraction with RNA ---", fraction))

  # Remove proteins without a mappable gene symbol
  prot_df <- filter(prot_df, !is.na(gene_symbol) & gene_symbol != "")
  message(sprintf("  Proteins with gene symbols: %d", nrow(prot_df)))

  # Inner join; many-to-many is expected (gene -> multiple protein groups or
  # multiple RNA entries); warn if the join is unexpectedly large.
  joined <- inner_join(prot_df, rna_df, by = "gene_symbol",
                       relationship = "many-to-many")
  message(sprintf("  After inner join: %d rows", nrow(joined)))

  # Deduplicate: per (gene_symbol, protein_groups) pair keep the RNA row with
  # the lowest rna_padj; ties broken by first occurrence.
  integrated <- joined %>%
    group_by(gene_symbol, protein_groups) %>%
    slice_min(order_by = rna_padj, n = 1, with_ties = FALSE) %>%
    ungroup()
  message(sprintf("  After RNA deduplication: %d rows", nrow(integrated)))

  # Flag proteins significant in >= 1 Spectronaut comparison
  integrated <- integrated %>%
    mutate(prot_any_significant = any_significant_protein(.))

  # Apply dual significance filter
  integrated_sig <- integrated %>%
    filter(prot_any_significant & !is.na(rna_padj) & rna_padj < 0.05)
  message(sprintf("  After significance filter (prot + rna_padj<0.05): %d rows",
                  nrow(integrated_sig)))

  integrated_sig
}

# =============================================================================
# SECTION 6: Run Integration for Both Fractions
# =============================================================================

integrated_soluble  <- integrate_omics(prot_soluble,  rna, "soluble")
integrated_cellular <- integrate_omics(prot_cellular, rna, "cellular")

write_csv(integrated_soluble,
          file.path(TABLES_DIR, "multi_omics_integrated_soluble.csv"),  na = "")
write_csv(integrated_cellular,
          file.path(TABLES_DIR, "multi_omics_integrated_cellular.csv"), na = "")
message("\n✔ Wrote multi_omics_integrated_soluble.csv and multi_omics_integrated_cellular.csv")

# =============================================================================
# SECTION 7: Normalization and Imputation
# =============================================================================

#' Median-centre and impute missing values in a proteomics data frame.
#'
#' Normalization: subtract per-sample median (log2 space) to correct for
#' global loading differences between samples.
#' Imputation: replace remaining NA values with (global minimum - 0.1).
#' This conservative "below-detection" strategy reflects the likely cause of
#' missingness (proteins below the detection limit) rather than random absence.
#'
#' @param prot_df   Cleaned proteomics data frame.
#' @param log2_cols Character vector of log2 abundance column names.
#' @return          Data frame with normalized and imputed abundance columns.
normalise_and_impute <- function(prot_df, log2_cols) {

  mat <- as.matrix(prot_df[, log2_cols])

  # --- 50%-detection filter: remove proteins missing in >50% of samples ------
  keep <- rowSums(!is.na(mat)) >= (ncol(mat) / 2)
  mat  <- mat[keep, ]
  prot_df <- prot_df[keep, ]
  message(sprintf("  50%% detection filter: %d proteins retained", nrow(prot_df)))

  # --- Median-centre per sample ----------------------------------------------
  sample_medians <- apply(mat, 2, median, na.rm = TRUE)
  mat_norm <- sweep(mat, 2, sample_medians, "-")

  # --- Impute remaining NAs with (global min - 0.1) -------------------------
  impute_val <- min(mat_norm, na.rm = TRUE) - 0.1
  mat_norm[is.na(mat_norm)] <- impute_val

  prot_df[, log2_cols] <- mat_norm
  prot_df
}

# Identify log2 columns in each fraction
log2_cols_sol  <- grep("\\.Quantity\\.log2$", names(prot_soluble),  value = TRUE)
log2_cols_cell <- grep("\\.Quantity\\.log2$", names(prot_cellular), value = TRUE)

message("\n--- Normalising soluble proteomics ---")
prot_sol_norm  <- normalise_and_impute(prot_soluble,  log2_cols_sol)

message("--- Normalising cellular proteomics ---")
prot_cell_norm <- normalise_and_impute(prot_cellular, log2_cols_cell)

# =============================================================================
# SECTION 8: PCA Visualisation Function
# =============================================================================

#' Generate and save a PCA plot from a proteomics abundance matrix.
#'
#' @param prot_df   Normalised proteomics data frame.
#' @param log2_cols Log2 abundance column names.
#' @param title     Plot title.
#' @param out_file  Output PDF path.
plot_pca <- function(prot_df, log2_cols, title, out_file) {

  mat <- as.matrix(prot_df[, log2_cols])

  # PCA on transposed matrix: rows = samples, columns = proteins
  pca_res <- prcomp(t(mat), scale. = TRUE)
  pct_var <- round(summary(pca_res)$importance[2, ] * 100, 1)

  pca_df <- as.data.frame(pca_res$x[, 1:2])
  pca_df$sample    <- rownames(pca_df)
  # Extract timepoint prefix (e.g. "D14" from "D14.R1.20747.Quantity.log2")
  pca_df$timepoint <- str_extract(pca_df$sample, "^[A-Z][A-Z0-9]+")

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = timepoint, label = sample)) +
    geom_point(size = 3, alpha = 0.85) +
    theme_bw(base_size = 12) +
    labs(
      title  = title,
      x      = paste0("PC1 (", pct_var[1], "% variance)"),
      y      = paste0("PC2 (", pct_var[2], "% variance)"),
      colour = "Timepoint"
    ) +
    theme(legend.position = "right")

  ggsave(out_file, plot = p, width = 7, height = 5, dpi = 300)
  message(sprintf("  ✔ Saved PCA: %s", out_file))
}

message("\n--- Generating PCA plots ---")
plot_pca(
  prot_sol_norm, log2_cols_sol,
  title    = "PCA: Soluble fraction (normalized + imputed)",
  out_file = file.path(PLOTS_DIR, "pca_samples.pdf")
)
plot_pca(
  prot_cell_norm, log2_cols_cell,
  title    = "PCA: Cellular fraction (filtered + normalized)",
  out_file = file.path(PLOTS_DIR, "qc_pca_filtered_normalized.pdf")
)

# =============================================================================
# SECTION 9: Heatmap Generation Function
# =============================================================================

#' Generate a Z-score heatmap of protein abundance across timepoints.
#'
#' @param prot_df      Data frame containing mean_log2_* columns and gene_symbol.
#' @param mean_cols    Character vector of mean_log2_* column names to plot.
#' @param top_n        Number of most-variable proteins to display.
#' @param title        Heatmap title.
#' @param out_file     Output PDF path.
#' @param col_labels   Optional vector of display labels for columns.
make_heatmap <- function(prot_df, mean_cols, top_n = 30,
                         title, out_file, col_labels = NULL) {

  # Aggregate duplicate gene symbols (take mean)
  mat_df <- prot_df %>%
    filter(!is.na(gene_symbol)) %>%
    group_by(gene_symbol) %>%
    summarise(across(all_of(mean_cols), function(x) mean(x, na.rm = TRUE)),
              .groups = "drop")

  mat <- mat_df %>%
    column_to_rownames("gene_symbol") %>%
    as.matrix()

  # Select top_n most variable proteins by row variance
  row_vars <- apply(mat, 1, var, na.rm = TRUE)
  top_idx  <- order(row_vars, decreasing = TRUE)[seq_len(min(top_n, nrow(mat)))]
  mat      <- mat[top_idx, , drop = FALSE]

  # Remove proteins with all-NA after subsetting
  mat <- mat[complete.cases(mat), , drop = FALSE]

  if (nrow(mat) == 0) {
    warning(sprintf("No complete rows for heatmap: %s", title))
    return(invisible(NULL))
  }

  # Row-scale (Z-score per protein) so colours reflect relative changes
  mat_scaled <- t(scale(t(mat)))
  mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]

  if (!is.null(col_labels)) colnames(mat_scaled) <- col_labels

  pdf(out_file, width = 8, height = max(6, ceiling(nrow(mat_scaled) / 4)))
  pheatmap(
    mat_scaled,
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,       # preserve temporal order
    scale          = "none",      # already Z-scored above
    color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks         = seq(-2.5, 2.5, length.out = 101),
    border_color   = NA,
    fontsize_row   = max(4, min(8, 120 / nrow(mat_scaled))),
    fontsize_col   = 11,
    main           = title,
    angle_col      = 0,
    treeheight_row = 30
  )
  dev.off()
  message(sprintf("  ✔ Saved heatmap: %s", out_file))
}

# =============================================================================
# SECTION 10: Heatmap 1 - Top 30 Soluble Proteins
# =============================================================================

message("\n--- Generating Heatmap 1: Top 30 soluble proteins ---")

mean_cols <- paste0("mean_log2_", TIMEPOINTS)

make_heatmap(
  prot_df    = prot_soluble,
  mean_cols  = mean_cols,
  top_n      = 30,
  title      = "Top 30 variable proteins (soluble fraction)\nZ-score log2 abundance",
  out_file   = file.path(PLOTS_DIR, "heatmap_top30_soluble.pdf"),
  col_labels = TIMEPOINTS
)

# =============================================================================
# SECTION 11: Heatmap 2 - Soluble vs. Cellular Comparison
# =============================================================================

message("--- Generating Heatmap 2: Soluble vs. cellular comparison ---")

#' Build a combined matrix of soluble and cellular mean log2 abundance,
#' select the top_n most variable genes, and produce a side-by-side heatmap.
#'
#' @param sol_df   Cleaned soluble proteomics data frame.
#' @param cell_df  Cleaned cellular proteomics data frame.
#' @param top_n    Number of genes to display.
#' @param out_file Output PDF path.
heatmap_sol_vs_cell <- function(sol_df, cell_df, top_n = 30, out_file) {

  mean_cols <- paste0("mean_log2_", TIMEPOINTS)

  # Aggregate per gene (mean across duplicate protein groups)
  agg <- function(df) {
    df %>%
      filter(!is.na(gene_symbol)) %>%
      group_by(gene_symbol) %>%
      summarise(across(all_of(mean_cols), function(x) mean(x, na.rm = TRUE)),
                .groups = "drop")
  }

  sol_agg  <- agg(sol_df)  %>% rename_with(~ paste0(.x, "_sol"),  all_of(mean_cols))
  cell_agg <- agg(cell_df) %>% rename_with(~ paste0(.x, "_cell"), all_of(mean_cols))

  combined <- inner_join(sol_agg, cell_agg, by = "gene_symbol")

  # Select top_n genes by variance across all 12 mean columns
  value_cols  <- setdiff(names(combined), "gene_symbol")
  combined$var_all <- apply(combined[, value_cols], 1, var, na.rm = TRUE)

  top_genes <- combined %>%
    arrange(desc(var_all)) %>%
    slice_head(n = top_n) %>%
    select(-var_all)

  mat <- top_genes %>%
    column_to_rownames("gene_symbol") %>%
    as.matrix()

  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) == 0) {
    warning("No complete rows for soluble vs cellular heatmap")
    return(invisible(NULL))
  }

  mat_scaled <- t(scale(t(mat)))
  mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]

  # Column annotation: which fraction each column belongs to
  annotation_col <- data.frame(
    Fraction = rep(c("Soluble", "Cellular"), each = length(TIMEPOINTS))
  )
  rownames(annotation_col) <- colnames(mat_scaled)

  # Attach GO Biological Process (first term only) as row annotation
  go_lookup <- bind_rows(
    select(sol_df,  gene_symbol, go_biological_process),
    select(cell_df, gene_symbol, go_biological_process)
  ) %>%
    filter(!is.na(go_biological_process) & go_biological_process != "") %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    mutate(
      # Extract first semicolon-delimited GO term; truncate to 40 chars for legibility
      GO_BP = str_trunc(str_extract(go_biological_process, "^[^;]+"), 40)
    )

  annotation_row <- data.frame(
    GO_BP = go_lookup$GO_BP[match(rownames(mat_scaled), go_lookup$gene_symbol)]
  )
  rownames(annotation_row) <- rownames(mat_scaled)
  # Replace NA annotations with empty string so pheatmap renders cleanly
  annotation_row$GO_BP[is.na(annotation_row$GO_BP)] <- ""

  pdf(out_file, width = 12, height = max(6, ceiling(nrow(mat_scaled) / 3)))
  pheatmap(
    mat_scaled,
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,
    scale          = "none",
    color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks         = seq(-2.5, 2.5, length.out = 101),
    border_color   = NA,
    fontsize_row   = max(4, min(8, 120 / nrow(mat_scaled))),
    fontsize_col   = 9,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    main           = "Top 30 proteins: Soluble vs. Cellular fractions\nZ-score log2 abundance",
    angle_col      = 0,
    treeheight_row = 30
  )
  dev.off()
  message(sprintf("  ✔ Saved heatmap: %s", out_file))
}

heatmap_sol_vs_cell(
  sol_df   = prot_soluble,
  cell_df  = prot_cellular,
  top_n    = 30,
  out_file = file.path(PLOTS_DIR, "heatmap_soluble_vs_cellular.pdf")
)

# =============================================================================
# SECTION 12: Integration Summary Report
# =============================================================================

message("\n--- Writing integration summary ---")

#' Collect summary statistics for one fraction.
fraction_summary <- function(prot_raw_df, prot_clean_df, integrated_df, fraction) {
  sig_cols  <- grep("^is\\.significant\\.", names(prot_clean_df), value = TRUE)
  n_sig_any <- if (length(sig_cols) > 0) {
    sig_mat <- prot_clean_df[, sig_cols, drop = FALSE] %>%
      mutate(across(everything(), ~ as.logical(.)))
    sum(rowSums(sig_mat, na.rm = TRUE) >= 1, na.rm = TRUE)
  } else NA_integer_

  list(
    fraction            = fraction,
    n_proteins_raw      = nrow(prot_raw_df),
    n_proteins_clean    = nrow(prot_clean_df),
    n_proteins_sig_any  = n_sig_any,
    # integrated_df is already filtered for both prot + RNA significance
    n_integrated_genes  = nrow(integrated_df)
  )
}

# Reload raw counts for summary (nrow before filtering)
prot_soluble_raw  <- read_tsv(SOLUBLE_PROT_PATH,  show_col_types = FALSE)
prot_cellular_raw <- read_tsv(CELLULAR_PROT_PATH, show_col_types = FALSE)

sol_summary  <- fraction_summary(prot_soluble_raw,  prot_soluble,  integrated_soluble,  "soluble")
cell_summary <- fraction_summary(prot_cellular_raw, prot_cellular, integrated_cellular, "cellular")

summary_lines <- c(
  "=============================================================",
  "  Multi-Omics Integration Summary",
  sprintf("  Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "=============================================================",
  "",
  "--- RNA-seq (DESeq2) ---",
  sprintf("  Total genes loaded      : %d", nrow(rna)),
  sprintf("  Genes with padj < 0.05  : %d",
          sum(!is.na(rna$rna_padj) & rna$rna_padj < 0.05)),
  "",
  "--- Soluble Fraction ---",
  sprintf("  Proteins (raw)          : %d", sol_summary$n_proteins_raw),
  sprintf("  Proteins (QC-filtered)  : %d", sol_summary$n_proteins_clean),
  sprintf("  Proteins (sig >= 1 cmp) : %s",
          ifelse(is.na(sol_summary$n_proteins_sig_any), "N/A",
                 sol_summary$n_proteins_sig_any)),
  sprintf("  Genes in integrated set : %d", sol_summary$n_integrated_genes),
  "",
  "--- Cellular Fraction ---",
  sprintf("  Proteins (raw)          : %d", cell_summary$n_proteins_raw),
  sprintf("  Proteins (QC-filtered)  : %d", cell_summary$n_proteins_clean),
  sprintf("  Proteins (sig >= 1 cmp) : %s",
          ifelse(is.na(cell_summary$n_proteins_sig_any), "N/A",
                 cell_summary$n_proteins_sig_any)),
  sprintf("  Genes in integrated set : %d", cell_summary$n_integrated_genes),
  "",
  "--- Timepoints ---",
  paste0("  ", paste(TIMEPOINTS, collapse = ", ")),
  "",
  "=============================================================",
  "  Output files",
  "=============================================================",
  paste0("  ", file.path(TABLES_DIR, "soluble_proteomics_clean.csv")),
  paste0("  ", file.path(TABLES_DIR, "cellular_proteomics_clean.csv")),
  paste0("  ", file.path(TABLES_DIR, "multi_omics_integrated_soluble.csv")),
  paste0("  ", file.path(TABLES_DIR, "multi_omics_integrated_cellular.csv")),
  paste0("  ", file.path(PLOTS_DIR,  "heatmap_top30_soluble.pdf")),
  paste0("  ", file.path(PLOTS_DIR,  "heatmap_soluble_vs_cellular.pdf")),
  paste0("  ", file.path(PLOTS_DIR,  "pca_samples.pdf")),
  paste0("  ", file.path(PLOTS_DIR,  "qc_pca_filtered_normalized.pdf")),
  ""
)

summary_path <- file.path(TABLES_DIR, "integration_summary.txt")
writeLines(summary_lines, summary_path)
message(sprintf("✔ Wrote integration summary: %s", summary_path))
message("\n=== Pipeline complete ===")
