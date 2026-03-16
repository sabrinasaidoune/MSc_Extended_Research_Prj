---
title: "GSEA and Time-Series Pathway Clustering"
output: html_document
date: "2026-03-16"
---

# ============================================================
# Script 11: GSEA + Time-Series Pathway Clustering Pipeline
# Mouse ovary proteomics (soluble & cellular) + transcriptomics
# Timepoints: D14, D21, M03, M06, M12, M18
# Comparison: young vs. old RNA-seq
# ============================================================

```{r setup}
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(ggridges)
  library(pheatmap)
  library(RColorBrewer)
  library(Mfuzz)
  library(Biobase)
})

set.seed(42)

# ── Paths ──────────────────────────────────────────────────
DATA_DIR    <- "results/tables"
OUT_TABLES  <- "results/gsea/tables"
OUT_PLOTS   <- "results/gsea/plots"

dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_PLOTS,  recursive = TRUE, showWarnings = FALSE)

TIMEPOINTS <- c("D14", "D21", "M03", "M06", "M12", "M18")
TP_PAIRS   <- list(
  c("D14", "D21"),
  c("D21", "M03"),
  c("M03", "M06"),
  c("M06", "M12"),
  c("M12", "M18")
)
```

# ============================================================
# HELPER FUNCTIONS
# ============================================================

```{r helpers}
# ── prepare_ranked_list ──────────────────────────────────────
# Convert gene symbols to Entrez IDs and create a named,
# sorted ranked vector for GSEA input.
#
# Args:
#   gene_symbols : character vector of gene symbols
#   fold_changes : numeric vector of log2 fold changes (same order)
#
# Returns: named numeric vector (Entrez ID → log2FC), sorted desc.
prepare_ranked_list <- function(gene_symbols, fold_changes) {
  stopifnot(length(gene_symbols) == length(fold_changes))

  df <- data.frame(
    symbol = gene_symbols,
    log2fc = fold_changes,
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(symbol), symbol != "", !is.na(log2fc), is.finite(log2fc))

  # Map gene symbols → Entrez IDs (mouse)
  mapping <- tryCatch(
    suppressMessages(
      bitr(df$symbol,
           fromType = "SYMBOL",
           toType   = "ENTREZID",
           OrgDb    = org.Mm.eg.db)
    ),
    error = function(e) {
      message("  [WARN] bitr failed: ", conditionMessage(e))
      data.frame(SYMBOL = character(), ENTREZID = character())
    }
  )

  if (nrow(mapping) == 0) {
    message("  [WARN] No Entrez mappings found; returning empty list")
    return(numeric(0))
  }

  df_mapped <- df %>%
    inner_join(mapping, by = c("symbol" = "SYMBOL")) %>%
    # Keep one entry per Entrez ID: highest |log2FC|
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) %>%
    ungroup()

  ranked <- setNames(df_mapped$log2fc, df_mapped$ENTREZID)
  sort(ranked, decreasing = TRUE)
}


# ── run_gsea_analysis ────────────────────────────────────────
# Run gseGO (BP) and gseKEGG on a ranked list.
#
# Args:
#   ranked_list : named numeric vector from prepare_ranked_list()
#   label       : string label used for messages
#   min_gs      : minimum gene set size (default 10)
#   max_gs      : maximum gene set size (default 500)
#   p_cutoff    : p-value cutoff (default 0.05)
#   n_perm      : number of permutations (default 1000 for speed)
#
# Returns: list with elements $go_bp and $kegg (or NULL if failed)
run_gsea_analysis <- function(ranked_list,
                              label    = "analysis",
                              min_gs   = 10,
                              max_gs   = 500,
                              p_cutoff = 0.05,
                              n_perm   = 1000) {

  if (length(ranked_list) < min_gs) {
    message("  [SKIP] ", label, ": ranked list too short (", length(ranked_list), ")")
    return(list(go_bp = NULL, kegg = NULL))
  }

  message("  Running GSEA for: ", label,
          "  (n=", length(ranked_list), " genes)")

  # ── GO Biological Process ──────────────────────────────────
  go_bp <- tryCatch({
    suppressMessages(
      gseGO(
        geneList     = ranked_list,
        OrgDb        = org.Mm.eg.db,
        ont          = "BP",
        keyType      = "ENTREZID",
        minGSSize    = min_gs,
        maxGSSize    = max_gs,
        pvalueCutoff = p_cutoff,
        pAdjustMethod = "BH",
        nPermSimple  = n_perm,
        verbose      = FALSE
      )
    )
  }, error = function(e) {
    message("  [WARN] gseGO failed for ", label, ": ", conditionMessage(e))
    NULL
  })

  n_go <- if (!is.null(go_bp)) nrow(as.data.frame(go_bp)) else 0
  message("    GO-BP significant pathways: ", n_go)

  # ── KEGG ───────────────────────────────────────────────────
  kegg <- tryCatch({
    suppressMessages(
      gseKEGG(
        geneList     = ranked_list,
        organism     = "mmu",
        minGSSize    = min_gs,
        maxGSSize    = max_gs,
        pvalueCutoff = p_cutoff,
        pAdjustMethod = "BH",
        nPermSimple  = n_perm,
        verbose      = FALSE
      )
    )
  }, error = function(e) {
    message("  [WARN] gseKEGG failed for ", label, ": ", conditionMessage(e))
    NULL
  })

  n_kegg <- if (!is.null(kegg)) nrow(as.data.frame(kegg)) else 0
  message("    KEGG significant pathways: ", n_kegg)

  list(go_bp = go_bp, kegg = kegg)
}


# ── create_temporal_matrix ───────────────────────────────────
# Build a pathway × timepoint NES matrix from a named list of
# GSEA result objects (one per timepoint comparison).
#
# Args:
#   gsea_list : named list; names are comparison labels
#               (e.g. "D14_D21"), values are run_gsea_analysis()
#               output lists (containing $go_bp and $kegg)
#   db        : which database to use, "go_bp" or "kegg"
#
# Returns: data.frame (pathways as rows, comparisons as columns,
#          NES values as cells; NA where pathway not significant)
create_temporal_matrix <- function(gsea_list, db = "go_bp") {

  all_nes <- lapply(names(gsea_list), function(nm) {
    res_obj <- gsea_list[[nm]][[db]]
    if (is.null(res_obj)) return(NULL)
    df <- as.data.frame(res_obj)
    if (nrow(df) == 0) return(NULL)
    df %>%
      select(ID, Description, NES) %>%
      rename(!!nm := NES)
  })

  all_nes <- Filter(Negate(is.null), all_nes)

  if (length(all_nes) == 0) {
    message("  [WARN] No GSEA results available for temporal matrix")
    return(data.frame())
  }

  # Full outer join so pathways not found in all comparisons get NA
  mat <- Reduce(function(a, b) full_join(a, b, by = c("ID", "Description")),
                all_nes)
  mat
}


# ── plot_temporal_clusters ───────────────────────────────────
# Visualise fuzzy cluster profiles as line plots.
#
# Args:
#   mfuzz_obj    : result of mfuzz() call (object of class ExpressionSet)
#   cl           : cluster result from mfuzz() (list with $cluster, $membership)
#   timepoint_labels : character vector for x-axis labels
#   title        : plot title
#   outfile      : path to save PDF (NULL = return ggplot object)
#
# Returns: ggplot object (and writes PDF if outfile is provided)
plot_temporal_clusters <- function(mfuzz_obj, cl,
                                   timepoint_labels = NULL,
                                   title  = "Temporal Pathway Clusters",
                                   outfile = NULL) {

  expr_mat  <- exprs(mfuzz_obj)
  n_cluster <- max(cl$cluster)

  if (is.null(timepoint_labels)) {
    timepoint_labels <- colnames(expr_mat)
  }

  # Build data frame for ggplot
  plot_data <- lapply(seq_len(nrow(expr_mat)), function(i) {
    pathway_name <- rownames(expr_mat)[i]
    cl_id        <- cl$cluster[i]
    membership   <- cl$membership[i, cl_id]
    data.frame(
      pathway    = pathway_name,
      cluster    = paste0("Cluster ", cl_id),
      membership = membership,
      timepoint  = timepoint_labels,
      NES        = as.numeric(expr_mat[i, ]),
      stringsAsFactors = FALSE
    )
  })

  plot_data <- bind_rows(plot_data)

  p <- ggplot(plot_data,
              aes(x = timepoint, y = NES, group = pathway,
                  colour = membership)) +
    geom_line(alpha = 0.6, linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~cluster, scales = "free_y") +
    scale_colour_gradient(low = "lightblue", high = "darkblue",
                          name = "Membership") +
    scale_x_discrete(limits = timepoint_labels) +
    labs(title = title, x = "Timepoint Comparison", y = "NES") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey90"))

  if (!is.null(outfile)) {
    ggsave(outfile, p,
           width = 12, height = ceiling(n_cluster / 2) * 4 + 1,
           limitsize = FALSE)
    message("  Saved cluster plot → ", outfile)
  }

  invisible(p)
}


# ── safe_save_csv ────────────────────────────────────────────
# Save a GSEA result object (or plain data frame) to CSV.
safe_save_csv <- function(obj, path) {
  if (is.null(obj)) {
    message("  [SKIP] NULL object; not saving ", basename(path))
    return(invisible(NULL))
  }
  df <- if (inherits(obj, "data.frame")) obj else as.data.frame(obj)
  write_csv(df, path, na = "")
  message("  Saved → ", path, "  (", nrow(df), " rows)")
}


# ── safe_dotplot ─────────────────────────────────────────────
# Save a dot plot for a single GSEA result object.
safe_dotplot <- function(gsea_obj, title, outfile, top_n = 15) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("  [SKIP] No results for dotplot: ", title)
    return(invisible(NULL))
  }
  p <- tryCatch(
    dotplot(gsea_obj, showCategory = top_n, title = title,
            font.size = 9) +
      theme(plot.title = element_text(size = 10)),
    error = function(e) {
      message("  [WARN] dotplot failed: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(p)) {
    ggsave(outfile, p, width = 8, height = 6)
    message("  Saved dotplot → ", outfile)
  }
  invisible(p)
}


# ── safe_ridgeplot ───────────────────────────────────────────
# Save a ridge plot for a single GSEA result object.
safe_ridgeplot <- function(gsea_obj, title, outfile, top_n = 20) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("  [SKIP] No results for ridgeplot: ", title)
    return(invisible(NULL))
  }
  p <- tryCatch(
    ridgeplot(gsea_obj, showCategory = top_n) +
      labs(title = title) +
      theme(plot.title = element_text(size = 10)),
    error = function(e) {
      message("  [WARN] ridgeplot failed: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(p)) {
    ht <- max(5, min(top_n * 0.4 + 2, 20))
    ggsave(outfile, p, width = 9, height = ht)
    message("  Saved ridgeplot → ", outfile)
  }
  invisible(p)
}


# ── safe_cnetplot ────────────────────────────────────────────
# Save a cnet enrichment map (gene-pathway network) for the top
# results of a single GSEA object.
safe_cnetplot <- function(gsea_obj, ranked_list, title, outfile,
                          top_n = 5) {
  if (is.null(gsea_obj) || nrow(as.data.frame(gsea_obj)) == 0) {
    message("  [SKIP] No results for cnetplot: ", title)
    return(invisible(NULL))
  }
  p <- tryCatch({
    top <- head(gsea_obj, top_n)
    cnetplot(top, foldChange = ranked_list,
             colorEdge = TRUE,
             cex_label_category = 0.7,
             cex_label_gene     = 0.5) +
      labs(title = title) +
      theme(plot.title = element_text(size = 10))
  }, error = function(e) {
    message("  [WARN] cnetplot failed: ", conditionMessage(e)); NULL
  })
  if (!is.null(p)) {
    ggsave(outfile, p, width = 12, height = 10)
    message("  Saved cnetplot → ", outfile)
  }
  invisible(p)
}
```

# ============================================================
# PART 1 — LOAD DATA
# ============================================================

```{r load_data}
message("=== Loading data ===")

soluble  <- read_csv(file.path(DATA_DIR, "multi_omics_integrated_soluble.csv"),
                     show_col_types = FALSE)
cellular <- read_csv(file.path(DATA_DIR, "multi_omics_integrated_cellular.csv"),
                     show_col_types = FALSE)
rna_full <- read_csv(file.path(DATA_DIR, "deseq2_full_results.csv"),
                     show_col_types = FALSE)

message("  Soluble proteins : ", nrow(soluble))
message("  Cellular proteins: ", nrow(cellular))
message("  RNA genes        : ", nrow(rna_full))

# Map Ensembl IDs → gene symbols for the full RNA dataset
rna_mapped <- tryCatch({
  suppressMessages(
    bitr(rna_full$gene,
         fromType = "ENSEMBL",
         toType   = c("SYMBOL", "ENTREZID"),
         OrgDb    = org.Mm.eg.db)
  )
}, error = function(e) {
  message("[WARN] Ensembl→Symbol mapping failed: ", conditionMessage(e))
  data.frame(ENSEMBL = character(), SYMBOL = character(),
             ENTREZID = character())
})

rna_annotated <- rna_full %>%
  inner_join(rna_mapped, by = c("gene" = "ENSEMBL")) %>%
  # Keep one row per gene (highest |log2FC| when duplicated)
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(log2FoldChange), is.finite(log2FoldChange))

message("  RNA genes after mapping: ", nrow(rna_annotated))
```

# ============================================================
# PART 1 — GSEA: TEMPORAL PROTEOMICS
# ============================================================

```{r gsea_temporal}
message("\n=== GSEA: Temporal Proteomics ===")

# Compute log2FC for each consecutive timepoint pair from mean_log2 columns.
# log2FC = mean_log2_later - mean_log2_earlier (difference of log2 values).

run_temporal_gsea <- function(prot_df, fraction_label) {
  results <- list()

  for (pair in TP_PAIRS) {
    tp_early <- pair[1]
    tp_late  <- pair[2]
    label    <- paste0(tp_early, "_", tp_late)

    col_early <- paste0("mean_log2_", tp_early)
    col_late  <- paste0("mean_log2_", tp_late)

    if (!col_early %in% names(prot_df) || !col_late %in% names(prot_df)) {
      message("  [SKIP] Missing columns for ", label)
      next
    }

    fc_df <- prot_df %>%
      filter(!is.na(gene_symbol), gene_symbol != "") %>%
      mutate(log2fc = !!sym(col_late) - !!sym(col_early)) %>%
      filter(!is.na(log2fc), is.finite(log2fc))

    ranked <- prepare_ranked_list(fc_df$gene_symbol, fc_df$log2fc)

    message(sprintf("  [%s] %s  n_genes=%d", fraction_label, label,
                    length(ranked)))

    gsea_res <- run_gsea_analysis(ranked, label = paste(fraction_label, label),
                                  n_perm = 1000)
    results[[label]] <- gsea_res
  }

  results
}

gsea_soluble  <- run_temporal_gsea(soluble,  "Soluble")
gsea_cellular <- run_temporal_gsea(cellular, "Cellular")
```

# ============================================================
# PART 1 — GSEA: RNA (young vs. old)
# ============================================================

```{r gsea_rna}
message("\n=== GSEA: RNA young vs. old ===")

ranked_rna <- setNames(rna_annotated$log2FoldChange, rna_annotated$ENTREZID)
ranked_rna <- sort(ranked_rna, decreasing = TRUE)

gsea_rna <- run_gsea_analysis(ranked_rna, label = "RNA_young_vs_old",
                              n_perm = 1000)
```

# ============================================================
# PART 1 — SAVE GSEA RESULTS AS CSV
# ============================================================

```{r save_gsea_csv}
message("\n=== Saving GSEA CSVs ===")

# ── Soluble temporal ──────────────────────────────────────────
sol_rows <- lapply(names(gsea_soluble), function(cmp) {
  go   <- gsea_soluble[[cmp]]$go_bp
  kegg <- gsea_soluble[[cmp]]$kegg
  bind_rows(
    if (!is.null(go))   as.data.frame(go)   %>% mutate(db = "GO_BP", comparison = cmp),
    if (!is.null(kegg)) as.data.frame(kegg) %>% mutate(db = "KEGG",  comparison = cmp)
  )
})
sol_df <- bind_rows(sol_rows)
safe_save_csv(sol_df, file.path(OUT_TABLES, "gsea_soluble_temporal.csv"))

# ── Cellular temporal ─────────────────────────────────────────
cell_rows <- lapply(names(gsea_cellular), function(cmp) {
  go   <- gsea_cellular[[cmp]]$go_bp
  kegg <- gsea_cellular[[cmp]]$kegg
  bind_rows(
    if (!is.null(go))   as.data.frame(go)   %>% mutate(db = "GO_BP", comparison = cmp),
    if (!is.null(kegg)) as.data.frame(kegg) %>% mutate(db = "KEGG",  comparison = cmp)
  )
})
cell_df <- bind_rows(cell_rows)
safe_save_csv(cell_df, file.path(OUT_TABLES, "gsea_cellular_temporal.csv"))

# ── RNA ───────────────────────────────────────────────────────
rna_rows <- bind_rows(
  if (!is.null(gsea_rna$go_bp))
    as.data.frame(gsea_rna$go_bp) %>% mutate(db = "GO_BP"),
  if (!is.null(gsea_rna$kegg))
    as.data.frame(gsea_rna$kegg)  %>% mutate(db = "KEGG")
)
safe_save_csv(rna_rows, file.path(OUT_TABLES, "gsea_rna_young_vs_old.csv"))
```

# ============================================================
# PART 1 — GSEA VISUALISATIONS
# ============================================================

```{r gsea_dotplots}
message("\n=== Generating dot plots ===")

# ── Soluble: one page per comparison (GO-BP) ──────────────────
sol_dp_pages <- lapply(names(gsea_soluble), function(cmp) {
  safe_dotplot(
    gsea_soluble[[cmp]]$go_bp,
    title   = paste0("Soluble GO-BP: ", cmp),
    outfile = NULL   # collect, will combine to PDF
  )
})
sol_dp_pages <- Filter(Negate(is.null), sol_dp_pages)

if (length(sol_dp_pages) > 0) {
  pdf(file.path(OUT_PLOTS, "gsea_dotplot_soluble_temporal.pdf"),
      width = 8, height = 6)
  for (p in sol_dp_pages) print(p)
  dev.off()
  message("  Saved → gsea_dotplot_soluble_temporal.pdf")
}

# ── Cellular: one page per comparison (GO-BP) ────────────────
cell_dp_pages <- lapply(names(gsea_cellular), function(cmp) {
  safe_dotplot(
    gsea_cellular[[cmp]]$go_bp,
    title   = paste0("Cellular GO-BP: ", cmp),
    outfile = NULL
  )
})
cell_dp_pages <- Filter(Negate(is.null), cell_dp_pages)

if (length(cell_dp_pages) > 0) {
  pdf(file.path(OUT_PLOTS, "gsea_dotplot_cellular_temporal.pdf"),
      width = 8, height = 6)
  for (p in cell_dp_pages) print(p)
  dev.off()
  message("  Saved → gsea_dotplot_cellular_temporal.pdf")
}

# ── RNA ───────────────────────────────────────────────────────
safe_dotplot(
  gsea_rna$go_bp,
  title   = "RNA GO-BP: young vs. old",
  outfile = file.path(OUT_PLOTS, "gsea_dotplot_rna_young_old.pdf")
)
```

```{r gsea_ridgeplots}
message("\n=== Generating ridge plots ===")

# Combine all soluble GO-BP results into one list for a combined ridge plot
sol_combined_go <- tryCatch({
  valid <- Filter(function(x) !is.null(x$go_bp) &&
                    nrow(as.data.frame(x$go_bp)) > 0,
                  gsea_soluble)
  if (length(valid) == 0) NULL else valid[[1]]$go_bp   # use first non-NULL
}, error = function(e) NULL)

safe_ridgeplot(
  sol_combined_go,
  title   = "Soluble GO-BP NES distribution",
  outfile = file.path(OUT_PLOTS, "gsea_ridgeplot_soluble.pdf")
)

cell_combined_go <- tryCatch({
  valid <- Filter(function(x) !is.null(x$go_bp) &&
                    nrow(as.data.frame(x$go_bp)) > 0,
                  gsea_cellular)
  if (length(valid) == 0) NULL else valid[[1]]$go_bp
}, error = function(e) NULL)

safe_ridgeplot(
  cell_combined_go,
  title   = "Cellular GO-BP NES distribution",
  outfile = file.path(OUT_PLOTS, "gsea_ridgeplot_cellular.pdf")
)
```

```{r gsea_cnetplots}
message("\n=== Generating enrichment maps (cnetplots) ===")

# Use the first comparison with results for each fraction
first_sol_cmp <- Filter(function(x) !is.null(x$go_bp) &&
                           nrow(as.data.frame(x$go_bp)) > 0,
                         gsea_soluble)

if (length(first_sol_cmp) > 0) {
  cmp_name <- names(first_sol_cmp)[1]
  # Rebuild the ranked list for this comparison to pass foldChange
  pair      <- strsplit(cmp_name, "_")[[1]]
  fc_df_s   <- soluble %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    mutate(log2fc = !!sym(paste0("mean_log2_", pair[2])) -
                    !!sym(paste0("mean_log2_", pair[1]))) %>%
    filter(!is.na(log2fc), is.finite(log2fc))
  ranked_s  <- prepare_ranked_list(fc_df_s$gene_symbol, fc_df_s$log2fc)

  safe_cnetplot(
    first_sol_cmp[[1]]$go_bp,
    ranked_list = ranked_s,
    title   = paste0("Soluble GO-BP enrichment map: ", cmp_name),
    outfile = file.path(OUT_PLOTS, "enrichment_maps_soluble.pdf")
  )
}
```

```{r gsea_summary_table}
message("\n=== Building GSEA cross-comparison summary table ===")

# Identify pathways significant in ≥2 comparisons (soluble)
if (nrow(sol_df) > 0) {
  sol_summary <- sol_df %>%
    filter(!is.na(NES)) %>%
    group_by(ID, Description, db) %>%
    summarise(
      n_comparisons = n_distinct(comparison),
      mean_NES      = mean(NES, na.rm = TRUE),
      min_p_adjust  = min(p.adjust, na.rm = TRUE),
      comparisons   = paste(sort(unique(comparison)), collapse = "; "),
      .groups = "drop"
    ) %>%
    filter(n_comparisons >= 2) %>%
    arrange(desc(n_comparisons), min_p_adjust)

  safe_save_csv(sol_summary,
                file.path(OUT_TABLES, "gsea_soluble_cross_comparison_summary.csv"))
  message("  Soluble pathways significant in ≥2 comparisons: ", nrow(sol_summary))
}
```

# ============================================================
# PART 2 — TIME-SERIES PATHWAY CLUSTERING
# ============================================================

```{r temporal_matrix}
message("\n=== Building temporal NES matrices ===")

# ── Soluble ────────────────────────────────────────────────────
nes_mat_sol <- create_temporal_matrix(gsea_soluble, db = "go_bp")

if (nrow(nes_mat_sol) > 0) {
  safe_save_csv(nes_mat_sol,
                file.path(OUT_TABLES, "pathway_timeseries_matrix_soluble.csv"))
  message("  Soluble NES matrix: ", nrow(nes_mat_sol), " pathways × ",
          ncol(nes_mat_sol) - 2, " comparisons")
}

# ── Cellular ───────────────────────────────────────────────────
nes_mat_cell <- create_temporal_matrix(gsea_cellular, db = "go_bp")

if (nrow(nes_mat_cell) > 0) {
  safe_save_csv(nes_mat_cell,
                file.path(OUT_TABLES, "pathway_timeseries_matrix_cellular.csv"))
  message("  Cellular NES matrix: ", nrow(nes_mat_cell), " pathways × ",
          ncol(nes_mat_cell) - 2, " comparisons")
}
```

```{r fuzzy_clustering}
message("\n=== Fuzzy c-means clustering (mfuzz) ===")

run_mfuzz_clustering <- function(nes_mat, fraction_label,
                                 n_clusters = 5, fuzz_m = 2) {

  if (is.null(nes_mat) || nrow(nes_mat) == 0) {
    message("  [SKIP] Empty NES matrix for ", fraction_label)
    return(list(eset = NULL, cl = NULL, clusters_df = NULL))
  }

  # ── Prepare numeric matrix ─────────────────────────────────
  comparison_cols <- setdiff(names(nes_mat), c("ID", "Description"))

  if (length(comparison_cols) < 2) {
    message("  [SKIP] Fewer than 2 comparisons for clustering")
    return(list(eset = NULL, cl = NULL, clusters_df = NULL))
  }

  mat <- nes_mat %>%
    column_to_rownames("Description") %>%
    select(all_of(comparison_cols)) %>%
    as.matrix()

  # Remove pathways with ≥50% NA
  na_frac <- rowMeans(is.na(mat))
  mat     <- mat[na_frac < 0.5, , drop = FALSE]

  if (nrow(mat) < n_clusters) {
    message("  [SKIP] Not enough pathways (", nrow(mat),
            ") for ", n_clusters, " clusters")
    return(list(eset = NULL, cl = NULL, clusters_df = NULL))
  }

  # Replace remaining NA with 0 (neutral enrichment)
  mat[is.na(mat)] <- 0

  # ── Standardise ────────────────────────────────────────────
  eset <- new("ExpressionSet",
              exprs = mat)
  eset <- standardise(eset)

  # ── Estimate optimal m ─────────────────────────────────────
  m_est <- mestimate(eset)
  message("  Estimated fuzziness m = ", round(m_est, 3),
          " (using m=", fuzz_m, ")")

  # ── Fuzzy c-means ──────────────────────────────────────────
  set.seed(42)
  cl <- mfuzz(eset, c = n_clusters, m = fuzz_m)

  message("  Cluster sizes (", fraction_label, "): ",
          paste(tabulate(cl$cluster), collapse = ", "))

  # ── Cluster assignments data frame ─────────────────────────
  clusters_df <- data.frame(
    Description = rownames(mat),
    cluster     = cl$cluster,
    stringsAsFactors = FALSE
  )

  # Add per-cluster membership scores
  membership_df <- as.data.frame(cl$membership)
  colnames(membership_df) <- paste0("membership_cl", seq_len(n_clusters))
  clusters_df <- bind_cols(clusters_df, membership_df)

  # Join back the ID column from original matrix
  clusters_df <- clusters_df %>%
    left_join(nes_mat %>% select(ID, Description) %>% distinct(),
              by = "Description")

  list(eset = eset, cl = cl, clusters_df = clusters_df)
}

mfuzz_sol  <- run_mfuzz_clustering(nes_mat_sol,  "Soluble",  n_clusters = 5)
mfuzz_cell <- run_mfuzz_clustering(nes_mat_cell, "Cellular", n_clusters = 5)

# ── Save cluster assignments ────────────────────────────────
if (!is.null(mfuzz_sol$clusters_df)) {
  safe_save_csv(mfuzz_sol$clusters_df,
                file.path(OUT_TABLES, "pathway_clusters_soluble.csv"))
}
if (!is.null(mfuzz_cell$clusters_df)) {
  safe_save_csv(mfuzz_cell$clusters_df,
                file.path(OUT_TABLES, "pathway_clusters_cellular.csv"))
}
```

```{r cluster_plots}
message("\n=== Generating cluster temporal profile plots ===")

# ── Soluble ────────────────────────────────────────────────────
if (!is.null(mfuzz_sol$eset) && !is.null(mfuzz_sol$cl)) {
  comparison_labels <- colnames(exprs(mfuzz_sol$eset))
  plot_temporal_clusters(
    mfuzz_obj        = mfuzz_sol$eset,
    cl               = mfuzz_sol$cl,
    timepoint_labels = comparison_labels,
    title            = "Soluble — Pathway Temporal Clusters",
    outfile          = file.path(OUT_PLOTS,
                                 "pathway_timeseries_clusters_soluble.pdf")
  )
}

# ── Cellular ───────────────────────────────────────────────────
if (!is.null(mfuzz_cell$eset) && !is.null(mfuzz_cell$cl)) {
  comparison_labels <- colnames(exprs(mfuzz_cell$eset))
  plot_temporal_clusters(
    mfuzz_obj        = mfuzz_cell$eset,
    cl               = mfuzz_cell$cl,
    timepoint_labels = comparison_labels,
    title            = "Cellular — Pathway Temporal Clusters",
    outfile          = file.path(OUT_PLOTS,
                                 "pathway_timeseries_clusters_cellular.pdf")
  )
}
```

```{r temporal_heatmaps}
message("\n=== Generating temporal pathway heatmaps ===")

save_nes_heatmap <- function(nes_mat, title, outfile) {
  if (is.null(nes_mat) || nrow(nes_mat) == 0) {
    message("  [SKIP] Empty matrix for heatmap: ", title)
    return(invisible(NULL))
  }

  comparison_cols <- setdiff(names(nes_mat), c("ID", "Description"))
  if (length(comparison_cols) == 0) return(invisible(NULL))

  mat <- nes_mat %>%
    column_to_rownames("Description") %>%
    select(all_of(comparison_cols)) %>%
    as.matrix()

  # Keep pathways with at least one non-NA value
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  if (nrow(mat) == 0) return(invisible(NULL))

  # Replace NAs with 0 for display
  mat[is.na(mat)] <- 0

  # Limit to top 50 by variance for readability
  if (nrow(mat) > 50) {
    rv   <- apply(mat, 1, var)
    mat  <- mat[order(rv, decreasing = TRUE)[seq_len(50)], , drop = FALSE]
  }

  ht <- max(6, min(nrow(mat) * 0.3 + 2, 24))

  tryCatch(
    pheatmap(
      mat,
      cluster_rows   = TRUE,
      cluster_cols   = FALSE,
      color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks         = seq(-3, 3, length.out = 101),
      border_color   = NA,
      fontsize_row   = 7,
      fontsize_col   = 9,
      main           = title,
      angle_col      = 45,
      na_col         = "grey80",
      filename       = outfile,
      width          = 10,
      height         = ht
    ),
    error = function(e)
      message("  [WARN] heatmap failed: ", conditionMessage(e))
  )
  message("  Saved heatmap → ", outfile)
}

save_nes_heatmap(
  nes_mat_sol,
  title   = "Soluble — Pathway NES across timepoint comparisons",
  outfile = file.path(OUT_PLOTS, "pathway_temporal_heatmap_soluble.pdf")
)

save_nes_heatmap(
  nes_mat_cell,
  title   = "Cellular — Pathway NES across timepoint comparisons",
  outfile = file.path(OUT_PLOTS, "pathway_temporal_heatmap_cellular.pdf")
)
```

# ============================================================
# PART 3 — INTEGRATED CROSS-OMICS ANALYSIS
# ============================================================

```{r cross_omics}
message("\n=== Cross-omics integration ===")

# Identify pathways significant in BOTH proteomics (sol or cell) and RNA.
# Build a comparison matrix: rows = shared pathways,
# columns = Soluble NES, Cellular NES, RNA NES, p-adjust values.

extract_gsea_summary <- function(gsea_result, fraction_label,
                                 db = "go_bp") {
  obj <- gsea_result[[db]]
  if (is.null(obj)) return(NULL)
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return(NULL)
  df %>%
    select(ID, Description, NES, p.adjust) %>%
    rename(
      !!paste0("NES_",      fraction_label) := NES,
      !!paste0("padj_",     fraction_label) := p.adjust
    )
}

# ── RNA summary ───────────────────────────────────────────────
rna_sum <- extract_gsea_summary(gsea_rna, "RNA")

# ── Soluble: pool across comparisons, keep best padj per pathway
sol_cross <- if (nrow(sol_df) > 0) {
  sol_df %>%
    filter(db == "GO_BP", !is.na(NES)) %>%
    group_by(ID, Description) %>%
    slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(ID, Description, NES, p.adjust) %>%
    rename(NES_Soluble = NES, padj_Soluble = p.adjust)
} else NULL

# ── Cellular: pool across comparisons, keep best padj per pathway
cell_cross <- if (nrow(cell_df) > 0) {
  cell_df %>%
    filter(db == "GO_BP", !is.na(NES)) %>%
    group_by(ID, Description) %>%
    slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(ID, Description, NES, p.adjust) %>%
    rename(NES_Cellular = NES, padj_Cellular = p.adjust)
} else NULL

# ── Merge all three ───────────────────────────────────────────
shared_pathways <- Reduce(
  function(a, b) {
    if (is.null(a) || is.null(b)) return(NULL)
    inner_join(a, b, by = c("ID", "Description"))
  },
  Filter(Negate(is.null), list(sol_cross, cell_cross, rna_sum))
)

if (!is.null(shared_pathways) && nrow(shared_pathways) > 0) {
  # Flag concordance (same direction = concordant)
  shared_pathways <- shared_pathways %>%
    mutate(
      concordant_sol_rna =
        case_when(
          !is.na(NES_Soluble) & !is.na(NES_RNA) ~
            sign(NES_Soluble) == sign(NES_RNA),
          TRUE ~ NA
        ),
      concordant_cell_rna =
        case_when(
          !is.na(NES_Cellular) & !is.na(NES_RNA) ~
            sign(NES_Cellular) == sign(NES_RNA),
          TRUE ~ NA
        )
    ) %>%
    arrange(padj_RNA)

  safe_save_csv(shared_pathways,
                file.path(OUT_TABLES, "shared_pathways_summary.csv"))
  message("  Shared pathways (proteomics + RNA): ", nrow(shared_pathways))
  message("  Concordant sol-RNA: ",
          sum(shared_pathways$concordant_sol_rna,  na.rm = TRUE))
  message("  Concordant cell-RNA: ",
          sum(shared_pathways$concordant_cell_rna, na.rm = TRUE))

  # ── Cross-omics comparison heatmap ─────────────────────────
  nes_cols <- intersect(c("NES_Soluble", "NES_Cellular", "NES_RNA"),
                        names(shared_pathways))
  if (length(nes_cols) >= 2) {
    hm_mat <- shared_pathways %>%
      slice_head(n = min(40, nrow(shared_pathways))) %>%
      column_to_rownames("Description") %>%
      select(all_of(nes_cols)) %>%
      as.matrix()

    hm_mat[is.na(hm_mat)] <- 0
    ht <- max(6, min(nrow(hm_mat) * 0.35 + 2, 20))

    tryCatch(
      pheatmap(
        hm_mat,
        cluster_rows   = TRUE,
        cluster_cols   = FALSE,
        color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
        breaks         = seq(-3, 3, length.out = 101),
        border_color   = NA,
        fontsize_row   = 7,
        fontsize_col   = 10,
        main           = "Cross-omics pathway NES comparison",
        filename       = file.path(OUT_PLOTS,
                                   "pathway_heatmap_comparison.pdf"),
        width  = 7,
        height = ht
      ),
      error = function(e)
        message("  [WARN] cross-omics heatmap failed: ", conditionMessage(e))
    )
    message("  Saved → pathway_heatmap_comparison.pdf")
  }
} else {
  message("  [INFO] No shared pathways found across all three omics layers")
  # Still create an empty summary CSV
  safe_save_csv(
    data.frame(
      ID = character(), Description = character(),
      NES_Soluble = numeric(), padj_Soluble = numeric(),
      NES_Cellular = numeric(), padj_Cellular = numeric(),
      NES_RNA = numeric(), padj_RNA = numeric()
    ),
    file.path(OUT_TABLES, "shared_pathways_summary.csv")
  )
}
```

# ============================================================
# PART 3 — SUMMARY REPORT
# ============================================================

```{r summary_report}
message("\n=== Writing summary report ===")

report_path <- file.path(OUT_TABLES, "gsea_timeseries_summary.txt")
sink(report_path)
cat("=============================================================\n")
cat("GSEA + Time-Series Pathway Clustering — Summary Report\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================\n\n")

# ── 1. Input data ─────────────────────────────────────────────
cat("INPUT DATA\n----------\n")
cat("Soluble proteins  :", nrow(soluble),  "\n")
cat("Cellular proteins :", nrow(cellular), "\n")
cat("RNA genes (mapped):", nrow(rna_annotated), "\n\n")

# ── 2. GSEA pathway counts ────────────────────────────────────
cat("SIGNIFICANT PATHWAYS PER COMPARISON\n------------------------------------\n")
cat("Soluble (GO-BP):\n")
for (cmp in names(gsea_soluble)) {
  n <- if (!is.null(gsea_soluble[[cmp]]$go_bp))
         nrow(as.data.frame(gsea_soluble[[cmp]]$go_bp)) else 0
  cat(sprintf("  %-12s %d\n", cmp, n))
}
cat("Soluble (KEGG):\n")
for (cmp in names(gsea_soluble)) {
  n <- if (!is.null(gsea_soluble[[cmp]]$kegg))
         nrow(as.data.frame(gsea_soluble[[cmp]]$kegg)) else 0
  cat(sprintf("  %-12s %d\n", cmp, n))
}
cat("Cellular (GO-BP):\n")
for (cmp in names(gsea_cellular)) {
  n <- if (!is.null(gsea_cellular[[cmp]]$go_bp))
         nrow(as.data.frame(gsea_cellular[[cmp]]$go_bp)) else 0
  cat(sprintf("  %-12s %d\n", cmp, n))
}
cat("RNA young vs. old:\n")
cat(sprintf("  GO-BP: %d\n",
  if (!is.null(gsea_rna$go_bp)) nrow(as.data.frame(gsea_rna$go_bp)) else 0))
cat(sprintf("  KEGG : %d\n\n",
  if (!is.null(gsea_rna$kegg)) nrow(as.data.frame(gsea_rna$kegg)) else 0))

# ── 3. Cluster characteristics ────────────────────────────────
cat("PATHWAY CLUSTER CHARACTERISTICS\n--------------------------------\n")

describe_clusters <- function(mfuzz_res, fraction_label, nes_mat) {
  if (is.null(mfuzz_res$cl)) {
    cat(fraction_label, ": clustering not performed\n")
    return(invisible(NULL))
  }
  cat(fraction_label, ":\n")
  cl_tbl <- table(mfuzz_res$cl$cluster)
  for (k in sort(as.integer(names(cl_tbl)))) {
    cat(sprintf("  Cluster %d: %d pathways\n", k, cl_tbl[k]))

    # Top 5 by membership in this cluster
    if (!is.null(mfuzz_res$clusters_df)) {
      mem_col <- paste0("membership_cl", k)
      top5    <- mfuzz_res$clusters_df %>%
        filter(cluster == k) %>%
        {if (mem_col %in% names(.)) arrange(., desc(!!sym(mem_col))) else .} %>%
        slice_head(n = 5) %>%
        pull(Description)
      if (length(top5) > 0)
        cat("    Top pathways:", paste(top5, collapse = "; "), "\n")
    }

    # Identify temporal pattern from mean NES profile
    if (!is.null(nes_mat) && nrow(nes_mat) > 0) {
      pathway_ids <- mfuzz_res$clusters_df %>%
        filter(cluster == k) %>% pull(Description)
      sub_mat <- nes_mat %>%
        filter(Description %in% pathway_ids) %>%
        select(-ID, -Description) %>%
        as.matrix()
      if (nrow(sub_mat) > 0) {
        mean_profile <- colMeans(sub_mat, na.rm = TRUE)
        direction    <- if (mean_profile[length(mean_profile)] >
                            mean_profile[1]) "increasing" else "decreasing"
        cat("    Temporal trend: ", direction, "\n")
      }
    }
  }
  cat("\n")
}

describe_clusters(mfuzz_sol,  "Soluble",  nes_mat_sol)
describe_clusters(mfuzz_cell, "Cellular", nes_mat_cell)

# ── 4. Cross-omics concordance ────────────────────────────────
cat("CROSS-OMICS CONCORDANCE\n-----------------------\n")
shared_path <- file.path(OUT_TABLES, "shared_pathways_summary.csv")
if (file.exists(shared_path)) {
  sp <- read_csv(shared_path, show_col_types = FALSE)
  cat("Shared pathways (sol + cell + RNA):", nrow(sp), "\n")
  if ("concordant_sol_rna" %in% names(sp))
    cat("Concordant Soluble–RNA:", sum(sp$concordant_sol_rna,  na.rm = TRUE), "\n")
  if ("concordant_cell_rna" %in% names(sp))
    cat("Concordant Cellular–RNA:", sum(sp$concordant_cell_rna, na.rm = TRUE), "\n")
} else {
  cat("Shared pathways file not found\n")
}

cat("\n=== End of report ===\n")
sink()

message("  Saved summary report → ", report_path)
message("\n=== Pipeline complete ===")
```
