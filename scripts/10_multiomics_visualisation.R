---
title: "multiomics_visualisation"
output: html_document
date: "2026-03-16"
---

Full genes Heatmap
```{r}

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

dat <- read_csv("multiomics_results/multi_omics_significant_soluble.csv", show_col_types = FALSE)

mean_cols <- c("mean_log2_D14", "mean_log2_D21", "mean_log2_M03",
               "mean_log2_M06", "mean_log2_M12", "mean_log2_M18")

mat <- dat %>%
  group_by(gene_symbol) %>%
  summarise(across(all_of(mean_cols), mean, na.rm = TRUE)) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

colnames(mat) <- c("D14", "D21", "M03", "M06", "M12", "M18")

# normalization and transposing

mat_scaled <- t(scale(t(mat)))

# Remove any genes with NA across all timepoints
mat_scaled <- mat_scaled[complete.cases(mat_scaled), ]

png("heatmap_multi_omics_soluble.png", width = 8, height = 14,
    units = "in", res = 300)

pheatmap(
  mat_scaled,
  cluster_rows   = TRUE,       
  cluster_cols   = FALSE,      # keep time points in order
  scale          = "none",    
  color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks         = seq(-2.5, 2.5, length.out = 101),
  border_color   = NA,
  fontsize_row   = 6,
  fontsize_col   = 11,
  main           = "Multi-omics significant genes\nProtein abundance (Z-score) across timepoints",
  angle_col      = 0,
  treeheight_row = 40
)

dev.off()
```

Multi_omics_integrated tables heatmap
```{r}
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# ------------------------
# 1️⃣ Read CSVs
# ------------------------
soluble <- read_csv("multiomics_results/multi_omics_significant_soluble.csv", show_col_types = FALSE)
cellular <- read_csv("multiomics_results/multi_omics_significant_cellular.csv", show_col_types = FALSE)


sol_cols <- c("mean_log2_D14", "mean_log2_D21", "mean_log2_M03",
              "mean_log2_M06", "mean_log2_M12", "mean_log2_M18")
cell_cols <- sol_cols

# ------------------------
# 2️⃣ Aggregate duplicates and rename
# ------------------------
soluble_agg <- soluble %>%
  group_by(gene_symbol) %>%
  summarise(across(all_of(sol_cols), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  rename_with(~sub("mean_log2_", "", .x), all_of(sol_cols))

cellular_agg <- cellular %>%
  group_by(gene_symbol) %>%
  summarise(across(all_of(cell_cols), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  rename_with(~sub("mean_log2_", "", .x), all_of(cell_cols))

# ------------------------
# 3️⃣ Merge fractions
# ------------------------
combined <- inner_join(soluble_agg, cellular_agg, by = "gene_symbol", suffix = c("_sol", "_cell"))

combined$variance <- apply(combined[,-1], 1, var, na.rm = TRUE)

top_genes <- combined %>%
  arrange(desc(variance)) %>%
  slice_head(n = 30) %>%
  dplyr::select(-variance) %>%
  as.data.frame()
# ------------------------
# 5️⃣ Attach GO Biological Process annotations
# ------------------------
top_genes_go <- top_genes %>%
  dplyr::left_join(
    as.data.frame(soluble) %>%
      dplyr::select(gene_symbol, go_biological_process) %>%
      dplyr::distinct(gene_symbol, .keep_all = TRUE),
    by = "gene_symbol"
  ) %>%
  as.data.frame()

# Simplify GO annotation: first term if multiple
top_genes_go$GO_BP_short <- sapply(strsplit(top_genes_go$go_biological_process, ";"), `[`, 1)

# ------------------------
# 6️⃣ Prepare matrix and scale
# ------------------------
mat <- top_genes_go %>%
  column_to_rownames("gene_symbol") %>%
  dplyr::select(dplyr::ends_with("_sol"), dplyr::ends_with("_cell")) %>%
  as.matrix()

# Scale and transpose (genes on x-axis, samples on y-axis)
mat_scaled <- scale(mat)
mat_scaled_t <- t(mat_scaled)

sample_order <- c("D14_sol", "D21_sol", "M03_sol", "M06_sol", "M12_sol", "M18_sol",
                  "D14_cell", "D21_cell", "M03_cell", "M06_cell", "M12_cell", "M18_cell")
mat_scaled_t <- mat_scaled_t[sample_order, ]
# ------------------------
# 7️⃣ Prepare annotations
# ------------------------
# Samples/timepoints annotation (rows now - on y-axis)
annotation_row <- data.frame(
  Fraction = rep(c("Soluble","Cellular"), each = 6)
)
rownames(annotation_row) <- rownames(mat_scaled_t)

# Genes annotation (columns now - on x-axis)
annotation_col <- data.frame(
  GO_BP = top_genes_go$GO_BP_short
)
rownames(annotation_col) <- colnames(mat_scaled_t)

# ------------------------
# 8️⃣ Heatmap color palette
# ------------------------
pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# ------------------------
# 9️⃣ Draw and export heatmap
# ------------------------
pheatmap(
  mat_scaled_t,
  cluster_rows = FALSE,    # ✅ Set to FALSE to preserve chronological order
  cluster_cols = TRUE,     # cluster genes (x-axis)
  color = pal,
  fontsize_row = 10,       # y-axis (samples)
  fontsize_col = 8,        # x-axis (genes)
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  main = "Top 30 proteins: Genes (x-axis) vs Samples (y-axis) with GO BP",
  filename = "top30_genes_xaxis_GO_heatmap_multiomics.pdf",
  width = 14,
  height = 10
)

```

PCA plots

```{r}

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
  out_file = file.path(PLOTS_DIR, "soluble_pca_samples.pdf")
)
plot_pca(
  prot_cell_norm, log2_cols_cell,
  title    = "PCA: Cellular fraction (filtered + normalized)",
  out_file = file.path(PLOTS_DIR, "cellular_pca_samples.pdff")
)
```


