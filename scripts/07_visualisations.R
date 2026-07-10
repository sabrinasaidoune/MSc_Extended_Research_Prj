---
title: "visualisation"
output: html_document
date: "2026-03-10"
---
PCA Plot

```{r}
# Load DESeq2 if not already loaded
library(DESeq2)

# Run DESeq to estimate size factors, dispersion, and fit the model
dds <- DESeq(dds)

# Get the results table 
res <- results(dds)

# shrink log2 fold changes for more stable estimates
resLFC <- lfcShrink(dds, coef=2, type="apeglm")  

# Convert results to a data frame for easier handling
res_df <- as.data.frame(res)

# Order results by adjusted p-value (most significant first)
res_df <- res_df[order(res_df$padj), ]

# Optional: extract only significant genes (adjusted p < 0.05)
sig_genes <- subset(res_df, padj < 0.05)

# Variance stabilizing transformation for visualization 
rld <- rlog(dds, blind=FALSE)

# PCA plot using sample metadata
plotPCA(rld, intgroup="age")

```

Volcano Plot


Load libraries
```{r}

#Load libraries

library(EnhancedVolcano)
library(dplyr)
library(biomaRt)
library(stringr)


# Connect to Ensembl (Mouse)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "useast"
)


# Files to plot

files <- c("deseq_results.csv")

make_volcano <- function(file) {
  df <- read.csv(file, row.names = 1)
  df <- df %>% filter(!is.na(padj))
  
  df$ensembl_gene_id <- gsub("\\..*", "", rownames(df))
  
  annotation <- getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol", "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = unique(df$ensembl_gene_id),
    mart       = mart
  )
  
  annotation <- annotation[!duplicated(annotation$ensembl_gene_id), ]
  
  df <- left_join(df, annotation, by = "ensembl_gene_id")
  
  df$mgi_symbol <- ifelse(
    !is.na(df$mgi_symbol) & df$mgi_symbol != "",
    df$mgi_symbol,
    ifelse(
      !is.na(df$external_gene_name) & df$external_gene_name != "",
      df$external_gene_name,
      df$ensembl_gene_id
    )
  )
  
  # Remove genes with no annotation (still showing as Ensembl IDs)
  DEG_annotOrdered <- df %>% 
    arrange(padj) %>%
    filter(!grepl("^ENSMUSG", mgi_symbol))
  
  DEG_annotOrdered <- df %>% arrange(padj)
  top_genes        <- head(DEG_annotOrdered$mgi_symbol, 15)
  
  newStudyTitle <- str_to_title(gsub("\\.csv", "", file))
  output_name   <- gsub("\\.csv", "_volcano.pdf", file)
  
  y_max <- ceiling(max(-log10(DEG_annotOrdered$padj), na.rm = TRUE)) + 1

  p <- EnhancedVolcano(
    DEG_annotOrdered,
    lab = DEG_annotOrdered$mgi_symbol,
    selectLab = top_genes,
    parseLabels = FALSE,
    col = c("#8491B4", "#00A087", "#6DBBE1", "red"),
    x = "log2FoldChange",
    y = "padj",
    FCcutoff = 2.0,
    pCutoff = 0.05,
    pointSize = 1,
    labSize = 2.0,
    legendIconSize = 4,
    colAlpha = 0.4,
    borderWidth = 0.2,
    xlim = c(-6, 6),
    ylim = c(0, y_max),
    title = newStudyTitle,
    subtitle = NULL,
    axisLabSize = 10,
    titleLabSize = 12,
    legendLabSize = 8,
    captionLabSize = 8,
    drawConnectors = TRUE,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    xlab = bquote(~Log[2] ~ "fold change"),
    ylab = bquote(Padj)
  )
  
  ggsave(filename = output_name, plot = p, width = 17, height = 20, units = "cm")
  cat("Saved:", output_name, "\n")
}

# Run over all files
for (file in files) {
  tryCatch(
    make_volcano(file),
    error = function(e) cat("Error processing", file, ":", conditionMessage(e), "\n")
  )
}
```


Heatmap (combined)

```{r}
library(DESeq2)
library(pheatmap)

# Combine top genes 
combined_genes <- unique(c(top15_old$ensembl_gene_id, top15_young$ensembl_gene_id))

# Subset VST expression matrix 
expr_mat <- assay(vsd)[combined_genes, , drop=FALSE]

gene_labels <- c(top15_old$gene_label, top15_young$gene_label)
names(gene_labels) <- combined_genes
rownames(expr_mat) <- gene_labels[rownames(expr_mat)]

# Sample annotation
sample_annot <- as.data.frame(colData(vsd)[, "age", drop=FALSE])
colnames(sample_annot) <- "Age"
rownames(sample_annot) <- colnames(expr_mat)

#  Rename columns to old_ovary / young_ovary 
colnames(expr_mat) <- gsub("_old$", "_old_ovary", colnames(expr_mat))
colnames(expr_mat) <- gsub("_young$", "_young_ovary", colnames(expr_mat))
rownames(sample_annot) <- colnames(expr_mat)

# Plot heatmap 
pheatmap(expr_mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = sample_annot,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         angle_col = 45,
         main = "Top 15 Genes: Old vs Young")

```


Heatmap (Old and Young Separate with top 15 genes showing expression)

```{r}
```{r}
library(DESeq2)
library(pheatmap)

# --- Step 1: Top 15 old genes only ---
old_genes <- top15_old$ensembl_gene_id

# --- Step 2: Subset VST expression matrix ---
expr_mat <- assay(vsd)[old_genes, , drop=FALSE]

gene_labels <- top15_old$gene_label
names(gene_labels) <- old_genes
rownames(expr_mat) <- gene_labels[rownames(expr_mat)]

# --- Step 3: Sample annotation ---
sample_annot <- as.data.frame(colData(vsd)[, "age", drop=FALSE])
colnames(sample_annot) <- "Age"
rownames(sample_annot) <- colnames(expr_mat)

# --- Step 4: Rename columns ---
colnames(expr_mat) <- gsub("_old$", "_old_ovary", colnames(expr_mat))
colnames(expr_mat) <- gsub("_young$", "_young_ovary", colnames(expr_mat))
rownames(sample_annot) <- colnames(expr_mat)

# --- Step 5: Colour palette ---
cols <- colorRampPalette(c("#F7F7F7", "#D1E5F0", "#92C5DE", "#6DBBE1",
                            "#8491B4", "#8485B4", "#9D6A40", "#987884",
                            "#D887D5", "#522371"))(100)

# --- Step 6: Save and plot ---
png("hmTopGenes_Old.jpeg", width = 70, height = 25, units = "cm", res = 300)

pheatmap::pheatmap(expr_mat,
                   scale          = "row",
                   cluster_rows   = TRUE,
                   cluster_cols   = TRUE,
                   annotation_col = sample_annot,
                   show_rownames  = TRUE,
                   show_colnames  = TRUE,
                   fontsize_row   = 8,
                   fontsize_col   = 7,
                   color          = cols,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   border_color   = "#8491B4",
                   cellwidth      = 10,
                   cellheight     = 12,
                   angle_col      = 45,
                   main           = "Top 15 Genes: Old Ovary")

graphics.off()
```


```{r}
library(DESeq2)
library(pheatmap)

# --- Step 1: Top 15 young genes only ---
young_genes <- top15_young$ensembl_gene_id

# --- Step 2: Subset VST expression matrix ---
expr_mat <- assay(vsd)[young_genes, , drop=FALSE]

gene_labels <- top15_young$gene_label
names(gene_labels) <- young_genes
rownames(expr_mat) <- gene_labels[rownames(expr_mat)]

# --- Step 3: Sample annotation ---
sample_annot <- as.data.frame(colData(vsd)[, "age", drop=FALSE])
colnames(sample_annot) <- "Age"
rownames(sample_annot) <- colnames(expr_mat)

# --- Step 4: Rename columns ---
colnames(expr_mat) <- gsub("_old$", "_old_ovary", colnames(expr_mat))
colnames(expr_mat) <- gsub("_young$", "_young_ovary", colnames(expr_mat))
rownames(sample_annot) <- colnames(expr_mat)

# --- Step 5: Colour palette ---
cols <- colorRampPalette(c("#F7F7F7", "#D1E5F0", "#92C5DE", "#6DBBE1",
                            "#8491B4", "#8485B4", "#9D6A40", "#987884",
                            "#D887D5", "#522371"))(100)

# --- Step 6: Save and plot ---
png("hmTopGenes_Young.jpeg", width = 70, height = 25, units = "cm", res = 300)

pheatmap::pheatmap(expr_mat,
                   scale          = "row",
                   cluster_rows   = TRUE,
                   cluster_cols   = TRUE,
                   annotation_col = sample_annot,
                   show_rownames  = TRUE,
                   show_colnames  = TRUE,
                   fontsize_row   = 8,
                   fontsize_col   = 7,
                   color          = cols,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   border_color   = "#8491B4",
                   cellwidth      = 10,
                   cellheight     = 12,
                   angle_col      = 45,
                   main           = "Top 15 Genes: Young Ovary")

graphics.off()
```

