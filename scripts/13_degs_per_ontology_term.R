---
title: "DEGs per Ontological Term"
output: html_document
date: "2026-03-26"
---

# ============================================================
# Script 13: DEGs per Ontological Term
# Count differentially expressed genes (DEGs) annotated to
# each Gene Ontology (GO) term for mouse ovary RNA-seq data.
# Input:  results/tables/deseq2_full_results.csv
# Output: results/tables/degs_per_ontology_term.csv
# ============================================================

```{r setup}
suppressPackageStartupMessages({
  library(tidyverse)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

# ── Parameters ─────────────────────────────────────────────
PADJ_THRESHOLD <- 0.05
ONTOLOGY       <- "BP"   # "BP" | "MF" | "CC" | "ALL"

# ── Paths ──────────────────────────────────────────────────
INPUT_FILE  <- "results/tables/deseq2_full_results.csv"
OUTPUT_FILE <- "results/tables/degs_per_ontology_term.csv"

dir.create(dirname(OUTPUT_FILE), recursive = TRUE, showWarnings = FALSE)
```

# ============================================================
# STEP 1 – Load DESeq2 results
# ============================================================

```{r load-results}
message("Reading DESeq2 results from: ", INPUT_FILE)

deseq2_res <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)

# The gene column contains Ensembl IDs (ENSMUSG format)
if (!"gene" %in% colnames(deseq2_res)) {
  stop("Expected a 'gene' column with Ensembl IDs in ", INPUT_FILE)
}

message("Total genes in results: ", nrow(deseq2_res))
```

# ============================================================
# STEP 2 – Filter significant DEGs
# ============================================================

```{r filter-degs}
degs <- deseq2_res %>%
  filter(!is.na(padj), padj < PADJ_THRESHOLD)

message("Significant DEGs (padj < ", PADJ_THRESHOLD, "): ", nrow(degs))

ensembl_all  <- deseq2_res$gene
ensembl_degs <- degs$gene
```

# ============================================================
# STEP 3 – Map Ensembl IDs to Entrez IDs (required by org.Mm.eg.db)
# ============================================================

```{r map-entrez}
message("Mapping Ensembl IDs to Entrez IDs via org.Mm.eg.db ...")

map_to_entrez <- function(ensembl_ids) {
  AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys      = ensembl_ids,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
}

entrez_all  <- map_to_entrez(ensembl_all)
entrez_degs <- map_to_entrez(ensembl_degs)

# Remove unmapped (NA) IDs
entrez_all  <- na.omit(entrez_all)
entrez_degs <- na.omit(entrez_degs)

message("Mapped background genes: ", length(entrez_all),
        " / ", length(ensembl_all))
message("Mapped DEGs:             ", length(entrez_degs),
        " / ", length(ensembl_degs))

if (length(entrez_degs) == 0) {
  stop("No DEGs could be mapped to Entrez IDs. ",
       "Check that 'gene' column contains valid ENSMUSG identifiers.")
}
```

# ============================================================
# STEP 4 – Retrieve GO annotations and count DEGs per term
# ============================================================

```{r count-degs-per-term}
message("Retrieving GO annotations (ontology = '", ONTOLOGY, "') ...")

# Helper: build a long-format data frame mapping Entrez IDs to GO terms,
# including sub-ontology label (BP / MF / CC).
build_go_map <- function(entrez_ids) {
  # AnnotationDbi::select returns one row per (gene × GO term) pair
  go_raw <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys    = as.character(entrez_ids),
    columns = c("GO", "ONTOLOGY"),
    keytype = "ENTREZID"
  )

  go_raw %>%
    as_tibble() %>%
    rename(entrez_id = ENTREZID, go_id = GO, ontology = ONTOLOGY) %>%
    filter(!is.na(go_id), !is.na(ontology))
}

go_map_all  <- build_go_map(entrez_all)
go_map_degs <- build_go_map(entrez_degs)

# Filter by requested sub-ontology
if (ONTOLOGY != "ALL") {
  go_map_all  <- go_map_all  %>% filter(ontology == ONTOLOGY)
  go_map_degs <- go_map_degs %>% filter(ontology == ONTOLOGY)
}

message("GO annotation pairs (background): ", nrow(go_map_all))
message("GO annotation pairs (DEGs):        ", nrow(go_map_degs))
```

# ============================================================
# STEP 5 – Build summary table
# ============================================================

```{r build-summary}
message("Building summary table ...")

# Count total unique genes per GO term (background)
total_per_term <- go_map_all %>%
  group_by(go_id, ontology) %>%
  summarise(total_genes = n_distinct(entrez_id), .groups = "drop")

# Count DEGs per GO term
degs_per_term <- go_map_degs %>%
  group_by(go_id, ontology) %>%
  summarise(n_degs = n_distinct(entrez_id), .groups = "drop")

# Retrieve GO term names from org.Mm.eg.db
go_term_names <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = unique(go_map_all$go_id),
  columns = "TERM",
  keytype = "GO"
) %>%
  as_tibble() %>%
  rename(go_id = GO, term_name = TERM) %>%
  distinct(go_id, term_name)

# Join everything together
summary_table <- total_per_term %>%
  left_join(degs_per_term, by = c("go_id", "ontology")) %>%
  mutate(n_degs = replace_na(n_degs, 0L)) %>%
  left_join(go_term_names, by = "go_id") %>%
  mutate(
    pct_degs = round(n_degs / total_genes * 100, 2)
  ) %>%
  filter(n_degs > 0) %>%
  arrange(desc(n_degs)) %>%
  select(go_id, term_name, ontology, n_degs, total_genes, pct_degs)

message("GO terms with at least one DEG: ", nrow(summary_table))
```

# ============================================================
# STEP 6 – Save results
# ============================================================

```{r save-results}
write.csv(summary_table, OUTPUT_FILE, row.names = FALSE)
message("Results saved to: ", OUTPUT_FILE)

# Print top 20 terms for quick inspection
message("\nTop 20 GO terms by number of DEGs:")
print(head(summary_table, 20))
```
