---
title: "proteomics_clean"
output: html_document
date: "2026-03-16"
---
Soluble/cellular Proteomics Clean

```{r}
library(tidyverse)

prot_path_cellular  <- "proteomics/Singh-Cellular-20713-4-V2/tables/Table-19.txt"
prot_path_soluble <- "proteomics/Singh-Soluble-20747-7/tables/Table-17.txt"

prot_raw_soluble <- read_tsv(prot_path_soluble, show_col_types = FALSE)
prot_raw_cellular <- read_tsv(prot_path_cellular, show_col_types = FALSE)

meta_cols <- c(
  "PG.ProteinGroups", "PG.Genes",
  "PG.ProteinDescriptions", "PG.UniProtIds", "PG.FastaFiles",
  "PG.Database", "PG.OrganismId", "PG.MolecularWeight",
  "PG.NrOfPrecursorsIdentified (Experiment-wide)",
  "GOMF name", "GOBP slim name", "KEGG name", "GOCC slim name"
)

# Define cleaning function to avoid repetition
clean_proteomics <- function(prot_raw, output_name) {
  
  # Per-sample log2 abundance columns
  log2_cols <- grep("\\.Quantity\\.log2$", names(prot_raw), value = TRUE)
  
  # Statistical result columns
  stat_cols <- grep(
    "^(p\\.value|adj\\.p\\.value|log2\\.fold\\.change|is\\.significant|t\\.statistic)\\.",
    names(prot_raw), value = TRUE
  )
  
  # Per-sample precursor count columns
  precursor_quant_cols <- grep("\\.NrOfPrecursorsUsedForQuantification$",
                               names(prot_raw), value = TRUE)
  
  prot_clean <- prot_raw %>%
    select(all_of(c(meta_cols, log2_cols, stat_cols))) %>%
    
    mutate(
      gene_symbol         = str_split(PG.Genes, ";", simplify = TRUE)[, 1],
      gene_symbol         = na_if(str_trim(gene_symbol), ""),
      uniprot_id          = str_split(PG.UniProtIds, ";", simplify = TRUE)[, 1],
      protein_description = str_split(PG.ProteinDescriptions, ";", simplify = TRUE)[, 1],
      mol_weight_da       = as.numeric(
                              str_split(PG.MolecularWeight, ";", simplify = TRUE)[, 1]
                            )
    ) %>%
    
    # Convert NaN -> NA
    mutate(across(all_of(log2_cols), ~ na_if(., NaN))) %>%
    
    # Per-row QC metrics
    mutate(
      n_samples_detected   = rowSums(!is.na(across(all_of(log2_cols)))),
      pct_samples_detected = round(100 * n_samples_detected / length(log2_cols), 1),
      
      median_precursors_used = if (length(precursor_quant_cols) > 0)
        apply(prot_raw[, precursor_quant_cols], 1,
              function(x) median(x[x > 0], na.rm = TRUE))
      else NA_real_
    ) %>%
    
    # Filter for reliability
    filter(is.na(median_precursors_used) | median_precursors_used >= 2) %>%
    
    # Per-timepoint mean log2 abundance
    {
      tp <- c("D14", "D21", "M03", "M06", "M12", "M18")
      df <- .
      for (t in tp) {
        cols_t <- grep(paste0("^", t, "\\..*\\.Quantity\\.log2$"),
                       names(df), value = TRUE)
        if (length(cols_t) > 0) {
          df[[paste0("mean_log2_", t)]] <- rowMeans(df[, cols_t], na.rm = TRUE)
        }
      }
      df
    } %>%
    
    rename(
      protein_groups        = PG.ProteinGroups,
      n_precursors_exp      = `PG.NrOfPrecursorsIdentified (Experiment-wide)`,
      go_molecular_function = `GOMF name`,
      go_biological_process = `GOBP slim name`,
      go_cellular_component = `GOCC slim name`,
      kegg_pathway          = `KEGG name`
    ) %>%
    
    select(
      # Identifiers
      gene_symbol, uniprot_id, protein_groups, protein_description,
      mol_weight_da, n_precursors_exp,
      # QC
      n_samples_detected, pct_samples_detected,
      median_precursors_used,
      # Per-timepoint summaries
      starts_with("mean_log2_"),
      # Per-sample log2 values
      all_of(log2_cols),
      # Statistics
      all_of(stat_cols),
      # Annotations
      go_molecular_function, go_biological_process,
      go_cellular_component, kegg_pathway
    )
  
  write_csv(prot_clean, output_name, na = "")
  message("✔ Wrote ", output_name, " (", nrow(prot_clean), " rows)")
  
  return(prot_clean)
}

# Apply to both datasets
prot_cellular_clean <- clean_proteomics(prot_raw_cellular, "cellular_proteomics_clean_v2.csv")
prot_soluble_clean <- clean_proteomics(prot_raw_soluble, "soluble_proteomics_clean_v2.csv")
```

