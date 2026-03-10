---
title: "tximport_salmon"
output: html_document
date: "2026-02-05"
---


```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("tximport")

install.packages("rlang")
```

```{r}
library(rtracklayer)

# Import GTF
gtf <- rtracklayer::import("Mus_musculus.GRCm39.109.gtf")

# Keep only rows that have a transcript_id
tx_gtf <- gtf[!is.na(mcols(gtf)$transcript_id)]

# Build tx2gene mapping
tx2gene <- data.frame(
  TXNAME = gsub("\\.\\d+$", "", mcols(tx_gtf)$transcript_id),  # strip version
  GENEID = mcols(tx_gtf)$gene_id
)

# Optional: remove duplicates (each transcript should map to one gene)
tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), ]

# Save to CSV
write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)

# Check first rows
head(tx2gene)

list.files("salmon_quant_results", pattern="quant.sf", full.names=TRUE)


```

```{r}
library(tximport)

# Make sure the path is correct
files <- list.files("salmon_quant_results", pattern="quant.sf", full.names=TRUE)

# Read your tx2gene mapping
tx2gene <- read.csv("tx2gene.csv")

# Strip transcript versions in tx2gene to match Salmon IDs
tx2gene$TXNAME <- gsub("\\.\\d+$", "", tx2gene$TXNAME)

# Import with tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = TRUE)

# Check
names(txi)

head(txi$countsFromAbundance)

```

```{r}
# Export gene counts
write.csv(txi$counts, file = "gene_counts.csv", row.names = TRUE)

# Export TPM values
write.csv(txi$abundance, file = "gene_abundance_tpm.csv", row.names = TRUE)

# Export effective lengths
write.csv(txi$length, file = "gene_effective_length.csv", row.names = TRUE)

```

