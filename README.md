## Using Omics Approaches to Understand Reproductive Transitions and Ovarian Ageing Across the Life Course

Sabrina Souhaya Saidoune

Supervisors: Dr Michelle Simon & Dr Kim Jonas

This repository contains scripts, data, and results for a multi-omics study of ovarian ageing.

### Graphical Abstract



### Structure

- `data/` - raw and processed data
- `scripts/` - analysis scripts
- `notebooks/` - R markdown exploratory analysis
- `results/` - QC, DESeq2 results, plots, enrichment analysis
- `proteomics/` - proteomics raw and processed data
- `figures/` - dissertation figures and graphical abstract
- `environment/` - R environment configuration

### Instructions

1. Set up R environment: 
2. Run scripts in order:
   - `01_fastqc.sh`
   - `02_trimming_cutadapt.sh`
   - `03_salmon_index.sh`
   - `04_salmon_quantification.sh`
   - `05_tximport.R`
   - `06_deseq2_analysis.R`
   - `07_visualisations.R`
   - `08_go_enrichment.R`

3. Check results in `results/` folder.
