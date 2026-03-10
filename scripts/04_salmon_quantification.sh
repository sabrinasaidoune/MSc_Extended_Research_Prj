#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --output=salmon_quant_%j.out
#SBATCH --error=salmon_quant_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

#------------------------
# Paths on HPC
#------------------------
FASTQ_DIR=/scratch/prj/wch_ovarian_ageing/research_project/fastqc_trimmed
INDEX_DIR=/scratch/prj/wch_ovarian_ageing/research_project/genome_mm10/Salmon_transcripts_index_v2
QUANT_DIR=/scratch/prj/wch_ovarian_ageing/research_project/salmon_quant_v2

# Create output directory if it doesn't exist
mkdir -p "$QUANT_DIR"

#------------------------
# Quantify all paired-end FASTQs
#------------------------
for R1 in "$FASTQ_DIR"/trimmed_*_R1_001.fastq.gz; do

    # Extract sample name
    SAMPLE=$(basename "$R1" _R1_001.fastq.gz)

    # Define R2
    R2="$FASTQ_DIR/${SAMPLE}_R2_001.fastq.gz"

    # Define output directory
    OUT_DIR="$QUANT_DIR/$SAMPLE"

    # Skip if R2 is missing
    if [ ! -f "$R2" ]; then
        echo "WARNING: Missing R2 for $SAMPLE — skipping"
        continue
    fi

    echo "======================================"
    echo "Running Salmon quant for sample: $SAMPLE"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "Output: $OUT_DIR"
    echo "======================================"

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        --validateMappings \
            -p 8 \
        -o "$OUT_DIR"

    echo "Finished: $SAMPLE"
    echo
done

echo "All Salmon quantifications completed."

