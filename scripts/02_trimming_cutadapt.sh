#!/bin/bash
#SBATCH --job-name=cutadapt_trim
#SBATCH --output=cutadapt_trim_%j.log
#SBATCH --error=cutadapt_trim_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=msc_appbio

module load py-cutadapt

RAW_DIR=/scratch/users/k25041394/research_project/fastq_raw
OUT_DIR=/scratch/users/k25041394/research_project/fastqc_trimmed

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

cd "$RAW_DIR" || exit 1

for R1 in *_R1_001.fastq.gz; do
    R2=${R1/_R1_/_R2_}

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: Missing R2 for $R1 — skipping"
        continue
    fi

    OUT_R1="$OUT_DIR/trimmed_$R1"
    OUT_R2="$OUT_DIR/trimmed_$R2"

    echo "Trimming:"
    echo "  Input R1: $R1"
    echo "  Input R2: $R2"
    echo "  Output R1: $OUT_R1"
    echo "  Output R2: $OUT_R2"

    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -q 20,20 \
        -m 30 \
        -o "$OUT_R1" \
        -p "$OUT_R2" \
        "$R1" "$R2"

    echo "Finished $R1"
    echo "----------------------------"
done

echo "All trimming completed safely."