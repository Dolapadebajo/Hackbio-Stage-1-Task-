#!/usr/bin/env bash
# Run FastQC + MultiQC on raw sequencing data (resume-safe, path-proof).
set -euo pipefail

# ---------- Config / Auto-paths ----------
# Project root = parent dir of this script
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Allow optional override: bash raw_qc.sh data/raw_data_50
RAW_DIR_REL="${1:-data/raw_data}"                  # relative to project root
RAW_DATA_DIR="$PROJECT_ROOT/$RAW_DIR_REL"

# Output directory
QC_OUTPUT_DIR="$PROJECT_ROOT/results/raw_fastqc_reports"

# Threads for FastQC
THREADS="${THREADS:-4}"

# ---------- Checks ----------
# Check tools exist
command -v fastqc   >/dev/null 2>&1 || { echo "ERROR: fastqc not found in PATH"; exit 1; }
command -v multiqc  >/dev/null 2>&1 || { echo "ERROR: multiqc not found in PATH"; exit 1; }

# Create output dir
mkdir -p "$QC_OUTPUT_DIR"

# Ensure raw data dir exists
if [[ ! -d "$RAW_DATA_DIR" ]]; then
  echo "ERROR: Raw data directory not found: $RAW_DATA_DIR"
  echo "Tip: Your project root is: $PROJECT_ROOT"
  exit 1
fi

# Ensure there are .fastq.gz files
shopt -s nullglob
files=( "$RAW_DATA_DIR"/*.fastq.gz )
if [[ ${#files[@]} -eq 0 ]]; then
  echo "ERROR: No .fastq.gz files in: $RAW_DATA_DIR"
  exit 1
fi
shopt -u nullglob

echo "== FastQC start =="
echo " Input : $RAW_DATA_DIR"
echo " Output: $QC_OUTPUT_DIR"
echo " Threads: $THREADS"

# ---------- FastQC ----------
fastqc --threads "$THREADS" --outdir "$QC_OUTPUT_DIR" "${files[@]}"

echo "FastQC complete. Reports written to: $QC_OUTPUT_DIR"

# ---------- MultiQC ----------
echo "Running MultiQC..."
multiqc "$QC_OUTPUT_DIR" \
  --outdir "$QC_OUTPUT_DIR" \
  --filename "multiqc_report_raw.html" \
  --quiet

echo "MultiQC done: $QC_OUTPUT_DIR/multiqc_report_raw.html"
echo "QC finished successfully."

