#!/usr/bin/env bash
# Trim low-quality bases/adapters with fastp, then QC with FastQC + MultiQC.
# Resume-safe and path-proof for the Dolapo project layout.
set -euo pipefail

# -------- Auto paths & settings --------
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Input dir (override with: bash trim_and_qc.sh data/raw_data)
SRC_REL="${1:-data/raw_data_50}"         # default to the 50-sample subset
SRC_DIR="$PROJECT_ROOT/$SRC_REL"

TRIM_DIR="$PROJECT_ROOT/data/trimmed"
FASTP_DIR="$PROJECT_ROOT/results/fastp_reports"
QC_TRIM_DIR="$PROJECT_ROOT/results/trimmed_fastqc_reports"

THREADS="${THREADS:-4}"

# Tools check
for t in fastp fastqc multiqc; do
  command -v "$t" >/dev/null 2>&1 || { echo "ERROR: $t not found in PATH"; exit 1; }
done

# Create output dirs
mkdir -p "$TRIM_DIR" "$FASTP_DIR" "$QC_TRIM_DIR" "$PROJECT_ROOT/meta"

# -------- Validate input --------
if [[ ! -d "$SRC_DIR" ]]; then
  echo "ERROR: Source directory not found: $SRC_DIR"
  echo "Tip: project root is: $PROJECT_ROOT"
  exit 1
fi

# Build list of basenames that have BOTH mates present
comm -12 \
  <(ls "$SRC_DIR"/*_1.fastq.gz 2>/dev/null | sed 's#.*/##; s/_1\.fastq\.gz$//' | sort) \
  <(ls "$SRC_DIR"/*_2.fastq.gz 2>/dev/null | sed 's#.*/##; s/_2\.fastq\.gz$//' | sort) \
  > "$PROJECT_ROOT/meta/paired_to_trim.txt"

N=$(wc -l < "$PROJECT_ROOT/meta/paired_to_trim.txt" || echo 0)
if [[ "$N" -eq 0 ]]; then
  echo "ERROR: No paired FASTQs found in $SRC_DIR"
  exit 1
fi

echo "== Trimming with fastp =="
echo " Source : $SRC_DIR"
echo " Trimmed: $TRIM_DIR"
echo " Reports: $FASTP_DIR"
echo " Pairs   : $N"
echo

# -------- Trim loop (resume-safe) --------
while read -r S; do
  in1="$SRC_DIR/${S}_1.fastq.gz"
  in2="$SRC_DIR/${S}_2.fastq.gz"

  out1="$TRIM_DIR/${S}_1.fastq.gz"
  out2="$TRIM_DIR/${S}_2.fastq.gz"

  html="$FASTP_DIR/${S}.fastp.html"
  json="$FASTP_DIR/${S}.fastp.json"

  # Skip if already trimmed
  if [[ -s "$out1" && -s "$out2" ]]; then
    echo "[SKIP] $S (trimmed files already exist)"
    continue
  fi

  echo "[fastp] $S"
  fastp \
    -i "$in1" -I "$in2" \
    -o "$out1" -O "$out2" \
    -h "$html" -j "$json" \
    --thread "$THREADS" \
    --detect_adapter_for_pe \
    --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
    --qualified_quality_phred 20 --unqualified_percent_limit 40 \
    --length_required 50 \
    --compression 6 \
    --verbose

done < "$PROJECT_ROOT/meta/paired_to_trim.txt"

echo
echo "== FastQC on trimmed reads =="
# gather all trimmed R1 + R2
shopt -s nullglob
trimmed_files=( "$TRIM_DIR"/*_1.fastq.gz "$TRIM_DIR"/*_2.fastq.gz )
shopt -u nullglob

if [[ ${#trimmed_files[@]} -eq 0 ]]; then
  echo "ERROR: No trimmed files found in $TRIM_DIR"
  exit 1
fi

fastqc --threads "$THREADS" --outdir "$QC_TRIM_DIR" "${trimmed_files[@]}"
echo "FastQC complete → $QC_TRIM_DIR"

echo "== MultiQC (trimmed QC) =="
multiqc "$QC_TRIM_DIR" \
  --outdir "$QC_TRIM_DIR" \
  --filename "multiqc_report_trimmed.html" \
  --quiet

echo
echo "All done."
echo "Trimmed FASTQs : $TRIM_DIR"
echo "fastp reports  : $FASTP_DIR"
echo "MultiQC (trim) : $QC_TRIM_DIR/multiqc_report_trimmed.html"
