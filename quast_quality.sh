#!/usr/bin/env bash
# Run QUAST on all assemblies (path-proof, resume-friendly).
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ASM_DIR="$PROJECT_ROOT/asm"                           # where per-sample SPAdes outputs live
QC_DIR="$PROJECT_ROOT/results/asm_quast"              # QUAST output
THREADS="${THREADS:-4}"

# Use min500 contigs if available; fall back to contigs.fasta
USE_MIN500="${USE_MIN500:-1}"   # set to 0 to force full contigs

command -v quast >/dev/null 2>&1 || { echo "ERROR: quast not found in PATH"; exit 1; }

mkdir -p "$QC_DIR"

# Collect assembly FASTAs
assemblies=()
while IFS= read -r -d '' d; do
  if [[ "$USE_MIN500" -eq 1 && -s "$d/contigs.min500.fasta" ]]; then
    assemblies+=( "$d/contigs.min500.fasta" )
  elif [[ -s "$d/contigs.fasta" ]]; then
    assemblies+=( "$d/contigs.fasta" )
  fi
done < <(find "$ASM_DIR" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

if [[ ${#assemblies[@]} -eq 0 ]]; then
  echo "ERROR: No assemblies found in $ASM_DIR"
  exit 1
fi

echo "== QUAST on ${#assemblies[@]} assemblies =="
echo " Output: $QC_DIR"
quast -t "$THREADS" -o "$QC_DIR" "${assemblies[@]}"

echo "QUAST finished. Open: $QC_DIR/report.html"
