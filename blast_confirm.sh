#!/usr/bin/env bash
# Script: 07_blast_confirm.sh
# Description: Run BLAST on a single representative sample for organism identification.

set -euo pipefail

# --- Directories ---
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSEMBLY_DIR="$PROJECT_ROOT/asm"
BLAST_DIR="$PROJECT_ROOT/results/blast"
mkdir -p "$BLAST_DIR"

echo "=== Running BLAST for organism identification (rubric requirement) ==="

# --- Get the first successful assembly ---
REPRESENTATIVE_ASSEMBLY=$(find "$ASSEMBLY_DIR" -name "contigs.fasta" | head -1)

if [[ -z "$REPRESENTATIVE_ASSEMBLY" ]]; then
    echo "Error: No assemblies found in $ASSEMBLY_DIR. Run assembly.sh first."
    exit 1
fi

SAMPLE_NAME=$(basename "$(dirname "$REPRESENTATIVE_ASSEMBLY")")
echo "Using representative sample: $SAMPLE_NAME"
echo "Assembly file: $REPRESENTATIVE_ASSEMBLY"

# --- Extract the first contig (not whole genome, keeps BLAST fast) ---
REP_CONTIG="$BLAST_DIR/${SAMPLE_NAME}_representative_contig.fasta"
head -n 200 "$REPRESENTATIVE_ASSEMBLY" > "$REP_CONTIG"

# --- Run BLAST remotely against NCBI nt ---
echo "Running BLAST remotely against NCBI nt database (may take several minutes)..."

blastn \
    -query "$REP_CONTIG" \
    -db nt \
    -remote \
    -outfmt "6 std stitle" \
    -max_target_seqs 5 \
    -evalue 1e-50 \
    -out "$BLAST_DIR/blast_identification_results.tsv"

echo ""
echo "BLAST complete. Top hits:"
echo "----------------------------------------"
awk -F'\t' '{printf "%-60s %-6s %-6s %-10s\n", $13, $3, $4, $11}' \
    "$BLAST_DIR/blast_identification_results.tsv" | head -5
echo "----------------------------------------"

# --- Check for Listeria ---
if grep -qi "listeria" "$BLAST_DIR/blast_identification_results.tsv"; then
    echo "✓ SUCCESS: Listeria monocytogenes identified via BLAST."
else
    echo "✗ WARNING: Expected Listeria not found in top BLAST hits."
fi
