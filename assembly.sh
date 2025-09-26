#!/usr/bin/env bash
# Assemble paired-end reads with SPAdes (resume-safe, path-proof).
set -euo pipefail

# -------- Auto paths & settings --------
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Input dir (override with: bash assembly.sh data/trimmed_alt )
TRIM_REL="${1:-data/trimmed}"
TRIM_DIR="$PROJECT_ROOT/$TRIM_REL"

ASM_DIR="$PROJECT_ROOT/asm"
LOG_DIR="$PROJECT_ROOT/logs"        # summary logs (per-sample logs live in asm/<sample>/spades.log)
THREADS="${THREADS:-4}"
MEM_GB="${MEM_GB:-8}"               # Set lower if RAM is tight; raise if you can
CARE="-k auto --careful --only-assembler"

# Tools check
command -v spades.py >/dev/null 2>&1 || { echo "ERROR: spades.py not found in PATH"; exit 1; }

mkdir -p "$ASM_DIR" "$LOG_DIR" "$PROJECT_ROOT/meta"

# -------- Validate input --------
if [[ ! -d "$TRIM_DIR" ]]; then
  echo "ERROR: Trimmed reads directory not found: $TRIM_DIR"
  echo "Tip: project root is: $PROJECT_ROOT"
  exit 1
fi

# Build list of basenames that have BOTH mates present
comm -12 \
  <(ls "$TRIM_DIR"/*_1.fastq.gz 2>/dev/null | sed 's#.*/##; s/_1\.fastq\.gz$//' | sort) \
  <(ls "$TRIM_DIR"/*_2.fastq.gz 2>/dev/null | sed 's#.*/##; s/_2\.fastq\.gz$//' | sort) \
  > "$PROJECT_ROOT/meta/paired_for_assembly.txt"

N=$(wc -l < "$PROJECT_ROOT/meta/paired_for_assembly.txt" || echo 0)
if [[ "$N" -eq 0 ]]; then
  echo "ERROR: No paired trimmed FASTQs found in $TRIM_DIR"
  exit 1
fi

echo "== SPAdes assembly =="
echo " Input : $TRIM_DIR"
echo " Output: $ASM_DIR"
echo " Pairs : $N"
echo " Threads: $THREADS  RAM: ${MEM_GB}G"
echo

# -------- Assembly loop (resume-safe) --------
while read -r S; do
  R1="$TRIM_DIR/${S}_1.fastq.gz"
  R2="$TRIM_DIR/${S}_2.fastq.gz"
  OUT="$ASM_DIR/$S"
  LOG="$OUT/spades.log"

  mkdir -p "$OUT"

  # Resume skip: if finished before, skip
  if [[ -s "$OUT/contigs.fasta" ]] && grep -qi "finished successfully" "$LOG" 2>/dev/null; then
    echo "[SKIP] $S (already finished)"
    continue
  fi

  # Clean partials from previous crash for this sample (safe)
  find "$OUT" -maxdepth 1 -type d -name "K??" -exec rm -rf {} + 2>/dev/null || true
  rm -rf "$OUT/tmp" 2>/dev/null || true

  echo "[SPAdes] $S"
  spades.py -1 "$R1" -2 "$R2" -o "$OUT" -t "$THREADS" -m "$MEM_GB" $CARE &>> "$LOG"

  # Keep a size-filtered contig set too (≥500 bp) using awk fallback
  if [[ -s "$OUT/contigs.fasta" ]]; then
    awk 'BEGIN{RS=">"; ORS=""}
         NR>1{
           header=$0; sub(/\n.*/,"",header);
           seq=$0; sub(/^[^\n]*\n/,"",seq); gsub(/\n/,"",seq);
           if (length(seq)>=500) print ">"header"\n"seq"\n";
         }' "$OUT/contigs.fasta" > "$OUT/contigs.min500.fasta" || true
  fi
done < "$PROJECT_ROOT/meta/paired_for_assembly.txt"

echo
echo "SPAdes step finished. Assemblies in: $ASM_DIR"
