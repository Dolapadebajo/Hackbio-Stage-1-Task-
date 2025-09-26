#!/usr/bin/env bash
set -euo pipefail

# Project root = parent of this script
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

PARENT_DIR="$PROJECT_ROOT/data/raw_data"
DEST_DIR="$PROJECT_ROOT/data/raw_data_50"
MANIFEST="$PROJECT_ROOT/meta/samples.txt"

mkdir -p "$DEST_DIR" "$PROJECT_ROOT/meta"

# Build list of basenames that have BOTH mates present
comm -12 \
  <(ls "$PARENT_DIR"/*_1.fastq.gz 2>/dev/null | sed 's#.*/##; s/_1\.fastq\.gz$//' | sort) \
  <(ls "$PARENT_DIR"/*_2.fastq.gz 2>/dev/null | sed 's#.*/##; s/_2\.fastq\.gz$//' | sort) \
  > "$PROJECT_ROOT/meta/paired_all.txt"

PAIRED_COUNT=$(wc -l < "$PROJECT_ROOT/meta/paired_all.txt")
if [[ "$PAIRED_COUNT" -lt 50 ]]; then
  echo "Only $PAIRED_COUNT complete pairs found in $PARENT_DIR (need 50)."
  echo "Add/finish downloads, then re-run."
  exit 1
fi

# Pick 50 at random
shuf "$PROJECT_ROOT/meta/paired_all.txt" | head -50 > "$MANIFEST"

# Copy the chosen pairs (use mv if you truly want to move them)
while read -r S; do
  cp -p "$PARENT_DIR/${S}_1.fastq.gz" "$DEST_DIR/"
  cp -p "$PARENT_DIR/${S}_2.fastq.gz" "$DEST_DIR/"
done < "$MANIFEST"

echo "Selected $(wc -l < "$MANIFEST") samples."
echo "Pairs copied to: $DEST_DIR"
echo "Manifest written to: $MANIFEST"
