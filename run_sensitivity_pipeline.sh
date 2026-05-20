#!/usr/bin/env bash
# Run the xGen UMI pipeline on the 5-sample sensitivity panel.
# Usage: ./run_sensitivity_pipeline.sh [threads] [jobs]
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh && conda activate ngs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTQ_DIR="$SCRIPT_DIR/fastq_sens"
OUT_DIR="$SCRIPT_DIR/out_sens"
THREADS="${1:-8}"
JOBS="${2:-3}"

mkdir -p "$OUT_DIR"

echo "Running xGen UMI sensitivity pipeline on pos1 pos2 pos3 pos4 neg5"
echo "  FASTQs: $FASTQ_DIR"
echo "  Output: $OUT_DIR"
echo "  Threads: $THREADS  Jobs: $JOBS"
echo ""

time parallel -j "$JOBS" \
    "$SCRIPT_DIR/run_sample.sh {} $THREADS $FASTQ_DIR $OUT_DIR > $OUT_DIR/{}.log 2>&1" \
    ::: pos1 pos2 pos3 pos4 neg5

echo ""
echo "Pipeline complete. Run: python3 benchmark_sensitivity.py"
