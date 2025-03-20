#!/bin/bash
RAW_DIR=$1
RESULTS_DIR=$2

# Run QC
fastqc "$RAW_DIR"/*.csv -o "$RESULTS_DIR/qc"
multiqc "$RESULTS_DIR/qc" -o "$RESULTS_DIR/qc"
