#!/bin/bash
DATA_DIR=$1

wget -P "$DATA_DIR" ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234245/suppl/GSE234245_FUS-counts_norm_vst.csv.gz
gunzip "$DATA_DIR"/*.gz
