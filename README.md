# 🧬 ALS Nextflow Pipeline

A comprehensive Nextflow pipeline for multiomic analysis of ALS (Amyotrophic Lateral Sclerosis) data. This pipeline covers everything from data acquisition to visualization and supports additional analyses such as differential expression, functional enrichment, multi-model comparisons, and optional multiomic integration.

---

## Table of Contents
1. [Introduction](#introduction)
2. [Features](#features)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Step-by-Step Workflow](#step-by-step-workflow)
    - [1. Data Acquisition](#1-data-acquisition)
    - [2. Quality Control (Optional)](#2-quality-control-optional)
    - [3. Read Alignment (Optional)](#3-read-alignment-optional)
    - [4. Count Matrix Generation (Optional)](#4-count-matrix-generation-optional)
    - [5. Differential Gene Expression Analysis](#5-differential-gene-expression-analysis)
    - [6. Functional Enrichment Analysis](#6-functional-enrichment-analysis)
    - [7. Multi-model and Sex-specific Comparison](#7-multi-model-and-sex-specific-comparison)
    - [8. Multiomic Integration (Optional)](#8-multiomic-integration-optional)
    - [9. Visualization and Reporting](#9-visualization-and-reporting)
6. [Configuration](#configuration)
7. [Troubleshooting](#troubleshooting)
8. [Contributing](#contributing)
9. [License](#license)

---

## Introduction

The **ALS Nextflow Pipeline** automates the analysis of sequencing data related to ALS. It is designed to be modular, reproducible, and scalable—allowing you to download raw data, perform quality control, execute differential expression analysis, carry out pathway enrichment, and integrate multiomic datasets.

---

## Features

- **Automated Data Acquisition:** Download pre-processed count matrices or raw SRA data using tools like wget, fastq-dump, and SRA Toolkit.
- **Quality Control:** Execute FastQC and MultiQC for assessing raw data quality.
- **Read Alignment:** (Optional) Align reads to a reference genome using STAR or HISAT2.
- **Count Matrix Generation:** (Optional) Generate read count matrices using HTSeq or featureCounts.
- **Differential Expression Analysis:** Identify differentially expressed genes with DESeq2/edgeR.
- **Functional Enrichment Analysis:** Perform GO and pathway enrichment using clusterProfiler and ReactomePA.
- **Multi-model and Sex-specific Analysis:** Compare signatures across different ALS models and between sexes.
- **Multiomic Integration:** (Optional) Integrate transcriptomic data with proteomics/metabolomics using WGCNA and OmicsIntegrator.
- **Visualization:** Generate heatmaps, volcano plots, and comprehensive reports.

---

## Requirements

- [Nextflow](https://www.nextflow.io/)
- [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/Downloads)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://github.com/ewels/MultiQC)

---

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/TripathiNoSekai/ALS_Pipeline.git
   cd ALS_Pipeline
   ```

2. **Install Nextflow (if not already installed):**
   ```bash
   curl -s https://get.nextflow.io | bash
   sudo mv nextflow /usr/local/bin/
   ```

3. **Install SRA Toolkit, FastQC, and MultiQC:**
   - Follow instructions on the respective websites.
   - For MultiQC, you can install via pip:
     ```bash
     pip install multiqc
     ```

4. **Make the scripts executable:**
   ```bash
   chmod +x scripts/*.sh
   ```

---

## Step-by-Step Workflow

### 1. Data Acquisition
**Tools:** wget, fastq-dump, SRA Toolkit

- **A) Download the Pre-processed Count Matrices:**  
  ```bash
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234245/suppl/GSE234245_C9orf72-counts_norm_vst.csv.gz  
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234245/suppl/GSE234245_FUS-counts_norm_vst.csv.gz  
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234245/suppl/GSE234245_SOD1-counts_norm_vst.csv.gz  
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234245/suppl/GSE234245_TDP43-counts_norm_vst.csv.gz  
  ```

- **B) Download Raw SRA Data (Optional):**  
  ```bash
  prefetch PRJNA980618  
  fastq-dump --split-files SRRXXXXXXX
  ```

---

### 2. Quality Control (Optional)
**Tools:** FastQC, MultiQC, Trimmomatic

- **Run FastQC:**
  ```bash
  fastqc *.fastq  
  multiqc .
  ```

- **Trim Low-quality Reads:**
  ```bash
  trimmomatic PE -phred33 input_R1.fastq input_R2.fastq \
    output_R1_paired.fastq output_R1_unpaired.fastq \
    output_R2_paired.fastq output_R2_unpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
  ```

---

### 3. Read Alignment (Optional)
**Tools:** STAR, HISAT2

- **Index the Genome:**
  ```bash
  STAR --runMode genomeGenerate --genomeDir ./genome --genomeFastaFiles mm10.fa --runThreadN 8
  ```

- **Align the Reads:**
  ```bash
  STAR --runThreadN 8 --genomeDir ./genome --readFilesIn input_R1.fastq input_R2.fastq --outFileNamePrefix sample1
  ```

---

### 4. Count Matrix Generation (Optional)
**Tools:** HTSeq, featureCounts

- **Count Reads:**
  ```bash
  htseq-count -f bam -r pos -s no -t exon -i gene_id sample1.bam genes.gtf > sample1_counts.txt
  ```

---

### 5. Differential Gene Expression Analysis
**Tools:** DESeq2, edgeR (R)

- **Example R Script:**
  ```r
  # Load Libraries
  library(DESeq2)

  # Import Data
  counts <- read.csv("GSE234245_FUS-counts_norm_vst.csv", row.names = 1)
  coldata <- data.frame(
    condition = factor(c("control", "ALS")),
    row.names = colnames(counts)
  )

  # Create DESeq2 Object and Analyze
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  ```

---

### 6. Functional Enrichment Analysis
**Tools:** clusterProfiler, ReactomePA, GOseq

- **GO Enrichment Example:**
  ```r
  library(clusterProfiler)
  
  ego <- enrichGO(gene = rownames(res)[res$padj < 0.05],
                  OrgDb = org.Mm.eg.db, 
                  keyType = "ENSEMBL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)
  
  dotplot(ego, showCategory = 10)
  ```

- **Pathway Analysis (Reactome):**
  ```r
  library(ReactomePA)
  
  rpa <- enrichPathway(gene = rownames(res)[res$padj < 0.05], 
                       organism = "mouse", 
                       pvalueCutoff = 0.05)
  
  dotplot(rpa, showCategory = 10)
  ```

---

### 7. Multi-model and Sex-specific Comparison

- **Combine DEGs from All Four Models:**  
  Identify shared and unique gene signatures using Venn diagrams.

- **Sex-stratified Analysis Example:**
  ```r
  dds_sex <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ sex + condition)
  dds_sex <- DESeq(dds_sex)
  ```

---

### 8. Multiomic Integration (Optional)
**Tools:** WGCNA, OmicsIntegrator

- **Correlation Analysis:**  
  Map RNA expression profiles to proteomics/metabolomics data.

- **WGCNA Network Example:**
  ```r
  library(WGCNA)
  net = blockwiseModules(datExpr, power = 6)
  ```

---

### 9. Visualization and Reporting
**Tools:** pheatmap, EnhancedVolcano (R)

- **Heatmap Example:**
  ```r
  library(pheatmap)
  pheatmap(assay(vst(dds)), cluster_rows = TRUE, cluster_cols = TRUE)
  ```

- **Volcano Plot Example:**
  ```r
  library(EnhancedVolcano)
  EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue')
  ```

---

## Configuration

You can customize the pipeline parameters by editing the `nextflow.config` file. For example:

```groovy
params {
    input   = 'samplesheet.csv'          // Input file containing SRA IDs
    outdir  = './results'                // Output directory
    threads = 8                          // Number of parallel threads
    mem     = '16GB'                     // Memory allocation
}
```

---

## Troubleshooting

**Nextflow not found:**  
Ensure that Nextflow is installed and added to your PATH:
```bash
which nextflow
```

**Permission errors:**  
Ensure that all scripts have execution permissions:
```bash
chmod +x scripts/*.sh
```

**Missing tools:**  
Verify that SRA Toolkit, FastQC, and MultiQC are installed:
```bash
which fastqc
which multiqc
```

---

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository:
   ```bash
   git fork https://github.com/TripathiNoSekai/ALS_Pipeline.git
   ```
2. Create a new branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```
3. Make your changes and commit:
   ```bash
   git add .
   git commit -m "Add new feature"
   ```
4. Push your changes:
   ```bash
   git push origin feature/your-feature-name
   ```
5. Open a pull request.

---

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.
