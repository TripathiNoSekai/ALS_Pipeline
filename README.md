# 🧬 ALS Pipeline

A Nextflow pipeline for the automated analysis of ALS (Amyotrophic Lateral Sclerosis) sequencing data. This pipeline downloads SRA datasets, performs quality control and preprocessing, and aggregates data into summary reports.

---

## Table of Contents
1. [Introduction](#introduction)
2. [Features](#features)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Workflow](#workflow)
7. [Output](#output)
8. [Configuration](#configuration)
9. [Troubleshooting](#troubleshooting)
10. [Contributing](#contributing)
11. [License](#license)

---

## Introduction

The **ALS Nextflow Pipeline** is designed for processing and analyzing sequencing data related to ALS. It automates tasks such as data downloading, quality control, preprocessing, and report generation, providing a reproducible and modular framework for bioinformatics analyses.

---

## Features

- **Automated SRA Data Download:** Uses [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/Downloads) to fetch raw sequencing data.
- **Quality Control:** Performs QC using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- **Data Aggregation:** Aggregates QC reports into a single report using [MultiQC](https://github.com/ewels/MultiQC).
- **Reproducibility:** Runs on [Nextflow](https://www.nextflow.io/) for parallelized and scalable workflows.
- **Modular Design:** Easily customizable and extendable pipeline.

---

## Requirements

Ensure you have the following tools installed:

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
   - **SRA Toolkit:** Follow the instructions on [SRA Toolkit Downloads](https://github.com/ncbi/sra-tools/wiki/Downloads).  
   - **FastQC:** Download from [FastQC Download Page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and follow installation instructions.  
   - **MultiQC:** Install via pip:
     ```bash
     pip install multiqc
     ```

4. **Make the scripts executable:**
   ```bash
   chmod +x scripts/*.sh
   ```

---

## Usage

1. **Run the pipeline:**
   ```bash
   nextflow run main.nf -profile standard
   ```

2. **Optional Parameters:**  
   You can specify custom input and output paths:
   ```bash
   nextflow run main.nf --input samplesheet.csv --outdir ./results
   ```

---

## Workflow

The pipeline consists of the following steps:

1. **Download SRA Data:**  
   - Uses SRA Toolkit to download raw sequencing data in FASTQ format.
   - Script: `scripts/download.sh`

2. **Quality Control (QC):**  
   - Executes FastQC on the raw reads and generates individual QC reports.
   - Script: `scripts/preprocessing.sh`

3. **Data Aggregation:**  
   - Aggregates all FastQC reports into a single HTML summary using MultiQC.
   - Script: `scripts/multiqc.sh` *(if provided; otherwise, refer to MultiQC documentation)*

---

## Output

After executing the pipeline, you will find:

- **Raw Data:**  
  - Located in `results/raw_data/` (FASTQ files).
- **QC Reports:**  
  - Individual FastQC reports in `results/fastqc/`.
- **Aggregated Report:**  
  - A consolidated MultiQC report at `results/multiqc_report.html`.
- **Logs:**  
  - Execution logs stored in `results/logs/`.

---

## Configuration

You can customize pipeline parameters by editing the `nextflow.config` file. For example:

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
Ensure Nextflow is installed and in your PATH:
```bash
which nextflow
```

**Permission errors:**  
Ensure scripts have execution permissions:
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
