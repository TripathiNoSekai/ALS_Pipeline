# Base Image
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install core dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    build-essential \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    software-properties-common \
    apt-utils \
    libharfbuzz-dev \
    libfribidi-dev \
    r-base \
    r-base-dev \
    python3-pip \
    python3-dev \
    python3-venv \
    && apt-get clean

# Install SRA Toolkit
RUN apt-get install -y sra-toolkit

# Install Bioinformatics Tools
RUN apt-get install -y \
    fastqc \
    multiqc \
    trimmomatic \
    star

# Install R and Bioconductor Packages
COPY R-packages.R /R-packages.R
RUN Rscript /R-packages.R
RUN R -e "install.packages('jsonlite', repos='http://cran.rstudio.com/')"


# Install Python Libraries
COPY requirements.txt /requirements.txt
RUN pip3 install --upgrade pip
RUN pip3 install -r /requirements.txt

# Install Nextflow
RUN wget -qO- https://get.nextflow.io | bash && mv nextflow /usr/local/bin/

# Set working directory
WORKDIR /app

# Copy scripts
COPY ./scripts /app/scripts
COPY ./main.nf /app/main.nf
COPY ./nextflow.config /app/nextflow.config

# Expose volumes for data and results
VOLUME ["/data", "/results"]

ENTRYPOINT ["nextflow", "run", "/app/main.nf"]
