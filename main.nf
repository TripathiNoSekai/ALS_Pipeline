#!/usr/bin/env nextflow

// Define channels
params.dataDir = "/data"
params.resultsDir = "/results"

process download_data {
    publishDir "$params.resultsDir/raw", mode: 'copy'

    """
    bash scripts/download.sh $params.dataDir
    """
}

process preprocessing {
    input:
    path rawData from download_data.out
    
    publishDir "$params.resultsDir/processed", mode: 'copy'
    
    """
    bash scripts/preprocessing.sh $rawData $params.resultsDir
    """
}

process de_analysis {
    input:
    path processedData from preprocessing.out

    publishDir "$params.resultsDir/DEG", mode: 'copy'

    """
    Rscript scripts/de_analysis.R $processedData $params.resultsDir
    """
}

process pathway_analysis {
    input:
    path deResults from de_analysis.out

    publishDir "$params.resultsDir/pathway", mode: 'copy'

    """
    Rscript scripts/pathway_analysis.R $deResults $params.resultsDir
    """
}

process visualization {
    input:
    path pathwayResults from pathway_analysis.out

    publishDir "$params.resultsDir/visuals", mode: 'copy'

    """
    Rscript scripts/visualization.R $pathwayResults $params.resultsDir
    """
}

// Workflow execution
workflow {
    download_data()
    preprocessing()
    de_analysis()
    pathway_analysis()
    visualization()
}
