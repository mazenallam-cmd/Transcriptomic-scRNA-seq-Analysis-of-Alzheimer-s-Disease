# Transcriptomic scRNA-seq Analysis of Alzheimer's Disease

## Overview

This repository contains a bioinformatics pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data to investigate transcriptomic alterations associated with Alzheimer's disease. The analysis focuses on identifying differential gene expression patterns between Alzheimer's disease (AD) samples and control samples across different age intervals.

The project aims to explore how gene expression changes within individual brain cells during neurodegeneration and aging. By analyzing gene expression at single-cell resolution, this workflow enables identification of cell populations, disease-associated genes, and biological pathways involved in Alzheimer's pathology.

## Objectives

The main goals of this analysis are to:

* Process and analyze single-cell RNA sequencing data
* Identify distinct cellular populations in brain tissue
* Detect differentially expressed genes between Alzheimer's disease and control samples
* Compare gene expression patterns across different age groups
* Visualize transcriptomic differences using heatmaps and clustering analysis

## Dataset

The dataset used in this project consists of single-cell RNA sequencing data derived from brain tissue samples of individuals diagnosed with Alzheimer's disease and healthy controls. The samples are stratified into different age intervals to examine age-associated transcriptomic changes.

The expression matrix contains gene expression counts for thousands of genes measured across thousands of individual cells.

## Analysis Workflow

The analysis pipeline was implemented using the R toolkit **Seurat** for single-cell transcriptomic analysis.

The workflow includes the following steps:

### 1. Data Loading

Raw gene expression matrices are imported into R and converted into a Seurat object for downstream analysis.

### 2. Quality Control

Low-quality cells are removed based on:

* Number of detected genes per cell
* Total UMI counts
* Percentage of mitochondrial gene expression

This step ensures that only high-quality cells are retained for analysis.

### 3. Data Normalization

Gene expression values are normalized to correct for differences in sequencing depth across cells. Log-normalization is applied to stabilize variance across the dataset.

### 4. Feature Selection

Highly variable genes are identified to capture genes with the most biological variation across cells.

### 5. Dimensionality Reduction

Principal Component Analysis (PCA) is performed to reduce dimensionality, followed by non-linear visualization techniques such as UMAP to visualize cellular heterogeneity.

### 6. Cell Clustering

Cells are grouped into clusters based on transcriptional similarity using graph-based clustering algorithms.

### 7. Differential Gene Expression Analysis

Differential expression analysis is conducted to identify genes that are significantly upregulated or downregulated between Alzheimer's disease samples and control samples within specific age intervals.

### 8. Data Visualization

Heatmaps are generated to visualize the expression patterns of differentially expressed genes across conditions and age groups.

## Repository Structure

```
Transcriptomic_ScRNA_analysis
│
├── DE_per_age_interval_AD_vs_CTRL
│   ├── Heatmaps_DE_genes_Age_60_70.pdf
│
├── scripts
│   └── analysis_pipeline.R
│
├── results
│   └── differential_expression_results
│
└── README.md
```

## Key Outputs

The main outputs of this project include:

* Differentially expressed gene lists
* Heatmap visualizations of gene expression
* Clustered cell populations
* Comparative analysis of gene expression between AD and control samples across age groups

## Biological Significance

Understanding transcriptomic differences between Alzheimer's disease and healthy brain tissue can help identify genes and pathways involved in neurodegeneration. Single-cell transcriptomic analysis enables researchers to study disease mechanisms at cellular resolution and identify potential therapeutic targets.

## Tools and Technologies

The analysis pipeline uses the following tools:

* **Seurat**
* R programming language
* Single-cell transcriptomics analysis workflows
* Heatmap visualization methods

## Future Work

Future improvements to this project may include:

* Cell type annotation of identified clusters
* Pathway enrichment analysis using biological pathway databases
* Protein-protein interaction network analysis
* Cell-cell communication analysis
* Trajectory analysis to study disease progression

## Author

Mazen Allam
Biotechnology Student – Cairo University
Junior Bioinformatician interested in transcriptomics, neuroscience, and data-driven biological research.


