# Differential Centrality-Ensemble Analysis Based on Gene Expression Profiles and Protein-Protein Interaction Network

This R package provides a comprehensive framework for identifying and prioritizing disease-related genes based on gene expression profiles and protein-protein interaction (PPI) networks.

## Function Overview

The `DiffCentEn_function` function takes a list of two elements:
- `data`: A data frame of gene expression data with gene symbols as column names and class labels (Tumor, Normal) in the last column.
- `topGenes`: A data frame of differentially expressed genes (DEGs) with columns `Gene.symbol`, `adj.P.Val`, and `logFC`.

It performs the following steps:
1. Preprocesses gene names to handle different formats.
2. Construct a candidate gene pool based on DEA with a loose cutoff.
3. Select top discriminative genes using the Information Gain (IG) filter approach.
4. Construct a PPI network using the STRING database and calculate edge weights based on gene expression correlations.
5. Analyzes the PPI network using centrality measures (Eigenvector and Betweenness) for Tumor and Normal samples.
6. Combines centrality rankings into an ensemble ranking to prioritize key genes.

## Installation

# Install required dependencies
install.packages(c("devtools","dplyr", "tibble", "FSelectorRcpp", "igraph", "data.table", "afc"))

```R
library(devtools)
devtools::install_github("nazpashaei/DiffCentEn")


## Example data preparation
data <- list(
  data = data.frame(
    gene1 = c(1, 2),
    gene2 = c(3, 4),
    class = c("Tumor", "Normal")
  ),
  topGenes = data.frame(
    Gene.symbol = c("gene1", "gene2"),
    adj.P.Val = c(0.01, 0.02),
    logFC = c(2, -2)
  )
)

## Run Differential Centrality-Ensemble Analysis
library(DiffCentEn)
KeyGenes <- DiffCentEn_function(data)
View(KeyGenes)

##Citation

