# Differential Centrality-Ensemble Analysis Based on Gene Expression Profiles and Protein-Protein Interaction Network

This R package provides a comprehensive framework for identifying and prioritizing disease-related genes based on gene expression profiles and protein-protein interaction (PPI) networks.

## Function Overview

The `DiCE_function` function takes a list of two elements:
- `data`: A data frame of gene expression data with gene symbols as column names and class labels (Tumor, Normal) in the last column.
- `topGenes`: A data frame of differentially expressed genes (DEGs) with columns `Gene.symbol`, `adj.P.Val`, and `logFC`.

It performs the following steps:
1. Preprocesses gene names to handle different formats.
2. Construct a candidate gene pool based on DEA with a loose cutoff.
3. Select top discriminative genes using the Information Gain (IG) filter approach.
4. Construct a PPI network using the STRING database and calculate edge weights based on gene expression correlations.
5. Analyzes the PPI network using centrality measures (Eigenvector and Betweenness) for Tumor and Normal samples.
6. Combines centrality rankings into an ensemble ranking to prioritize key genes.

```R
##Installation

Step 1: Install Required Dependencies
First, ensure that all necessary packages are installed:
install.packages(c("devtools", "dplyr", "tibble", "FSelectorRcpp", "igraph", "data.table", "NetWeaver"))

Step 2: Install Bioconductor Package
Ensure you have the BiocManager package installed, then use it to install the STRINGdb package from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("STRINGdb")

Step 3: Install Package from GitHub
Use devtools to install the DiffCentEn package from GitHub:
devtools::install_github("nazpashaei/DiCE")

## Example Data Preparation

This section provides an example of how to prepare and save the data required for the `DiCE_function`.

### Data Structure
- `data`: A list of two data frames:
  1. `data`: Gene expression data with class labels. The last column should be the class label ("Tumor" and "Normal").
  2. `topGenes`: Information about differentially expressed genes (DEGs), including gene symbols, adjusted p-values, and log fold changes.


data <- list(
  data = data.frame( 
    gene1 = c(0.938, 1.203),
    gene2 = c(2.107, 2.057),
    class = c("Tumor", "Normal")
  ),
  topGenes = data.frame(
    Gene.symbol = c("gene1", "gene2"),
    adj.P.Val = c(0.01, 0.02),
    logFC = c(2, -2)
  )
)
saveRDS(object = data,
        file = "data.RDS",
        compress = FALSE)

#########Parameters
--data: A list of two data frames:
*data – A data frame containing gene expression values. Gene symbols should be column names, and the last column must  contain class labels (e.g., "Tumor", "Normal").
*topGenes – A data frame with results from differential expression analysis. Expected columns: Gene.symbol, adj.P.Val, and logFC.

--regulation_status: A character string indicating which genes to consider based on their regulation status. Must be one of:
"Up": Upregulated genes
"Down": Downregulated genes
"Both": Both upregulated and downregulated genes

--species: A character string specifying the organism. Must be one of:
"human"
"mouse"
"rat"

## Run Differential Centrality-Ensemble Analysis
library(DiCE)
data <- readRDS("~/Ovarian_cancer.RDS");#Downloading and Reading an RDS File
KeyGenes <- DiCE_function(data,regulation_status = "Up");
View(KeyGenes)

##Citation
DiCE: differential centrality-ensemble analysis based on gene expression profiles and protein-protein interaction network
Elnaz Pashaei, Sheng Liu, Kailing Li, Yong Zang, Lei Yang, Tim Lautenschlaeger, Jun Huang, Xin Lu, Jun Wan
bioRxiv 2025.03.14.638654; doi: https://doi.org/10.1101/2025.03.14.638654
