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
install.packages(c("devtools", "dplyr", "tibble", "FSelectorRcpp", "igraph", "data.table", "afc"))

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

This section provides an example of how to prepare and save the data required for the `DiffCentEn_function`.

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

#### Parameters
- `data`: A list of two data frames containing gene expression data with a class label at the last column and a list of information about DEGs analysis.
- `regulation_status`: A character Value indicating the regulation status of genes. It must be one of the following:
  - `"Up"`: For upregulated genes.
  - `"Down"`: For downregulated genes.
  - `"Both"`: For both upregulated and downregulated genes.


## Run Differential Centrality-Ensemble Analysis
library(DiCE)
data <- readRDS("~/Ovarian_cancer.RDS");#Downloading and Reading an RDS File
KeyGenes <- DiCE_function(data,regulation_status = "Up");
View(KeyGenes)

##Citation

