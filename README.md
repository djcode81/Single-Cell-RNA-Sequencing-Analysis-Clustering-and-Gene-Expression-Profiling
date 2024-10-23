# Single-Cell RNA Sequencing Analysis: Clustering and Gene Expression Profiling

This project focuses on analyzing single-cell RNA sequencing data. It includes clustering and gene expression profiling using PCA, UMAP, and Louvain clustering techniques. The analysis is performed on the 10k human diseased PBMC (ALL) dataset, which was processed through Scanpy.

## Project Overview

1. **Data Preprocessing**: The raw HDF5 data was preprocessed, normalized, and filtered to remove low-quality cells.
2. **Clustering**: Cells were clustered using the Louvain algorithm to identify distinct cell populations.
3. **Dimensionality Reduction**: Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP) were used to reduce the dimensionality of the data for visualization.
4. **Differential Gene Expression**: Gene expression profiling was performed to identify marker genes across clusters.
5. **Visualization**: Heatmaps and UMAP plots were used to visualize clusters and differentially expressed genes.

## Files

- `analysis.py`: This is the python script for the analysis.
- `processed_data.h5ad`: The processed dataset used for the analysis.
- `results/`: Contains output files such as heatmaps, PCA, and UMAP plots.
- `10k_5p_Human_diseased_PBMC_ALL_Fix_count_filtered_feature_bc_matrix.h5`: This the dataset used for the project, source : https://www.10xgenomics.com/datasets/10k_5p_Human_diseased_PBMC_ALL_Fresh

## Analysis Pipeline

The analysis workflow is as follows:

1. **Data Loading and Preprocessing**:
   - The scRNA-seq data is loaded from the provided `.h5` file.
   - Quality control steps, including filtering cells and genes, are applied.
   - The data is normalized and log-transformed for further analysis.
   
2. **Dimensionality Reduction and Clustering**:
   - Highly variable genes are identified.
   - Principal Component Analysis (PCA) is performed for dimensionality reduction.
   - UMAP is used to visualize the data in two dimensions.
   - Clustering is conducted using the Louvain algorithm to identify distinct groups of cells.

3. **Differential Gene Expression**:
   - The top-ranked differentially expressed genes for each cluster are calculated using a t-test.
   - The ranked genes are saved in a CSV file and used for further downstream analysis.

4. **Visualization**:
   - UMAP plots are generated to visualize the clustering results.
   - Heatmaps and violin plots are created to visualize the expression of top genes across clusters.


## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/djcode81/Single-Cell-RNA-Sequencing-Analysis-Clustering-and-Gene-Expression-Profiling.git

2. After setting up your environment, run the following command:

   ```bash
   python analysis.py

