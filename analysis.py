import scanpy as sc
import h5py
import matplotlib.pyplot as plt

import pandas as pd


adata = sc.read_10x_h5("/Users/dheerajpv/Documents/Intern/Single-Cell RNA-Seq Analysis of Immune Cell Populations/10k_5p_Human_diseased_PBMC_ALL_Fix_count_filtered_feature_bc_matrix.h5")

adata.var_names_make_unique()

print(adata)

adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filtering cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)



sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var['highly_variable']]



sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')  # Example of coloring by a specific gene

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

sc.tl.louvain(adata)
sc.pl.umap(adata, color=['louvain', 'CST3'])  # Color by clusters and gene


sc.pl.umap(adata, color=['louvain'])


adata.write("processed_scRNAseq_data.h5ad")

import scanpy as sc

adata = sc.read("processed_scRNAseq_data.h5ad")

sc.pl.umap(adata, color='louvain')  # Example: visualize UMAP with clusters

print(adata.obs['louvain'].unique())

sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')

sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

adata.X = adata.X.astype('float32')
adata.obs['louvain'] = adata.obs['louvain'].astype('category')

def improved_rank_genes_groups(adata, groupby):
    names_list = []
    scores_list = []
    pvals_list = []
    pvals_adj_list = []
    logfoldchanges_list = []


    result = sc.tl.rank_genes_groups(adata, groupby, method='t-test', copy=True)

    for group in result.uns['rank_genes_groups']['names'].dtype.names:
        names_list.append(result.uns['rank_genes_groups']['names'][group])
        scores_list.append(result.uns['rank_genes_groups']['scores'][group])
        pvals_list.append(result.uns['rank_genes_groups']['pvals'][group])
        pvals_adj_list.append(result.uns['rank_genes_groups']['pvals_adj'][group])
        logfoldchanges_list.append(result.uns['rank_genes_groups']['logfoldchanges'][group])


    df = pd.DataFrame({
        'names': pd.concat(names_list, axis=1),
        'scores': pd.concat(scores_list, axis=1),
        'pvals': pd.concat(pvals_list, axis=1),
        'pvals_adj': pd.concat(pvals_adj_list, axis=1),
        'logfoldchanges': pd.concat(logfoldchanges_list, axis=1)
    })

    return df

top_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
print(top_genes)


# Ensuring that the data is in float format for plotting
adata.X = adata.X.astype(float)


sc.pl.heatmap(adata, var_names=top_genes.values.flatten(), groupby='louvain', use_raw=False)


sc.pl.violin(adata, keys=top_genes.values.flatten()[:10], groupby='louvain', use_raw=False)


sc.get.rank_genes_groups_df(adata, group='0')  # Replace '0' with a cluster number

# Saving top-ranked genes
top_genes.to_csv("top_ranked_genes.csv")

# Saving differential expression results
adata.write("processed_data_with_ranked_genes.h5ad")
