# major codes refer to https://github.com/wanglab-broad/mCNS-atlas/tree/main/Imputation

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

import warnings 
warnings.filterwarnings('ignore')

if not os.path.exists('output'):
    os.makedirs('output')

def _get_final_mapping_data(sample = 'JS40', n_neighbors = 100):
    with open(f'../v5_seurat_all_genes/output/distances_{sample}_{n_neighbors}.pkl', 'rb') as f:
        distances = pickle.load(f) # n_cells * n_neighbors
    with open(f'../v5_seurat_all_genes/output/indices_{sample}_{n_neighbors}.pkl', 'rb') as f:
        indices = pickle.load(f) # n_cells * n_neighbors

    adata = sc.read(f'../v5_seurat_all_genes/output/adata_imput_expr_{sample}_{n_neighbors}.h5ad')
    df_dorc = pd.read_csv(
        '/home/jinmr2/sample_integration/three_samples/imputation/v6_map_ATAC/single-cell_atac_data_reference/p2g_500Kbp_all_dorc_thresh10.tsv', 
        sep = '\t', index_col = 0, header = 0
    )
    adata_sc_imputation_X = df_dorc.values

    for idx, i in enumerate(indices):
        if idx % 500 == 0:
            print(sample, idx)
        test = 1/distances[idx]
        test = test/test.sum()
        test1 = np.dot(test, adata_sc_imputation_X[i, :])
        if idx == 0:
            st_imputation = test1
        else:
            st_imputation = np.vstack((st_imputation, test1))

    adata_imput_expr = AnnData(
        X = st_imputation,
        obs = adata.obs.copy(),
        var = pd.DataFrame(index = df_dorc.columns)
    )
    adata_imput_expr.write(f'output/adata_imput_atac_dorc_{sample}_{n_neighbors}_final.h5ad')
    
if __name__ == '__main__':
    samples = ['JS34', 'JS35', 'JS36', 'JS40']

    for sample in samples:
        _get_final_mapping_data(sample = sample, n_neighbors = 100)