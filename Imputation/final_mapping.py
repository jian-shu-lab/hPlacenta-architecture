# major codes refer to https://github.com/wanglab-broad/mCNS-atlas/tree/main/Imputation

import pickle
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

import warnings 
warnings.filterwarnings('ignore')

def _get_final_mapping_data(sample = 'JS40', n_neighbors = 100):
    with open(f'output/distances_{sample}_{n_neighbors}.pkl', 'rb') as f:
        distances = pickle.load(f) # n_cells * n_neighbors
    with open(f'output/indices_{sample}_{n_neighbors}.pkl', 'rb') as f:
        indices = pickle.load(f) # n_cells * n_neighbors

    adata = sc.read(f'output/adata_imput_expr_{sample}_{n_neighbors}.h5ad')

    adata_sc = sc.read_h5ad('/home/jinmr2/20230814_iss_hplac_r1-r6/data/single_cell_reference/data/hplacenta_gene_matrix.h5ad')
    adata_sc_imputation_X = adata_sc.X.toarray()

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
        var = pd.DataFrame(index = adata_sc.var.index)
    )
    adata_imput_expr.write(f'output/adata_imput_expr_{sample}_{n_neighbors}_final.h5ad')
    
if __name__ == '__main__':
    samples = ['JS34', 'JS35', 'JS36', 'JS40']

    for sample in samples:
        _get_final_mapping_data(sample = sample, n_neighbors = 100)