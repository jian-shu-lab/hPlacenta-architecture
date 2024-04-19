# major codes refer to https://github.com/wanglab-broad/mCNS-atlas/tree/main/Imputation

import os
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

from anndata import AnnData
from scipy.stats import pearsonr

import warnings
warnings.filterwarnings('ignore')

if not os.path.exists('output'):
    os.makedirs('output')

def load_data(min_counts = 80):
    ### read single cell data
    adata_sc=sc.read_h5ad('/home/jinmr2/20230814_iss_hplac_r1-r6/data/single_cell_reference/data/hplacenta_gene_matrix.h5ad')
    adata_sc.var.index=adata_sc.var.index.str.upper()
    adata_sc.var_names_make_unique()

    ###read STARmap data
    adata = sc.read_h5ad('/home/jinmr2/sample_integration/three_samples/raw_clustermap_all.h5ad')
    bdata = sc.read_h5ad('/home/jinmr2/sample_integration/three_samples/processed_clustermap_all_v2.h5ad')
    adata = adata[bdata.obs.index.values, :].copy()
    adata.obs = bdata.obs.copy()

    sc.pp.filter_cells(adata, min_counts = min_counts)
    predict_genelist_whole = [i for i in adata.var.index.values if i in adata_sc.var.index.values] 
    predict_genelist = predict_genelist_whole

    adata_predict = adata.copy()
    adata_org=adata.copy()
    adata=adata_predict
    adata.var.index=adata.var.index.str.upper()

    list_of_variable_names=adata.var.index.intersection(adata_sc.var.index)
    adata_subset = adata[:, list_of_variable_names]

    list_of_variable_names=adata_sc.var.index.intersection(adata.var.index)
    adata_sc_subset = adata_sc[:, list_of_variable_names]
    return adata_subset, adata_sc_subset, predict_genelist

def integrate_query_ref(adata_subset, adata_sc_subset):
    ### using previous seurat integration results (instead of harmony)
    adata_subset.obs['dataset']='st'
    adata_sc_subset.obs['dataset']='scrna'
    combine_adata=adata_subset.concatenate(adata_sc_subset, index_unique=None)

    umaps_query = pd.read_csv('../../umaps_query.csv', header = 0, index_col = 0)
    umaps_query.columns = ['UMAP_1', 'UMAP_2']
    umaps_ref = pd.read_csv('../../umaps_ref.csv', header = 0, index_col = 0)
    umaps_combined = pd.concat([umaps_query, umaps_ref], axis = 0)

    print(f'number of cells in umaps_combined and combine_adata: {umaps_combined.shape[0]}, {combine_adata.shape[0]}')
    print(f'number of overlapping cells: {len(set(umaps_combined.index.values).intersection(combine_adata.obs.index.values))}')

    umaps_combined = umaps_combined.loc[combine_adata.obs.index.values, :]
    print(umaps_combined.head())

    combine_adata.obsm['X_umap'] = umaps_combined.values

    return combine_adata

def draw_integration(combine_adata):
    
    if not os.path.exists('figures'):
        os.makedirs('figures')
    
    cluster_pl = sns.color_palette("tab20",4)
    cluster_pl[0]=np.array([252,40,3])/255
    cluster_pl[1]=np.array([3, 169, 252])/255
    cluster_pl[2]=np.array([4,217,61])/255
    sns.palplot(cluster_pl)

    ### save umap plot
    fig,axes = plt.subplots(figsize=(17,8), ncols = 2, nrows = 1)
    ax=sc.pl.umap(combine_adata, color='dataset',legend_fontsize=13,ax=axes[0],show=False)
    ax=sc.pl.umap(combine_adata, color='Clusters',legend_fontsize=13,ax=axes[1],show=False)
    ax.set_title('Harmony')
    ax.title.set_fontsize(20)
    fig.savefig(f'figures/umap_integration_datasets_all_samples.png', bbox_inches='tight', dpi=300)

    # umap for st cells only
    combind_adata_st = combine_adata[combine_adata.obs['dataset'] == 'st', :]
    adata_processed = sc.read('/home/jinmr2/sample_integration/three_samples/processed_clustermap_all_v2.h5ad')
    adata_processed = adata_processed[combind_adata_st.obs.index.values, :]
    combind_adata_st.obs['leiden'] = adata_processed.obs['leiden'].values
    fig,axes = plt.subplots(figsize=(17,8), ncols = 2, nrows = 1)
    sc.pl.umap(combind_adata_st, color='sample',legend_fontsize=13,ax=axes[0],show=False)
    sc.pl.umap(combind_adata_st, color='leiden',legend_fontsize=13,ax=axes[1],show=False)
    ax.set_title('Harmony')
    ax.title.set_fontsize(20)
    fig.savefig(f'figures/umap_query_samples_all_samples.png', bbox_inches='tight', dpi=300)

def compute_imputation(combine_adata, adata_sc, predict_genelist, sample):
    ### load data
    adata_sc=sc.read_h5ad('/home/jinmr2/20230814_iss_hplac_r1-r6/data/single_cell_reference/data/hplacenta_gene_matrix.h5ad')
    adata_sc.var.index=adata_sc.var.index.str.upper()
    adata_sc.var_names_make_unique()

    ########prepare for imputation
    combine_adata_st=combine_adata[combine_adata.obs['dataset']=='st',:].copy()
    combine_adata_st = combine_adata_st[combine_adata_st.obs['sample'] == sample, :]
    combine_adata_sc=combine_adata[combine_adata.obs['dataset']=='scrna',:].copy()

    print(f'dimension of combine_adata_st: {combine_adata_st.shape}')
    print(f'dimension of combine_adata_sc: {combine_adata_sc.shape}')
    print(combine_adata_sc.obsm['X_umap'][:,:10])
    print(combine_adata_st.obsm['X_umap'][:,:10])

    adata_sc_imputation=adata_sc.copy()
    adata_sc_imputation_subset=adata_sc_imputation[:,[list(adata_sc_imputation.var.index).index(i) for i in predict_genelist]]
    singlecell_rawexpr=adata_sc_imputation_subset.X
    shape1=singlecell_rawexpr.shape[1]

    adata = sc.read_h5ad('/home/jinmr2/sample_integration/three_samples/raw_clustermap_all.h5ad')
    adata = adata[adata.obs['sample'] == sample, :]
    sc.pp.filter_cells(adata, min_counts = 80)
    adata_real=adata[:,[list(adata.var.index).index(i) for i in predict_genelist]]
    adata_realexpr=adata_real.X

    pearsonrlist={'gene':[],'value':[],'pvalue':[],'n_neighbors':[]}

    for n_neighbors in tqdm([5, 10, 20, 25, 40, 50, 100, 200, 400]):
        neigh = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto').fit(combine_adata_sc.obsm['X_umap'])
        distances, indices = neigh.kneighbors(combine_adata_st.obsm['X_umap'])

        indices_all=np.concatenate(indices)

        adata_st_imputation=[]

        shape1=singlecell_rawexpr.shape[1]
        for i in range(0,indices_all.shape[0],100000):

            test1=singlecell_rawexpr[indices_all[i:i+100000]].toarray()
            test2=test1.reshape(int(test1.shape[0]/n_neighbors),n_neighbors,shape1)
            test3=np.mean(test2,axis=1)
            if i==0:
                adata_st_imputation=test3
            else:
                adata_st_imputation=np.concatenate((adata_st_imputation,test3),axis=0)

        adata_imput_expr = AnnData(
            X = adata_st_imputation,
            obs = combine_adata_st.obs,
            var = pd.DataFrame(index = predict_genelist)
        )
        adata_imput_expr.var_names_make_unique()
        adata_imput_expr.write(f'output/adata_imput_expr_{sample}_{n_neighbors}.h5ad')

        # compute pearsonr
        focus_gene=0
        gene_name=predict_genelist
        dis1=adata_imput_expr.X
        dis2=adata_realexpr.toarray()

        # total correlation
        pearsonrlist['gene'].append('overlapping_genes')
        pearsonrlist['value'].append(pearsonr(dis1.reshape(-1),dis2.reshape(-1))[0])
        pearsonrlist['pvalue'].append(pearsonr(dis1.reshape(-1),dis2.reshape(-1))[1])
        pearsonrlist['n_neighbors'].append(n_neighbors)

        # per-gene correlation
        df_pearsonrGenes = pd.DataFrame(index = predict_genelist)
        pearsonrGenes = [pearsonr(dis1[:,i],dis2[:,i]) for i in range(dis1.shape[1])]
        Genes = predict_genelist
        df_pearsonrGenes['genes'] = 'overlapping_genes'
        df_pearsonrGenes['pearsonr'] = [i[0] for i in pearsonrGenes]
        df_pearsonrGenes['pvalue'] = [i[1] for i in pearsonrGenes]
        df_pearsonrGenes.to_csv(f'output/pearsonrGenes_{sample}_{n_neighbors}.csv')
        
    pd.DataFrame(pearsonrlist).to_csv(f'output/pearsonrlist_{sample}.csv')

if __name__ == '__main__':
    samples = ['JS34', 'JS35', 'JS36', 'JS40']
    
    adata_subset, adata_sc_subset, predict_genelist = load_data()
    print(f'dimensions of adata_subset: {adata_subset.shape}, dimensions of adata_sc_subset: {adata_sc_subset.shape}, number of genes: {len(predict_genelist)}')
    combine_adata = integrate_query_ref(adata_subset, adata_sc_subset)
    combine_adata.write_h5ad(f'output/adata_seurat_integrate_all_samples.h5ad')

    draw_integration(combine_adata)
    for sample in samples:
        compute_imputation(combine_adata, adata_sc_subset, predict_genelist, sample)