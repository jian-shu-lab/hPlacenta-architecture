import os
import numpy as np
import scanpy as sc
import tifffile as tiff
import matplotlib as mpl
import multiprocessing as mp
import matplotlib.pyplot as plt


### params
version = ''
scale = 0.0007
s = 5.0
vmax_pct = 95
vmin_pct = 0
min_expr_imputation = 0.5
min_expr_observed = 0.05
n_neighbors = 100

# scale bar metadata
pixel_size = 141.98 / 1000 # um
scale_factor = 100 / pixel_size
scale_x = 600
fig_folder = 'figures_imputation_comparison_v6'

plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'font.family': 'Arial'})

if not os.path.exists(fig_folder):
    os.makedirs(fig_folder)

def _get_expression_ax(x, y, expr, dapi, ax, tile = 'Imputed'):
    ax.imshow(dapi, cmap = 'Greys', alpha = 1.0)

    ## remove low expression genes
    # if tile == 'Imputed':
    #     expr[expr < min_expr_imputation] = 0
    # elif tile == 'Observed':
    #     expr[expr < min_expr_observed] = 0
    
    # x, y, expr = x[expr > 0], y[expr > 0], expr[expr > 0]
    
    vmin = -0.01

    vmin, vmax = np.percentile(expr, vmin_pct), np.percentile(expr, vmax_pct)
    y = dapi.shape[0] - y
    ax.scatter(x, y, c=expr, cmap='Reds', s = s, vmax=vmax, vmin=vmin)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=expr.max())    # colorbar
    sm = plt.cm.ScalarMappable(cmap='Reds', norm=norm)

    scale_y = np.max(y) - 200
    ax.plot([scale_x, scale_x + scale_factor], [scale_y, scale_y], color='black', linewidth=5)
    ax.annotate(f'100µm', (scale_x + scale_factor / 2, scale_y - 150), color='black', ha='center', fontsize = 14)

    sm.set_array([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(tile)
    return ax, sm

def draw_expression_two_subplots(adata, adata_observed, sample, gene):
    # imputed
    x, y = adata.obs['x_pixel'].values, adata.obs['y_pixel'].values
    expr = adata[:, gene].X.toarray().flatten()

    x_range = x.max() - x.min()
    y_range = y.max() - y.min()

    dapi = tiff.imread(f'/home/jinmr2/sample_integration/data/dapi_{sample}.tif')

    if gene in adata_observed.var.index.values:
        fig, axes = plt.subplots(figsize=(x_range * scale * 2, y_range * scale), nrows = 1, ncols = 2)
        axes[0], sm = _get_expression_ax(x, y, expr, dapi, axes[0], tile = 'Imputed')
        fig.colorbar(sm, ax=axes[0], pad=0.01, shrink=0.8)
        
        # observed
        x_observed, y_observed = adata_observed.obs['x_pixel'].values, adata_observed.obs['y_pixel'].values
        expr_observed = adata_observed[:, gene].X.toarray().flatten()
        axes[1], sm = _get_expression_ax(x_observed, y_observed, expr_observed, dapi, axes[1], tile = 'Observed')
        fig.colorbar(sm, ax=axes[1], pad=0.01, shrink=0.8)
        fig.subplots_adjust(wspace=0.05, hspace=0.1)
    else:
        fig, axes = plt.subplots(figsize=(x_range * scale, y_range * scale), nrows = 1, ncols = 1)
        axes, sm = _get_expression_ax(x, y, expr, dapi, axes, tile = 'Imputed')
        fig.colorbar(sm, ax=axes, pad=0.01, shrink=0.8)    
    
    fig.savefig(f'{fig_folder}/{gene}_{sample}{version}.png', dpi=300)
    plt.close(fig)

def _draw_imputation_samples(adata, adata_observed, sample, genes):
    # multiprocessing
    pool = mp.Pool(processes=mp.cpu_count())
    results = [pool.apply_async(draw_expression_two_subplots, args=(adata, adata_observed, sample, gene)) for gene in genes]
    output = [p.get() for p in results]
    pool.close()

if __name__ == '__main__':
    # 1. load data
    samples = ['JS34', 'JS35', 'JS36', 'JS40']
    adatas = []
    adatas_observed = []
    dapis = []
    for sample in samples:
        adatas.append(sc.read(f'/home/jinmr2/sample_integration/three_samples/imputation/v5_seurat_all_genes/output/adata_imput_expr_{sample}_{n_neighbors}_final.h5ad'))
        adatas_observed.append(sc.read(f'/home/jinmr2/sample_integration/three_samples/raw_clustermap_{sample}.h5ad'))
    
    # 2. select genes
    genes_selection = ['COL3A1', 'CYP19A1', 'HAPLN3', 'AOC1', 'MSI2', 'ERVH48-1', 'BACH2', 'MYCNUT', 'KANK1', 'EFEMP1', 'TFPI',
                       'TFPI2', 'QSOX1', 'GDF15', 'MSX1', 'PXDN', 'BMP2', 'CGA', 'PAGE4', 'NOTUM', 'VIM', 'KDR', 'CD163', 'PARP1',
                       'IDH1', 'ABL1', 'SMARCA4', 'CEBPA', 'HRAS', 'CPS1', 'SPINT1', 'MKI67', 'CDK1', 'CDK7', 'HSPG2', 'ADAM19', 
                       'ITGB4', 'PAPPA2', 'KRT7', 'KRT19', 'KRT23', 'ADAM12', 'EBI3', 'ADAMTS20', 'ANK2', 'COL27A1', 'COL4A1']
    genes_all = genes_selection

    # 3. draw side-by-side comparison
    for sample, adata, adata_observed in zip(samples, adatas, adatas_observed):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        sc.pp.filter_cells(adata_observed, min_counts=80)
        sc.pp.normalize_total(adata_observed, target_sum=1e3)
        sc.pp.log1p(adata_observed)

        _draw_imputation_samples(adata, adata_observed, sample, genes_all)