import os
import numpy as np
import pandas as pd 
import scanpy as sc
import tifffile as tiff

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

import warnings
warnings.filterwarnings('ignore')
from tqdm import tqdm

if not os.path.exists('slide_results'):
    os.mkdir('slide_results')

pixel_size = 141.98 / 1000; scale_factor = 200 / pixel_size; scale_x = 100 # scale bar
def _draw_whole_slice_multichannels(sample, num_rows, num_cols):
    # load spots
    df_spots = pd.read_csv(f'/home/jinmr2/sample_integration/data/spots_global_{sample}.txt',sep='\t',header=0,index_col=0, keep_default_na=False)
    x, y, features = df_spots.x.values, df_spots.y.values, df_spots.target.values 

    ### multiple color
    fig, ax = plt.subplots(figsize = (2*num_cols, 2*num_rows))
    genes = ['PAGE4', 'CGA', 'VIM', 'NOTUM']
    # colors = ['orange', 'dark magenta', 'turquoise', 'dark blue']
    colors = ['#FFA500', '#8B008B', '#40E0D0', '#00008B']
    for gene, color in zip(genes, colors):
        xt, yt = x[features == gene], y[features == gene]
        ax.scatter(xt, yt, s = 0.35, color = color, label = gene)

    xf, yf = x[~np.isin(features, genes)], y[~np.isin(features, genes)]
    ax.scatter(xf, yf, s = 0.025, c = '#B0C4DE', alpha = 0.05)

    scale_y = 100
    plt.plot([scale_x, scale_x + scale_factor], [scale_y, scale_y], color='black', linewidth=6)
    plt.annotate(f'200 µm', (scale_x + scale_factor / 2, scale_y + 200), color='black', ha='center', fontsize = 30)

    ax.legend(fontsize = 20, markerscale = 20, loc='center left', bbox_to_anchor=(1.05, 1.0))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xticklabels([]); ax.set_yticklabels([])
    genes_str = '_'.join(genes)
    fig.savefig(f'slide_results/sample_{sample}_whole_slice_{genes_str}.png', dpi = 300, bbox_inches = 'tight')

if __name__ == '__main__':
    samples = ['JS34', 'JS35', 'JS36', 'JS40']
    nrows_list = [6, 9, 7, 9]
    ncols_list = [10, 11, 10, 7]

    for sample, nrows, ncols in zip(samples, nrows_list, ncols_list):
        _draw_whole_slice_multichannels(sample, nrows, ncols)