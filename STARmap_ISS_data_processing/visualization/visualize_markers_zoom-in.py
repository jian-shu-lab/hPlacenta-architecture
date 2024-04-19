import os
import numpy as np
import pandas as pd 
from tqdm import tqdm
import tifffile as tiff
import multiprocessing as mp
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.colors import LinearSegmentedColormap

import warnings
warnings.filterwarnings('ignore')
folder = 'tiles_results/List-5'
if not os.path.exists(folder):
    os.mkdir(folder)

### load data from excel
def _extract_info_from_excel(
        file_name, 
        num_genes = 4, gene_colors = ['#FF0000', '#40E0D0', '#FFD700', '#C0C0C0'], tiles_columns = ['tiles','tiles.1','tiles.2'],
        cell_type = 'unknown', sample = 'JS35', annotation = 'list-3b', important_order = 'left-to-right',
):
    df_genes = pd.read_excel(file_name)
    df_genes = df_genes.dropna(axis = 0, how = 'all')
    df_genes = df_genes[df_genes['sample'] == sample]
    df_genes = df_genes.iloc[:, 1:]

    n_max_regions = len(tiles_columns)
    genes_all = []; colors_all = []; overview_tiles_all = []; cell_types_all = []; annotations_all = []

    for i in range(df_genes.shape[0]):
        genes = df_genes.iloc[i, :num_genes].values
        genes = [i for i in genes if not pd.isna(i)]
        genes = [i.replace('ERV48-1', 'ERVH48-1') for i in genes]
        genes = [i.replace('BACH', 'BACH2') for i in genes]
        colors = gene_colors[:len(genes)]

        if important_order == 'left-to-right':
            genes = genes[::-1]
            colors = colors[::-1]

        regions = df_genes.iloc[i, -n_max_regions:].values
        regions = [i for i in regions if not pd.isna(i)]
        num_regions = len(regions)
        for j in range(num_regions):
            genes_all.append(genes)
            colors_all.append(colors)
            
            region_vector = list(str(regions[j]).split('/'))
            region_vector = [int(k) for k in region_vector]
            overview_tiles_all.append(region_vector)

            cell_types_all.append(cell_type)
            annotations_all.append(annotation)

    print(f'number of files in the end: {len(genes_all)}')
    return genes_all, colors_all, overview_tiles_all, cell_types_all, annotations_all

### parameters
tile_size = 2048
updated_tile_size = 2048 - 2*102
color_dict = {
    '#FF8C00': 'orange', '#8B008B': 'magenta', 
    '#00008B': 'blue', '#FF0000': 'red', 
    '#40E0D0': 'turquoise', '#228B22': 'forest_green', 
    '#FFA500': 'orange', '#FFD700': 'gold',
    '#C0C0C0': 'gray'}
font = font_manager.FontProperties(family='Arial')
# scale bar
pixel_size = 141.98 / 1000 # um
scale_factor = 20 / pixel_size  # Adjust this value as needed
scale_x = 100    # X-coordinate for the scale bar
s_bg = 0.040
s_gene_list = [2.0]
#  create colormap
colors = [(0, 0, 0), (0, 0, 1)]  # Black to blue
custom_cmap = LinearSegmentedColormap.from_list('custom_gray_blue', colors, N=256)
cmaps = ['gray']

### draw plots for local regions
def draw_gene_local_region(df_spots, dapi, gene, color, overview_tile, cell_type, annotation, s_gene, cmap = custom_cmap):
    unique_genes = df_spots.target.unique()
    gene_str = '_'.join(gene)
    color_str = '_'.join([color_dict[i] for i in color])
    cell_type_str = cell_type.replace('/', '-')
    unit_scale = 32 / 16596 * 2
    overview_tile = [int(i) for i in overview_tile]
    overview_tile_str = '_'.join([str(i) for i in overview_tile])
    print(f'running plotting for gene {gene_str} in tiles {overview_tile}')

    # overview
    df_overview = df_spots[df_spots['tile'].isin(overview_tile)].copy()
    xmin, xmax, ymin, ymax = df_overview.x.min(), df_overview.x.max(), df_overview.y.min(), df_overview.y.max()
    x_size = unit_scale * (xmax - xmin)
    y_size = unit_scale * (ymax - ymin)
    fig, ax = plt.subplots(figsize = (x_size, y_size), nrows = 1, ncols = 1)
    x, y, features = df_overview.x.values, df_overview.y.values, df_overview.target.values
    x = x - xmin 
    # y = y - ymin
    y = ymax - y

    for gene_i, color_i in zip(gene, color):
        if gene_i not in unique_genes:
            raise Exception(f'{gene_i} is not available in list')
        xt, yt = x[features == gene_i], y[features == gene_i]
        ax.scatter(xt, yt, s = s_gene**2, color = color_i, label = gene_i, marker = 'o')

    xf, yf = x[~np.isin(features, gene)], y[~np.isin(features, gene)]
    ax.scatter(xf, yf, s = s_bg**2, c = '#B0C4DE', alpha = 0.05, marker = 'o')

    # scale bar
    scale_y = np.max(y) - 100
    plt.plot([scale_x, scale_x + scale_factor], [scale_y, scale_y], color='white', linewidth=3)
    plt.annotate(f'20 µm', (scale_x + scale_factor / 2, scale_y - 30), color='white', ha='center')

    ax.axis('off')
    ax.set_facecolor('white')

    dapi_tile = dapi[ymin:ymax, xmin:xmax]; dapi_tile = dapi_tile[::-1, :] # 0 (dark) ~255 (white / bright)
    ax.imshow(dapi_tile, cmap = cmap, alpha = 1, vmin = 0, vmax = 255)
    
    stain_scale = 'Gray' if cmap == 'gray' else 'Blue'
    plt.savefig(
        f'{folder}/marker_{annotation}_{gene_str}_{color_str}_cell_type_{cell_type_str}_overview_tiles_{overview_tile_str}_dot_size={s_gene}_stain_{stain_scale}.png', 
        bbox_inches = 'tight', transparent = True, dpi = 300
    )


if __name__ == '__main__':
    samples = ['JS34', 'JS35', 'JS40']

    for sample in samples:
        # 1. extract genes to plot
        genes_all, colors_all, overview_tiles_all, cell_types_all, annotations_all = _extract_info_from_excel(
            file_name = 'List-5.xlsx',
            num_genes = 4, gene_colors = ['#FF0000', '#40E0D0', '#FFD700', '#C0C0C0'], tiles_columns = ['tiles', 'tiles.1'],
            cell_type = 'unknown', sample = sample, annotation = f'list-5_{sample}', important_order='left-to-right',
        )

        # 2. load data
        df_spots = pd.read_csv(f'/home/jinmr2/sample_integration/data/spots_global_{sample}.txt',sep='\t',header=0,index_col=0)
        xmax, ymax = df_spots.x.max(), df_spots.y.max()
        df_spots = df_spots[df_spots['x'] != xmax]
        df_spots = df_spots[df_spots['y'] != ymax]
        dapi = tiff.imread(f'/home/jinmr2/sample_integration/data/dapi_{sample}.tif')
        dapi = dapi[::-1, :]

        # 3. plotting
        unique_genes = df_spots.target.unique()
        for cmap in cmaps:
            for s_gene in s_gene_list:
                for gene, color, overview_tile, cell_type, annotation in zip(genes_all, colors_all, overview_tiles_all, cell_types_all, annotations_all):
                    draw_gene_local_region(df_spots, dapi, gene, color, overview_tile, cell_type, annotation, s_gene, cmap = cmap)