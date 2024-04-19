import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

# 0. params
pixel_size = 141.98 / 1000; scale_factor = 200 / pixel_size; scale_x = 100 # scale bar
color_dict = {
    'EVT2': '#208A42', #
    'EVT3': '#89288F', #
    'Endothelial_cells': '#D51F26', #
    'Fibroblast1': '#F47D2B', #
    'Fibroblast2': '#FEE500', #
    'Hofbauer cells': '#8A9FD1', #
    'STB': '#90D5E4', #
    'maternal macrophages ': '#C06CAB', #
    'unknown2': '#F37B70', #
    'unknown3': '#9983BD', #
    'vCTB2': '#3BBCA9', #
    'vCTBp': '#0C727C',
}

def _draw_segmentation(sample, nrows, ncols):
    # 1. load data
    df_spots = pd.read_csv(f'../data/spots_global_{sample}.txt',sep='\t',header=0,index_col=0, keep_default_na=False)
    xmax, ymax = df_spots.x.max(), df_spots.y.max()
    df_spots = df_spots[df_spots['x'] != xmax]
    df_spots = df_spots[df_spots['y'] != ymax]

    # 2. add annotations
    df_anno = pd.read_csv('label_transferred_annotations_all.csv', index_col = 0)
    df_anno = df_anno[df_anno['sample'] == sample].copy()
    df_anno['index'] = [i.split(f'{sample}_')[1].split('-')[0] for i in df_anno.index.values]
    df_anno.set_index('index', inplace = True)

    df_spots = df_spots.loc[df_spots['cell'].isin(df_anno.index.values), :].copy()
    print(f'sample {sample}: number of cells in the spots file: {len(np.unique(df_spots.cell.values))}')
    mapped_columns = ['predicted.celltype.score', 'predicted.celltype']
    for col in mapped_columns:
        df_spots[col] = df_anno.loc[df_spots.cell.values, col].values 

    # 3. scatter plots
    fig, ax = plt.subplots(figsize = (ncols*2, nrows*2))
    unique_cell_types = np.unique(df_spots['predicted.celltype'])
    x, y, cells, cell_types = df_spots.x.values, df_spots.y.values, df_spots['cell'].values, df_spots['predicted.celltype'].values
    for cell_type in unique_cell_types:
        xc, yc = x[cell_types == cell_type], y[cell_types == cell_type]
        ax.scatter(xc, yc, s = 0.5, color = color_dict[cell_type], alpha = 1.0, label = cell_type)

    # 4. draw convex hull for cell boundaries
    for cell in tqdm(np.unique(cells)):
        x_cell, y_cell = x[cells == cell], y[cells == cell]
        points = np.array([x_cell, y_cell]).T
        try:
            hull = ConvexHull(points)
            for simplex in hull.simplices:
                ax.plot(points[simplex, 0], points[simplex, 1], 'k-', linewidth = 0.5)
        except:
            continue

    # 5. add scale bar
    scale_y = 100
    plt.plot([scale_x, scale_x + scale_factor], [scale_y, scale_y], color='black', linewidth=6)
    plt.annotate(f'200 µm', (scale_x + scale_factor / 2, scale_y + 200), color='black', ha='center', fontsize = 30)

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())    
    ax.axis('off')

    # add legend to the right of the figure
    ax.legend(markerscale = 30, fontsize = 30, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.savefig(f'starmap_{sample}_transcripts_cells_annotations.png', dpi = 300, bbox_inches = 'tight')

if __name__ == '__main__':
    samples = ['JS34', 'JS35', 'JS36', 'JS40']
    nrows_list = [6, 9, 7, 9]
    ncols_list = [10, 11, 10, 7]
    n_tiles_list = [42, 48, 44, 43]

    for sample, nrows, ncols in zip(samples, nrows_list, ncols_list):
        _draw_segmentation(sample, nrows, ncols)