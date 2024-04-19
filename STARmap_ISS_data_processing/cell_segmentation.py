import os
import numpy as np
import pandas as pd
import tifffile as tiff
import multiprocessing as mp
import matplotlib.pyplot as plt
from ClusterMap.clustermap import *

version = 'v1'
if not os.path.exists('tiles_segmentation_results'):
    os.mkdir('tiles_segmentation_results')

def load_data(tile):
    # 1. load spots
    folder_spots = '/home/jinmr2/20230925_iss_hplacenta_jd34/data/spots_analysis/tiles/'
    folder_dapi = '/home/jinmr2/20230925_iss_hplacenta_jd34/data/dapi_images/'
    folder_save = 'tiles_segmentation_results/'
    df_spots = pd.read_csv(f'{folder_spots}tile_{tile}.txt', sep = '\t', header = 0, index_col = 0, keep_default_na = False)
    df_spots = df_spots[df_spots['target'] != 'nan'].copy()
    df_spots['spot_indices'] = df_spots.index.values
    df_spots['gene_name'] = df_spots['target'].values # dummy column
    df_spots['x'] = df_spots['x'].astype(np.int16)
    df_spots['y'] = df_spots['y'].astype(np.int16)
    df_spots.rename(columns={'x':'spot_location_1','y':'spot_location_2'}, inplace = True)

    # 2. load image data
    img = tiff.imread(f'{folder_dapi}dapi_tile_whole_{tile}_v3.tif')

    # 3. get gene list
    gene_name_dict = dict(zip(df_spots['gene_name'].unique(), np.arange(1, 1+len(df_spots['gene_name'].unique()))))
    df_spots['gene'] = [gene_name_dict[i] for i in df_spots['gene_name'].values]
    df_spots['index'] = np.arange(0, df_spots.shape[0])
    df_spots.set_index('index', inplace = True)
    gene_list = df_spots['gene'].unique()

    # 4. draw plot 
    fig, ax = plt.subplots(figsize = (8,8))
    x, y = df_spots.spot_location_1.values, df_spots.spot_location_2.values 
    ax.scatter(x, y, s = 0.05, color = 'blue', alpha = 0.5)
    plt.imshow(img)
    fig.savefig(f'{folder_save}tile_{tile}_dapi_spots.png', bbox_inches='tight')

    return df_spots, img, gene_list

def define_model(df_spots, img, gene_list, xy_radius = 40, sigma = 1):
    model = ClusterMap(
        spots = df_spots,
        dapi = img,
        gene_list = gene_list,
        num_dims = 2,
        xy_radius = xy_radius,
        z_radius = 0,
        fast_preprocess = True,
        gauss_blur = True,
        sigma = sigma,
    )
    return model

def preprocess(model, pct_filter = 0):
    range_1 = model.spots.spot_location_2.max() - model.spots.spot_location_2.min()
    range_2 = model.spots.spot_location_1.max() - model.spots.spot_location_1.min()
    # model.dapi_binary = np.zeros((range_1, range_2)) # create an empty binary dapi
    
    model.preprocess(LOF = True, pct_filter = pct_filter)
    
    # model.dapi_binary = None

def run_segmentation(
    model, cell_num_threshold = 0.01, dapi_grid_interval = 5, add_dapi = True, use_genedis = True,
    savename = None,
):
    model.segmentation(
        cell_num_threshold = cell_num_threshold, 
        dapi_grid_interval = dapi_grid_interval, 
        add_dapi=add_dapi,
        use_genedis=True,
    )
    if savename is not None:
        model.spots.to_csv(savename, sep = '\t')
    else:
        model.spots.to_csv('whole_clustermap_noimage_result.txt', sep = '\t')

def run_pipeline_tile(
    tile,
    use_dapi = False,
    xy_radius = 55,
    cell_num_threshold = 0.0001,
    dapi_grid_interval = 3,
    pct_filter = 1.0
):
    folder = 'tiles_segmentation_results/'
    df_spots, img, gene_list = load_data(tile = tile)
    tiff.imwrite(f'{folder}dapi_clustermap_cropped_tile_{tile}.tif', img)
    print(df_spots.spot_location_1.min(), df_spots.spot_location_1.max())
    print(df_spots.spot_location_2.min(), df_spots.spot_location_2.max())
    print(img.shape)
    print(f'setting up the model for tile {tile}')
    model = define_model(df_spots, img, gene_list, xy_radius = xy_radius, sigma = 1)
    
    if use_dapi:
        print(f'preprocessing for tile {tile}')
        preprocess(model, pct_filter = pct_filter)
    print(f'running segmentation for tile {tile}')
    run_segmentation(
        model, cell_num_threshold = cell_num_threshold, dapi_grid_interval = dapi_grid_interval, add_dapi = True, use_genedis = True,
        savename = f'{folder}spots_decoded_all_stitched_globalcoords_{version}_tile_{tile}.txt',
    )
    print(f'segmentation completed for tile {tile}')

if __name__ == '__main__':

    tiles = np.arange(42)
    n_process = 14
    with mp.Pool(processes=n_process) as pool:
        items = list(tiles)
        results = pool.map(run_pipeline_tile, items)
