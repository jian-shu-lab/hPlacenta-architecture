import os 
import pickle
import warnings
from tqdm import tqdm
 
import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
if not os.path.exists('tiles'):
    os.mkdir('tiles')

# params
row_dict = {
    1: 5, 2: 5, 3: 5, 
    4: 4, 5: 4, 6: 4, 7: 4, 8: 4, 9: 4, 
    10: 3, 11: 3, 12: 3, 13: 3, 14: 3, 15: 3, 16: 3, 17: 3, 18: 3, 
    19: 2, 20: 2, 21: 2, 22: 2, 23: 2, 24: 2, 25: 2, 26: 2, 27: 2, 28: 2, 
    29: 1, 30: 1, 31: 1, 32: 1, 33: 1, 34: 1, 35: 1, 36: 1, 
    37: 0, 38: 0, 39: 0, 40: 0, 41: 0, 42: 0,
}

col_dict = {
    10:0,28:0,
    11:1,27:1,
    9:2,12:2,26:2,29:2,
    1:3,8:3,13:3,25:3,30:3,42:3,
    2:4,7:4,14:4,24:4,31:4,41:4,
    3:5,6:5,15:5,23:5,32:5,40:5,
    5:6,16:6,22:6,33:6,39:6,
    4:7,17:7,21:7,34:7,38:7,
    18:8,20:8,35:8,37:8,
    19:9,36:9,
}

### parameters
tiles = np.arange(42) # 0~41
num_z = 89
num_rows = 6
num_cols = 10
tile_size = size_x = size_y = 2048
overlap_x = overlap_y = 102
updated_tile_size = size_x - overlap_x * 2
version = 'v1'

df_spots_new = pd.DataFrame()
for tile in tqdm(tiles):
    print(f'Running tile {tile}')
    df_spots_tile = pd.read_csv(f'tiles_segmentation_results/spots_decoded_all_stitched_globalcoords_{version}_tile_{tile}.txt', sep = '\t', header = 0, index_col = 0, keep_default_na = False)
    df_spots_tile.rename(columns={'spot_location_1':'x', 'spot_location_2':'y'}, inplace = True)
    df_spots_tile['clustermap'] = [f'tile_{tile}_cell_{i+1}' if i!=-1 else 'unsegmented' for i in df_spots_tile.clustermap.values]
    print(f'For tile {tile}, xmin - {df_spots_tile.x.min()}, xmax - {df_spots_tile.x.max()}, ymin - {df_spots_tile.y.min()}, ymax - {df_spots_tile.y.max()}')  
    row = row_dict[tile+1]
    col = col_dict[tile+1]

    df_spots_tile = df_spots_tile[(df_spots_tile.y >= overlap_y) & (df_spots_tile.y <= tile_size - overlap_y)]
    df_spots_tile = df_spots_tile[(df_spots_tile.x >= overlap_x) & (df_spots_tile.x <= tile_size - overlap_x)]
    df_spots_tile.x = df_spots_tile.x.values - overlap_x
    df_spots_tile.y = df_spots_tile.y.values - overlap_y
    df_spots_tile.y = updated_tile_size - df_spots_tile.y

    df_spots_tile.y = df_spots_tile.y + row * updated_tile_size
    df_spots_tile.x = df_spots_tile.x + col * updated_tile_size

    df_spots_new = pd.concat([df_spots_new, df_spots_tile], axis = 0)

df_spots_new['index'] = np.arange(len(df_spots_new))
df_spots_new = df_spots_new.set_index('index')

df_spots_new['cell'] = df_spots_new.clustermap.values
df_spots_new.to_csv(f'spots_decoded_all_stitched_globalcoords_clustermap_{version}.txt', sep = '\t')