import os 
import warnings
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
import tifffile as tiff
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
if not os.path.exists('tiles'):
    os.mkdir('tiles')

### 0. parameters (tile location)
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

### 0. parameters (others)
tiles = np.arange(42) # 0~41
num_z = 89
num_rows = 6
num_cols = 10
tile_size = size_x = size_y = 2048
overlap_x = overlap_y = 102
updated_tile_size = size_x - overlap_x * 2
version = 'v1'

keys = ['spot_id','x','y','z','features','target','radius','intensity', 'xc','yc','zc','x_min','x_max','y_min','y_max','z_min','z_max'] # the information to extract from Starfish result
keys_xy = ['x','y','z','xc','yc','zc','x_min','x_max','y_min','y_max','z_min','z_max']
root_folder = '/home/jinmr2/20230925_iss_hplacenta_jd34/data/tile_decoded/'
df_spots_all = pd.DataFrame()

### 1. get individual tiles 
for tile in tqdm(tiles):
    tile_folder = f'{root_folder}tile_{tile}'
    try:
        with open(f'{tile_folder}_decoded_v2.pkl', 'rb') as f:
            decoded = pickle.load(f)
    except:
        print(f'No decoded.pkl in tile_{tile}')
        continue
    for k in keys:
        if hasattr(decoded, k):
            values = decoded[k].values
        else:
             continue
        if k == 'spot_id':
            df_spots = pd.DataFrame(index = np.arange(len(values)))
            df_spots[k] = [f'tile{tile}_{v}' for v in values]
        else:            
            df_spots[k] = values
            
    df_spots['tile'] = tile
    df_spots.to_csv(f'tiles/tile_{tile}.txt', sep = '\t')
    df_spots_all = pd.concat([df_spots_all, df_spots], axis = 0)

### 2. stitch spots across tiles
df_spots_new = pd.DataFrame()
for tile in tqdm(tiles):
    print(f'Running tile {tile}')
    df_spots_tile = pd.read_csv(f'tiles/tile_{tile}.txt', sep = '\t', header = 0, index_col = 0, keep_default_na = False)
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

df_spots_new.to_csv(f'spots_decoded_all_stitched_globalcoords_{version}.txt', sep = '\t')
df_spots_new = pd.read_csv(f'spots_decoded_all_stitched_globalcoords_{version}.txt', sep = '\t', header = 0, index_col = 0, keep_default_na = False)

### 3. plot stitched spots   
x, y = df_spots_new.x.values, df_spots_new.y.values
fig, ax = plt.subplots(figsize = (num_cols*2,num_rows*2))
plt.scatter(x, y, s = 0.0005)
plt.savefig(f'spots_decoded_all_manual_stitched_{version}.png', dpi = 300)

df_spots_new = df_spots_new[df_spots_new.target != ''].copy()
x, y = df_spots_new.x.values, df_spots_new.y.values
fig, ax = plt.subplots(figsize = (num_cols*2,num_rows*2))
plt.scatter(x, y, s = 0.0005)
plt.savefig(f'spots_decoded_all_manual_stitched_validGenes_{version}.png', dpi = 300)

### 4. stitch dapi
img_whole = tiff.imread('/home/jinmr2/20230925_iss_hplacenta_jd34/data/20230918_1Kgenes_hPlacenta_JS34_R1_RAW_ch00.tif')
img_merged = np.zeros(
    (num_rows * updated_tile_size, num_cols * updated_tile_size),
    dtype = np.uint8
)

### 5. plot stitched dapi
for tile in tqdm(tiles):
    img_tile = img_whole[num_z * tile : num_z * (tile + 1)].max(axis = 0)    
    row, col = row_dict[tile+1], col_dict[tile+1]

    img_tile = img_tile[overlap_y:(size_y - overlap_y), :]
    img_tile = img_tile[:, overlap_x:(size_x - overlap_x)]
    img_tile = img_tile[::-1,:]
    img_merged[(row * updated_tile_size):((row+1) * updated_tile_size), (col * updated_tile_size) : ((col+1) * updated_tile_size)] = img_tile

img_merged = img_merged[::-1, :]
tiff.imwrite(f'dapi_all_{version}.tif', img_merged)
print(img_merged.shape)
plt.figure(figsize = (22, 18))
plt.imshow(img_merged, cmap = 'gray')
plt.axis('off')
plt.savefig(f'dapi_all_{version}.png', dpi = 300, bbox_inches = 'tight')