import os
import numpy as np
from tqdm import tqdm
import tifffile as tiff
import matplotlib as mpl
import matplotlib.pyplot as plt

# tile -> row
row_dict = {
    1: 5, 2: 5, 3: 5, 
    4: 4, 5: 4, 6: 4, 7: 4, 8: 4, 9: 4, 
    10: 3, 11: 3, 12: 3, 13: 3, 14: 3, 15: 3, 16: 3, 17: 3, 18: 3, 
    19: 2, 20: 2, 21: 2, 22: 2, 23: 2, 24: 2, 25: 2, 26: 2, 27: 2, 28: 2, 
    29: 1, 30: 1, 31: 1, 32: 1, 33: 1, 34: 1, 35: 1, 36: 1, 
    37: 0, 38: 0, 39: 0, 40: 0, 41: 0, 42: 0,
}

# tile -> column
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

dapi_source_file = '/home/jinmr2/20230925_iss_hplacenta_jd34/data/20230918_1Kgenes_hPlacenta_JS34_R1_RAW_ch00.tif'

### parameters
tiles = np.arange(42) # tiles
num_z = 89 # number of z slices
num_rows = 6 # number of rows
num_cols = 10 # number of columns
size_x = size_y = 2048 # size of each tile (pixel)
overlap_x = overlap_y = 102 # overlap between tiles (pixel)
updated_tile_size = size_x - overlap_x * 2 # updated tile size (pixel)


### read data and get individual dapi
img_whole = tiff.imread(dapi_source_file)
if not os.path.exists('dapi_images'):
    os.mkdir('dapi_images')
for tile in tiles:
    print(tile)
    img_tile = img_whole[num_z * tile : num_z * (tile + 1)].max(axis = 0)
    tiff.imwrite(f'dapi_images/dapi_tile_{tile}.tif', img_tile)
    plt.figure(figsize = (2*num_cols, 2*num_rows))
    plt.imshow(img_tile, cmap = 'gray')
    plt.axis('off')
    plt.savefig(f'dapi_images/dapi_tile_{tile}.png', dpi = 300, bbox_inches = 'tight')

### get stitched image for dapi
img_whole = tiff.imread(dapi_source_file)
img_merged = np.zeros(
    (num_rows * updated_tile_size, num_cols * updated_tile_size),
    dtype = np.uint8
)

### merge data tile by tile
for tile in tqdm(tiles):
    img_tile = img_whole[num_z * tile : num_z * (tile + 1)].max(axis = 0)    
    row, col = row_dict[tile+1], col_dict[tile+1]

    img_tile = img_tile[overlap_y:(size_y - overlap_y), :]
    img_tile = img_tile[:, overlap_x:(size_x - overlap_x)]
    img_tile = img_tile[::-1,:]
    img_merged[(row * updated_tile_size):((row+1) * updated_tile_size), (col * updated_tile_size) : ((col+1) * updated_tile_size)] = img_tile
    tiff.imwrite(f'dapi_images/dapi_tile_whole_{tile}_v3.tif', img_tile)

### save
img_merged = img_merged[::-1, :]
tiff.imwrite('dapi_all_v3.tif', img_merged)
print(img_merged.shape)
plt.figure(figsize = (2*num_cols, 2*num_rows))
plt.imshow(img_merged, cmap = 'gray')
plt.axis('off')
plt.savefig('dapi_all_v3.png', dpi = 300, bbox_inches = 'tight')