import os
import numpy as np
from tqdm import tqdm
import tifffile as tiff
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

if not os.path.exists('spots_images'):
    os.mkdir('spots_images')

def _get_tile_row_col(root, num_tiles):
    tiles = np.arange(num_tiles)
    x = []; y = []
    for tile in tiles:
        x.append(float(root.find(f".//Tile[@FieldX='{tile}']").get('PosX')))
        y.append(float(root.find(f".//Tile[@FieldX='{tile}']").get('PosY')))

    xpos = np.sort(np.unique(x))
    xpos_dict = dict(zip(xpos, np.arange(len(xpos))))

    ypos = np.sort(np.unique(y))
    ypos_dict = dict(zip(ypos, np.arange(len(ypos))))

    row_dict = {}
    col_dict = {}
    for tile in tiles:
        row_dict[tile+1] = ypos_dict[y[tile]]
        col_dict[tile+1] = xpos_dict[x[tile]]
    return row_dict, col_dict

def _stitch_spots_images(sample, batch, r, ch, vmax = 30, overlap_pixels = 27):
    root_folder = '/home/jinmr2/ISH_20230928/20230928_ish_js36-js39-js40_r1-r2/'
    if r == 'Round1':
        root_folder = root_folder + '20230922 hPLAC_ISH_JS36_JS39_JS40_R1_OVERVIEW_20X/'
    elif r == 'Round2':
        root_folder = root_folder + '20230920 hPLAC_ISH_JS36_JS39_JS40_R2_OVERVIEW_20X/'
    root = ET.parse(f'{root_folder}{sample}_{batch}_{r}/MetaData/{sample}_{batch}_{r}_Properties.xml').getroot()
    num_tiles = len(root.findall(".//Tile"))

    spots_source_file = f'{root_folder}{sample}_{batch}_{r}/{sample}_{batch}_{r}_RAW_{ch}.tif'
    img_whole = tiff.imread(spots_source_file)

    row_dict, col_dict = _get_tile_row_col(root, num_tiles=num_tiles)

    ### parameters
    tiles = np.arange(num_tiles)
    num_z = int(img_whole.shape[0] / num_tiles)
    num_rows = max(row_dict.values()) + 1
    num_cols = max(col_dict.values()) + 1
    size_x = size_y = 512
    overlap_x = overlap_y = overlap_pixels
    updated_tile_size = size_x - overlap_x * 2

    ### get stitched image for spots_images
    img_merged = np.zeros(
        (num_rows * updated_tile_size, num_cols * updated_tile_size),
        dtype = np.uint8
    )

    # append data tile by tile
    for tile in tqdm(tiles):
        img_tile = img_whole[num_z * tile : num_z * (tile + 1)].max(axis = 0)    
        row, col = row_dict[tile+1], col_dict[tile+1]

        img_tile = img_tile[overlap_y:(size_y - overlap_y), :]
        img_tile = img_tile[:, overlap_x:(size_x - overlap_x)]
        img_merged[(row * updated_tile_size):((row+1) * updated_tile_size), (col * updated_tile_size) : ((col+1) * updated_tile_size)] = img_tile

    img_merged = img_merged[::-1, :]
    tiff.imwrite(f'spots_images/{sample}_{batch}_{r}_{ch}_all.tif', img_merged)

    plt.figure(figsize = (2*num_cols, 2*num_rows))
    plt.imshow(img_merged, cmap = 'gray', vmax = vmax)
    plt.axis('off')
    plt.savefig(f'spots_images/{sample}_{batch}_{r}_{ch}_all.png', dpi = 300, bbox_inches = 'tight')

if __name__ == '__main__':
    samples = ['JS36', 'JS39', 'JS40']
    batches = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']
    rs = ['Round1', 'Round2']
    channels = ['ch01', 'ch02', 'ch03', 'ch04']

    for sample in samples:
        for batch in batches:
            for r in rs:
                for ch in channels:
                    print(f'{sample}_{batch}_{r}_{ch}')
                    _stitch_spots_images(sample, batch, r, ch, vmax = 30, overlap_pixels = 27)