import os
import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

def _find_local_maxima(
    image,
    neighborhood_size = 5,
    threshold = 50,
):
    data_max = filters.maximum_filter(image, neighborhood_size)
    maxima = (image == data_max)
    data_min = filters.minimum_filter(image, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    x, y = [], []

    for dy,dx in slices:
        x_center = (dx.start + dx.stop - 1)/2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1)/2    
        y.append(y_center)

    return num_objects, x, y

def spots_calling(
    image_ch, 
    channel_name,
    local_maxima_neighbors = 10, 
    local_maxima_diff_threshold = 50,
):
    num_objects, x, y = _find_local_maxima(
        image_ch,
        neighborhood_size = local_maxima_neighbors,
        threshold = local_maxima_diff_threshold,
    )

    return x, y, num_objects

### params
samples = ['JS36', 'JS39', 'JS40']
Groups = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']
Rs = ['Round1', 'Round2']
channels = ['ch01', 'ch02', 'ch03', 'ch04']
colors = ['#4B4C76', '#687DA5', '#A8CED1', '#91B894', '#CED87F', '#E2CB6B', '#D1A289', '#9E5B52']
root_folder = '/home/jinmr2/ISH_20230928/spots_images/'
save_folder = 'figures_v3'
scale = 1/300
dapi_amplification = 20
# scale bar
pixel_size = 1.137 # um
bar_size = 500 # um
scale_factor = bar_size / pixel_size  # Adjust this value as needed
scale_x = 600    # X-coordinate for the scale bar

if not os.path.exists(save_folder):
    os.mkdir(save_folder)

df_genes = pd.read_excel('/home/jinmr2/ISH_hPLAC_8genes_6groups_Group-Round-Channel_(13).xlsx', sheet_name='Sheet3')
df_genes.set_index('Gene', inplace = True)

for sample in samples:
    for G in Groups:
        dapi = tiff.imread(f'/home/jinmr2/ISH_20230928/dapi_images/{sample}_{G}_Round1_dapi_all_v3.tif')
        fig = plt.figure(figsize = (dapi.shape[1] * scale, dapi.shape[0] * scale))
        dapi = dapi * dapi_amplification
        dapi[dapi > 255] = 255
        dapi = dapi.astype(np.uint8)
        plt.imshow(dapi, cmap = 'gray', vmax = 255, alpha = 1.0)

        cell_types = []
        for R in Rs:
            for channel in channels:
                if R == 'Round1':
                    img = tiff.imread(f'{root_folder}{sample}_{G}_{R}_{channel}_all.tif')
                elif R == 'Round2':
                    img = tiff.imread(f'{root_folder}{sample}_{G}_{R}_{channel}_all_aligned.tif')

                G2 = 'Group' + G.split('G')[1]
                channel2 = int(channel.split('ch0')[1])
                gene_name = df_genes[(df_genes['Group'] == G2) & (df_genes['Round'] == R) & (df_genes['Channel'] == channel2)].index[0]

                xo, yo = np.where(img>0)
                x, y, num_objects = spots_calling(
                    image_ch = img,
                    channel_name = f'{G}_{R}_{channel}_{gene_name}',
                    local_maxima_neighbors = 15, 
                    local_maxima_diff_threshold = 100,
                )
                print(f'{sample}_{G}_{R}_{channel}: {gene_name}, {len(x) / 1e6} million spots out of {len(xo) / 1e6} million foreground pixels')
                color = df_genes.loc[gene_name, 'Color']
                cell_type = df_genes.loc[gene_name, 'Cell type']
                cell_types.append(cell_type)
                label = f'{gene_name} ({cell_type})' if not pd.isna(cell_type) else gene_name
                plt.scatter(x, y, color = color, s = 1.5, label = label, alpha = 0.8)

        # sort them by legend labels
        handles, labels = plt.gca().get_legend_handles_labels()
        sorted_indices = np.argsort(cell_types)
        handles = [handles[i] for i in sorted_indices]
        labels = [labels[i] for i in sorted_indices]

        scale_y = np.max(y) - 200
        plt.plot([scale_x, scale_x + scale_factor], [scale_y, scale_y], color='white', linewidth=5)
        plt.annotate(f'{bar_size}µm', (scale_x + scale_factor / 2, scale_y - 75), color='white', ha='center', fontsize = 24)

        plt.legend(handles, labels, markerscale = 10, fontsize = 20)
        plt.axis('off')
        plt.savefig(f'{save_folder}/{sample}_{G}.png', dpi = 300, bbox_inches = 'tight')
