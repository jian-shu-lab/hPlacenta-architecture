# total 18 pairs of r1 and r2 (6 groups * 3 samples)
# first, find the shift of dapi between r1 and r2
# apply the shift to spots images (4 channels) in r2 to r1 (4 channels)
# the final output for each concat pair should be 512 * 512 * 8

import os
import cv2
import sys
import numpy as np
from tqdm import tqdm
import tifffile as tiff
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.ndimage import fourier_shift


if not os.path.exists('concat_2_rounds'):
    os.makedirs('concat_2_rounds')

def _get_fft_2d(img1, img2):
    if img2.shape != img1.shape:
        img2 = cv2.resize(img2, (img1.shape[1], img1.shape[0]))
    
    plt.imshow(img2, cmap = 'gray', vmax = 75)
    plt.axis('off')
    plt.savefig('resized_image.png', dpi = 300, bbox_inches = 'tight')
    print(f'shape of image 1 and image 2: {img1.shape}, {img2.shape}')

    # step 1
    fft_img1 = np.fft.fftn(img1)
    fft_img2 = np.fft.fftn(img2)

    # step 2
    cross_correlation = np.fft.ifftn(fft_img1 * np.conj(fft_img2))

    # step 3
    shifts = np.unravel_index(np.argmax(np.abs(cross_correlation)), cross_correlation.shape)
    print(f'shifts: {shifts}')
    translation = []
    for shift, dim in zip(shifts, img1.shape):
        if shift > dim // 2:
            shift -= dim
        translation.append(shift)
    print(f'translation: {translation}')

    return translation

def _apply_shift(img2, translation):
    # step 4
    aligned_img2 = fourier_shift(np.fft.fftn(img2), translation)
    aligned_img2 = np.fft.ifftn(aligned_img2).real.astype(img2.dtype)

    shift_x, shift_y = translation
    shift_x = int(shift_x); shift_y = int(shift_y)
    print(shift_x, shift_y)
    if shift_x > 0:
        aligned_img2[:shift_x, :] = 0
    elif shift_x < 0:
        aligned_img2[shift_x:, :] = 0
    if shift_y > 0:
        aligned_img2[:, :shift_y] = 0
    elif shift_y < 0:
        aligned_img2[:, shift_y:] = 0
    
    plt.imshow(aligned_img2, cmap = 'gray', vmax = 75)
    plt.axis('off')
    plt.savefig('registered_image.png', dpi = 300, bbox_inches = 'tight')

    return aligned_img2

if __name__ == '__main__':
    R1_folder = R2_folder = 'dapi_images/'
    
    samples = ['JS36', 'JS39', 'JS40']
    Gs = ['G1', 'G2', 'G3', 'G4', 'G5', 'G6']
    channels = ['ch01', 'ch02', 'ch03', 'ch04']
    
    for sample in samples:
        for G in Gs:
            print(f'processing {sample}_{G}')
            dapi_r1 = tiff.imread(f'{R1_folder}/{sample}_{G}_Round1_dapi_all_v3.tif')
            dapi_r2 = tiff.imread(f'{R2_folder}/{sample}_{G}_Round2_dapi_all_v3.tif')
            shift = _get_fft_2d(dapi_r1, dapi_r2)
            dapi_r2_aligned = _apply_shift(dapi_r2, shift)
            tiff.imsave(f'{R2_folder}{sample}_{G}_Round2_dapi_all_v3_aligned.tif', dapi_r2_aligned)

            for channel in tqdm(channels):
                img_r2 = tiff.imread(f'spots_images/{sample}_{G}_Round2_{channel}_all.tif')
                img_aligned = _apply_shift(img_r2, shift)
                tiff.imsave(f'spots_images/{sample}_{G}_Round2_{channel}_all_aligned.tif', img_aligned)
