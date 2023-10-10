These code snippets were used to preprocess In-situ hybridization (ISH) data, call local maxima signals and draw multiplexed images.

### Stitch tiles
- Stitch tiles from raw ISH data into concatenate images for both DAPI and gene staining signals `stitch_nuclei_images.py` and `stitch_gene_channels.py`

### Align images across rounds
- Images for same sample across different rounds are aligned fast fourier transform (FFT) `round_alignment.py`

### Local maxima detection and visualization
- Signal detection was conducted by finding local maxima in the image `visualize_multiplexed_ISH.py`