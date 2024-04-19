These code snippets were used to preprocess STARmap data, segment cells, annotate cell types and visualize results. Detailed steps include:

### Image registration
- After image deconvolution (see method section in paper), images for same tile across different rounds are aligned to the first round 3D fast fourier transform (FFT) `image_registration.py`

### Prepare Starfish input
- We used [Starfish](https://spacetx-starfish.readthedocs.io/en/latest/) package to process registered STARmap images. It requires specific format as the input for spot calling and gene decoding. We prepared the input based on their tutorial: `create_starfish_format_data.py`

### Spot calling and gene decoding
- After preparation of Starfish input, we called valid spots using BlobDetector method, and then decoded spots with PerRoundMaxChannel method. Other steps, including noise removal and signal normalization, can also be used based on the data quality: `run_starfish_pipeline.py`

### Stitch tiles
- Stitch tiles from processed data into concatenate images for both DAPI and gene staining signals `stitch_nuclei_image.py` and `stitch_decoded_spots.py`

### Cell segmentation
- We ran [ClusterMap](https://github.com/wanglab-broad/ClusterMap) to segment cells from decoded transcripts tile by tile. We didn't use DAPI staining since it may generate too conservative segmentation masks: `cell_segmentation.py`

### Cell annotation from reference single cell data
- We used the single-cell multiomics data we generated in this study to annotate STARmap segmented cells. We utilized Seurat v4 and transferred labels (embeddings) to the STARmap data: `integration_label_transfer`

### Cell-cell interaction inferferce
- We inferred cell-cell interactions using CellChat with imputed gene expression profiles. More details can be found in: `infer_cell_interaction.R`

### Visualization
- Results were visualized in various ways, which can be found here: `visualization`

