The imputation strategy refers to [Shi et al.](https://www.nature.com/articles/s41586-023-06569-5). This code was inspired by their method and [scripts](https://github.com/wanglab-broad/mCNS-atlas/tree/main/Imputation), and made some modification based on our own integration approach.

### Intermediate mapping
- Specifically, for each of all 1,002 genes in the STARmap data, we performed an intermediate mapping to align each STARmap cell with the nearest neighbors in reference single-cell dataset. This step is used to find the best neighborhood parameter. `intermediate_mapping.py`

### Final imputation
- For every cell, we calculated each gene’s imputed expression level as the weighted average of the gene’s expression across the associated set of scRNA-seq atlas cells, where weights were proportional to the reciprocol of distances: `final_mapping.py`

### Visualization
- Visualize imputation results and observed STARmap data for user-defined genes side-by-side: `visualize_imputation_results.py`