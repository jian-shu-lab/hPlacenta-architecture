custom_lib_path <- "/home/jinmr2/R/x86_64-pc-linux-gnu-library/4.0"
.libPaths(custom_lib_path)

### load environment
library(Matrix)
library(anndata)
library(Seurat)
library(reticulate)
use_python('/home/kjin2/pyenv/run_starfish/bin/python3.8')
setwd('/home/jinmr2/sample_integration/three_samples/')

### load query data
data <- read_h5ad("raw_clustermap_all.h5ad")
saveRDS(data,'raw_clustermap_all.rds')
query <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
query <- subset(query,nCount_RNA >= 80)
query <- NormalizeData(query, verbose = FALSE)
query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = 1000,
                                        verbose = FALSE)

### load reference data
ref = readRDS('/home/jinmr2/20230814_iss_hplac_r1-r6/data/single_cell_reference/data/integrate_overlap_v3.rds')

### transfer embedding
transfer.anchors <- FindTransferAnchors(
  reference = ref, 
  query = query,
  dims = 1:30, 
  reference.reduction = "pca"
)
saveRDS(transfer.anchors, 'anchors_transfer_all.rds')
predictions <- TransferData(
  anchorset = transfer.anchors,
  anchors, 
  refdata = ref$Clusters,
  dims = 1:30
)

### label transfer
query <- AddMetaData(query, metadata = predictions)
saveRDS(query, 'transfer_clustermap_all_v1.rds')

### umap integration
ref <- RunUMAP(
  ref, 
  dims = 1:30, 
  reduction = "pca", 
  return.model = TRUE
)
query <- MapQuery(
  anchorset = transfer.anchors, 
  reference = ref, 
  query = query,
  refdata = list(celltype = "Clusters"), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)
saveRDS(query, paste0('transfer_clustermap_all_v2.rds'))
write.csv(query@meta.data, paste0('label_transferred_annotations_all.csv'))

# save umaps
write.csv(ref@reductions$umap@cell.embeddings, 'umaps_ref.csv')
query = readRDS('transfer_clustermap_all_v2.rds')
write.csv(query@reductions$ref.umap@cell.embeddings, 'umaps_query.csv')

### visualizations
# ref.umap
colnames(query@meta.data)
DimPlot(query, reduction = 'ref.umap', group.by = 'predicted.celltype',label=T)
# DimPlot(query, reduction = 'ref.umap', group.by = 'leiden',label=T)
DimPlot(ref, reduction = 'umap', group.by = 'Clusters',label = T)
