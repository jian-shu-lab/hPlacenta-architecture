library(Seurat)
library(SeuratDisk)

setwd('/home/jinmr2/20230814_iss_hplac_r1-r6/data/single_cell_reference/data/')
Convert('raw.h5ad', dest = 'h5seurat', overwrite = T)
raw = LoadH5Seurat('raw.h5seurat', meta.data = F, misc = F)

df_anno = read.csv('ctype.txt', sep = '\t', header = 1, row.names = 1, na.strings = 'NA')
raw@meta.data = df_anno

# #### use overlapping genes
# df_genes = read.csv('/home/jinmr2/20230814_iss_hplac_r1-r6/data/tile_data_registered/tile_1/genes.csv',header=F)
# overlap_genes = intersect(df_genes$V1, rownames(raw))
# raw = raw[overlap_genes,]
# hplac.list <- SplitObject(raw, split.by = "Sample")
# for (i in 1:length(hplac.list)) {
#   hplac.list[[i]] <- NormalizeData(hplac.list[[i]], verbose = FALSE)
#   hplac.list[[i]] <- FindVariableFeatures(hplac.list[[i]], selection.method = "vst", nfeatures = 500,
#                                              verbose = FALSE)
# }
# hplac.anchors <- FindIntegrationAnchors(object.list = hplac.list, dims = 1:50)
# saveRDS(hplac.anchors,'anchors_integrate_overlap_v2.rds')
# hplac.integrated <- IntegrateData(anchorset = hplac.anchors, dims = 1:50)
# hplac.integrated <- ScaleData(hplac.integrated, verbose = FALSE)
# hplac.integrated <- RunPCA(hplac.integrated, npcs = 50, verbose = FALSE)
# hplac.integrated <- RunUMAP(hplac.integrated, reduction = "pca", dims = 1:50, verbose = FALSE)

# saveRDS(hplac.integrated,'integrate_overlap_v2.rds')

### use overlapping genes (rPCA)
df_genes = read.csv('/home/jinmr2/20230814_iss_hplac_r1-r6/data/tile_data_registered/tile_1/genes.csv',header=F)
overlap_genes = intersect(df_genes$V1, rownames(raw))
raw = raw[overlap_genes,]
hplac.list <- SplitObject(raw, split.by = "Sample")
for (i in 1:length(hplac.list)) {
  hplac.list[[i]] <- NormalizeData(hplac.list[[i]], verbose = FALSE)
  hplac.list[[i]] <- FindVariableFeatures(hplac.list[[i]], selection.method = "vst", nfeatures = 1000,
                                             verbose = FALSE)
}
features = SelectIntegrationFeatures(object.list = hplac.list)
hplac.list = lapply(hplac.list, function(x) {
  x = ScaleData(x, features = features, verbose = FALSE)
  x = RunPCA(x, features = features, npcs = 30, verbose = FALSE)
})
hplac.anchors <- FindIntegrationAnchors(
  object.list = hplac.list, dims = 1:30, anchor.features = features, reduction = "rpca"
)
saveRDS(hplac.anchors,'anchors_integrate_overlap_v3.rds')
hplac.integrated <- IntegrateData(anchorset = hplac.anchors, dims = 1:30)
hplac.integrated <- ScaleData(hplac.integrated, verbose = FALSE)
hplac.integrated <- RunPCA(hplac.integrated, npcs = 30, verbose = FALSE)
hplac.integrated <- RunUMAP(
  hplac.integrated, reduction = "pca", dims = 1:30, verbose = FALSE, return_model = T
)
saveRDS(hplac.integrated,'integrate_overlap_v3.rds')
