## Code for processing ATAC and RNA Slide-tags placenta data

############ 0. Read in data and libraries ############
library(Seurat)
library(qs)
library(tidyverse)
library(reticulate)
library(harmony)
library(ArchR)
library(tidyverse)
library(Signac)
library(scatterpie)
library(spdep)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ArchRtoSignac)
library(ggplot2)
library(harmony)
library(parallel)
library(future)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

set.seed(42)
addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores()/2))

# structure of this analysis
1. Process RNA data
2. Process ATAC data
3. WNN for combination of RNA and ATAC
4. CellChat

############ 1. Process RNA data ############
# process single-cell reference
sc <- import('scanpy')
adata <- sc$read_h5ad('../data/sc_hplacenta_gene_matrix.h5ad')
expr_matrix <- Matrix::t(adata$X)
colnames(expr_matrix) <- rownames(adata$obs)
rownames(expr_matrix) <- rownames(adata$var)
meta <- adata$obs
ref <- CreateSeuratObject(counts = expr_matrix, meta.data = meta)
ref <- NormalizeData(ref) %>% FindVariableFeatures() 
ref@assays$RNA@data <- ref@assays$RNA@counts
ref <- ScaleData(ref)
ref <- RunPCA(ref)
ref <- FindNeighbors(ref, dims = 1:20)
ref <- FindClusters(ref, resolution = 1.5)
ref <- RunUMAP(ref, reduction = "pca", dims = 1:30)
saveRDS(ref, '../data/sc_reference.rds')
ref <- readRDS('../data/sc_reference.rds')

## specify runs to analyse:
run <- c("JS35", "JS36_3", "JS40")
run_name <- c("JS35", "JS36", "JS40")
runs_list_full <- c("JS35", "JS36_3", "JS40")
runs_list_full_rnaonly <- c("JS35", "JS36_3", "JS40")
run_week <- c()
## Mouse OR Human:
species <- "Human" 

# integrate and process slide-tag RNA
list_seurats <- list()
for(i in seq_along(run)){
  cellranger_output <- paste0("../data/matrices/", run[i], "RNA_only/raw_feature_bc_matrix")
  data <- Read10X(cellranger_output,
                  gene.column = 2,
                  cell.column = 1,
                  unique.features = TRUE,
                  strip.suffix = TRUE
                  )
  list_seurats[[i]] <- CreateSeuratObject(counts = data)
  # only take forward mapped cells for each run
    # read bc files
  bc_multi <- read.table(paste0("../data/matrices/", runs_list_full[i], "/MULTI/filtered_feature_bc_matrix/barcodes.tsv.gz"))
  bc_RNA_only <- read.table(paste0("../data/matrices/", runs_list_full_rnaonly[i], "RNA_only/filtered_feature_bc_matrix/barcodes.tsv.gz"))
  # calculate intersect
  mapped_cells <- sub("-1", "", union(bc_multi$V1, bc_RNA_only$V1))
  list_seurats[[i]] <-subset(list_seurats[[i]], cells = mapped_cells)
  # add metadata
  list_seurats[[i]] <- RenameCells(list_seurats[[i]], add.cell.id = run_name[i])
  list_seurats[[i]]$group <- run_name[i]
}

int <- Reduce(function(x, y) merge(x, y), list(list_seurats[[1]], list_seurats[[2]], list_seurats[[3]]))
int <- SCTransform(int)
int <- NormalizeData(int) %>% FindVariableFeatures() %>% ScaleData()
int <- RunPCA(int)
int <- RunUMAP(int, reduction = "pca", dims = 1:30)

## add mitochondrial %
## add the proportion gene expression that comes from ribosomal proteins 
## Percentage hemoglobin genes - includes all genes starting with HB except HBP.
## Percentage platlet
int[["percent.mt"]] <- PercentageFeatureSet(int, pattern = "^MT-")
int[["percent.ribo"]] <- PercentageFeatureSet(int, "^RP[SL]")
int[["percent.hb"]] <- PercentageFeatureSet(int, "^HB[^(P)]")
int[["percent.plat"]] <- PercentageFeatureSet(int, "PECAM1|PF4")

# use harmony to remove batch effects
int <- RunHarmony(int, "group")
int <- RunUMAP(int, reduction = "harmony", dims = 1:20)
int <- FindNeighbors(int, reduction = "harmony", dims = 1:20)
int <- FindClusters(int, resolution = 1.5)
p1 <- DimPlot(int, split.by = "group", label = T) + coord_fixed()
print(p1)
p2 <- DimPlot(int, label = T) + coord_fixed()
print(p2)
ggsave(paste0("../output/int_GEX_harmony_clusters_facet.pdf"),plot=p1,width=12,height=6)
ggsave(paste0("../output/int_GEX_harmony_clusters.pdf"),plot=p2,width=8,height=8)

# predict clusters based on single-cell reference
anchors <- FindTransferAnchors(ref, int, dims = 1:30, reference.reduction = "pca", reference.assay = "RNA", query.assay = "RNA")
#anchors <- FindTransferAnchors(ref, int, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchors, refdata = ref$Clusters, dims = 1:30)
int <- AddMetaData(int, metadata = predictions[,'predicted.id', drop=F], col.name = 'predicted.id_int')
#int <- FindNeighbors(int, dims = 1:20)
#int <- FindClusters(int, resolution = 1)
p1 <- DimPlot(int, label = T, group.by = 'predicted.id_int')
print(p1)
ggsave(paste0("../output/int_GEX_predictedid_UMAP.pdf"),plot=p1,width=8,height=8)

# plot expression of marker genes
marker_list <- list(
  vCTB = c('PEG10', 'PAGE4', 'CDH1'),
  vCTBp = c('MKI67', 'TOP2A', 'PCNA'),
  vCTB1 = c('ARL15', 'TBL1X', 'SIK3'),
  vCTB2 = c('IFI6', 'EFEMP1', 'SMAGP'),
  vCTB3 = c('MSI2', 'LRP5', 'AHCYL2'),
  EVT = c('ITGA5', 'HLA-G', 'DIO2'),
  EVT1 = c('EGFR', 'FBN2', 'UTRN'),
  EVT2 = c('CCNE1', 'HAPLN3', 'LY6E'),
  EVT3 = c('AOC1', 'PAPP2', 'PRG'),
  STB = c('ERVFRD-1', 'CYP19A1', 'IGSF5'),
  Endo = c('PECAM1', 'EGFL7', 'KDR'),
  Macrophages = c('SPP1', 'CD14', 'CD83'),
  Maternal_MAC = c('CD74', 'HLA-DRA', 'LYZ'),
  HBC = c('F13A1', 'LYVE1', 'ADAMTS17'),
  Myeloid_unk = c('FGF13', 'ALPK3'),
  FIB = c('COL3A1', 'SOX5', 'LAMA2'),
  FIB1 = c('PDGFRB', 'PLCB1', 'AGTR1'),
  FIB2 = c('PDGFRA', 'COL15A1', 'CXCL14'),
  Unknown1 = c('HGF', 'DCN', 'CDH11'),
  Unknown2 = c('FAM155A', 'ALDH1A2', 'RORB'),
  Unknown3 = c('RSPO3', 'PODN', 'BIN1')
)

DefaultAssay(int) <- "SCT"
int <- AddModuleScore(int, features = marker_list, name = "_score")
names(int@meta.data)[grep("_score", names(int@meta.data))] <- names(marker_list)

p1 <- FeaturePlot(int, features = c("vCTB1", "vCTBp", "vCTB2", "vCTB3"), order = T, pt.size = .1) + coord_fixed()
p2 <- FeaturePlot(int, features = c("EVT", "EVT1", "EVT2", "EVT3"), order = T, pt.size = .1) + coord_fixed()
p3 <- FeaturePlot(int, features = c("Macrophages", "Maternal_MAC", "Myeloid_unk", "HBC"), order = T, pt.size = .1) + coord_fixed()
print(p1)
print(p2)
print(p3)
ggsave(paste0("../output/int_GEX_markers1_UMAP.pdf"),plot=p1,width=8,height=8)
ggsave(paste0("../output/int_GEX_markers2_UMAP.pdf"),plot=p2,width=8,height=8)
ggsave(paste0("../output/int_GEX_markers3_UMAP.pdf"),plot=p3,width=8,height=8)

p1 <- VlnPlot(int, features = c("percent.mt", "percent.ribo", "percent.hb", "percent.plat"), pt.size = .1)
p2 <- FeaturePlot(int, features = c("percent.mt", "percent.ribo", "percent.hb", "percent.plat"), order = T, pt.size = .1) + coord_fixed()
p3 <- DimPlot(int, label = T) + coord_fixed()
p4 <- DotPlot(int, features = marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)
print(p2)
print(p3)
print(p4)
ggsave(paste0("../output/int_GEX_qcmetrics_violin.pdf"),plot=p1,width=10,height=10)
ggsave(paste0("../output/int_GEX_qcmetrics_UMAP.pdf"),plot=p2,width=8,height=8)
ggsave(paste0("../output/int_GEX_clusters_UMAP.pdf"),plot=p3,width=8,height=8)
ggsave(paste0("../output/int_GEX_markers_dotplot.pdf"),plot=p4,width=17,height=10)

# violin plots for markers
p1 <- VlnPlot(int, features = unlist(marker_list))
print(p1)
ggsave(paste0("../output/int_GEX_markers_violin.pdf"),plot=p1,width=15,height=40)

# create a heatmap plot for proportions of preditcted ids for each seurat cluster
proportions <- table(int@meta.data$seurat_clusters, int@meta.data$predicted.id_int)
proportions <- prop.table(proportions, margin = 1)

# Convert table to data frame
proportions <- as.data.frame.matrix(proportions)
proportions$seurat_clusters <- rownames(proportions)

# Reshape data frame to long format for plotting
proportions_long <- gather(proportions, key = "predicted_id_int", value = "proportion", -seurat_clusters)

# Create the heatmap plot
p1 <- ggplot(proportions_long, aes(x = predicted_id_int, y = seurat_clusters, fill = proportion)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Predicted ID Int", y = "Seurat Clusters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave(paste0("../output/int_GEX_clusterbypredictedid_heatmap.pdf"),plot=p1,width=10,height=10)

# Take spatiall mapped cells

mapped_cells <- c()
for(i in seq_along(run)){
  df_coords <- read.delim(paste0("../output/", run[i], "/coords_", run[i], ".csv"), sep = " ")
  df_coords <- df_coords[which(df_coords$clusters == 1), ]
  mapped_cells <- c(mapped_cells, paste0(run_name[i], "_", df_coords$cb))
}

int_mapped <- subset(int, cells = mapped_cells)
int_mapped

df_difference <- data.frame(all = table(int@meta.data$seurat_clusters), mapped = table(int_mapped@meta.data$seurat_clusters))
rownames(df_difference) <- df_difference$mapped.Var1
df_difference <- df_difference[,-1]
colnames(df_difference) <- c("all", "cluster", "mapped")
df_difference$pc_mapped <- (df_difference$mapped/df_difference$all)*100
df_difference$pc_lost <- (1-(df_difference$mapped/df_difference$all))*100

p1 <- ggplot(df_difference, aes(x = "x", y = cluster, fill = pc_lost)) +
  geom_tile() +
  labs(x = "", y = "Seurat Clusters") +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  geom_text(aes(label = round(pc_lost, 2)), color = "black", size = 3, vjust = 0.5)

print(p1)
ggsave(paste0("../output/int_GEX_pcmappedcellsbycluster_heatmap.pdf"),plot=p1,width=8,height=8)


proportions <- table(int_mapped@meta.data$seurat_clusters, int_mapped@meta.data$predicted.id_int)
proportions <- prop.table(proportions, margin = 1)

# Convert table to data frame
proportions <- as.data.frame.matrix(proportions)
proportions$seurat_clusters <- rownames(proportions)

# Reshape data frame to long format for plotting
proportions_long <- gather(proportions, key = "predicted_id_int", value = "proportion", -seurat_clusters)

# Create the heatmap plot
p1 <- ggplot(proportions_long, aes(x = predicted_id_int, y = seurat_clusters, fill = proportion)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Predicted ID Int", y = "Seurat Clusters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave(paste0("../output/int_GEX_clusterbypredictedid_mappedonly_heatmap.pdf"),plot=p1,width=10,height=10)

# Find markers for each cluster
markers_all <- FindAllMarkers(int)
write.table(markers_all, "../output/int_all_markers.txt", sep = "\t", row.names = TRUE)

top_markers <- markers_all %>%
  group_by(cluster) %>%
  top_n(50, wt = avg_log2FC) %>%
  ungroup()

write.table(top_markers, "../output/int_top50_markers.txt", sep = "\t", row.names = TRUE)

# remove low-quality clusters and assign labels
int_mapped <- int_mapped[,!int_mapped$seurat_clusters %in% c(0,5,7,8)]
int_mapped$Clusters.2 <- as.character(int_mapped$seurat_clusters) %>% {
  #.[. %in% c('')] <- 'NK/T cells'
  .[. %in% c('3', '9', '19', '23')] <- 'STB'
  .[. %in% c('1', '17', '18')] <- 'FIB2'
  .[. %in% c('2', '6', '12', '14')] <- 'vCTB'
  .[. %in% c('13')] <- 'EVT'
  .[. %in% c('11')] <- 'EVT-progenitor'
  #.[. %in% c('')] <- 'vCTB-proliferative'
  .[. %in% c('10', '20')] <- 'FIB1'
  .[. %in% c('22')] <- 'Erythroblasts'
  .[. %in% c('15')] <- 'Endothelial'
  .[. %in% c('21')] <- 'STB-progenitor'
  .[. %in% c('4', '16')] <- 'Hofbauer cells'
  ;.}

int_mapped$stage <- as.character(int_mapped$group) %>% {
  .[. %in% c('JS35')] <- '9w_1d'
  .[. %in% c('JS36')] <- '11w_1d'
  .[. %in% c('JS40')] <- '8w_2d'
  ;.}
int_mapped$stage <- factor(int_mapped$stage, levels = c("8w_2d", "9w_1d", "11w_1d"))

# Add spatial coordinates
df_coords_all <- data.frame()
for(i in seq_along(run)){
  df_coords <- read.table(paste0("../output/", run[i], "/coords_", run[i], ".csv"), header = T)
  df_coords <- df_coords[which(df_coords$clusters == 1), ]
  rownames(df_coords) <- paste0(run_name[i], "_", df_coords$cb)
  df_coords <- df_coords[,c("clusters", "x", "y", "UMI_cluster", "UMI_noise", "cb")]
  df_coords$SB_UMI_total <- df_coords$UMI_cluster + df_coords$UMI_noise
  df_coords$SNR <- df_coords$UMI_cluster / (df_coords$UMI_noise + df_coords$UMI_cluster)
  df_coords_all <- rbind(df_coords_all, df_coords)
}

int_mapped <- AddMetaData(int_mapped, metadata = df_coords_all[,!names(df_coords_all) %in% "cb", drop=F], )

cells_mito_pass <- rownames(int_mapped@meta.data[which(int_mapped@meta.data$percent.mt < 50), ])

paste0("number of cells filtered out due to mito >50 % = ", length(colnames(int_mapped)) - length(cells_mito_pass))

int_mapped <- subset(int_mapped, cells = cells_mito_pass)

saveRDS(int, file = '../data/Integrated_RNA_3samples_unfiltered.rds')
saveRDS(int_mapped, file = '../data/Integrated_RNA_3samples.rds')

## plot a dotplot to explain how markers were chosen
int_toplotmarkers <- int
int_toplotmarkers$Clusters.2 <- as.character(int_toplotmarkers$seurat_clusters) %>% {
  #.[. %in% c('')] <- 'NK/T cells'
  .[. %in% c('3', '9', '19', '23')] <- 'STB'
  .[. %in% c('1', '17', '18')] <- 'FIB2'
  .[. %in% c('2', '6', '12', '14')] <- 'vCTB'
  .[. %in% c('13')] <- 'EVT'
  .[. %in% c('11')] <- 'EVT-progenitor'
  .[. %in% c('0', '5', '7', '8')] <- 'Unassigned'
  #.[. %in% c('')] <- 'vCTB-proliferative'
  .[. %in% c('10', '20')] <- 'FIB1'
  .[. %in% c('22')] <- 'Erythroblasts'
  .[. %in% c('15')] <- 'Endothelial'
  .[. %in% c('21')] <- 'STB-progenitor'
  .[. %in% c('4', '16')] <- 'Hofbauer cells'
  ;.}

int_toplotmarkers@meta.data$Clusters.assignment <- paste(int_toplotmarkers$seurat_clusters, paste0("(",int_toplotmarkers$Clusters.2,")"), sep = " - ")

## new markers
marker_list2 <- list(
  vCTB = c("PEG10", "PAGE4", "CDH1"), 
  'STB-progenitor' = c("OVOL1"), 
  STB = c("ERVFRD-1", "CYP19A1", "IGSF5"), 
  'EVT-progenitor' = c(), 
  #'EVT-progenitor' = c("EGFR", "UTRN", "ITGA2"), 
  EVT = c("HAPLN3", "LY6E", "CCNE1"), 
  FIB1 = c("PDGFRB", "PLCB1", "AGTR1"), 
  FIB2 = c("PDGFRA", "COL15A1", "CXCL14"), 
  'Hofbauer cells' = c("F13A1", "LYVE1", "ADAMTS17"), 
  Endothelial = c("PECAM1", "EGFL7", "KDR"), 
  Erythroblasts = c("HBA1", "HBA2", "HBG1")
)

# Get the cluster names from the marker list
cluster_order <- names(marker_list2)

# Extract the cluster assignments from the metadata
cluster_assignments <- unique(int_toplotmarkers@meta.data$Clusters.assignment)

# Remove the cluster number and hyphen from the cluster assignments
cluster_names <- gsub("^[0-9]+", "", cluster_assignments)
cluster_names <- gsub(" - ", "", cluster_names)
# Remove the parentheses from the cluster names
cluster_names <- gsub("\\((.*)\\)", "\\1", cluster_names)
cluster_names <- trimws(cluster_names)

new_order <- rev(order(match(cluster_names, cluster_order)))
#cluster_names <- cluster_names[new_order]

# Order the cluster names to match the marker list
int_toplotmarkers@meta.data$Clusters.assignment <- factor(int_toplotmarkers@meta.data$Clusters.assignment, levels = cluster_assignments[new_order])

celltypes <- unique(int_toplotmarkers@meta.data$Clusters.2)
int_toplotmarkers@meta.data$Clusters.2 <- factor(int_toplotmarkers@meta.data$Clusters.2, levels = celltypes[rev(order(match(celltypes, cluster_order)))])

p <- DotPlot(int_toplotmarkers, features = marker_list2, group.by = 'Clusters.assignment') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.x = element_text(angle = 90, hjust = 0, size = 8)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
print(p)
ggsave(paste0("../output/int_GEX_markers_celltypeassignment_merge_dotplot.pdf"),plot=p,width=17,height=10)

p2 <- DotPlot(int_toplotmarkers, features = marker_list2, group.by = 'Clusters.2') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.x = element_text(angle = 90, hjust = 0, size = 8)) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
print(p2)
ggsave(paste0("../output/int_GEX_markers_celltypeassignment_celltype_dotplot.pdf"),plot=p2,width=17,height=10)

# print UMAP of seurat clusters
pdf('../output/figs/RNA/RNA_UMAP_seurat_clusters.pdf', width = 7.5, height = 6)
DimPlot(int, label = T, group.by = 'seurat_clusters') +
    labs(title = "UMAP of RNA", x = "UMAP1", y = "UMAP2", color = "Cluster") +
  coord_fixed() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(), 
          legend.key.size = unit(.5, "cm"),
          panel.border = element_rect(linetype = "solid", color = "black", size = 1.2),
          plot.title = element_text(hjust = 0.5, size = 18)) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    coord_fixed()
dev.off()

int_all <- readRDS("../data/Integrated_RNA_3samples_unfiltered.rds")
df_coords_all <- data.frame()
for(i in seq_along(run)){
  df_coords <- read.table(paste0("../output/", run[i], "/coords_", run[i], ".csv"), header = T)
  rownames(df_coords) <- paste0(run_name[i], "_", df_coords$cb)
  df_coords <- df_coords[,c("clusters", "x", "y", "UMI_cluster", "UMI_noise", "cb")]
  df_coords$SB_UMI_total <- df_coords$UMI_cluster + df_coords$UMI_noise
  df_coords$SNR <- df_coords$UMI_cluster / (df_coords$UMI_noise + df_coords$UMI_cluster)
  df_coords_all <- rbind(df_coords_all, df_coords)
}
int_all <- AddMetaData(int_all, df_coords_all)
table(int_all@meta.data$clusters, int_all@meta.data$seurat_clusters)

############ 2. Process ATAC data ############

# process combined atac_fragments_sorted.bed.gz 
# add index 1,2,3 for sample JS35,JS36,JS40 
inputFiles <- '../data/combined_atac_fragments_sorted.bed.gz'
names(inputFiles) <- 'combined'
ArrowFiles <- createArrowFiles(inputFiles, force = TRUE, addTileMat = TRUE,
                               addGeneScoreMat = TRUE, minFrags = 0, minTSS = 0,
                               maxFrags = 5e+07, threads = 4, outputNames = names(inputFiles))
proj <- ArchRProject(ArrowFiles)
saveArchRProject(ArchRProj = proj, outputDirectory = "../output/ArchR", load = FALSE)

# add label from RNA analysis
seu <- readRDS('../data/Integrated_RNA_3samples.rds')
meta <- data.frame(cells = colnames(seu), 
                   clusters = seu@meta.data$Clusters.2,
                   sample = seu$group)
proj$barcode = gsub(paste0('combined', "#"), "", proj$cellNames)
updated_barcodes <- sapply(proj$barcode, function(bc) {
  if (grepl("-1$", bc)) {
    return(paste0("JS35_", sub("-1$", "", bc)))
  } else if (grepl("-2$", bc)) {
    return(paste0("JS36_", sub("-2$", "", bc)))
  } else if (grepl("-3$", bc)) {
    return(paste0("JS40_", sub("-3$", "", bc)))
  } else {
    return(bc)
  }
})
names(updated_barcodes) <- NULL
proj$barcode <- updated_barcodes
proj <- proj[proj$barcode %in% meta$cells, ]
meta <- meta[match(proj$barcode, meta$cells), ]
proj$cluster <- meta$clusters
proj$group <- meta$sample
saveRDS(proj, file = '../data/Integrated_ATAC_Proj_unfiltered.rds')

# quality control filtering for fragments
proj_filter <- proj[proj@cellColData$nFrags > 500, ]

# Dimension Reduction and Clustering
proj_filter <- addIterativeLSI(proj_filter, useMatrix = "TileMatrix", name = "IterativeLSI", 
                        outlierQuantiles = c(0, 1), force = TRUE)
proj_filter <- addHarmony(proj_filter, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "group")
proj_filter <- addClusters(proj_filter, reducedDims = "Harmony", force = TRUE)
proj_filter <- addUMAP(proj_filter, reducedDims = "Harmony", force = TRUE, minDist = 0.2, nNeighbors = 10)
plotEmbedding(ArchRProj = proj_filter, colorBy = "cellColData", name = "cluster", embedding = "UMAP")
p <- plotEmbedding(proj_filter, colorBy = "GeneScoreMatrix", name = unlist(marker_list)[!(unlist(marker_list) %in% c("PAPP2","PRG"))], 
                   embedding = "UMAP", imputeWeights = getImputeWeights(proj_filter))

# motifs
pathToMacs2 = findMacs2()
proj_filter <- addGroupCoverages(proj_filter, groupBy = "cluster", sampleLabels = "group", minCells = 1, force = TRUE) 
proj_filter <- addReproduciblePeakSet(proj_filter, groupBy = "cluster", pathToMacs2 = pathToMacs2, sampleLabels = "group", force = T, threads = 1)
proj_filter <- addPeakMatrix(proj_filter)
proj_filter <- addMotifAnnotations(proj_filter, motifSet = "cisbp", name = "Motif", force = TRUE)
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_filter, 
  useMatrix = "PeakMatrix", 
  groupBy = "cluster",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
motifsUp <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj_filter,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
motifsDo <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj_filter,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
saveRDS(proj_filter, file = '../data/Integrated_ATAC_Proj.rds')
write.csv(motifsUp@assays@data@listData[["mlog10Padj"]], file = '../output/motifsup.csv')
write.csv(motifsDo@assays@data@listData[["mlog10Padj"]], file = '../output/motifsdown.csv')

# spatial visualization of motifs of coustom TFs
ctfs <- read.csv('../data/motifs.csv')[,1][1:63]
proj_filter <- readRDS('../data/Integrated_ATAC_Proj.rds')
if("Motif" %ni% names(proj_filter@peakAnnotation)){proj_filter <- addMotifAnnotations(proj_filter, motifSet = "cisbp", name = "Motif")}
proj_filter <- addBgdPeaks(proj_filter)
proj_filter <- addDeviationsMatrix(proj_filter, peakAnnotation = "Motif", force = TRUE)

# extract data for custom plots
# "MotifMatrix" | "GeneScoreMatrix"
Matrix <- getMatrixFromProject(proj_filter, useMatrix = "GeneScoreMatrix",
                               useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads())
# data <- as.data.frame(Matrix@assays@data$z)
data <- as.data.frame(Matrix@assays@data$GeneScoreMatrix);rownames(data) <- Matrix@elementMetadata$name
colnames(data) <- gsub('combined#','',colnames(data))
updated_barcodes <- sapply(colnames(data), function(bc) {
  if (grepl("-1$", bc)) {
    return(paste0("JS35_", sub("-1$", "", bc)))
  } else if (grepl("-2$", bc)) {
    return(paste0("JS36_", sub("-2$", "", bc)))
  } else if (grepl("-3$", bc)) {
    return(paste0("JS40_", sub("-3$", "", bc)))
  } else {
    return(bc)
  }
})
names(updated_barcodes) <- NULL
colnames(data) <- updated_barcodes
# rownames(data) <- sub("_\\d+$", "", rownames(data))
data <- as.data.frame(t(data[ctfs,]))
coord <- data.frame('x'=int_mapped$x, 'y'=int_mapped$y)[rownames(data),]
data <- cbind(coord, data)
saveRDS(data, file = '../data/int_GeneScoreMatrix.rds')


Matrix <- getMatrixFromProject(proj_filter, useMatrix = "MotifMatrix",
                               useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads())
data <- as.data.frame(Matrix@assays@data$z)
colnames(data) <- gsub('combined#','',colnames(data))
updated_barcodes <- sapply(colnames(data), function(bc) {
  if (grepl("-1$", bc)) {
    return(paste0("JS35_", sub("-1$", "", bc)))
  } else if (grepl("-2$", bc)) {
    return(paste0("JS36_", sub("-2$", "", bc)))
  } else if (grepl("-3$", bc)) {
    return(paste0("JS40_", sub("-3$", "", bc)))
  } else {
    return(bc)
  }
})
names(updated_barcodes) <- NULL
colnames(data) <- updated_barcodes
rownames(data) <- sub("_\\d+$", "", rownames(data))
data <- as.data.frame(t(data[ctfs,]))
coord <- data.frame('x'=int_mapped$x, 'y'=int_mapped$y)[rownames(data),]
data <- cbind(coord, data)
saveRDS(data, file = '../data/int_MotifMatrix_Z.rds')


############ Spatial autocorrelation

# plot clusters or expressions of genes spatially
# based on Seurat function
plot_cluster_feature_moran <- function(seu, reduction, label = FALSE, group_by = NULL,
                                 split_by = NULL, features = NULL, stat, ...) {
  if (!is.null(group_by)){
    col_length = length(unique(eval(parse(text=paste0("seu$", group_by)))))
    my_palette = rep(c(Seurat:::DiscretePalette(26, palette = "alphabet"),
                       RColorBrewer::brewer.pal(9, "Set1"),
                       RColorBrewer::brewer.pal(8, "Dark2"),
                       RColorBrewer::brewer.pal(12, "Paired")), 
                     length.out = col_length)
    plot <- DimPlot(seu, reduction = reduction, label = label, 
                    group.by = group_by, split.by = split_by, 
                    cols = my_palette, ...) + 
      scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey") +
      coord_fixed(ratio = 1) +
      theme(strip.text = element_text(size = 8, face = "bold"),
            line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      labs(title = past0(sample), subtitle = paste0("Moran's I = ", stat))
    if (!is.null(split_by)) { plot <- plot + NoLegend() }
  }
  if (!is.null(features)) {
    plot <- FeaturePlot(seu, reduction = reduction, order = TRUE, features = features, ...) 
    plot <- plot 
    #& scale_colour_gradientn(colours = inlmisc::GetColors(256, scheme = "sunset")[45:250])
    plot <- plot & theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_blank(),
                         axis.text = element_blank(),
                         axis.title = element_blank(),
                         axis.line = element_blank())
  }
  if (is.null(group_by) & is.null(features)) {
    message('Please provide clusters or features.')
  }
  if (reduction == 'spatial') {
    plot <- plot & labs(x = 'x_um', y = 'y_um')
    plot <- plot & theme(panel.border = element_rect(linetype = "solid", color = "black", size = 1.2))

  }
  return(plot)
}


# RNA
seu <- readRDS('../data/Integrated_RNA_3samples.rds')

seu_list <- SplitObject(seu, split.by = 'stage')
for (n in names(seu_list)) {
  print(n)
  seu <- seu_list[[n]]
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
  hvg <- VariableFeatures(seu)
  coords <- data.frame(seu$x, seu$y)
  coords_nb <- knn2nb(knearneigh(coords, k = 5))
  gene_expression <- as.data.frame(as.matrix((seu[hvg, ]@assays$SCT@counts)))
  gene_expression <- gene_expression[rowSums(gene_expression) > 0, ]
  
  res_list <- list()
  for (gene in rownames(gene_expression)) {
    gene_expression_vector <- as.numeric(gene_expression[gene, ])
    listw <- nb2listw(coords_nb)
    moran_result <- moran.test(gene_expression_vector , listw)
    res_list[[gene]] <- moran_result
  }
  
  original_p_values <- sapply(res_list, function(x) x$p.value)
  adjusted_p_values <- p.adjust(original_p_values, method = "fdr")
  names(adjusted_p_values) <- names(res_list)
  significant_adjusted <- adjusted_p_values[adjusted_p_values < 0.05]
  
  moran_values <- as.data.frame(t(as.data.frame(sapply(res_list[names(significant_adjusted)], function(x) x$estimate))))
  moran_values <- moran_values[order(moran_values$`Moran I statistic`, decreasing = TRUE), ]
  top_genes <- rownames(moran_values[moran_values$`Moran I statistic` > 0.3, ])
  
  seu <- add_spatial_metadata(seu, coords, assay = "SCT") 
  pdf(paste0('../output/spatial_RNA_HVG_moranI_', n, '_0.3.pdf'), width = 5, height = 4)
  for (i in seq_along(top_genes)) {print(plot_cluster_feature_moran(seu, 'spatial', features = top_genes[i], stat = moran_values$`Moran I statistic`[i], pt.size = 1.2))}
  dev.off()
  
  results_list <- lapply(names(res_list), function(gene) {
    moranI <- res_list[[gene]]$estimate['Moran I statistic']
    pvalue <- original_p_values[gene]
    adjP <- adjusted_p_values[gene]
    return(data.frame(Gene = gene, MoranI = moranI, PValue = pvalue, AdjP = adjP))
  })
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  results_df <- results_df[order(results_df$AdjP), ]
  write.csv(results_df, file = paste0('../output/spatial_RNA_HVG_moranI_', n, '.csv'), row.names = FALSE)
}


# RNA - per cell type - top 500 hvg in each cell type

#' Function which prints a message using shell echo; useful for printing messages from inside mclapply when running in Rstudio # ref: https://stackoverflow.com/questions/17345837/printing-from-mclapply-in-r-studio
message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}

dir.create("../output/spatial_autocorrelation")
seu <- readRDS('../data/Integrated_RNA_3samples.rds')
seu_list <- SplitObject(seu, split.by = 'stage')

for (n in names(seu_list)) {
  print(n)
  seu_all <- seu_list[[n]]
  
  seu_list_by_cluster <- SplitObject(seu_all, split.by = 'Clusters.2')
  top_genes_df <- data.frame(Sample = character(), Cluster = character(), TopGenesCount = numeric(), stringsAsFactors = FALSE)
  
  # results <- mclapply(names(seu_list_by_cluster), function(m) {
  results <- lapply(names(seu_list_by_cluster), function(m) {
    print(paste0(m, "n cells = ", ncol(seu)))
    seu <- seu_list_by_cluster[[m]]
    if (ncol(seu) > 20) {
      seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 500)
      hvg <- VariableFeatures(seu)
      print(paste0("number of HVG found for ", n, ", ", m, ": ", length(hvg)))
      coords <- data.frame(seu$x, seu$y)
      
      coords_nb <- knn2nb(knearneigh(coords, k = 5))
      gene_expression <- as.data.frame(as.matrix((seu[hvg, ]@assays$SCT@counts)))
      gene_expression <- gene_expression[rowSums(gene_expression) > 0, ]
      
      res_list <- list()
      for (gene in rownames(gene_expression)) {
        gene_expression_vector <- as.numeric(gene_expression[gene, ])
        listw <- nb2listw(coords_nb)
        moran_result <- moran.test(gene_expression_vector , listw)
        res_list[[gene]] <- moran_result
      }
      
      original_p_values <- sapply(res_list, function(x) x$p.value)
      adjusted_p_values <- p.adjust(original_p_values, method = "fdr")
      names(adjusted_p_values) <- names(res_list)
      significant_adjusted <- adjusted_p_values[adjusted_p_values < 0.05]
      
      print(length(significant_adjusted))
      
      if (length(significant_adjusted) > 0 & length(moran_values$`Moran I statistic`) > 0) {
        
        moran_values <- as.data.frame(t(as.data.frame(sapply(res_list[names(significant_adjusted)], function(x) x$estimate))))
        moran_values <- moran_values[order(moran_values$`Moran I statistic`, decreasing = TRUE), ]
        
        top_genes <- rownames(moran_values[moran_values$`Moran I statistic` > 0.3, ])
        
        seu <- add_spatial_metadata(seu, coords, assay = "SCT")
        
        if (length(top_genes) > 0) {
          pdf(paste0('../output/spatial_autocorrelation/spatial_RNA_HVG_moranI_', n, "_", m, '_0.3.pdf'), width = 5, height = 4)
          for (i in seq_along(top_genes)) {
            print(plot_cluster_feature_moran(seu, 'spatial', features = top_genes[i], stat = moran_values$`Moran I statistic`[i], pt.size = 1.2))
          }
          dev.off()
        }
        
        results_list <- lapply(names(res_list), function(gene) {
          moranI <- res_list[[gene]]$estimate['Moran I statistic']
          pvalue <- original_p_values[gene]
          adjP <- adjusted_p_values[gene]
          return(data.frame(Gene = gene, MoranI = moranI, PValue = pvalue, AdjP = adjP))
        })
        results_df <- do.call(rbind, results_list)
        rownames(results_df) <- NULL
        
        results_df <- results_df[order(results_df$AdjP), ]
        write.csv(results_df, file = paste0('../output/spatial_autocorrelation/spatial_RNA_HVG_moranI_', n,'_', m, '.csv'), row.names = FALSE)
        
        return(data.frame(Sample = n, Cluster = m, TopGenesCount = length(top_genes)))
        # return(results_df)
      }
    } else {
      print(paste0("skipping ", m, ", not enough cells"))
      return(data.frame(Sample = n, Cluster = m, TopGenesCount = NA))
    }
    
  }
  # , mc.cores = 10
  )
  top_genes_df <- do.call(rbind, results)
  write.csv(top_genes_df, file = paste0('../output/spatial_autocorrelation/spatial_RNA_HVG_summary_', n, '.csv'), row.names = FALSE)
}



# ATAC - all cells - marker motifs

# ATAC
cluster_levels <- c("Endothelial", "Erythroblasts", "EVT", "EVT-progenitor", "FIB1", "FIB2", 
                    "Hofbauer cells", "STB", "STB-progenitor", "vCTB")
gex_cluster <- data.frame('x_um'=int_mapped$x, 'y_um'=int_mapped$y,'sample'=int_mapped$stage, 'stage' = int_mapped$stage,
                         'cluster'=factor(int_mapped$Clusters.2, levels=cluster_levels))
rownames(gex_cluster) <- rownames(int_mapped@meta.data)

dt = readRDS('../data/int_MotifMatrix_Z.rds')
dt = cbind(gex_cluster[rownames(gex_cluster), c('sample', 'stage'), drop=F], dt[rownames(gex_cluster),])
dt <- dt[!is.na(dt$x) , ]
dt_list <- split(dt, dt$stage)
for (n in names(dt_list)) {
  print(n)
  dt <- dt_list[[n]]
  #hvg <- colnames(dt)[5:ncol(dt)]
  coords <- data.frame(dt$x, dt$y)
  coords_nb <- knn2nb(knearneigh(coords, k = 5))
  gene_expression <- t(dt[,5:ncol(dt)])
  #gene_expression <- gene_expression[rowSums(gene_expression) > 0, ]
  
  res_list <- list()
  for (gene in rownames(gene_expression)) {
    gene_expression_vector <- as.numeric(gene_expression[gene, ])
    listw <- nb2listw(coords_nb)
    moran_result <- moran.test(gene_expression_vector , listw)
    res_list[[gene]] <- moran_result
  }
  
  original_p_values <- sapply(res_list, function(x) x$p.value)
  adjusted_p_values <- p.adjust(original_p_values, method = "fdr")
  names(adjusted_p_values) <- names(res_list)
  significant_adjusted <- adjusted_p_values[adjusted_p_values < 0.05]
  
  moran_values <- as.data.frame(t(as.data.frame(sapply(res_list[names(significant_adjusted)], function(x) x$estimate))))
  moran_values <- moran_values[order(moran_values$`Moran I statistic`, decreasing = TRUE), ]
  top_genes <- rownames(moran_values[moran_values$`Moran I statistic` > 0.3, ])
  
  # seu <- add_spatial_metadata(seu, coords, assay = "SCT") 
  # pdf(paste0('../output/spatial_RNA_HVG_moranI_', n, '_0.3.pdf'), width = 5, height = 4)
  # for (i in seq_along(top_genes)) {print(plot_cluster_feature_moran(seu, 'spatial', features = top_genes[i], stat = moran_values$`Moran I statistic`[i], pt.size = 1.2))}
  # dev.off()
  
  results_list <- lapply(names(res_list), function(gene) {
    moranI <- res_list[[gene]]$estimate['Moran I statistic']
    pvalue <- original_p_values[gene]
    adjP <- adjusted_p_values[gene]
    return(data.frame(Gene = gene, MoranI = moranI, PValue = pvalue, AdjP = adjP))
  })
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  results_df <- results_df[order(results_df$AdjP), ]
  write.csv(results_df, file = paste0('../output/spatial_MOTIF_moranI_', n, '.csv'), row.names = FALSE)
}


# ATAC - All cells - All motifs
Matrix <- getMatrixFromProject(proj_filter, useMatrix = "MotifMatrix",
                               useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads())
data <- as.data.frame(Matrix@assays@data$z)
colnames(data) <- gsub('combined#','',colnames(data))
updated_barcodes <- sapply(colnames(data), function(bc) {
  if (grepl("-1$", bc)) {
    return(paste0("JS35_", sub("-1$", "", bc)))
  } else if (grepl("-2$", bc)) {
    return(paste0("JS36_", sub("-2$", "", bc)))
  } else if (grepl("-3$", bc)) {
    return(paste0("JS40_", sub("-3$", "", bc)))
  } else {
    return(bc)
  }
})
names(updated_barcodes) <- NULL
colnames(data) <- updated_barcodes
rownames(data) <- sub("_\\d+$", "", rownames(data))
dt <- data
dt <- t(dt)
dt = merge(gex_cluster[rownames(gex_cluster), c('sample', 'stage','x_um','y_um'), drop=F], dt[which(rownames(dt) %in% rownames(gex_cluster)),], by = 0)
#dt <- dt[!is.na(dt$x_um) , ]
dt_list <- split(dt, dt$stage)
for (n in names(dt_list)) {
  print(n)
  dt <- dt_list[[n]]
  #hvg <- colnames(dt)[5:ncol(dt)]
  coords <- data.frame(dt$x_um, dt$y_um)
  coords_nb <- knn2nb(knearneigh(coords, k = 5))
  gene_expression <- t(dt[,6:ncol(dt)])
  #gene_expression <- gene_expression[rowSums(gene_expression) > 0, ]
  
res_list <- list()
skipped_genes <- c()
for (gene in rownames(gene_expression)) {
  gene_expression_vector <- as.numeric(gene_expression[gene, ])
  
  if (any(!is.finite(gene_expression_vector))) {
    # Skip genes with non-finite values
    skipped_genes <- c(skipped_genes, gene)
    next
  }
  
  listw <- nb2listw(coords_nb)
  moran_result <- moran.test(gene_expression_vector, listw)
  res_list[[gene]] <- moran_result
}

if (length(skipped_genes) > 0) {
  cat("The following genes were skipped due to non-finite values:", paste(skipped_genes, collapse = ", "), "\n")
}
  
  original_p_values <- sapply(res_list, function(x) x$p.value)
  adjusted_p_values <- p.adjust(original_p_values, method = "fdr")
  names(adjusted_p_values) <- names(res_list)
  significant_adjusted <- adjusted_p_values[adjusted_p_values < 0.05]
  
  moran_values <- as.data.frame(t(as.data.frame(sapply(res_list[names(significant_adjusted)], function(x) x$estimate))))
  moran_values <- moran_values[order(moran_values$`Moran I statistic`, decreasing = TRUE), ]
  top_genes <- rownames(moran_values[moran_values$`Moran I statistic` > 0.3, ])
  
  # seu <- add_spatial_metadata(seu, coords, assay = "SCT") 
  # pdf(paste0('../output/spatial_RNA_HVG_moranI_', n, '_0.3.pdf'), width = 5, height = 4)
  # for (i in seq_along(top_genes)) {print(plot_cluster_feature_moran(seu, 'spatial', features = top_genes[i], stat = moran_values$`Moran I statistic`[i], pt.size = 1.2))}
  # dev.off()
  
  results_list <- lapply(names(res_list), function(gene) {
    moranI <- res_list[[gene]]$estimate['Moran I statistic']
    pvalue <- original_p_values[gene]
    adjP <- adjusted_p_values[gene]
    return(data.frame(Gene = gene, MoranI = moranI, PValue = pvalue, AdjP = adjP))
  })
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  results_df <- results_df[order(results_df$AdjP), ]
  write.csv(results_df, file = paste0('../output/spatial_MOTIF_ALL_moranI_', n, '.csv'), row.names = FALSE)
}

# ATAC - per cell type - all motifs
dir.create("../output/spatial_autocorrelation_ATAC")

Matrix <- getMatrixFromProject(proj_filter, useMatrix = "MotifMatrix",
                               useSeqnames = NULL, verbose = TRUE, binarize = FALSE, threads = getArchRThreads())
data <- as.data.frame(Matrix@assays@data$z)
colnames(data) <- gsub('combined#','',colnames(data))
updated_barcodes <- sapply(colnames(data), function(bc) {
  if (grepl("-1$", bc)) {
    return(paste0("JS35_", sub("-1$", "", bc)))
  } else if (grepl("-2$", bc)) {
    return(paste0("JS36_", sub("-2$", "", bc)))
  } else if (grepl("-3$", bc)) {
    return(paste0("JS40_", sub("-3$", "", bc)))
  } else {
    return(bc)
  }
})
names(updated_barcodes) <- NULL
colnames(data) <- updated_barcodes
rownames(data) <- sub("_\\d+$", "", rownames(data))
dt <- data
dt <- t(dt)
dt = merge(gex_cluster[rownames(gex_cluster), c('sample', 'stage','x_um','y_um','cluster'), drop=F], dt[which(rownames(dt) %in% rownames(gex_cluster)),], by = 0)
#dt <- dt[!is.na(dt$x_um) , ]
dt_list <- split(dt, dt$stage)

for (n in names(dt_list)) {
  print(n)
  dt_all <- dt_list[[n]]
  dt_sublist <- split(dt_all, dt_all$cluster)
  top_genes_df <- data.frame(Sample = character(), Cluster = character(), TopGenesCount = numeric(), stringsAsFactors = FALSE)
  
  results <- lapply(names(dt_sublist), function(m) {
    
    dt <- dt_sublist[[m]]
    print(paste0(m, ", n cells = ", nrow(dt)))
    
  if (nrow(dt) > 20) {
  coords <- data.frame(dt$x_um, dt$y_um)
  coords_nb <- knn2nb(knearneigh(coords, k = 5))
  gene_expression <- t(dt[,6:ncol(dt)])
  #gene_expression <- gene_expression[rowSums(gene_expression) > 0, ]
  
res_list <- list()
skipped_genes <- c()
for (gene in rownames(gene_expression)) {
  gene_expression_vector <- as.numeric(gene_expression[gene, ])
  
  if (any(!is.finite(gene_expression_vector))) {
    # Skip genes with non-finite values
    skipped_genes <- c(skipped_genes, gene)
    next
  }
  
  listw <- nb2listw(coords_nb)
  moran_result <- moran.test(gene_expression_vector, listw)
  res_list[[gene]] <- moran_result
}

if (length(skipped_genes) > 0) {
  cat("The following genes were skipped due to non-finite values:", paste(skipped_genes, collapse = ", "), "\n")
}
  
  original_p_values <- sapply(res_list, function(x) x$p.value)
  adjusted_p_values <- p.adjust(original_p_values, method = "fdr")
  names(adjusted_p_values) <- names(res_list)
  significant_adjusted <- adjusted_p_values[adjusted_p_values < 0.05]
  
  if (length(significant_adjusted) > 0 & length(moran_values$`Moran I statistic`) > 0) {
    
  moran_values <- as.data.frame(t(as.data.frame(sapply(res_list[names(significant_adjusted)], function(x) x$estimate))))
  moran_values <- moran_values[order(moran_values$`Moran I statistic`, decreasing = TRUE), ]
  top_genes <- rownames(moran_values[moran_values$`Moran I statistic` > 0.3, ])
  
  
  # pdf(paste0('../output/spatial_RNA_HVG_moranI_', n, '_0.3.pdf'), width = 5, height = 4)
  # for (i in seq_along(top_genes)) {print(plot_cluster_feature_moran(seu, 'spatial', features = top_genes[i], stat = moran_values$`Moran I statistic`[i], pt.size = 1.2))}
  # dev.off()
if (length(top_genes) > 0) {  
df_melt <- dt_sublist[[m]][,-which(names(dt_sublist[[m]]) %in% c("Row.names", "sample", "stage", "cluster"))]  
df_melt <- reshape2::melt(df_melt, id.vars = c("x_um", "y_um"), variable.name = "Motif", value.name = "Gene_Score")
pdf(paste0('../output/spatial_autocorrelation_ATAC/spatial_ATAC_HVG_moranI_', n, "_", m, '_0.3.pdf'), width = 6, height = 5)
for(motif in top_genes){
    motif_data <- subset(df_melt, Motif == motif)
    motif_data <- motif_data[order(motif_data$Gene_Score),]
    motif_data <- motif_data %>% arrange(desc(is.na(motif_data)))
  p <- ggplot(motif_data, aes(x = x_um, y = y_um, fill = Gene_Score)) + 
    geom_point(size = 2, shape = 21, color = "grey") + 
      scale_fill_continuous_diverging("Blue-Red 3", na.value = "grey") +
    theme_bw() + 
    labs(title = motif) +
    theme(legend.position = "right") +
    guides(color = guide_colorbar(title = "Deviation \nScore")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          legend.position = "right",
          panel.border = element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 14, 
                                    face = "bold", margin = margin(r = 10))) +
    coord_fixed()
  print(p)
}
dev.off()
}  
  
  
  results_list <- lapply(names(res_list), function(gene) {
    moranI <- res_list[[gene]]$estimate['Moran I statistic']
    pvalue <- original_p_values[gene]
    adjP <- adjusted_p_values[gene]
    return(data.frame(Gene = gene, MoranI = moranI, PValue = pvalue, AdjP = adjP))
  })
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  results_df <- results_df[order(results_df$AdjP), ]
  write.csv(results_df, file = paste0('../output/spatial_autocorrelation_ATAC/spatial_MOTIF_ALL_moranI_', n,'_', m, '.csv'), row.names = FALSE)

       return(data.frame(Sample = n, Cluster = m, TopGenesCount = length(top_genes)))
        # return(results_df)
  }  
    } else {
      print(paste0("skipping ", m, ", not enough cells"))
      return(data.frame(Sample = n, Cluster = m, TopGenesCount = NA))
    }
    
  }
  # , mc.cores = 10
  )
  top_genes_df <- do.call(rbind, results)
  write.csv(top_genes_df, file = paste0('../output/spatial_autocorrelation_ATAC/spatial_ATAC_HVG_summary_', n, '.csv'), row.names = FALSE)
}



############ 3. WNN for combination of RNA and ATAC ############
addArchRGenome("hg38")
addArchRThreads(threads = floor(parallel::detectCores()/2))
pathToMacs2 = '/Users/arussell/opt/anaconda3/bin/macs2'
folder = '../data'
annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38") 
atac_peaks <- lapply(c('js35','js36','js40'), function(x){
  index = list('js35'='-1','js36'='-2','js40'='-3')[[x]]
  path <- paste0(folder, '/', x, '_filtered_feature_bc_matrix.h5')
  atac_count <- Read10X_h5(path)
  peaks <- atac_count$Peaks
  colnames(peaks) <- gsub('-1',index,colnames(peaks))
  return(peaks)
})
all_rows <- Reduce(union, lapply(atac_peaks, rownames))
complete_matrix <- function(mat, all_rows) {
  missing_rows <- setdiff(all_rows, rownames(mat))
  missing_mat <- Matrix(0, nrow = length(missing_rows), ncol = ncol(mat),
                        dimnames = list(missing_rows, colnames(mat)))
  rbind(mat, missing_mat)
}
atac_peaks_complete <- lapply(atac_peaks, complete_matrix, all_rows = all_rows)
merged_matrix <- do.call(cbind, atac_peaks_complete)
chrom_assay <- CreateChromatinAssay(
  counts = merged_matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = inputFiles,
  min.cells = 1,
  annotation = annotations
)
updated_barcodes <- sapply(colnames(chrom_assay), function(bc) {
  if (grepl("-1$", bc)) {
    return(paste0("JS35_", sub("-1$", "", bc)))
  } else if (grepl("-2$", bc)) {
    return(paste0("JS36_", sub("-2$", "", bc)))
  } else if (grepl("-3$", bc)) {
    return(paste0("JS40_", sub("-3$", "", bc)))
  } else {
    return(bc)
  }
})
names(updated_barcodes) <- NULL
colnames(chrom_assay) <- updated_barcodes
chrom_assay <- subset(chrom_assay, cells =intersect(updated_barcodes, colnames(int)))
int[["ATAC"]] <- chrom_assay

# Combine RNA and ATAC; Use harmony to remove batch effects
DefaultAssay(int) <- "ATAC"
int <- RunTFIDF(int)
int <- FindTopFeatures(int, min.cutoff = 'q0')
int <- RunSVD(int)
int <- RunHarmony(int, group.by.vars = 'stage', reduction.use = "lsi", reduction.save = "harmony_lsi", project.dim = FALSE)
int <- RunUMAP(int, reduction = 'harmony_lsi', dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
int <- FindMultiModalNeighbors(int, reduction.list = list("harmony", "harmony_lsi"), dims.list = list(1:20, 2:20))
int <- FindClusters(int, graph.name = "wsnn", algorithm = 1, resolution = 0.8)
int <- RunUMAP(int, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", min.dist = 0.1)
DimPlot(int, reduction = "wnn.umap", label = T)
DimPlot(int, reduction = "wnn.umap", label = T, group.by = 'Clusters.2')
saveRDS(int, file = '../data/Integrated_WNN_harmony.rds')







############ 4. CellChat ############

seu <- readRDS('../data/Integrated_RNA_3samples.rds')
seulist = SplitObject(seu, split.by = 'group')
for (i in names(seulist)) {
  saveRDS(seulist[[i]], file = paste0('../data/', i, '.rds'))
} 

save_dir = '../output/cci/'
dir.create(save_dir)
input_JS_list = c('40', '35', '36')
for (input_JS in input_JS_list) {
  seu = readRDS(paste0('../data/JS', input_JS,'.rds'))
  assay_use = 'SCT'
  idents_col = 'Clusters.2'
  data.input = GetAssayData(seu, slot = "data", assay = assay_use) 
  #data.input = as.matrix(seu@assays[[assay_use]]$data)
  meta = seu@meta.data[,c(idents_col), drop = F] 
  colnames(meta) <- "labels"
  spatial.locs = data.frame(seu$x, seu$y)
  scale.factors = list(spot.diameter = 1, spot = 1)
  cellchat <- createCellChat(object = data.input, 
                             meta = meta, 
                             group.by = "labels",
                             datatype = "spatial", 
                             coordinates = spatial.locs,
                             scale.factors = scale.factors
                             )
  cellchat@DB <- CellChatDB.human
  # # try to filter the pathways based on scRNA results,
  # # but failed due to poor enrichemnts cause by poor expression.
  # cellchat <- subsetData(cellchat, features = subset_pathway_name) 
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 4) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.1, 
                                distance.use = TRUE, interaction.length = 300, 
                                scale.distance = 1)
  cellchat <- filterCommunication(cellchat, min.cells = 1)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  future::plan("sequential")
  saveRDS(cellchat, file = paste0(save_dir, 'cellchat_300um_', input_JS, '.rds'))
}


# figs
cellchat.35 <- readRDS(paste0(save_dir, 'cellchat_300um_35.rds'))
cellchat.36 <- readRDS(paste0(save_dir, 'cellchat_300um_36.rds'))
cellchat.40 <- readRDS(paste0(save_dir, 'cellchat_300um_40.rds'))
cellchat.sc <- readRDS('../data/2023-08-25_cellchat.rds')
cellchat_list = list('js35'=cellchat.35, 'js36'=cellchat.36, 'js40'=cellchat.40)

# plot cell interacting of integrated counts for 3 samples
cluster_levels <- c("Endothelial", "Erythroblasts", "EVT", "EVT-progenitor", "FIB1", "FIB2", 
                    "Hofbauer cells", "STB", "STB-progenitor", "vCTB")
sample_palette <- setNames(c('#af8dc3','#d8b365','#5ab4ac'), c("8w_2d", "9w_1d", "11w_1d"))
sample_ns <- c('js40','js35','js36')
stage <- c("8w_2d","9w_1d","11w_1d")
data_mat <- lapply(1:3, function(i){
  cellchat = cellchat_list[[sample_ns[[i]]]]
  data = cellchat@net$count
  summary_data <- reshape2::melt(data, value.name = "interactions")
  colnames(summary_data)[1:2] <- c("source", "target")
  summary_data$sample <- stage[[i]]
  summary_data$lr <- paste0(summary_data$source,'_',summary_data$target)
  return(summary_data)
}) %>% do.call('rbind',.)
summary_data <- data_mat %>%
  dplyr::group_by(source, target) %>%
  dplyr::mutate(value = sum(interactions)) %>% 
  ungroup() 
summary_data <- summary_data %>%
  dplyr::group_by(source, target, sample) %>%
  dplyr::mutate(r = sqrt(value / pi)) %>% 
  ungroup()
summary_data <- summary_data %>%
  pivot_wider(names_from = sample, values_from = interactions, values_fill = 0) %>% 
  mutate(source_num = as.numeric(factor(source, levels = cluster_levels)), 
         target_num = as.numeric(factor(target, levels = cluster_levels)))
summary_data$r_1 <- (summary_data$r - min(summary_data$r)) / (max(summary_data$r) - min(summary_data$r)) * 0.1

pdf(paste0(save_dir, 'cellchat_interaction_counts_int.pdf'), width = 8, height = 6)
ggplot() + 
  geom_scatterpie(data=summary_data, 
                  aes(x=source_num, y=target_num, r=r*0.05), 
                  cols=c("8w_2d","9w_1d","11w_1d"),
                  legend_name = 'sample',) + 
  scale_x_continuous(breaks=c(1:10), labels=cluster_levels) + 
  scale_y_continuous(breaks=c(1:10), labels=cluster_levels) + 
  labs(x="source", y="target") + 
  scale_fill_manual(values = sample_palette) +
  coord_fixed() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# dotplot of filtered LRs
# 1.shared LRs between each slidetag sample with scRNA data
cellchat_list = lapply(names(cellchat_list), function(x){
  sample = cellchat_list[[x]]
  lr_sig = sample@LR[["LRsig"]]
  sc_lr_sig = cellchat.sc@LR[["LRsig"]]$interaction_name
  com_lr = intersect(lr_sig$interaction_name, sc_lr_sig)
  lr_sig = lr_sig[lr_sig$interaction_name %in% com_lr,]
  sample@LR[["LRsig"]] = lr_sig
  return(sample)
}) %>% {names(.) <- c('JS35','JS36','JS40');.}

pdf(paste0(save_dir, '1.cellchat_LR_3sample_individual.pdf'), height = 20, width = 8)
for (n in names(cellchat_list)) {
  cellchat = cellchat_list[[n]]
  p <- netVisual_bubble(cellchat, sources.use = 1:11, 
                        targets.use = 1:11, 
                        signaling = unique(cellchat@LR$LRsig$pathway_name), 
                        remove.isolate = TRUE) +
    labs(title = n)
  print(p)
}
dev.off()

# 2. shared LRs between all slidetag samples and scRNA data
com_st_lr = lapply(names(cellchat_list), function(x){
  sample = cellchat_list[[x]]
  lr_sig = sample@LR[["LRsig"]]
  sc_lr_sig = cellchat.sc@LR[["LRsig"]]$interaction_name
  com_lr = intersect(lr_sig$interaction_name, sc_lr_sig)
  lr_sig = lr_sig[lr_sig$interaction_name %in% com_lr,]
  lr_sig$sample <- x
  return(lr_sig)
}) %>% do.call('rbind', .)
common_interactions <- com_st_lr %>%
  group_by(interaction_name) %>%
  tally() %>%
  filter(n == length(unique(com_st_lr$sample))) %>%
  select(interaction_name)
result <- com_st_lr %>%
  filter(interaction_name %in% common_interactions$interaction_name)

pdf(paste0(save_dir, '2.cellchat_LR_3sample_shared.pdf'), height = 18, width = 8)
for (n in names(cellchat_list)) {
  cellchat = cellchat_list[[n]]
  p <- netVisual_bubble(cellchat, sources.use = 1:11, 
                        targets.use = 1:11, 
                        signaling = unique(result$pathway_name), 
                        remove.isolate = TRUE) +
    labs(title = n)
  print(p)
}
dev.off()

