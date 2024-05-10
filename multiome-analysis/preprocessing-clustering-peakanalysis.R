library(ArchR)
library(Seurat)
library(tibble)
library(dplyr)
library(ggplot2)
set.seed(1)
addArchRThreads(threads = 5)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## --------------------------------------------------------------------------------------------------------------------------
inputs = paste0("../data/Combined-Out-v1/", 
                c("3P/atac_fragments.tsv.gz",
                  "4P/atac_fragments.tsv.gz",
                  "6P/atac_fragments.tsv.gz",
                  "15P/atac_fragments.tsv.gz",
                  "JS34/atac_fragments.tsv.gz",
                  "JS35/atac_fragments.tsv.gz",
                  "JS36/atac_fragments.tsv.gz",
                  "JS40/atac_fragments.tsv.gz"))

names(inputs) = c("3P", "4P", "6P", "15P","JS34","JS35","JS36","JS40")


## --------------------------------------------------------------------------------------------------------------------------
addArchRGenome("hg38")


## --------------------------------------------------------------------------------------------------------------------------
#ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputs,
  sampleNames = names(inputs),
  addGeneScoreMat = TRUE
)


#ArchRProject
proj <- ArchRProject(ArrowFiles)

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "../ArchROutput",
  overwrite = FALSE,
  load = TRUE,
  dropCells = FALSE,
  logFile = createLogFile("saveArchRProject"),
  threads = getArchRThreads()
)

# Import each object separately
rna_3P <- import10xFeatureMatrix(input = "../data/Combined-Out-v1/3P/filtered_feature_bc_matrix.h5", names = "3P")
rna_4P <- import10xFeatureMatrix(input = "../data/Combined-Out-v1/4P/filtered_feature_bc_matrix.h5", names = "4P")
rna_6P <- import10xFeatureMatrix(input = "../data/Combined-Out-v1/6P/filtered_feature_bc_matrix.h5", names = "6P")
rna_15P <- import10xFeatureMatrix(input = "../data/Combined-Out-v1/15P/filtered_feature_bc_matrix.h5", names = "15P")
rna_JS34 <- import10xFeatureMatrix(input = "../data/Combined-Out-v1/JS34/filtered_feature_bc_matrix.h5", names = "JS34")
rna_JS35 <- import10xFeatureMatrix(input = "../data/Combined-Out-v1/JS35/filtered_feature_bc_matrix.h5", names = "JS35")
rna_JS36 <- import10xFeatureMatrix(input = "../data/Combined-Out-v1/JS36/filtered_feature_bc_matrix.h5", names = "JS36")
rna_JS40 <- import10xFeatureMatrix(input = "../data/Combined-Out-v1/JS40/filtered_feature_bc_matrix.h5", names = "JS40")

# Reorder rna_6P to match the order of rna_3P
rna_6P_reordered <- rna_6P[rownames(rna_3P), ]

seRNA <- list("3P" = rna_3P, "4P" = rna_4P, "6P" = rna_6P_reordered, "15P" = rna_15P, "JS34" = rna_JS34, "JS35" = rna_JS35, "JS36" = rna_JS36, "JS40" = rna_JS40)

rows1 <- rownames(seRNA[[1]])
rows2 <- rownames(seRNA[[2]])
rows3 <- rownames(seRNA[[3]])
rows4 <- rownames(seRNA[[4]])
rows5 <- rownames(seRNA[[5]])
rows6 <- rownames(seRNA[[6]])
rows7 <- rownames(seRNA[[7]])
rows8 <- rownames(seRNA[[8]])

rows <- intersect(rows1, rows2)
rows <- intersect(rows, rows3)
rows <- intersect(rows, rows4)
rows <- intersect(rows, rows5)
rows <- intersect(rows, rows6)
rows <- intersect(rows, rows7)
rows <- intersect(rows, rows8)

seRNA[[1]] <- seRNA[[1]][rows,]
seRNA[[2]] <- seRNA[[2]][rows,]
seRNA[[3]] <- seRNA[[3]][rows,]
seRNA[[4]] <- seRNA[[4]][rows,]
seRNA[[5]] <- seRNA[[5]][rows,]
seRNA[[6]] <- seRNA[[6]][rows,]
seRNA[[7]] <- seRNA[[7]][rows,]
seRNA[[8]] <- seRNA[[8]][rows,]

mask1_3 <- (seRNA[[1]]@rowRanges@ranges %in% seRNA[[3]]@rowRanges@ranges)
mask1_4 <- (seRNA[[1]]@rowRanges@ranges %in% seRNA[[4]]@rowRanges@ranges)
mask1_5 <- (seRNA[[1]]@rowRanges@ranges %in% seRNA[[5]]@rowRanges@ranges)
mask1_6 <- (seRNA[[1]]@rowRanges@ranges %in% seRNA[[6]]@rowRanges@ranges)
mask1_7 <- (seRNA[[1]]@rowRanges@ranges %in% seRNA[[7]]@rowRanges@ranges)
mask1_8 <- (seRNA[[1]]@rowRanges@ranges %in% seRNA[[8]]@rowRanges@ranges)
mask <- mask1_3 & mask1_4 & mask1_5 & mask1_6 & mask1_7 & mask1_8

seRNA[[1]] <- seRNA[[1]][mask,]
seRNA[[2]] <- seRNA[[2]][mask,]
seRNA[[3]] <- seRNA[[3]][mask,]
seRNA[[4]] <- seRNA[[4]][mask,]
seRNA[[5]] <- seRNA[[5]][mask,]
seRNA[[6]] <- seRNA[[6]][mask,]
seRNA[[7]] <- seRNA[[7]][mask,]
seRNA[[8]] <- seRNA[[8]][mask,]

seRNA <- Reduce("cbind", seRNA)


## --------------------------------------------------------------------------------------------------------------------------
# remove low/high gene cells
detected = colSums(assay(seRNA)>0)
mask1 = detected>500
mask2 = detected<6000

mask = mask1&mask2


plot(density(detected[mask]))

## --------------------------------------------------------------------------------------------------------------------------
seRNA_hqc = seRNA[,mask]

## --------------------------------------------------------------------------------------------------------------------------



## --------------------------------------------------------------------------------------------------------------------------
#Add scRNA
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA_hqc, force = TRUE, excludeChr = c("chrM"))



## --------------------------------------------------------------------------------------------------------------------------
p1 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

## --------------------------------------------------------------------------------------------------------------------------
p2 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2


## --------------------------------------------------------------------------------------------------------------------------
p3 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
p3


## --------------------------------------------------------------------------------------------------------------------------
p4 <- plotGroups(
  ArchRProj = proj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p4

## --------------------------------------------------------------------------------------------------------------------------
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)


## --------------------------------------------------------------------------------------------------------------------------
p1 <- plotFragmentSizes(ArchRProj = proj)
p1

## --------------------------------------------------------------------------------------------------------------------------

proj <- proj[proj$TSSEnrichment > 6 & proj$nFrags > 2500 & !is.na(proj$Gex_nUMI)]

## --------------------------------------------------------------------------------------------------------------------------

saveArchRProject(proj, outputDirectory = dirname(rstudioapi::getActiveDocumentContext()$path))

## --------------------------------------------------------------------------------------------------------------------------


#Remove doublets
doubScores <- addDoubletScores(input = proj)
proj <- filterDoublets(doubScores)


#LSI-ATAC
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC",
  force = T
)

#LSI-RNA
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  force = T
)




## --------------------------------------------------------------------------------------------------------------------------


# create a new column for gestational week
proj$gest_week <- NA

# manually input the gestational week for each sample
proj$gest_week[proj$Sample == "JS34"] <- "gw_early"
proj$gest_week[proj$Sample == "JS35"] <- "gw_late"
proj$gest_week[proj$Sample == "JS36"] <- "gw_late"
proj$gest_week[proj$Sample == "JS40"] <- "gw_early"
proj$gest_week[proj$Sample == "3P"] <- "gw_early"
proj$gest_week[proj$Sample == "4P"] <- "gw_late"
proj$gest_week[proj$Sample == "6P"] <- "gw_early"
proj$gest_week[proj$Sample == "15P"] <- "gw_late"


## --------------------------------------------------------------------------------------------------------------------------
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_RNA",
  name = "Harmony_RNA",
  groupBy = "Sample"
)

## --------------------------------------------------------------------------------------------------------------------------
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_ATAC",
  name = "Harmony_ATAC",
  groupBy = "Sample"
)


## --------------------------------------------------------------------------------------------------------------------------
mean(proj@cellColData@listData$nFrags)
mean(proj@cellColData@listData$Gex_nUMI)
mean(proj@cellColData@listData$Gex_nGenes)

median(proj@cellColData@listData$nFrags)
median(proj@cellColData@listData$Gex_nUMI)
median(proj@cellColData@listData$Gex_nGenes)



## --------------------------------------------------------------------------------------------------------------------------
#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("Harmony_ATAC", "Harmony_RNA"), name =  "LSI_Combined")

r = 0.3
proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = r, force = TRUE)

p = plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p, name = "UMAP-Combined-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p = plotEmbedding(proj, name = "Sample", embedding = "UMAP_Combined", size = 0.5, labelAsFactors=F, labelMeans=F)
plotPDF(p, name = "UMAP-Combined-Sample.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5 )


## Find DEGs for every cluster
#Identifying Marker Genes
markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneExpressionMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markers, cutOff = "FDR <= 0.1 & Log2FC >= 1.25")

# Get the top 50 genes per cluster based on FDR
topGenesList <- lapply(names(markerList), function(cluster) {
  df <- markerList[[cluster]]
  df <- df[order(df$FDR),]
  top50 <- head(df, 50)
  top50$cluster <- cluster
  return(top50)
})

# Combine the top 20 genes from each cluster into a single data frame
topGenes <- do.call(rbind, topGenesList)%>%as.data.frame()%>%dplyr::select(cluster,name, Log2FC, FDR, MeanDiff)

clabels = read.table('cluster_labels_2023_07_res03.tsv', sep='\t',row.names = 1)
clabels

for (c in row.names(clabels)) {
  print(clabels[c,1])
  tmp = proj$Clusters
  mask = tmp==c
  tmp[mask] = clabels[c,1]
  # tmp = str_replace(tmp, c, clabels[c,1])
  # print(unique(tmp))
  proj$Clusters=tmp
}

#UMAPs
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "Harmony_RNA", name = "UMAP_Harmony_RNA", minDist = 0.8, force = TRUE)
proj <- addUMAP(proj, reducedDims = "Harmony_ATAC", name = "UMAP_Harmony_ATAC", minDist = 0.8, force = TRUE)

# proj <- addUMAP(proj, reducedDims = "Harmony_RNA", name = "UMAP_Harmony_RNA", minDist = 0.8, force = TRUE)
p = plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Harmony_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p, name = "UMAP-Harmony-RNA-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
# proj <- addUMAP(proj, reducedDims = "Harmony_ATAC", name = "UMAP_Harmony_ATAC", minDist = 0.8, force = TRUE)
p = plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Harmony_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p, name = "UMAP-Harmony-ATAC-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


## --------------------------------------------------------------------------------------------------------------------------


saveArchRProject(ArchRProj = proj, outputDirectory = ".", load = FALSE)

## --------------------------------------------------------------------------------------------------------------------------
p = plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p, name = "UMAP-Combined-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

idxSample <- BiocGenerics::which(proj$Clusters %in% c("Myeloid_cells"))
cellsSample <- proj$cellNames[idxSample]
proj_Mac=proj[cellsSample, ]

proj_Mac <- addClusters(proj_Mac, reducedDims = "LSI_Combined", name = "Mac_Clusters", resolution = 0.1, force = TRUE)

proj_Mac <- addUMAP(proj_Mac, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)


p46 <- plotEmbedding(ArchRProj = proj_Mac, colorBy = "cellColData", name = "Mac_Clusters", embedding = "UMAP_Combined", baseSize = 10, plotAs = "points", size = 3)
p47 <- plotEmbedding(ArchRProj = proj_Mac, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined", baseSize = 10, plotAs = "points", size = 3)

ggAlignPlots(p46, p47, type = "h")

p47 = p47 + theme(legend.text = element_text(size = 8), legend.key.size = unit(0.2, "cm"))
plotPDF(p47, name = "UMAP-Combined-Mac_Sample.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p46 = p46 + theme(legend.text = element_text(size = 8))
plotPDF(p46, name = "UMAP-Combined-Mac_Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

## --------------------------------------------------------------------------------------------------------------------------

#Identifying Marker Genes
markers_Mac <- getMarkerFeatures(
  ArchRProj = proj_Mac, 
  useMatrix = "GeneExpressionMatrix", 
  groupBy = "Mac_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList_Mac <- getMarkers(markers_Mac, cutOff = "FDR <= 0.1 & Log2FC >= 1.25")

# Get the top 20 genes per cluster based on FDR
topGenesList_Mac <- lapply(names(markerList_Mac), function(cluster) {
  df <- markerList_Mac[[cluster]]
  df <- df[order(df$FDR),]
  top20 <- head(df, 20)
  top20$cluster <- cluster
  return(top20)
})

# Combine the top 20 genes from each cluster into a single data frame
topGenes_Mac <- do.call(rbind, topGenesList_Mac)%>%as.data.frame()%>%dplyr::select(cluster,name, Log2FC, FDR, MeanDiff)
write.table(topGenes_Mac, paste0("./ArchROutput/Plots/marker_gene_subclusters-Mac-", r, ".txt"), append = FALSE, sep = "\t", col.names = T, row.names = F)



## --------------------------------------------------------------------------------------------------------------------------
plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)

idxSample <- BiocGenerics::which(proj$Clusters %in% c("EVT"))
cellsSample <- proj$cellNames[idxSample]
proj_EVT=proj[cellsSample, ]

proj_EVT <- addClusters(proj_EVT, reducedDims = "LSI_Combined", name = "EVT_Clusters", resolution = 0.1, force = TRUE)

proj_EVT <- addUMAP(proj_EVT, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)


p46 <- plotEmbedding(ArchRProj = proj_EVT, colorBy = "cellColData", name = "EVT_Clusters", embedding = "UMAP_Combined", baseSize = 10, plotAs = "points", size = 3)
p47 <- plotEmbedding(ArchRProj = proj_EVT, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Combined", baseSize = 10, plotAs = "points", size = 3)

ggAlignPlots(p46, p47, type = "h")

p47 = p47 + theme(legend.text = element_text(size = 8), legend.key.size = unit(0.2, "cm"))
plotPDF(p47, name = "UMAP-Combined-EVT_Sample.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p46 = p46 + theme(legend.text = element_text(size = 8))
plotPDF(p46, name = "UMAP-Combined-EVT_Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

## --------------------------------------------------------------------------------------------------------------------------

#Identifying Marker Genes
markers_EVT <- getMarkerFeatures(
  ArchRProj = proj_EVT, 
  useMatrix = "GeneExpressionMatrix", 
  groupBy = "EVT_Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList_EVT <- getMarkers(markers_EVT, cutOff = "FDR <= 0.1 & Log2FC >= 1.25")

# Get the top 20 genes per cluster based on FDR
topGenesList_EVT <- lapply(names(markerList_EVT), function(cluster) {
  df <- markerList_EVT[[cluster]]
  df <- df[order(df$FDR),]
  top20 <- head(df, 20)
  top20$cluster <- cluster
  return(top20)
})

# Combine the top 20 genes from each cluster into a single data frame
topGenes_EVT <- do.call(rbind, topGenesList_EVT)%>%as.data.frame()%>%dplyr::select(cluster,name, Log2FC, FDR, MeanDiff)
write.table(topGenes_EVT, paste0("./ArchROutput/Plots/marker_gene_subclusters-EVT-", r, ".txt"), append = FALSE, sep = "\t", col.names = T, row.names = F)


## update cluster info of mac--------------------------------------------------------------------------------------------------------------------------

clabels = read.table('cluster_labels_2023_7_res03_mac.tsv', sep='\t',row.names = 1)
clabels

for (c in row.names(clabels)) {
  print(clabels[c,1])
  tmp = proj_Mac$Mac_Clusters
  mask = tmp==c
  tmp[mask] = clabels[c,1]
  proj_Mac$Mac_Clusters=tmp
}

## --------------------------------------------------------------------------------------------------------------------------

p = plotEmbedding(ArchRProj = proj_Mac, colorBy = "cellColData", name = "Mac_Clusters", embedding = "UMAP_Combined", baseSize = 10, plotAs = "points", size = 3)

plotPDF(p, name = "UMAP-Combined-Mac_Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
## --------------------------------------------------------------------------------------------------------------------------

# Update the clustering annotation to the original proj
new_clusters = proj_Mac$Mac_Clusters


idxSample <- BiocGenerics::which(proj$Clusters %in% c("Myeloid_cells"))
cellsSample <- proj$cellNames[idxSample]

proj$Clusters[match(cellsSample, proj$cellNames)] <- new_clusters

p = plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)

plotPDF(p, name = "UMAP-Combined-Clusters-detailed.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


## update cluster info of EVT --------------------------------------------------------------------------------------------------------------------------

clabels = read.table('cluster_labels_2023_7_res03_EVT.tsv', sep='\t',row.names = 1)
clabels

for (c in row.names(clabels)) {
  print(clabels[c,1])
  tmp = proj_EVT$EVT_Clusters
  mask = tmp==c
  tmp[mask] = clabels[c,1]
  # tmp = str_replace(tmp, c, clabels[c,1])
  # print(unique(tmp))
  proj_EVT$EVT_Clusters=tmp
}


## --------------------------------------------------------------------------------------------------------------------------

p = plotEmbedding(ArchRProj = proj_EVT, colorBy = "cellColData", name = "EVT_Clusters", embedding = "UMAP_Combined", baseSize = 10, plotAs = "points", size = 3)
plotPDF(p, name = "UMAP-Combined-EVT_Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


## --------------------------------------------------------------------------------------------------------------------------

# Update the clustering annotation to the original proj
new_clusters = proj_EVT$EVT_Clusters

idxSample <- BiocGenerics::which(proj$Clusters %in% c("EVT"))
cellsSample <- proj$cellNames[idxSample]

proj$Clusters[match(cellsSample, proj$cellNames)] <- new_clusters

p = plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined")#, size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p, name = "UMAP-Combined-Clusters-detailed.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


## --------------------------------------------------------------------------------------------------------------------------
saveArchRProject(ArchRProj = proj, outputDirectory = ".", load = FALSE)

##---------------- EVT marker gene UMAPs

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneExpressionMatrix", 
  name = topGenes_EVT[,2],
  embedding = "UMAP_Combined",
  log2Norm = T,
)

plotPDF(plotList = p, 
        name = paste0("UMAP-EVT_subclusterDEGs-GEM-log.pdf"), 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)


## ----------------------------------- Male gene UMAPs
##
male_genes = c("DDX3Y","EIF1AY","PRKY","ZFY","XIST")

p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneExpressionMatrix",
  name = male_genes,
  embedding = "UMAP_Combined",
  # highlightCells = cellNames,
  plotAs = "points"
)


plotPDF(plotList = p,
        name = "UMAP-MarkerGenes-GEM-MaleGenes.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)


## --------------------------------------------------------------------------------------------------------------------------
gmat = getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
feats = getFeatures(proj, useMatrix = "GeneExpressionMatrix")
rownames(gmat) = feats
male_genes = c("DDX3Y","EIF1AY","PRKY","ZFY","XIST")

for (g in male_genes){
  g_exist = as(assay(gmat[g,])>0, "integer")
  proj = addCellColData(proj, name = g, data = g_exist, cells = colnames(gmat), force = T)
}

cellCol = getCellColData(proj)

p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = male_genes,
  embedding = "UMAP_Combined",
  # highlightCells = cellNames,
  plotAs = "points"
)

plotPDF(plotList = p,
        name = "UMAP-MarkerGenes-GEM-MaleGenes-exists.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)


# show as table
aggregate(x=list(ncells_DDX3Y=cellCol$DDX3Y), by=list(Sample=cellCol$Sample), FUN=sum)


# Downstream ATAC multiome analysis
## --------------------------------------------------------------------------------------------------------------------------
#Identifying Marker Genes
proj$Clusters = as.character(proj$Clusters)
markersGE <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneExpressionMatrix", 
  groupBy = 'Clusters',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGE, cutOff = "FDR <= 0.1 & Log2FC >= 1.25")
markerList

markerGenes = c()
clusterNames = c()
n = 10
for (i in seq(length(markerList))){
  if (nrow(markerList[[i]])){
    markerGenes = c(markerGenes, markerList[[i]][1:n,'name'])
    clusterNames = c(clusterNames, rep(names(markerList)[[i]], n))
  } else {
    next
  }
}

## ----include=FALSE---------------------------------------------------------------------------------------------------------

#Track Plotting with ArchRBrowser
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
);

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5);


## ----include=FALSE---------------------------------------------------------------------------------------------------------
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters", force=T)
library(reticulate)
use_condaenv("base", required=T)
py_config()
pathToMacs2 <- findMacs2()


## --------------------------------------------------------------------------------------------------------------------------
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  pathToMacs2 =  pathToMacs2
)

getPeakSet(proj)


## --------------------------------------------------------------------------------------------------------------------------
proj <- addPeakMatrix(proj, force = T)

getAvailableMatrices(proj)



## --------------------------------------------------------------------------------------------------------------------------
library(Matrix)

pm = getMatrixFromProject(proj, useMatrix = "PeakMatrix")
pm_features = getPeakSet(proj)

var_names=list()
for (i in seq(nrow(pm))){
  seqnames = pm_features@seqnames
  ranges = pm_features@ranges
  
  var_names[i] = paste0(seqnames[1], ':',ranges@start[1], '-',ranges@start[1]+ranges@start[1]-1)
}

rownames(pm) = var_names

dir.create('./ArchROutput/filtered_peak_bc_matrix')
writeMM(assay(pm), "./ArchROutput/filtered_peak_bc_matrix/matrix.mtx")

cellCol = getCellColData(proj)

write.table(cellCol, './ArchROutput/filtered_peak_bc_matrix/cellCol.tsv', quote=FALSE, row.names = T, col.names = T, sep='\t')
write.table(data.frame(pm_features), './ArchROutput/filtered_peak_bc_matrix/peakAnnotations.tsv', quote=FALSE, row.names = T, col.names = T, sep='\t')


write.table(colnames(pm), './ArchROutput/filtered_peak_bc_matrix/barcodes.tsv', quote=FALSE, row.names = F, col.names = F)

####### export gene matrix
library(Matrix)
gmat = getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
gnames = getFeatures(proj, useMatrix = "GeneExpressionMatrix")

# save sparse matrix
sparse.gmat <- Matrix(assay(gmat) , sparse = T )
head(sparse.gmat)
dir.create('./ArchROutput/filtered_gene_bc_matrix')
writeMM(obj = sparse.gmat, file="./ArchROutput/filtered_gene_bc_matrix/matrix.mtx")

# save genes and cells names
write(x = gnames, file = "./ArchROutput/filtered_gene_bc_matrix/genes.tsv")
write(x = colnames(gmat), file = "./ArchROutput/filtered_gene_bc_matrix/barcodes.tsv")

cellCol = getCellColData(proj)

write.table(cellCol, file = "./ArchROutput/filtered_gene_bc_matrix/cellcol.tsv", sep = "\t", row.names = T, col.names = T, quote=F)

# export lsi

lsi = getReducedDims(proj, reducedDims = "LSI_Combined")

write.table(lsi, file = "./ArchROutput/filtered_gene_bc_matrix/lsi.tsv", sep = "\t", row.names = T, col.names = T, quote = FALSE)

# export rna lsi
lsi_rna = getReducedDims(proj, reducedDims = "LSI_RNA")
write.table(lsi_rna, file = "./ArchROutput/filtered_gene_bc_matrix/lsi_rna.tsv", sep = "\t", row.names = T, col.names = T, quote = FALSE)

# export atac lsi
lsi_atac = getReducedDims(proj, reducedDims = "LSI_ATAC")
write.table(lsi_atac, file = "./ArchROutput/filtered_gene_bc_matrix/lsi_atac.tsv", sep = "\t", row.names = T, col.names = T, quote = FALSE)

harmony_atac = getReducedDims(proj, reducedDims = "Harmony_ATAC")
write.table(lsi_atac, file = "./ArchROutput/filtered_gene_bc_matrix/harmony_atac.tsv", sep = "\t", row.names = T, col.names = T, quote = FALSE)

# export umap

umap = getEmbedding(proj, embedding = 'UMAP_Combined')

write.table(umap, file = "./ArchROutput/filtered_gene_bc_matrix/umap.tsv", sep = "\t", row.names = T, col.names = T, quote=F)

####### export gene score matrix
library(Matrix)
gmat = getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
gnames = getFeatures(proj, useMatrix = "GeneScoreMatrix")

# save sparse matrix
sparse.gmat <- Matrix(assay(gmat) , sparse = T )
head(sparse.gmat)
dir.create('./ArchROutput/filtered_genescore_bc_matrix')
writeMM(obj = sparse.gmat, file="./ArchROutput/filtered_genescore_bc_matrix/matrix.mtx")

# save genes and cells names
write(x = gnames, file = "./ArchROutput/filtered_genescore_bc_matrix/genes.tsv")
write(x = colnames(gmat), file = "./ArchROutput/filtered_genescore_bc_matrix/barcodes.tsv")

cellCol = getCellColData(proj)

write.table(cellCol, file = "./ArchROutput/filtered_genescore_bc_matrix/cellcol.tsv", sep = "\t", row.names = T, col.names = T, quote=F)

# export lsi

lsi = getReducedDims(proj, reducedDims = "LSI_Combined")

write.table(lsi, file = "./ArchROutput/filtered_genescore_bc_matrix/lsi.tsv", sep = "\t", row.names = T, col.names = T, quote = FALSE)

# export umap

umap = getEmbedding(proj, embedding = 'UMAP_Combined')

write.table(umap, file = "./ArchROutput/filtered_genescore_bc_matrix/umap.tsv", sep = "\t", row.names = T, col.names = T, quote=F)

## --------------------------------------------------------------------------------------------------------------------------
table(proj$Clusters)


## --------------------------------------------------------------------------------------------------------------------------
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks


## --------------------------------------------------------------------------------------------------------------------------
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

## --------------------------------------------------------------------------------------------------------------------------
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes,
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
  # features = markerList,
  upstream = 50000,
  downstream = 50000
)

plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)



proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=T)

## --------------------------------------------------------------------------------------------------------------------------
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichMotifs

## --------------------------------------------------------------------------------------------------------------------------
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")


## --------------------------------------------------------------------------------------------------------------------------
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)



## --------------------------------------------------------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg38)
if("Motif" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}

proj <- addBgdPeaks(proj, force=T)

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  matrixName = 'MotifMatrix',
  force = TRUE
)

## --------------------------------------------------------------------------------------------------------------------------
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


## --------------------------------------------------------------------------------------------------------------------------
markerMotifs <- getFeatures(proj, useMatrix = "MotifMatrix")
markerMotifs

## --------------------------------------------------------------------------------------------------------------------------
markerMotifs <- grep("deviations:", markerMotifs, value = TRUE)
markerMotifs

## ----include=FALSE---------------------------------------------------------------------------------------------------------
p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


## ----include=FALSE---------------------------------------------------------------------------------------------------------
p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP_Combined",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "UMAP-Combined-MotifMatrix", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


######### motifs on DORC genes

motifs <- c("NFKB1", "RUNX1", "SPI1", "ETV6", "RELB", "ESRRG")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("deviations:", markerMotifs, value = TRUE)
# markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
p3 = do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

plotPDF(p3, name = "Cowplot-Groups-Deviations-1", width = 10, height = 10, ArchRProj = proj, addDOC = FALSE)

# TP73, SNAI1, MESP2, SMARCC1
motifs <- c("TP73", "SNAI1", "MESP2", "SMARCC1")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("deviations:", markerMotifs, value = TRUE)
# markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
p3 = do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

plotPDF(p3, name = "Cowplot-Groups-Deviations-2", width = 10, height = 10, ArchRProj = proj, addDOC = FALSE)


# RUNX1, NFKB1, SPI1, ETV6, RELB, PITX1, STAT3, TCF21
motifs <- c("RUNX1", "NFKB1", "SPI1", "ETV6", "RELB", "PITX1", "STAT3", "TCF21")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("deviations:", markerMotifs, value = TRUE)
markerMotifs

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
p3 = do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

plotPDF(p3, name = "Cowplot-Groups-Deviations-3", width = 10, height = 10, ArchRProj = proj, addDOC = FALSE)

# ID3, ZEB1, BACH1, ESRRG
motifs <- c("ID3", "ZEB1", "BACH1", "ESRRG")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("deviations:", markerMotifs, value = TRUE)
markerMotifs

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
p3 = do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

plotPDF(p3, name = "Cowplot-Groups-Deviations-4", width = 10, height = 10, ArchRProj = proj, addDOC = FALSE)



## --------------------------------------------------------------------------------------------------------------------------
motifPositions <- getPositions(proj)
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions, 
  groupBy = "Clusters"
)

## --------------------------------------------------------------------------------------------------------------------------
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

## --------------------------------------------------------------------------------------------------------------------------
seTSS <- getFootprints(
  ArchRProj = proj, 
  positions = GRangesList(TSS = getTSS(proj)), 
  groupBy = "Clusters",
  flank = 2000
)


## --------------------------------------------------------------------------------------------------------------------------
plotFootprints(
  seFoot = seTSS,
  ArchRProj = proj, 
  normMethod = "None",
  plotName = "TSS-No-Normalization",
  addDOC = FALSE,
  flank = 2000,
  flankNorm = 100
)


## --------------------------------------------------------------------------------------------------------------------------
proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "LSI_Combined"
)

## --------------------------------------------------------------------------------------------------------------------------
cA <- getCoAccessibility(
  ArchRProj = proj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA


## ----include=FALSE---------------------------------------------------------------------------------------------------------
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(proj),
  verbose = FALSE
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)



## peak2gene--------------------------------------------------------------------------------------------------------------------------
kbp = 1000000
kbp_str = '1M'
proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  useMatrix = "GeneExpressionMatrix",
  maxDist = kbp
)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = F
)

# p2g

dir.create('./ArchROutput/p2g')
write.table(p2g, paste0('./ArchROutput/p2g/p2g_', kbp_str,'bp_all.tsv'), quote=FALSE, row.names = F, col.names = T, sep = '\t')

k=25
nPlot=5000
p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "Clusters", k=k, nPlot=nPlot)
plotPDF(p, 
        name = paste0("Plot-Peak2GeneHeatmap-k",k,"-nPlot",nPlot,".pdf"), 
        ArchRProj = proj, 
        addDOC = FALSE, width = 10, height = 15)

## --------------------------------------------------------------------------------------------------------------------------
kbp = 1000000
kbp_str = '1M'

# iterate through all celltypes
for (c in unique(proj$Clusters)) {
  cells = proj$cellNames[proj$Clusters == c]
  
  tryCatch({ 
    proj <- addPeak2GeneLinks(
      ArchRProj = proj,
      reducedDims = "LSI_Combined",
      useMatrix = "GeneExpressionMatrix",
      maxDist = kbp,
      cellsToUse = cells
    )
  }, error = function(cond) {
    next
  })
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = F
  )
  
  # p2g
  
  write.table(p2g, paste0('./ArchROutput/p2g/p2g_', kbp_str,'bp_', c, '.tsv'), quote=FALSE, row.names = F, col.names = T, sep = '\t')
}


## --------------------------------------------------------------------------------------------------------------------------
kbp = 500000
kbp_str = '500K'

# iterate through all celltypes
for (c in unique(proj$Clusters)) {
  cells = proj$cellNames[proj$Clusters == c]
  tryCatch({ 
    proj <- addPeak2GeneLinks(
      ArchRProj = proj,
      reducedDims = "LSI_Combined",
      useMatrix = "GeneExpressionMatrix",
      maxDist = kbp,
      cellsToUse = cells
    )
  }, error = function(cond){
    next
    print(c)
  })
  
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = F
  )
  
  # p2g
  
  write.table(p2g, paste0('./ArchROutput/p2g/p2g_', kbp_str,'bp_', c, '.tsv'), quote=FALSE, row.names = F, col.names = T, sep = '\t')
}


kbp = 500000
kbp_str = '500K'

proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  useMatrix = "GeneExpressionMatrix",
  maxDist = kbp,
)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = F
)
# p2g

write.table(p2g, paste0('./ArchROutput/p2g/p2g_', kbp_str,'bp_all.tsv'), quote=FALSE, row.names = F, col.names = T, sep = '\t')



kbp = 50000
kbp_str = '50K'

proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  useMatrix = "GeneExpressionMatrix",
  maxDist = kbp,
)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = F
)
# p2g

write.table(p2g, paste0('./ArchROutput/p2g/p2g_', kbp_str,'bp_all.tsv'), quote=FALSE, row.names = F, col.names = T, sep = '\t')


## ----include=FALSE---------------------------------------------------------------------------------------------------------
kbp = 1000000
kbp_str = '1M'
proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  useMatrix = "GeneExpressionMatrix",
  maxDist = kbp
)

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj),
  verbose = FALSE
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks-50K.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

####### p2g and UMAP for DORC genes
library(openxlsx)
library(stringr)
dorcGenes = read.xlsx("DORC Genes.xlsx",colNames = F)

mask = dorcGenes$X1 %in% rownames(seRNA_hqc)
dorcGenes = dorcGenes$X1[mask]


p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = dorcGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj),
  verbose = FALSE
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-DORC-Genes-with-Peak2GeneLinks-50K.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)


## UMAPs
p = plotEmbedding(proj, colorBy = "GeneExpressionMatrix", name = dorcGenes, embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p, name = "UMAP-Combined-DORCGenes.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


## --------------------------------------------------------------------------------------------------------------------------
seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Clusters")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGEM_MM <- correlateMatrices(
  ArchRProj = proj,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "LSI_Combined"
)

corGEM_MM$maxDelta <- rowData(seZ)[match(corGEM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGEM_MM <- corGEM_MM[order(abs(corGEM_MM$cor), decreasing = TRUE), ]
corGEM_MM <- corGEM_MM[which(!duplicated(gsub("\\-.*","",corGEM_MM[,"MotifMatrix_name"]))), ]
corGEM_MM$TFRegulator <- "NO"
corGEM_MM$TFRegulator[which(corGEM_MM$cor > 0.5 & corGEM_MM$padj < 0.01 & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.75))] <- "YES"
sort(corGEM_MM[corGEM_MM$TFRegulator=="YES",1])


p <- ggplot(data.frame(corGEM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGEM_MM$maxDelta)*1.05)
  )

# p

plotPDF(plotList = p,
        name = "Positive-TF-with-GeneExpressionMatrix.pdf",
        ArchRProj = proj,
        addDOC = FALSE, width = 5, height = 5)





#7.3 Identifying marker genes
markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1.25")

# Load the package
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Add dataframes to the workbook
for(i in seq_along(markerList)){
  # Add a new sheet to the workbook
  addWorksheet(wb, sheetName = names(markerList)[i])
  
  # Write data to the sheet
  writeData(wb, sheet = names(markerList)[i], x = markerList[[i]])
}

# Save the workbook to a file
saveWorkbook(wb, "GeneScore_Markers.xlsx", overwrite = TRUE)

top_30_list <- lapply(markerList, function(x) head(x, 30))

names_array <- unlist(lapply(top_30_list, function(x) x$name))

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.25", 
  labelMarkers = NULL,
  transpose = TRUE,
  nLabel= 5
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")


##HEATMAP TOP GENESCORE MARKERS
#EXTENDED DATA FIG 2F

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 12, height = 8, ArchRProj = proj, addDOC = FALSE)


##UMAPs TOP GENESCORE MARKERS
#7.5 Visualizing Marker Genes on an Embedding

#Rearrange for grid plotting
p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = names_array, 
  embedding = "UMAP_Combined",
  imputeWeights = getImputeWeights(proj)
)


for(x in names(top_30_list)){
  
  # Get the current set of markers
  current_markers <- top_30_list[[x]]$name
  if (length(current_markers)==0) {
    next
  }
  
  # Generate the plot
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = current_markers, 
    embedding = "UMAP_Combined",
    imputeWeights = getImputeWeights(proj)
  )
  
  # Define the name of the PDF file
  pdf_name <- paste0("Plot-UMAP-Marker-Genes-W-Imputation_", x, ".pdf")
  
  # Save the plot as a PDF
  plotPDF(plotList = p, 
          name = pdf_name, 
          ArchRProj = proj, 
          addDOC = FALSE, width = 5, height = 5)
}

##SELECTED GENES GENESCORE

# Define the markers for each cell type
markerList <- getMarkers(markersGS)

markers <- list(
  general_markers = c("NOTCH1", "CDX2", "ELF5", "BCAM", "TCF7", "YAP1", "OVOL1", "ITGA2", "WWTR1"),
  vCTB = c("PEG10", "PAGE4", "CDH1"),
  vCTBp = c("MKI67", "TOP2A", "PCNA"),
  vCTB1 = c("ARL15", "TBL1X", "SIK3"),
  vCTB2 = c("IFI6","EFEMP1", "SMAGP"),
  vCTB3 = c("MSI2", "LRP5", "AHCYL2"),
  EVT = c("ITGA5", "HLA-G", "DIO2"),
  EVT1 = c("EGFR", "FBN2", "UTRN"),
  EVT2 = c("CCNE1", "HAPLN3", "LY6E"),
  EVT3 = c("AOC1", "PAPPA2", "PRG2"),
  STB = c("ERVFRD-1", "CYP19A1", "IGSF5"),
  Endo = c("PECAM1", "EGFL7", "KDR"),
  Macrophages = c("SPP1", "CD14", "CD83"),
  Maternal_MAC = c("CD74", "HLA-DRA", "LYZ"),
  HBC = c("F13A1", "LYVE1", "ADAMTS17"),
  Myeloid_unk = c("FGF13", "ALPK3"),
  FIB = c("COL3A1", "SOX5", "LAMA2"),
  FIB1 = c("PDGFRB", "PLCB1", "AGTR1"),
  FIB2 = c("PDGFRA", "COL15A1", "CXCL14"),
  Unknown1 = c("HGF", "DCN", "CDH11"),
  Unknown2 = c("FAM155A", "ALDH1A2", "RORB"),
  Unknown3 = c("RSPO3", "PODN", "BIN1"),
  general_markers = c("NOTCH1", "CDX2", "ELF5", "BCAM", "TCF7", "YAP1", "OVOL1", "ITGA2", "WWTR1"),
  additional_markers = c("PLXNB2", "RUNX1", "c12orf75", "QSOX1", "RASGRF2", "JAK1", "FGFR1", "MYCN", "EBI3", "BRAF", "TBX3", "PDCD1LG2", "CD276")
)


# Iterate through each marker set
for(x in names(markers)){
  
  # Get the current set of markers
  current_markers <- markers[[x]]
  if (length(current_markers)==0){
    next
  }
  
  # Generate the plot
  p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = current_markers, 
    embedding = "UMAP_Combined",
    imputeWeights = getImputeWeights(proj)
  )
  
  # Define the name of the PDF file
  pdf_name <- paste0("Plot-UMAP-Selected-Genes-W-Imputation_", x, ".pdf")
  
  # Save the plot as a PDF
  plotPDF(plotList = p, 
          name = pdf_name, 
          ArchRProj = proj, 
          addDOC = FALSE, width = 5, height = 5)
}



## -----------------------------------------
##13.1 Motif Deviations
if("Motif" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}
# 
proj <- addBgdPeaks(proj, force=T)

proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  matrixName = 'MotifMatrix',
  force = TRUE
)


##FIG 2E
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotVarDev

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 10, height = 10, ArchRProj = proj, addDOC = FALSE)

motifs <- plotVarDev$data$name[1:25]

motifs <- sub("_\\d+$", "", motifs)

print(motifs)

write.xlsx(plotVarDev$data, "ChromVAR-motifs-ranked-13.1.xlsx")

markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Clusters", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)


##FIG 2F
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP_Combined",
  imputeWeights = getImputeWeights(proj)
)


plotPDF(p, name = "Plot-UMAP-z-scores-Motif-Matrix", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)




markerRNA <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP_Combined",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-z-scores-GeneScoreMatrix", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

markerRNA <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "GeneExpressionMatrix")
markerRNA

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneExpressionMatrix", 
  name = sort(markerRNA), 
  embedding = "UMAP_Combined",
  continuousSet = "blueYellow",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-z-scores-GeneExpressionMatrix", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


##----------------------
##15.4 Identification of Positive TF-Regulators

seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Clusters")

seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


corGSM_MM <- correlateMatrices(
  ArchRProj = proj,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "LSI_Combined"
)

corGSM_MM

corGIM_MM <- correlateMatrices(
  ArchRProj = proj,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "LSI_Combined"
)

corGIM_MM

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

p1 <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

p1

corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

p2 <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

p2


##FIG 2G
plotPDF(p1, name = "Plot-positive-TF-regulators-15.4-GeneScore", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
plotPDF(p2, name = "Plot-positive-TF-regulators-15.4-GeneExpression", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

corGSM_MM

write.xlsx(corGSM_MM, "corGSM_MM-positive-TF-regulators-15.4.xlsx")
write.xlsx(corGIM_MM, "corGEM_MM-positive-TF-regulators-15.4.xlsx")


##----------------------------------------------
##FIG 3B
##15.3
kbp = 1000000
kbp_str = '1M'
proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "LSI_Combined",
  useMatrix = "GeneExpressionMatrix",
  maxDist = kbp
)

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 10,
  returnLoops = F
)

p2g

p2g[[1]]

write.table(p2g, paste0('p2g_', kbp_str,'bp.tsv'), quote=FALSE, row.names = F, col.names = T, sep = '\t')

markerGenesp2g = c("NOTCH1", "CDX2", "ELF5", "BCAM", "TCF7", "YAP1", "OVOL1", "ITGA2", "WWTR1",
                   "PLXNB2", "RUNX1", "c12orf75", "QSOX1", "RASGRF2", "JAK1", "FGFR1", "MYCN", "EBI3", "BRAF", "TBX3", "PDCD1LG2", "CD276")

p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenesp2g, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj)
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)




##---------------------------
#plot extra markers 
genes <- c('ZSCAN16', 'LELP1',  'ANXA1', 'SNHG32', 
           'MANSC4', 'TBX3', 'TACC2', 'COBLL1', 'DIO2-AS1', 'SEMA6A-AS1', 'ATP11A', 
           'PXDN', 'PFKP', 'KIF12', 'HPCAL1', 'MYO7B',  'SLC7A1', 
           'AGTR1', 'PLA2G5', 'HEYL', 'MEG3', 'CAMTA1', 'MCC', 'DOCK10', 'CIITA')



p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneExpressionMatrix", 
  name = sort(genes), 
  embedding = "UMAP_Combined",
  #continuousSet = "blueRed",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-GeneExpressionMatrix-extra", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

##---------------------------
#plot extra markers TF regulators
genes <- c("TFAP2A", "GCM1", "STAT3", "TCF21", "NFIX", "NFIA", "TP63", "PITX1", "NFIB", "ESRRG", "NFKB1", "TFAP2C", "GRHL1", "RELB", "HOXA13",
           "BACH1", "SPI1", "FLI1", "BCL11A", "ETV6", "HMGA1", "NR2F6", "STAT6", "ELF2", "ERG", "RUNX1", "TFEC", "ESRRA", "ETS1", "MITF", "ENO1",
           "CEBPA", "NR4A1", "IRF8", "EBF1", "HOXC10", "REL", "ZNF148", "FOSL2", "JUNB")


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneExpressionMatrix", 
  name = sort(genes), 
  embedding = "UMAP_Combined",
  #continuousSet = "blueRed",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-GeneExpressionMatrix-pos-TFregulators", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)



##---------------------------
#plot extra markers TF regulators
genes <- c("TFAP2A", "GCM1", "STAT3", "TCF21", "NFIX", "NFIA", "TP63", "PITX1", "NFIB", "ESRRG", "NFKB1", "TFAP2C", "GRHL1", "RELB", "HOXA13",
           "BACH1", "SPI1", "FLI1", "BCL11A", "ETV6", "HMGA1", "NR2F6", "STAT6", "ELF2", "ERG", "RUNX1", "TFEC", "ESRRA", "ETS1", "MITF", "ENO1",
           "CEBPA", "NR4A1", "IRF8", "EBF1", "HOXC10", "REL", "ZNF148", "FOSL2", "JUNB")


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = sort(genes), 
  embedding = "UMAP_Combined",
  continuousSet = "blueRed",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-GeneScoreMatrix-pos-TFregulators", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


##---------------------------
#plot extra markers TF regulators
markerMotifs <- c(
  "z:TP73_705", "z:TP63_704", "z:SMARCC1_651", "z:SNAI1_199", "z:SNAI2_161", "z:ZEB1_157",
  "z:JUN_143", "z:FOSL1_142", "z:JUNB_139", "z:FOS_137", "z:BACH1_130", "z:JDP2_125",
  "z:JUND_124", "z:FOSB_121", "z:BACH2_113", "z:FOSL2_105", "z:TCF4_97", "z:MESP2_94",
  "z:ID4_75", "z:MESP1_69", "z:ID3_38", "z:TFAP2A_5", "z:TFAP2E_4", "z:TFAP2C_3",
  "z:TFAP2B_1", "deviations:TP73_705", "deviations:TP63_704", "deviations:SMARCC1_651",
  "deviations:SNAI1_199", "deviations:SNAI2_161", "deviations:ZEB1_157", "deviations:JUN_143",
  "deviations:FOSL1_142", "deviations:JUNB_139", "deviations:FOS_137", "deviations:BACH1_130",
  "deviations:JDP2_125", "deviations:JUND_124", "deviations:FOSB_121", "deviations:BACH2_113",
  "deviations:FOSL2_105", "deviations:TCF4_97", "deviations:MESP2_94", "deviations:ID4_75",
  "deviations:MESP1_69", "deviations:ID3_38", "deviations:TFAP2A_5", "deviations:TFAP2E_4",
  "deviations:TFAP2C_3", "deviations:TFAP2B_1"
)


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP_Combined",
  imputeWeights = getImputeWeights(proj)
)


plotPDF(p, name = "Plot-UMAP-MotifMatrix-pos-TFregulators", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


##---------------------------
#plot extra markers TF regulators
genes <- c("MGAT5")


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneExpressionMatrix", 
  name = sort(genes), 
  embedding = "UMAP_Combined",
  #continuousSet = "blueRed",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-GeneExpressionMatrix-MGAT5", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)






##---------------------------
#plot extra markers TF regulators
genes <- c("TFAP2C", "TCF21", "TFAP2A", "NFIX", "PITX1", "TP63", "NFIA-AS1", "TFAP2B", "HOXA13", "HOXB13", "TP73", "GATA6", "FOXS1", "RUNX1",
           "NFIB", "TEAD1", "FOXK1", "NFKB1", "GRHL1", "IRF2", "MESP2", "GCM1", "ETV6", "IRF4", "TEAD4", "ESRRA", "BCL11A", "RELA", "ASCL1",
           "NFYA", "ZBTB3", "CEBPB", "CRX", "HOXC10", "NR2F6", "GSC", "PPARD", "KLF16", "FLI1", "ESR1", "SP2", "RUNX3", "ERF", "ERG", "RELB",
           "HOXD12", "CEBPA-DT", "BACH2", "ETS1-AS1", "HMGA1", "HOXD13", "TBX1", "TEAD3", "NR2E1", "FOXJ1")


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = sort(genes), 
  embedding = "UMAP_Combined",
  #continuousSet = "blueRed",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-GeneScoreMatrix-pos-TFregulators-GASlist", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


##---------------------------
#plot extra markers TF regulators
markerMotifs <- getFeatures(proj, useMatrix = "MotifMatrix")
markerMotifs


markerMotifs_subset <- c(
  "z:TFAP2C_3", "z:TCF21_39", "z:TFAP2A_5", "z:NFIX_738", "z:PITX1_404", "z:TP63_704", "z:NFIA_742", "z:TFAP2B_1", "z:HOXA13_418", "z:HOXB13_494",
  "z:TP73_705", "z:GATA6_387", "z:FOXS1_372", "z:RUNX1_733", "z:NFIB_741", "z:TEAD1_796", "z:FOXK1_364", "z:NFKB1_719", "z:GRHL1_391", "z:IRF2_634",
  "z:MESP2_94", "z:GCM1_390", "z:ETV6_335", "z:IRF4_632", "z:TEAD4_797", "z:ESRRA_691", "z:BCL11A_194", "z:RELA_722", "z:ASCL1_843", "z:NFYA_288",
  "z:ZBTB3_269", "z:CEBPB_140", "z:CRX_410", "z:HOXC10_548", "z:NR2F6_687", "z:GSC_470", "z:PPARD_664", "z:KLF16_205", "z:FLI1_337", "z:ESR1_661",
  "z:SP2_232", "z:RUNX3_731", "z:ERF_325", "z:ERG_339", "z:RELB_718", "z:HOXD12_524", "z:CEBPA_155", "z:BACH2_113", "z:ETS1_332", "z:HMGA1_12",
  "z:HOXD13_465", "z:TBX1_792", "z:TEAD3_795", "z:NR2E1_665", "z:FOXJ1_853",
  "deviations:TFAP2C_3", "deviations:TCF21_39", "deviations:TFAP2A_5", "deviations:NFIX_738", "deviations:PITX1_404", "deviations:TP63_704", "deviations:NFIA_742", "deviations:TFAP2B_1", "deviations:HOXA13_418", "deviations:HOXB13_494",
  "deviations:TP73_705", "deviations:GATA6_387", "deviations:FOXS1_372", "deviations:RUNX1_733", "deviations:NFIB_741", "deviations:TEAD1_796", "deviations:FOXK1_364", "deviations:NFKB1_719", "deviations:GRHL1_391", "deviations:IRF2_634",
  "deviations:MESP2_94", "deviations:GCM1_390", "deviations:ETV6_335", "deviations:IRF4_632", "deviations:TEAD4_797", "deviations:ESRRA_691", "deviations:BCL11A_194", "deviations:RELA_722", "deviations:ASCL1_843", "deviations:NFYA_288",
  "deviations:ZBTB3_269", "deviations:CEBPB_140", "deviations:CRX_410", "deviations:HOXC10_548", "deviations:NR2F6_687", "deviations:GSC_470", "deviations:PPARD_664", "deviations:KLF16_205", "deviations:FLI1_337", "deviations:ESR1_661",
  "deviations:SP2_232", "deviations:RUNX3_731", "deviations:ERF_325", "deviations:ERG_339", "deviations:RELB_718", "deviations:HOXD12_524", "deviations:CEBPA_155", "deviations:BACH2_113", "deviations:ETS1_332", "deviations:HMGA1_12",
  "deviations:HOXD13_465", "deviations:TBX1_792", "deviations:TEAD3_795", "deviations:NR2E1_665", "deviations:FOXJ1_853"
)

# Assuming that you have already defined your markerMotifs and markerMotifs_subset lists
not_found_features <- setdiff(markerMotifs_subset, markerMotifs)

# Assuming that you have already defined your markerMotifs list
nfia_features <- grep("nfia", markerMotifs, value = TRUE, ignore.case = TRUE)

# Print the features containing the string "nfia"
print(nfia_features)




p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = markerMotifs_subset, 
  embedding = "UMAP_Combined",
  imputeWeights = getImputeWeights(proj)
)


plotPDF(p, name = "Plot-UMAP-MotifMatrix-pos-TFregulators-GASlist", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)


##---------------------------
#plot extra markers TF regulators
genes <- c("TEAD4", "TEAD1", "TEAD3", "NFE2", "CEBPB", "CEBPA", "ATF4", "EVX2", "EVX1", 
           "MEOX1", "LHX4", "KLF6", "KLF7", "KLF2", "KLF14", "KLF11", "SP5", "SP6", 
           "PBX3", "NFYA", "NFYB", "NFYC", "CTCF", "CTCFL", "ASCL2", "RBPJ", "IKZF1", "IKZF4")


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneExpressionMatrix", 
  name = sort(genes), 
  embedding = "UMAP_Combined",
  #continuousSet = "blueRed",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-GeneExpressionMatrix-pos-TFregulators-moreTFs", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

##---------------------------
#plot extra markers TF regulators
genes <- c("TEAD4", "TEAD1", "TEAD3", "NFE2", "CEBPB", "CEBPA", "ATF4", "EVX2", "EVX1", 
           "MEOX1", "LHX4", "KLF6", "KLF7", "KLF2", "KLF14", "KLF11", "SP5", "SP6", 
           "PBX3", "NFYA", "NFYB", "NFYC", "CTCF", "CTCFL", "ASCL2", "RBPJ", "IKZF1", "IKZF4")


p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = sort(genes), 
  embedding = "UMAP_Combined",
  #continuousSet = "blueRed",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p, name = "Plot-UMAP-GeneScoreMatrix-pos-TFregulators-moreTFs", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)



##------------------------------------
##11.3 Pairwise Testing Between Groups

library(openxlsx)

compareClusters <- function(proj, cluster1, cluster2) {
  # Marker feature analysis
  markerTest <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = cluster1,
    bgdGroups = cluster2
  )
  
  # Motif annotation
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)
  
  # Motif enrichment analysis (upregulated)
  motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  
  # Dataframe creation (upregulated)
  df_up <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  df_up <- df_up[order(df_up$mlog10Padj, decreasing = TRUE),]
  df_up$rank <- seq_len(nrow(df_up))
  
  # Export (upregulated)
  write.xlsx(df_up, paste0(cluster1, "_", cluster2, ".xlsx"))
  
  # Plot (upregulated)
  ggUp <- ggplot(df_up, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df_up[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(P-adj) Motif Enrichment") + 
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  # Motif enrichment analysis (downregulated)
  motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )
  
  # Dataframe creation (downregulated)
  df_down <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
  df_down <- df_down[order(df_down$mlog10Padj, decreasing = TRUE),]
  df_down$rank <- seq_len(nrow(df_down))
  
  # Export (downregulated)
  write.xlsx(df_down, paste0(cluster2, "_", cluster1, ".xlsx"))
  
  # Plot (downregulated)
  ggDo <- ggplot(df_down, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df_down[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  # Save plots to PDF
  plotPDF(ggUp, ggDo, name = paste0(cluster1, "-vs-", cluster2, "-Markers-Motifs-Enriched"), width = 15, height = 15, ArchRProj = proj, addDOC = FALSE)
}


compare = c("EVT2","EVT3")
compareClusters(proj, "EVT1", compare)

compare = c("EVT1","EVT2","EVT3")
compareClusters(proj, "STB", compare)

compare1 = c("EVT1","EVT2","EVT3")
compare2 = c("vCTB1","vCTB2","vCTB3")
compareClusters(proj, compare1, compare2)

compareClusters(proj, "STB", compare1)

compareClusters(proj, "STB", compare2)

compareClusters(proj, "vCTB1", "vCTB2")

compareClusters(proj, "vCTB2", "vCTB3")

compareClusters(proj, "vCTB1", "vCTB3")

compare = c("EVT1","EVT3")
compareClusters(proj, "EVT2", compare)

compare = c("EVT1","EVT2")
compareClusters(proj, "EVT3", compare)




## --------------------------------------------------------------------------------------------------------------------------
