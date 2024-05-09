## Code for generating figures for Slide-tags placenta data

############ 0. Read in data, libraries, and functions ############
library(ArchR)
library(tidyverse)
library(Seurat)
library(wesanderson)
library(patchwork)
library(caret)
library(scatterpie)
library(CellChat)
set.seed(42)

dir.create('../output/figs')
dir.create('../output/figs/RNA')
dir.create('../output/figs/ATAC')
dir.create('../output/figs/Cellchat')

# color patterns and levels
cluster_palette_old <- setNames(c("#00A08A","#F2AD00","#F98400","#046C9A",
                              "#66C2A5","#FC8D62","#d7301f",
                              "#ECCBAE","#D69C4E","#E6A0C4"), 
                            c("Endothelial", "Erythroblasts", "EVT", "EVT-progenitor",
                              "FIB1", "FIB2", "Hofbauer cells", 
                              "STB", "STB-progenitor", "vCTB"))

cluster_palette <- setNames(c("#D42F2E",
                               "#F2AD00",
                               "darkgreen",
                               "darkseagreen2",
                               "#F28239",
                               "#FBE735",
                               "#8C9ECF",
                               "#92D3E3",
                               "deepskyblue4",
                               "#E6A0C4"), 
                            c("Endothelial", 
                              "Erythroblasts", 
                              "EVT", 
                              "EVT-progenitor",
                              "FIB1", 
                              "FIB2", 
                              "Hofbauer cells", 
                              "STB", 
                              "STB-progenitor", 
                              "vCTB"))

atac_palette <- setNames(c('#66c2a5','#fc8d62','#8da0cb','#e78ac3',
                           '#a6d854','#ffd92f','#e5c494','#b3b3b3'),c(paste0('C',1:8)))
wnn_palette <- setNames(c(RColorBrewer::brewer.pal(9, 'Set1'), 
                          RColorBrewer::brewer.pal(3, 'Set2')), 0:11)
sample_palette <- setNames(c('#af8dc3','#d8b365','#5ab4ac'), c("8w_2d", "9w_1d", "11w_1d"))
#cluster_levels <- c("Endothelial", "Erythroblasts", "EVT", "EVT-progenitor", "FIB1", "FIB2", "Hofbauer cells", "NK/T cells", "STB", "STB-precursor", "vCTB")
cluster_levels <- c("Endothelial", "Erythroblasts", "EVT", "EVT-progenitor", "FIB1", "FIB2", 
                    "Hofbauer cells", "STB", "STB-progenitor", "vCTB")
sample_qc_palette <- setNames(c('#C06DAA','#FBE735','#8C9ECF'), c("8w_2d", "9w_1d", "11w_1d"))
stage_levels <- c("8w_2d","9w_1d","11w_1d")

# functions for plots
# add coordinates to seurat objects as dims
add_spatial_metadata <- function(seu, assay = "RNA", suffix = NULL) {
  emb <- seu@meta.data[,c("x", "y")]
  colnames(emb) <- c("s_1","s_2")
  seu[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), 
                                           key = "s_", assay = assay)
  
  return(seu)
}

# plot clusters or expressions of genes spatially
# based on Seurat function
plot_cluster_feature <- function(seu, reduction, label = FALSE, group_by = NULL,
                                 split_by = NULL, features = NULL, ...) {
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
      coord_fixed(ratio = 1) +
      theme(strip.text = element_text(size = 8, face = "bold"),
            line = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      labs(title = sample)
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
# based on custom ggplots themes
cluster_distribution <- function(data, title = "UMAP of RNA", size = 0.8, alpha = 1,
                                 Dim1 = 'Dim1', Dim2 = 'Dim2', cluster = 'cluster',
                                 x_labs = "UMAP1", y_labs = "UMAP2",
                                 color_palette = cluster_palette, scale_x1 = "x1", scale_x2 = "x2"){
  ggplot(data = data, aes_string(x = Dim1, y = Dim2, color = cluster)) +
    geom_point(alpha = alpha, size = size) +
    labs(title = title, x = x_labs, y = y_labs, color = "Cluster") +
    scale_color_manual(values = color_palette) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(), 
          legend.key.size = unit(.5, "cm"),
          panel.border = element_rect(linetype = "solid", color = "black", size = 1.2),
          plot.title = element_text(hjust = 0.5, size = 18)) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    coord_fixed() +
    geom_segment(aes(x = scale_x1, xend = scale_x2, y = 100, yend = 100), color = "black", size = 1)
}

cluster_umap <- function(data, title = "UMAP of RNA", size = 0.8, alpha = 1,
                                 Dim1 = 'Dim1', Dim2 = 'Dim2', cluster = 'cluster',
                                 x_labs = "UMAP1", y_labs = "UMAP2",
                                 color_palette = cluster_palette){
  ggplot(data = data, aes_string(x = Dim1, y = Dim2, color = cluster)) +
    geom_point(alpha = alpha, size = size) +
    labs(title = title, x = x_labs, y = y_labs, color = "Cluster") +
    scale_color_manual(values = color_palette) +
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
}
# plot Ligand-Receptor genes and cells spatially
plot_LR <- function(seu, sender_cell, receiver_cell, sender_gene, 
                             receiver_gene, idents_col, assay_use) {
  # Helper function to create data.frame for a given cell and gene
  create_df <- function(cell, gene) {
    indices <- which(seu@meta.data[[idents_col]] == cell)
    data.frame(X = seu@meta.data$x_um[indices],
               Y = seu@meta.data$y_um[indices],
               gene_expr = seu@assays[[assay_use]]@data[gene, indices])
  }
  
  # Helper function to create a ggplot for a given data.frame
  create_plot <- function(df, title) {
    ggplot(df, aes(X, Y, color = gene_expr)) +   
      geom_point(size = 1) + labs(title = title) +
      theme_bw() +
      scale_color_gradientn(colours = inlmisc::GetColors(256, scheme = "sunset")) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  }
  
  sender <- create_df(sender_cell, sender_gene)
  receiver <- create_df(receiver_cell, receiver_gene)
  plot_sender <- create_plot(sender, paste0(sender_cell,':',sender_gene))
  plot_receiver <- create_plot(receiver, paste0(receiver_cell,':',receiver_gene))
  cowplot::plot_grid(plot_sender, plot_receiver)
}



# contents of figs
# 1. Plot clusters and genes based on RNA features (UMAP and Coordinates)
# 2. Plot clusters and TFs based on ATAC features (UMAP and Coordinates)
# 3. Plot distribution of clusters based on WNN features (UMAP and Coordinates)
# 4. Plot Cell interaction counts and L-R expression spatially


############ 1. Plot clusters and genes based on RNA features (UMAP and Coordinates)

int <- readRDS('../data/Integrated_RNA_3samples.rds')
rna_umap <- as.data.frame(int@reductions$umap@cell.embeddings)
colnames(rna_umap) <- c('Dim1','Dim2')
rna_umap$cluster <- factor(int$Clusters.2, levels = cluster_levels)
rna_umap$sample <- int$stage
rna_umap$seurat_cluster <- int$seurat_clusters

# 1.1. cluster proportion bar plots

cluster_proportion <- table(int$stage, int$Clusters.2)
cluster_proportion <- prop.table(cluster_proportion, margin = 1)
cp_melted <- as.data.frame(as.table(cluster_proportion))
colnames(cp_melted) <- c("Sample", "Cluster", "Proportion")
cp_melted$Sample <- factor(cp_melted$Sample, levels = stage_levels)
cp_melted$Cluster <- factor(cp_melted$Cluster, levels = cluster_levels)
pdf('../output/figs/RNA/Cluster_Proportion.pdf', width = 8, height = 4)
ggplot(cp_melted, aes(x = Sample, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_palette) + 
  theme_bw() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = "right",
        panel.border = element_rect(linetype = "solid", color = "black", size = 1.2))
dev.off()

# 1.2. UMAP of clusters
pdf('../output/figs/RNA/RNA_UMAP_Cluster.pdf', width = 7.5, height = 6)
cluster_umap(rna_umap, title = "Clusters (RNA)")
dev.off()

# 1.3. UMAP of samples
pdf('../output/figs/RNA/RNA_UMAP_sample.pdf', width = 7.5, height = 6)
cluster_umap(rna_umap, title = "Samples (RNA)", cluster = 'sample', alpha = 0.95,
                     color_palette = sample_palette)
dev.off()

# 1.4. UMAP of seurat clusters
pdf('../output/figs/RNA/RNA_UMAP_seurat_clusters.pdf', width = 7.5, height = 6)
cluster_umap(rna_umap, title = "Samples (RNA)", cluster = 'seurat_cluster', alpha = 0.95,
                     color_palette = unname(alphabet(n = length(unique(rna_umap$seurat_cluster)))))
dev.off()

scale_bar_xa <- c(5000,5000,5500)
scale_bar_xb <- c(5500,5500,6000)

# 1.5. spatial distribution of Clusters
st_cluster <- data.frame('x_um'=int$x, 'y_um'=int$y,'sample'=int$stage,
                         'cluster'=factor(int$Clusters.2, levels=cluster_levels))
rownames(st_cluster) <- rownames(int@meta.data)
spatial_list <- split(st_cluster, f=st_cluster$sample)
pdf('../output/figs/RNA/Spatial_Cluster.pdf', width = 7.5, height = 6)
for (n in seq_along(names(spatial_list))) {
  p <- cluster_distribution(spatial_list[[names(spatial_list)[n]]], title = paste0("Spatial Clusters of ",names(spatial_list)[n]),
                            Dim1 = 'x_um', Dim2 = 'y_um', size = 2,
                            x_labs = 'x_um', y_labs = 'y_um', scale_x1 = scale_bar_xa[n], scale_x2 = scale_bar_xb[n]) ; print(p)
}
dev.off()

# 1.6. plot expression of marker genes based on coordinats and UMAP
int <- add_spatial_metadata(int)
marker_list <- list(
  vCTB = c('PEG10', 'PAGE4', 'CDH1'),
  vCTBp = c('MKI67', 'TOP2A', 'PCNA'),
  vCTB1 = c('ARL15', 'TBL1X', 'SIK3'),
  vCTB2 = c('IFI6', 'EFEMP1', 'SMAGP'),
  vCTB3 = c('MSI2', 'LRP5', 'AHCYL2'),
  EVT = c('ITGA5', 'HLA-G', 'DIO2'),
  EVT1 = c('EGFR', 'FBN2', 'UTRN'),
  EVT2 = c('CCNE1', 'HAPLN3', 'LY6E'),
  STB = c('ERVFRD-1', 'CYP19A1', 'IGSF5'),
  Endo = c('PECAM1', 'EGFL7', 'KDR'),
  Macrophages = c('SPP1', 'CD14', 'CD83'),
  Maternal_MAC = c('CD74', 'HLA-DRA', 'LYZ'),
  HBC = c('F13A1', 'LYVE1', 'ADAMTS17'),
  Myeloid_unk = c('FGF13', 'ALPK3'),
  FIB = c('COL3A1', 'SOX5', 'LAMA2'),
  FIB1 = c('PDGFRB', 'PLCB1', 'AGTR1'),
  FIB2 = c('PDGFRA', 'COL15A1', 'CXCL14')
)
genes <- unlist(marker_list)
int_list <- SplitObject(int, split.by = 'stage')
for (sample in names(int_list)) {
  seu <- int_list[[sample]]
  pdf(paste0('../output/figs/RNA/Spatial_Genes_RNA_',sample,'.pdf'), width = 7.5, height = 6)
  for (i in genes) {print(plot_cluster_feature(seu, 'spatial', features = i, pt.size = 1.5))}
  dev.off()
}
pdf('../output/figs/RNA/UMAP_Genes_RNA.pdf', width = 7.5, height = 6)
for (i in genes) {print(plot_cluster_feature(int, 'umap', features = i, pt.size = 0.8))}
dev.off()

# selected marker genes
# AGTR1 		FIB1
# F13A1 		Hofbauer cells 
# ERVFRD−1 	STB
# HAPLN3   	EVT2
# DIO2   		EVT
# MKI67 		EVT−progenitor 
# PAGE4   	vCTB
# SLC4A1 		Erythroblasts 
# KDR   		Endothelial 
# PDGFRA 	  FIB2
genes = c('AGTR1','F13A1','ERVFRD-1','HAPLN3','DIO2','MKI67','PAGE4','SLC4A1','KDR','PDGFRA')
int_list <- SplitObject(int, split.by = 'stage')
for (sample in names(int_list)) {
  seu <- int_list[[sample]]
  plot_list <- lapply(genes, function(i) plot_cluster_feature(seu, 'spatial', features = i, pt.size = 1.5))
  combined_plot <- wrap_plots(plot_list, ncol = 5)
  pdf(paste0('../output/figs/RNA/Spatial_Selected_Genes_RNA_',sample,'.pdf'), width = 15, height = 6)
  print(combined_plot)
  dev.off()
}
plot_list <- lapply(genes, function(i) plot_cluster_feature(int, 'umap', features = i, pt.size = .5))
combined_plot <- wrap_plots(plot_list, ncol = 5)
pdf('../output/figs/RNA/UMAP_Selected_Genes_RNA.pdf', width = 15, height = 6)
print(combined_plot)
dev.off()





############ 2. Plot clusters and TFs based on ATAC features (UMAP and Coordinates)

proj <- readRDS('../data/Integrated_ATAC_Proj.rds')
atac_umap <- as.data.frame(proj@embeddings$UMAP@listData$df)
atac_meta <- data.frame('cellname'=proj$cellNames,'cluster'=proj$Clusters)
rownames(atac_meta) <- atac_meta$cellname
atac_meta <- atac_meta[rownames(atac_umap),]
atac_umap$cluster <- atac_meta$cluster
colnames(atac_umap)[1:2] <- c('Dim1','Dim2')
rownames(atac_umap) <- gsub("combined#", "", rownames(atac_umap))
updated_barcodes <- sapply(rownames(atac_umap), function(bc) {
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
atac_umap$cell <- updated_barcodes
meta <- data.frame(cells = colnames(int), sample = int$stage, clusters = int@meta.data$Clusters.2)
atac_umap$sample <- meta[atac_umap$cell, ]$sample

# 1. UMAP of clusters 
atac_umap$cluster <- factor(atac_umap$cluster, levels = paste0('C',1:8))
pdf('../output/figs/ATAC/ATAC_UMAP_Cluster_ATAC_based.pdf', width = 7.5, height = 6)
cluster_umap(atac_umap, title = "Clusters (ATAC)", 
                     color_palette = atac_palette)
dev.off()

# 2. UMAP of samples
pdf('../output/figs/ATAC/ATAC_UMAP_sample.pdf', width = 7.5, height = 6)
cluster_umap(atac_umap, title = "Samples (ATAC)", 
                     cluster = 'sample', alpha = 0.95,
                     color_palette = sample_palette)
dev.off()

# 3. spatial
## ATAC-only cells
coord = data.frame('x_um'=int$x, 'y_um'=int$y)
rownames(atac_umap) <- atac_umap$cell
atac_umap <- cbind(atac_umap, coord[rownames(atac_umap),])
spatial_list <- split(atac_umap, f=atac_umap$sample)
pdf('../output/figs/ATAC/ATAC_Spatial_Cluster.pdf', width = 7.5, height = 6)
for (n in names(spatial_list)) {
  p <- cluster_distribution(spatial_list[[n]], title = paste0("Spatial Clusters of ",n),
                            Dim1 = 'x_um', Dim2 = 'y_um', size = 2, 
                            color_palette = atac_palette,
                            x_labs = 'x_um', y_labs = 'y_um', scale_x1 = scale_bar_xa[n], scale_x2 = scale_bar_xb[n]) ; print(p)
}
dev.off()


# 5. spatial visualization of custom motifs using deviations score (from chromVAR) 
dt = readRDS('../data/int_MotifMatrix_Z.rds')
dt = cbind(st_cluster[rownames(st_cluster), c('sample', 'x_um', 'y_um'), drop=F], dt[rownames(st_cluster),])
sample_list <- split(dt, f=dt$sample)
i <- 0
scale_bar_x <- c(0,0,500)
for (sample in names(sample_list)) {
  i <- i + 1
  data = sample_list[[sample]]
  data$sample <- NULL
  # change to x and y if you want only ATAC cells
  data_melt <- reshape2::melt(data, id.vars = c("x_um", "y_um"), variable.name = "Motif", value.name = "Deviation_Score")
  motifs <- unique(data_melt$Motif)
  pdf(paste0('../output/figs/ATAC/int_motifs_deviation_score_spatial_',sample,'.pdf'), width = 6, height = 5)
  for(motif in motifs){
    motif_data <- subset(data_melt, Motif == motif)
    motif_data <- motif_data[order(motif_data$Deviation_Score),]
    motif_data <- motif_data %>% arrange(desc(is.na(motif_data)))
    p <- ggplot((motif_data), aes(x = x_um, y = y_um, fill = Deviation_Score)) + 
      # also need to add x and y instead of x_um and y_um to get only ATAC cells
      geom_point(size = 2, shape = 21, color = "grey") + 
      scale_fill_continuous_diverging("Blue-Red 3", na.value = "grey") +
      #scale_color_gradientn(colors = inlmisc::GetColors(256, scheme = "sunset")) + 
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
    coord_fixed()  +
    geom_segment(aes(x = 5000+scale_bar_x[i], xend = 5500+scale_bar_x[i], y = 100, yend = 100), color = "black", size = 1)
    print(p)
  }
  dev.off()
}


# 6. spatial visualization of custom motifs using GeneScoreMatrix
dt = readRDS('../data/int_GeneScoreMatrix.rds')
dt = cbind(st_cluster[rownames(st_cluster), c('sample', 'x_um', 'y_um'), drop=F], dt[rownames(st_cluster),])

sample_list <- split(dt, f=dt$sample)
for (sample in names(sample_list)) {
  data = sample_list[[sample]]
  data$sample <- NULL
  # to get ATAC-only, change this to xy and merge with rna_umap instead
  data_melt <- reshape2::melt(data, id.vars = c("x_um", "y_um"), variable.name = "Motif", value.name = "Gene_Score")
  motifs <- unique(data_melt$Motif)
  pdf(paste0('../output/figs/ATAC/int_motifs_gene_score_spatial_',sample,'.pdf'), width = 6, height = 5)
  for(motif in motifs){
    motif_data <- subset(data_melt, Motif == motif)
    motif_data <- motif_data[order(motif_data$Gene_Score),]
    motif_data <- motif_data %>% arrange(desc(is.na(motif_data)))
    p <- ggplot((motif_data), aes(x = x_um, y = y_um, fill = Gene_Score)) + 
      # also need to add x and y instead of x_um and y_um to get only ATAC cells
      geom_point(size = 2, shape = 21, color = "grey") + 
      scale_fill_viridis_c(option = "plasma", direction = -1, na.value = "grey") +
      #scale_color_gradientn(colors = inlmisc::GetColors(256, scheme = "sunset")) + 
      theme_bw() + 
      labs(title = motif) +
      theme(legend.position = "right") +
      guides(color = guide_colorbar(title = "Gene \nScore")) +
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
    coord_fixed()  +
    geom_segment(aes(x = 5000+scale_bar_x[i], xend = 5500+scale_bar_x[i], y = 100, yend = 100), color = "black", size = 1)
    print(p)
  }
  dev.off()
}



########### Quality control plots

plot_violin_qc <- function(data, y_variables, new_names = NULL) {
  
  # A function factory for minor log breaks
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
  
  for (y_var in y_variables) {
    # Convert factor variables to numeric if needed
    if (is.factor(data[[y_var]])) {
      data[[y_var]] <- as.numeric(as.character(data[[y_var]]))
    }
    
    # Check if the conversion was successful
    if (is.numeric(data[[y_var]])) {
      # Check for NA values in the variable
      if (any(is.na(data[[y_var]]))) {
        message(paste("Skipping variable:", y_var, " - contains NA values"))
        next
      }
      
      # Adding a small constant to prevent log(0)
      #data[[y_var]] <- data[[y_var]] + 1
      
      # Calculate the median for each group
      median_values <- aggregate(data[[y_var]], by = list(data$group), FUN = median)
      names(median_values) <- c("group", "median")
      
      # Decide whether to use log scale
      use_log_scale <- all(data[[y_var]] > 1) && max(data[[y_var]]) / min(data[[y_var]]) > 3  # Adjust the threshold as needed
      
    # Get new variable name if provided, otherwise use the original name
    if (!is.null(new_names) && y_var %in% names(new_names)) {
      y_var_name <- new_names[[y_var]]
    } else {
      message(paste(y_var, " - not in new_names"))
      y_var_name <- y_var
    }
      
      gg <- ggplot(data, aes_string(x = "group", y = y_var, fill = "stage")) +
        geom_violin(colour = "black") +
      #geom_text(data = median_values, aes(label = paste("Median:", round(median, 2)), x = group, y = -Inf), vjust = -1.5) +
        theme_minimal() +
        scale_fill_manual(values = sample_qc_palette, name = NULL, labels = c("W8-2", "W9", "W11")) +
        scale_color_manual(values = c("black")) +
        theme(
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
        ggtitle(paste(y_var_name)) + 
        #scale_x_discrete(labels = median_values$median) +
        guides(color = "none")
      
      # Conditionally add log scale
      if (use_log_scale) {
        gg <- gg + scale_y_log10(breaks = breaks, minor_breaks = minor_breaks, labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        annotation_logticks(sides = "l")
      }
      
      # Assign median values as x-axis labels
if (!is.null(median_values$median)) {
  gg <- gg + scale_x_discrete(labels = round(median_values$median, 2))
}
      
      # save plot
      ggsave(paste0("../output/figs/qc_metrics_violin_", y_var, ".pdf"),
             plot = gg, device = "pdf", dpi = 320, width = 3, height = 4)
    } else {
      message(paste("Skipping variable:", y_var, " - unable to convert to numeric"))
      next
    }
  }
}

meta <- (int@meta.data)
meta$group <- factor(meta$group, levels = c("JS40", "JS35", "JS36"))  
new_names <- list("nCount_RNA" = "nUMI (gene expression)",
                  "nFeature_RNA" = "nGenes",
                  "percent.mt" = "% Mitochondrial transcripts",
                  "SB_UMI_total" = "nUMI (spatial barcode)",
                  "SNR" = "proportion signal spatial barcodes")
y_variables_list <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "SB_UMI_total", "SNR")
plot_violin_qc(meta, y_variables_list, new_names)

df_coords_all <- data.frame()
for(i in seq_along(run)){
  df_coords <- read.table(paste0("../output/", run[i], "/coords_", run[i], ".csv"), header = T)
  rownames(df_coords) <- paste0(run_name[i], "_", df_coords$cb)
  df_coords <- df_coords[,c("clusters", "x", "y", "UMI_cluster", "UMI_noise", "cb")]
  df_coords_all <- rbind(df_coords_all, df_coords)
}
df_coords_all$group <- sub("_.*", "", rownames(df_coords_all))
 df_coords_all$stage <- NA
 df_coords_all$stage[which(df_coords_all$group == "JS40")] <- "8w_2d"
 df_coords_all$stage[which(df_coords_all$group == "JS35")] <- "9w_1d"
 df_coords_all$stage[which(df_coords_all$group == "JS36")] <- "11w_1d"

meta_mapped <- df_coords_all %>%
               group_by(group) %>%
              summarize(proportion_clusters_1 = mean(clusters == 1, na.rm = TRUE))
meta_mapped$stage <- c("9w_1d", "11w_1d", "8w_2d")
meta_mapped$pc <- meta_mapped$proportion_clusters_1 * 100
meta_mapped$stage <- factor(meta_mapped$stage, levels = c("8w_2d", "9w_1d", "11w_1d"))
meta_mapped <- meta_mapped %>% arrange(stage)
  

# Convert 'stage' to a factor with custom order
meta_mapped$stage <- factor(meta_mapped$stage, levels = c("8w_2d", "9w_1d", "11w_1d"))

p1 <- ggplot(meta_mapped, aes(fill = stage, y = pc, x = stage)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
        scale_fill_manual(values = sample_qc_palette, name = NULL, labels = c("W8-2", "W9", "W11")) +
        theme(
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) + 
  scale_x_discrete(labels = round(meta_mapped$pc, 2)) +
        ggtitle("% cells spatially mapped") +
        guides(color = "none")

print(p1)

ggsave(paste0("../output/figs/qc_metrics_mapped_cells.pdf"),
             plot = p1, device = "pdf", dpi = 320, width = 3, height = 4)

plot_violin_qc <- function(data, y_variables, new_names = NULL) {
  
  # A function factory for minor log breaks
  breaks <- 10^(-10:10)
  minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
  
  # Reorder levels of group variable
  data$group <- factor(data$group, levels = c("JS40", "JS35", "JS36"))
  
  for (y_var in y_variables) {
    # Convert factor variables to numeric if needed
    if (is.factor(data[[y_var]])) {
      data[[y_var]] <- as.numeric(as.character(data[[y_var]]))
    }
    
    # Check if the conversion was successful
    if (is.numeric(data[[y_var]])) {
      # Check for NA values in the variable
      if (any(is.na(data[[y_var]]))) {
        message(paste("Skipping variable:", y_var, " - contains NA values"))
        next
      }
      
      # Adding a small constant to prevent log(0)
      #data[[y_var]] <- data[[y_var]] + 1
      
      # Calculate the median for each group
      median_values <- aggregate(data[[y_var]], by = list(data$group), FUN = median)
      names(median_values) <- c("group", "median")
      
      # Decide whether to use log scale
      use_log_scale <- all(data[[y_var]] > 1) && max(data[[y_var]]) / min(data[[y_var]]) > 3  # Adjust the threshold as needed
      
      # Get new variable name if provided, otherwise use the original name
      if (!is.null(new_names) && y_var %in% names(new_names)) {
        y_var_name <- new_names[[y_var]]
      } else {
        message(paste(y_var, " - not in new_names"))
        y_var_name <- y_var
      }
      
      #data$group <- factor(data$group, levels = c("JS40", "JS35", "JS36"))
      
      gg <- ggplot(data, aes_string(x = "group", y = y_var, fill = "stage")) +
        geom_violin(colour = "black") +
        theme_minimal() +
        scale_fill_manual(values = sample_qc_palette, name = NULL, labels = c("W8-2", "W9", "W11")) +
        scale_color_manual(values = c("black")) +
        theme(
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          axis.line = element_line(color = "black"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
        ggtitle(paste(y_var_name)) + 
        guides(color = "none")
      
      # Conditionally add log scale
      if (use_log_scale) {
        gg <- gg + scale_y_log10(breaks = breaks, minor_breaks = minor_breaks, labels = scales::trans_format("log10", scales::math_format(10^.x))) +
          annotation_logticks(sides = "l")
      }
      
      # Assign median values as x-axis labels
      if (!is.null(median_values$median)) {
        gg <- gg + scale_x_discrete(labels = round(median_values$median, 2))
      }
      
      # save plot
      ggsave(paste0("../output/figs/qc_metrics_violin_", y_var, ".pdf"),
             plot = gg, device = "pdf", dpi = 320, width = 3, height = 4)
    } else {
      message(paste("Skipping variable:", y_var, " - unable to convert to numeric"))
      next
    }
  }
}


# plot QC metrics from ATAC
meta <- as.data.frame(proj@cellColData)
meta$stage <- NA
meta$stage[which(meta$group == "JS40")] <- "8w_2d"
meta$stage[which(meta$group == "JS35")] <- "9w_1d"
meta$stage[which(meta$group == "JS36")] <- "11w_1d"
meta$group <- factor(meta$group, levels = c("JS40", "JS35", "JS36"))
new_names <- list("TSSEnrichment" = "TSS Enrichment",
                  "nFrags" = "nFragments",
                  "ReadsInPeaks" = "Reads in Peaks",
                  "FRIP" = "Fraction of reads in peaks")
y_variables_list <- c("TSSEnrichment", "nFrags", "ReadsInPeaks", "FRIP")
meta$group <- factor(meta$group, levels = c("JS40", "JS35", "JS36"))
meta$stage <- factor(meta$stage, levels = c("8w_2d", "9w_1d", "11w_1d"))
plot_violin_qc(meta, y_variables_list, new_names)

print("cells in ArchR:")
table(proj@cellColData$group)
print("cells in Seurat:")
table(int@meta.data$group)

# Make QC table

df_qc_summary <- merge(as.data.frame(table(int@meta.data$group)), as.data.frame(table(proj@cellColData$group)), by = "Var1")
colnames(df_qc_summary) <- c("ID", "Number of Cells Total", "Number of Cells Pass QC ATAC")
df_qc_summary$Donor <- as.character(df_qc_summary$ID) %>% {
  .[. %in% c('JS35')] <- 'W9'
  .[. %in% c('JS36')] <- 'W11'
  .[. %in% c('JS40')] <- 'W8-2'
  ;.}
df_qc_summary
write.csv(df_qc_summary, "../output/figs/qc_table_summary.csv")

qc_rna <- as.data.frame(int@meta.data[,c("nCount_RNA", "nFeature_RNA", "group", "percent.mt", "SB_UMI_total", "Clusters.2")])
qc_atac <- as.data.frame(proj@cellColData[ ,c("TSSEnrichment", "nFrags", "FRIP", "group")])
rownames(qc_atac) <- gsub("combined#", "", rownames(qc_atac))
updated_barcodes <- sapply(rownames(qc_atac), function(bc) {
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
qc_atac$cell <- updated_barcodes

df_qc_full <- merge(qc_rna, qc_atac, by = "group", by.x = 0, by.y = "cell", all.x = T)
df_qc_full <- df_qc_full %>%
  group_by(group.x, Clusters.2) %>%
  summarize(
    count_RNA = sum(!is.na(nCount_RNA)),
    count_ATAC = sum(!is.na(TSSEnrichment)),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    median_SB_UMI_total = median(SB_UMI_total, na.rm = TRUE),
    median_TSSEnrichment = median(TSSEnrichment, na.rm = TRUE),
    median_nFrags = median(nFrags, na.rm = TRUE),
    median_FRIP = median(FRIP, na.rm = TRUE)
  ) %>%
  ungroup()
df_qc_full
write.csv(df_qc_full, "../output/figs/qc_table_full.csv")

# for reviewers, a table comparing number of cells in dataset
df_comparison <- merge(as.data.frame(table(int_old@meta.data$group)), as.data.frame(table(int@meta.data$group)), by = "Var1")
colnames(df_comparison) <- c("ID", "old count", "post-revisions count")
totals <- data.frame("ID" = "total", "old count" = sum(df_comparison$`old count`), "post-revisions count" = sum(df_comparison$`post-revisions count`))
colnames(totals) <- c("ID", "old count", "post-revisions count")
df_comparison <- rbind(df_comparison, totals)
df_comparison$pc_change <- ((df_comparison$`post-revisions count`/df_comparison$`old count`)*100)
df_comparison$pc_improvment <- df_comparison$pc_change-100
df_comparison








############ 4. Plot Cell interaction counts and L-R expression spatially

# 4.1. cell interaction counts
cellchat.35 <- readRDS('../output/cci/cellchat_300um_35.rds')
cellchat.36 <- readRDS('../output/cci/cellchat_300um_36.rds')
cellchat.40 <- readRDS('../output/cci/cellchat_300um_40.rds')
cellchat.sc <- readRDS('../data/2023-08-25_cellchat.rds')
cellchat_list = list('35'=cellchat.35, '36'=cellchat.36, '40'=cellchat.40)

sample_ns <- c('40','35','36')
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

pdf('../output/figs/Cellchat/cellchat_interaction_counts_int.pdf', width = 8, height = 6)
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

cellchat_list = lapply(names(cellchat_list), function(x){
  sample = cellchat_list[[x]]
  lr_sig = sample@LR[["LRsig"]]
  sc_lr_sig = cellchat.sc@LR[["LRsig"]]$interaction_name
  com_lr = intersect(lr_sig$interaction_name, sc_lr_sig)
  lr_sig = lr_sig[lr_sig$interaction_name %in% com_lr,]
  sample@LR[["LRsig"]] = lr_sig
  return(sample)
}) %>% {names(.) <- c('35','36','40');.}

com_st_lr = lapply(names(cellchat_list), function(x){
  sample = cellchat_list[[x]]
  lr_sig = sample@LR[["LRsig"]]
  sc_lr_sig = cellchat.sc@LR[["LRsig"]]$interaction_name
  com_lr = intersect(lr_sig$interaction_name, sc_lr_sig)
  lr_sig = lr_sig[lr_sig$interaction_name %in% com_lr,]
  lr_sig$sample <- x
  return(lr_sig)
}) %>% do.call('rbind', .)

pdf('../output/figs/Cellchat/cellchat_LR_3sample.pdf', height = 5, width = 8)
for (n in names(cellchat_list)) {
  cellchat = cellchat_list[[n]]
  p <- netVisual_bubble(cellchat, sources.use = 1:10, 
                        targets.use = 1:10, 
                        signaling = com_st_lr_filter$pathway_name, 
                        remove.isolate = TRUE) +
    labs(title = n)
  print(p)
}
dev.off()
