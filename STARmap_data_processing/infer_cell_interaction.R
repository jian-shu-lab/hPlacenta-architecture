custom_lib_path <- "/home/jinmr2/R/x86_64-pc-linux-gnu-library/4.0"
.libPaths(custom_lib_path)

library(future)
library(CellChat)
library(patchwork)
library(Seurat)
library(anndata)
library(reticulate)
library(scatterpie)
library(ggplot2)
library(tidyr)
use_python('/home/kjin2/pyenv/run_starfish/bin/python3.8')
options(stringsAsFactors = FALSE)

# params
setwd('/home/jinmr2/sample_integration/three_samples/cell_chat_v2/')
save_dir = '/home/jinmr2/sample_integration/three_samples/cell_chat_v2/'

# sc data
cellchat.sc <- readRDS('/home/jinmr2/sample_integration/three_samples/cell_chat/2023-08-25_cellchat.rds')
subset_pathway_name <- cellchat.sc@LR$LRsig$pathway_name # didn't use
subset_lr <- cellchat.sc@LR$LRsig$interaction_name


seu.anno = read_h5ad('../processed_clustermap_all_v2.h5ad')
df_anno = seu.anno$obs
cluster_levels = levels(df_anno$predicted.celltype)

distance = 10
input_JS_list = c('JS40', 'JS35', 'JS36','JS34')
for (input_JS in input_JS_list) {
  root_folder = '/home/jinmr2/sample_integration/three_samples/imputation/v5_seurat_all_genes/output/'
  seu.st = read_h5ad(paste0(root_folder, 'adata_imput_expr_', input_JS, '_100_final_integer.h5ad'))
  seu.st = CreateSeuratObject(counts = t(seu.st$X), meta.data = seu.st$obs)
  df_anno_sample = df_anno[colnames(seu.st),]
  seu.st@meta.data = df_anno_sample
  seu.st <- NormalizeData(seu.st, verbose = FALSE)

  # seu = seu.st[,seu.st$sample == input_JS]
  seu = seu.st
  assay_use = 'RNA'
  idents_col = 'predicted.celltype'
  data.input = as.matrix(seu@assays[[assay_use]]@data)
  meta = seu@meta.data[,c(idents_col), drop = F] 
  
  colnames(meta) <- "labels"
  spatial.locs = data.frame(seu$x, seu$y)
  scale.factors = list(spot.diameter = 1, spot = 1)
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                             datatype = "spatial", coordinates = spatial.locs,
                             scale.factors = scale.factors)
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat@meta$labels = droplevels(cellchat@meta$labels, exclude = setdiff(levels(cellchat@meta$labels),unique(cellchat@meta$labels)))
  cellchat@idents = factor(cellchat@idents, levels = levels(cellchat@meta$labels))
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0.1, 
                                distance.use = TRUE, interaction.length = distance, 
                                scale.distance = 1)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, file = paste0(save_dir, '/cellchat_net_', distance, 'um_', input_JS, '.csv'))
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  future::plan("sequential")
  saveRDS(cellchat, file = paste0(save_dir, '/cellchat_', distance, 'um_', input_JS, '.rds'))
}

# figs
cellchat.34 <- readRDS(paste0('cellchat_', distance, 'um_JS34.rds'))
cellchat.35 <- readRDS(paste0('cellchat_', distance, 'um_JS35.rds'))
cellchat.36 <- readRDS(paste0('cellchat_', distance, 'um_JS36.rds'))
cellchat.40 <- readRDS(paste0('cellchat_', distance, 'um_JS40.rds'))
cellchat.sc <- readRDS('/home/jinmr2/sample_integration/three_samples/cell_chat/2023-08-25_cellchat.rds')
cellchat_list = list('js34'=cellchat.34, 'js35'=cellchat.35, 'js36'=cellchat.36, 'js40'=cellchat.40)


# plot cell interacting of integrated counts for 3 samples
sample_ns <- c('js34','js35','js36','js40')
stage <- c("7w_6d","9w_1d","11w_1d","8w_2d")
data_mat <- lapply(1:4, function(i){
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
write.csv(summary_data, paste0('summary_interaction_counts_',distance,'_um.csv'))

pdf(paste0('figs/0.cellchat_LR_4sample_counts_', distance,'um.pdf'), height = 8, width = 8)
ggplot() + 
  geom_scatterpie(data=summary_data, 
                  aes(x=source_num, y=target_num, r=r*0.05), 
                  cols=c("7w_6d","9w_1d","11w_1d","8w_2d"),
                  legend_name = 'sample',) + 
  scale_x_continuous(breaks=c(1:12), labels=cluster_levels) + 
  scale_y_continuous(breaks=c(1:12), labels=cluster_levels) + 
  labs(x="source", y="target") + 
  # scale_fill_manual(values = sample_palette) +
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
}) %>% {names(.) <- c('JS34','JS35','JS36','JS40');.}

pdf(paste0('figs/1.cellchat_LR_4sample_individual_', distance,'um.pdf'), height = 20, width = 14)
for (n in names(cellchat_list)) {
  cellchat = cellchat_list[[n]]
  p <- netVisual_bubble(cellchat, sources.use = 1:12, 
                        targets.use = 1:12, 
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

pdf(paste0('figs/2.cellchat_LR_4sample_shared_', distance,'um.pdf'), height = 20, width = 14)
for (n in names(cellchat_list)) {
  cellchat = cellchat_list[[n]]
  p <- netVisual_bubble(cellchat, sources.use = 1:12, 
                        targets.use = 1:12, 
                        signaling = unique(result$pathway_name), 
                        remove.isolate = TRUE) +
    labs(title = n)
  print(p)
}
dev.off()

