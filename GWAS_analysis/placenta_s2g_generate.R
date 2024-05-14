library(reticulate)
library(anndata)
library(SparseM)
library(rhdf5)
library(testit)
library(data.table)


h5ad_file_RNA = "~/data/hplacenta_gene_matrix.h5ad"
h5ad_file_ATAC = "~/data/hplacenta_peak_matrix.h5ad"
featuredata_RNA = h5read(h5ad_file_RNA, name = "var")[['_index']]
featuredata_ATAC = h5read(h5ad_file_ATAC, name = "var")[['_index']]
h5read(h5ad_file_ATAC, name = "var")[['__categories']]
h5ls(h5ad_file_ATAC)

metadata = h5read(h5ad_file_ATAC, name = "obs")

chrvals_ATAC = as.character(sapply(featuredata_ATAC, function(x) return(strsplit(x, ":")[[1]][1])))
pos_start_ATAC = as.numeric(sapply(sapply(featuredata_ATAC, function(x) return(strsplit(x, ":")[[1]][2])),
                         function(y) return(strsplit(y, "-")[[1]][1])))
pos_end_ATAC = as.numeric(sapply(sapply(featuredata_ATAC, function(x) return(strsplit(x, ":")[[1]][2])),
                                   function(y) return(strsplit(y, "-")[[1]][2])))

ll = list.files("~/data/p2g")
ll = ll[c(1, 2, 3:5, 10, 11, 12, 13, 14, 15, 16:19,
          21:34)]
ll2 = as.character(sapply(ll, function(x) return(strsplit(x, ".tsv")[[1]][1])))
for(numl in 1:length(ll)){
  p2g_Celltype = read.table(paste0("~/data/p2g/", ll[numl]),  header=T)
  idxATAC = p2g_Celltype$idxATAC+1
  idxRNA = p2g_Celltype$idxRNA+1
  outdf = cbind.data.frame(chrvals_ATAC[idxATAC], pos_start_ATAC[idxATAC], pos_end_ATAC[idxATAC], featuredata_RNA[idxRNA])
  colnames(outdf) = c("CHR", "START", "END", "GENE")
  write.table(outdf, file = paste0("~/data/p2g_links/", ll2[numl], ".S2G"),
              row.names = F, col.names = T, sep = "\t", quote=F)
  cat("We are at tissue:", numl, "\n")
}




