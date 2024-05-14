
tt = list.files("~/data/Placenta_programs_2023/placenta")[-1]

human_mouse = read.delim("~/data/HMD_HumanPhenotype.rpt", header=F)

ncbi_tabb = read.table("~/data/NCBI37.3.ensembl.gene.loc", header=T)

celltypes = as.character(sapply(tt, function(x) return(strsplit(x, ".txt")[[1]][1])))

ll = list()
for(numl in 1:length(celltypes)){
  tabb = read.table(paste0("~/data/Placenta_programs_2023/placenta", "/", tt[numl]), header=F)
  genes1 = tabb[which(tabb[,2] > 0.8), 1]
  ncbi_genes1_human= ncbi_tabb$NCBI[match(genes1, ncbi_tabb$HGNC)]
  ncbi_genes1_human = ncbi_genes1_human[!is.na(ncbi_genes1_human)]
  ll[[numl]] = ncbi_genes1_human
  cat("We are at celltype:", numl, "\n")
}
names(ll) = celltypes

ll2 = lapply(ll, function(x) paste0(x, collapse = " "))
outdf = cbind(names(ll2), ll2)

write.table(outdf, file = "~/data/placenta_celltypes.set",
            row.names = F, col.names = F, sep = "\t", quote=F)


#traitnames = list.files("/data/deyk/kushal/Placenta/data/GENE_CELL_0", pattern = ".gene_trait.genes.out")
#traitnames2 = as.character(sapply(traitnames, function(x) return(strsplit(x, ".gene_trait")[[1]][1])))
#write.table(traitnames2, file = "/data/deyk/kushal/Placenta/data/MAGMA/traits_placenta.txt",
#            row.names = F, col.names = F, sep = "\t", quote=F)



