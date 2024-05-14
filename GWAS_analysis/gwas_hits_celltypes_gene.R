snp_tabb = read.table("~/data/GWAS_hits_pregnancy_traits.txt", header=T)
library(GenomicRanges)
gr_gwas  = GRanges(
  seqnames = Rle(paste0(snp_tabb$CHR)),
  ranges = IRanges(snp_tabb$BP_hg38-5, end = snp_tabb$BP_hg38+5),
  SNP = snp_tabb$SNP,
  Trait = snp_tabb$Trait,
  P = snp_tabb$P)

celltypes = c("all", "Endothelial_cells", "EVT1", "EVT2", "EVT3", "Fibroblast1", "Fibroblast2",
              "Hofbauer_cells", "maternal_macrophages", "STB", "vCTB1", "vCTB2", "vCTB3", "vCTBp")

merged_df = c()
for(numl in 1:length(celltypes)){
  tabb = read.table(paste0("~/data/p2g_links//",
                         "p2g_1Mbp_", celltypes[numl], ".S2G"), header=T)
  tabb = tabb[which(!is.na(tabb$START)), ]
  gr_tabb = GRanges(seqnames = Rle(tabb$CHR),
                    ranges = IRanges(tabb$START, end = tabb$END),
                    Gene = tabb$GENE)
  overlaps = data.frame(findOverlaps(gr_gwas, gr_tabb))
  if(nrow(overlaps) > 0){
    merged_df = rbind(merged_df,
                      cbind.data.frame(data.frame(gr_gwas)[overlaps[,1], ],
                                       data.frame(gr_tabb)[overlaps[,2], ],
                                       celltypes[numl]))
  }
}



write.table(merged_df, file = "~/data/GWAS_hits_celltype_gene_placenta_1Mbp.txt",
            row.names = F, col.names = T, sep = "\t", quote=F)
