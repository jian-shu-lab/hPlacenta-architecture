
options(scipen = 999)
celltypes = c("EVT1", "EVT2",
              "EVT3", "Endothelial",
              "Fibroblast1", "Fibroblast2", "Hofbauer",
              "STB", "Unknown3", "Maternal.Macrophage",
              "vCTB1", "vCTB2", "vCTB3",
              "vCTBp")

celltypes2 = c( "EVT1", "EVT2",
                "EVT3", "Endothelial_cells",
                "Fibroblast1", "Fibroblast2", "Hofbauer_cells",
                "STB", "unknown3", "maternal_macrophages",
                "vCTB1", "vCTB2", "vCTB3",
                "vCTBp")


for(numl in 1:length(celltypes)){
  genescore_tabb = read.table(paste0("~/data/Placenta_programs_2023/placenta", "/",
                                     celltypes[numl], ".txt"))
  s2g_tabb = read.table(paste0("~/data/p2g_links/",
                               "p2g_1Mbp_", celltypes2[numl], ".S2G"), header=T)
  scores = genescore_tabb[match(s2g_tabb$GENE, genescore_tabb[,1]), 2]
  outdf = cbind.data.frame(s2g_tabb[, 1], as.numeric(s2g_tabb[,2]), as.numeric(s2g_tabb[,3]), scores)
  outdf2 = outdf[which(outdf[,1] %in% paste0("chr", 1:22)), ]
  write.table(outdf2, file = paste0("~/data/BEDFILES/SHAREseq_S2G/placenta_S2G/",
                                   celltypes[numl], "_1Mbp.bed"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  cat("We are at celltype:", numl, "\n")
}


s2g_all_tabb = c()
for(numl in 1:length(celltypes)){
  s2g_tabb = read.table(paste0("~/data/p2g_links/",
                               "p2g_1Mbp_", celltypes2[numl], ".S2G"), header=T)
  s2g_all_tabb = rbind(s2g_all_tabb, s2g_tabb)
}
outdf = cbind.data.frame(s2g_all_tabb[,1], as.numeric(s2g_all_tabb[,2]), as.numeric(s2g_all_tabb[,3]), 1)
outdf2 = outdf[which(outdf[,1] %in% paste0("chr", 1:22)), ]
write.table(outdf2, file = paste0("~/data/BEDFILES/SHAREseq_S2G/placenta_S2G/",
                                 "ALL", "_1Mbp.bed"),
            row.names = F, col.names = F, sep = "\t", quote=F)






options(scipen = 999)
celltypes = c("EVT1", "EVT2",
              "EVT3", "Endothelial",
              "Fibroblast1", "Fibroblast2", "Hofbauer",
              "STB", "Unknown3", "Maternal.Macrophage",
              "vCTB1", "vCTB2", "vCTB3",
              "vCTBp")

celltypes2 = c( "EVT1", "EVT2",
                "EVT3", "Endothelial_cells",
                "Fibroblast1", "Fibroblast2", "Hofbauer_cells",
                "STB", "unknown3", "maternal_macrophages",
                "vCTB1", "vCTB2", "vCTB3",
                "vCTBp")


for(numl in 1:length(celltypes)){
  s2g_tabb = read.table(paste0("~/data/p2g_links/",
                               "p2g_1Mbp_", celltypes2[numl], ".S2G"), header=T)
  outdf = cbind.data.frame(s2g_tabb[,1], as.numeric(s2g_tabb[,2]), as.numeric(s2g_tabb[,3]), 1)
  outdf2 = outdf[which(outdf[,1] %in% paste0("chr", 1:22)), ]
  write.table(outdf2, file = paste0("~/data/BEDFILES/SHAREseq_S2G/placenta_S2G/",
                                   celltypes[numl], "_1Mbp_onlyS2G.bed"),
              row.names = F, col.names = F, sep = "\t", quote=F)
  cat("we are at celltype:", celltypes[numl], "\n")
}




