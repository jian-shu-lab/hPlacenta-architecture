
all_gwas_hits = read.table("~/data/GWAS_hits_pregnancy_traits.txt", header=T)
merged_bimtabb_list = list()
for(numchr in 1:22){
  bimtabb = read.table(paste0("~/data/1000G_BIMS_hg38/", "1000G.EUR.QC.", numchr, ".bim"))
  merged_bimtabb_list[[numchr]] = bimtabb
}
merged_bimtabb = do.call(rbind, merged_bimtabb_list)

length(intersect(all_gwas_hits$SNP, merged_bimtabb[,2]))

ll = list.files("~/data/ANNOTATIONS_hg38/SHAREseq_S2G/placenta_S2G/")
enrvec = c()
sd_enrvec = c()
for(numl in 1:length(ll)){
  ss1 = c()
  ss2 = c()
  ss3 = c()
  ss4 = c()
  for(numchr in 1:22){
    tabb = data.frame(fread(paste0("~/data/ANNOTATIONS_hg38/SHAREseq_S2G/placenta_S2G/", ll[numl], "/",
                                   ll[numl], ".", numchr, ".annot.gz")))
    common_snps = intersect(all_gwas_hits$SNP, tabb$SNP)
    ss1 = c(ss1, sum(tabb[match(common_snps, tabb$SNP), 5]))
    ss2 = c(ss2, sum(tabb[, 5]))
    ss3 = c(ss3, nrow(tabb))
    ss4 = c(ss4, length(common_snps))
  }
  enrvec = c(enrvec, (sum(ss1)/sum(ss4))/(sum(ss2)/sum(ss3)))
  boot_enrvec = c()
  for(nboot in 1:100){
    idx = sample(1:length(ss1), length(ss1), replace=T)
    bss1 = ss1[idx]
    bss2 = ss2[idx]
    bss3 = ss3[idx]
    bss4 = ss4[idx]
    boot_enrvec = c(boot_enrvec, (sum(bss1)/sum(bss4))/(sum(bss2)/sum(bss3)))
  }
  sd_enrvec = c(sd_enrvec, sd(boot_enrvec))
  cat("We are at cell type:", numl, "\n")
}

merged_df = cbind.data.frame(enrvec, sd_enrvec)
colnames(merged_df) = c("ENR", "s.ENR")
rownames(merged_df) = ll

write.table(merged_df, file = "~/data/ENR_Placenta_celltypeS2G_pregnancy_traits.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)

celltypes2 = list.files("~/data/ANNOTATIONS_hg38/placenta")
enrvec = c()
sd_enrvec = c()
for(numl in 1:length(celltypes2)){
  for(numchr in 1:22){
    tabb = data.frame(fread(paste0("~/data/ANNOTATIONS_hg38/placenta/",
                                   celltypes2[numl], "/",
                                   "ABC_Road_GI_PLCNT", "/",
                                   "ABC_Road_GI_PLCNT", ".", numchr, ".annot.gz")))
    common_snps = intersect(all_gwas_hits$SNP, tabb$SNP)
    ss1 = c(ss1, sum(tabb[match(common_snps, tabb$SNP), 5]))
    ss2 = c(ss2, sum(tabb[, 5]))
    ss3 = c(ss3, nrow(tabb))
    ss4 = c(ss4, length(common_snps))
  }
  enrvec = c(enrvec, (sum(ss1)/sum(ss4))/(sum(ss2)/sum(ss3)))
  boot_enrvec = c()
  for(nboot in 1:100){
    idx = sample(1:length(ss1), length(ss1), replace=T)
    bss1 = ss1[idx]
    bss2 = ss2[idx]
    bss3 = ss3[idx]
    bss4 = ss4[idx]
    boot_enrvec = c(boot_enrvec, (sum(bss1)/sum(bss4))/(sum(bss2)/sum(bss3)))
  }
  sd_enrvec = c(sd_enrvec, sd(boot_enrvec))
  cat("We are at cell type:", numl, "\n")
}


merged_df = cbind.data.frame(enrvec, sd_enrvec)
colnames(merged_df) = c("ENR", "s.ENR")
rownames(merged_df) = ll

write.table(merged_df, file = "~/data/ENR_Placenta_sclinker_pregnancy_traits.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)
