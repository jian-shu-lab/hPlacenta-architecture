
library(data.table)
ll = list.files("~/data/ANNOTATIONS_hg38/SHAREseq_S2G/placenta_S2G/")
#ll = ll[-grep("onlyS2G", ll)]
ll = ll[grep("onlyS2G", ll)]

all_gwas_hits = read.table("~/data/GWAS_hits_pregnancy_traits.txt", header=T)
traitnames=names(table(all_gwas_hits$Trait))[which(table(all_gwas_hits$Trait) > 50)]

enrmat = c()
sd_enrmat = c()

for(traitnum in 1:length(traitnames)){
  all_gwas_hits2 = all_gwas_hits[which(all_gwas_hits$Trait == traitnames[traitnum]), ]
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
      common_snps = intersect(all_gwas_hits2$SNP, tabb$SNP)
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
  enrmat = cbind(enrmat, enrvec)
  sd_enrmat = cbind(sd_enrmat, sd_enrvec)
  cat("We are at trait:", traitnum, "\n")
}


rownames(enrmat) = ll
colnames(enrmat) = traitnames

rownames(sd_enrmat) = ll
colnames(sd_enrmat) = traitnames

ll = list("ENR" = enrmat, "sENR" = sd_enrmat)
save(ll, file = "~/data/ENR_eachtrait_celltype_S2G_onlyS2G.rda")
