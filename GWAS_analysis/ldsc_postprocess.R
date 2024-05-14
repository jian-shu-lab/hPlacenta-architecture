library(data.table)
#library(rmeta)
library(R.utils)
source("/data/deyk/kushal/Placenta/codes/rmeta/R/meta.R")

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
annot_cell = toString(args[1])
results_cell = toString(args[2])
annot_modules = list.files(results_cell)
trait = toString(args[3])
trait2 = strsplit(trait, ".sumstats")[[1]][1]
flag = as.numeric(as.character(toString(args[4])))
index_in_results = as.numeric(as.character(toString(args[5])))


get_sd_annot = function(cell_path, annot_index = 1, flag=0){
  if(flag == 0){
    if(file.exists(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))){
      sd_annot = get(load(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda")))
      return(sd_annot)
    }else{
      flag = 1
    }}

  if(flag == 1){
    num = 0
    den = 0
    ll <- list.files(cell_path, pattern = ".annot.gz")
    for(m in 1:length(ll)){
      dat <- data.frame(fread(paste0("zcat ", cell_path, "/", ll[m])))
      num = num  + (nrow(dat)-1) * var(dat[,4+annot_index])
      den = den + (nrow(dat)-1)
      rm(dat)
    }
  }

  estd_sd_annot = sqrt(num/den)
  save(estd_sd_annot, file = paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))
  return(estd_sd_annot)
}

base_index = index_in_results

for(annot_id in 1:length(annot_modules)){
  final_df = c()
  annot_names = list.dirs(paste0(results_cell, "/", annot_modules[annot_id]), full.names = F)[-1]
  for(aa in 1:length(annot_names)){
    cell_path = paste0(annot_cell, "/", annot_modules[annot_id], "/", annot_names[aa])
    sd_annot1=get_sd_annot(cell_path, annot_index=base_index, flag = flag)
    result.file=paste0(results_cell, "/", annot_modules[annot_id], "/", annot_names[aa], "/",
                       trait2, ".sumstats.part_delete")
    new_table=read.table(result.file,header=F)
    logfile = paste(results_cell, "/", annot_modules[annot_id], "/", annot_names[aa], "/",
                    trait2,".sumstats.log", sep="")
    log = read.table(logfile,h=F,fill=T)
    h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
    Mref = 5961159
    coef1=sd_annot1*Mref/h2g
    sc=c()
    for(i in 1:dim(new_table)[1]){
      tau1=as.numeric(new_table[i,base_index])
      taus1=tau1*coef1
      sc=c(sc,taus1)
    }
    mean_sc=mean(sc)
    se_sc=sqrt(199**2/200*var(sc))
    p_sc=pnorm(abs(mean_sc/se_sc), 0, 1, lower.tail=F)*2

    result2.file=paste0(results_cell, "/", annot_modules[annot_id], "/", annot_names[aa], "/",
                        trait2, ".sumstats.results")
    res2=read.table(result2.file,header=T)
    myenr    =  res2$Enrichment[base_index] #step1
    myenr_sd = res2$Enrichment_std_error[base_index] #step3
    myenr_p = res2$Enrichment_p[base_index]
    tau_z = res2$Coefficient_z.score[base_index]
    ptau_z = pnorm(res2$Coefficient_z.score[base_index], 0, 1, lower.tail=F)
    final_df = rbind(final_df, c(mean_sc, se_sc, p_sc, myenr, myenr_sd, myenr_p, tau_z, ptau_z))
  }

  rownames(final_df) = annot_names
  colnames(final_df) = c("tau-star", "se-tau-star", "p-tau-star", "E", "se(E)", "p(E)", "tau_z", "ptau_z")
  write.table(final_df, file = paste0(results_cell, "/", annot_modules[annot_id], "/",
                                      trait2, "_ldsc_postprocess.txt"), quote=F, sep = "\t")
}

