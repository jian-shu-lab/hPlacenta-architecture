library(data.table)
#library(rmeta)
library(R.utils)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
annot_cell = toString(args[1])
index_in_results=as.numeric(as.character(toString(args[2])))

source("/data/deyk/kushal/Placenta/codes/rmeta/R/meta.R")

annot_modules = list.files(annot_cell)

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
  annot_names = list.dirs(paste0(annot_cell, "/", annot_modules[annot_id]), full.names = F)[-1]
  for(aa in 1:length(annot_names)){
    cell_path = paste0(annot_cell, "/", annot_modules[annot_id], "/", annot_names[aa])
    sd_annot1=get_sd_annot(cell_path, annot_index=base_index, flag = 0)
  }
}

