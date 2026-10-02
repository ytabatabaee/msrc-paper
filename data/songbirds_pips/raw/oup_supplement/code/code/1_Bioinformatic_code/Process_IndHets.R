args <- commandArgs(trailingOnly = TRUE)
library(reshape2)

inds <- read.table(args[1])
regn <- args[3]
reg <- args[4]


colnames(inds) <- "sample"

inds$PHet <- NA
inds$CHet <- NA

for(i in inds$sample){

Ppath=paste0(i,"_", regn,"_PIP.sfs")
PSFS <- scan(Ppath, skip=1)

if(length(PSFS)>0){  
PSFSlength <- sum(PSFS)

inds$PHet[inds$sample==i] <- PSFS[2]/PSFSlength
}else{NULL}



  
Cpath=paste0(i,"_", reg,"_COMP.sfs")
  
CSFS <- scan(Cpath, skip=1)


if(length(CSFS)>0){  
CSFSlength <- sum(CSFS)

inds$CHet[inds$sample==i] <- CSFS[2]/CSFSlength
}else{NULL}


}

inds$regn <- regn

write.csv(inds, file=args[2])