library("IRanges") # v 2.36.0
library("dplyr") # v 1.1.4


args <- commandArgs(trailingOnly = TRUE)

outname <- args[1]
spip <- args[2]
pstart <- args[3]
pstop <- args[4]
outcsv <- args[5]
reg <- args[6]
invmeta <- read.table(args[7])
genotypes <- read.csv(args[8])
allsamps <- read.csv(args[9])

# to extract the chromosome names from the pip names
genotypes$reg <- gsub("A", "", genotypes$regn)
genotypes$reg <- gsub("B", "", genotypes$reg)
genotypes$reg <- gsub(".1C", ".1", genotypes$reg)



genotypes$u.record <- paste(genotypes$species, genotypes$regn)
invmeta$u.record <- paste(invmeta$species, invmeta$regn)


genotypes$shared.pip <- NA
for(i in unique(invmeta$u.record[!is.na(invmeta$shared.pip)])){

  genotypes$shared.pip[genotypes$u.record==i] <- invmeta$shared.pip[invmeta$u.record==i]

  }
  
  
  
  
spp <- unique(invmeta$species)
  

picklist <-  matrix(list(), nrow=length(spp), ncol=1)
counter <- 1


for(sp in spp){



if(length(genotypes$sample[genotypes$species==sp&genotypes$shared.pip==spip&!is.na(genotypes$shared.pip)])>0){

genos <- unique(genotypes$genotype[genotypes$shared.pip==spip&genotypes$species==sp&!is.na(genotypes$shared.pip)])
genos <- genos[genos!="AB"]
  
  newpicks <- "dummy"
  
  for(geno in genos){
  pick <- sample(genotypes$sample[genotypes$species==sp&genotypes$genotype==geno
                                  &genotypes$shared.pip==spip&!is.na(genotypes$shared.pip)], 1)
  newpicks <- c(newpicks, pick)
    }
  newpicks <- newpicks[newpicks!="dummy"]

  } else {

  newpicks <- sample(allsamps$sample[allsamps$species==sp], 1)

  }

picksdf <- allsamps[allsamps$sample%in%newpicks,]
picksdf <- left_join(picksdf, genotypes[genotypes$shared.pip==spip,], by="sample")

picklist[[counter]] <- picksdf
                           
counter <- counter+1
  }



df <- do.call(rbind,picklist[,1])




write.csv(df, file=outcsv)













