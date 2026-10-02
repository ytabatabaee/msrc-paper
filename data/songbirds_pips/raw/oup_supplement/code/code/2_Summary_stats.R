#####################################################################################################################
# load libraries
#####################################################################################################################
library(tidyverse) # v 2.0.0
library(geosphere) # v 1.5-18
library(reshape2) # v 1.4.4
library(HardyWeinberg) # v 1.7.8
library(emmeans) # v 1.10.0

#####################################################################################################################
# load relevant datasets
#####################################################################################################################

# Be sure to set your working directory to where the relevant csv files are kept. We provide these
# SI datasets within the "code" folder containing this script.

# metadata about each PIP. This is Dataset S1.

invmeta <- read.csv("DatasetS1.csv")

# metadata about each individual, including their genotypes. This is Dataset S4.

indmeta <- read.csv("DatasetS4.csv")

# metadata about pips that might be shared across species boundaries. This is Dataset S5.
# note that these are named with a different naming system because a given potentially overlapping PIP 
# could be named with an "A" in some species and not others

sharemeta <- read.csv("DatasetS5.csv")


# for some of our analyses, we want invmeta in the individual level data too, so we left_join them

indmeta <- left_join(indmeta, invmeta, by=c("Species", "PIP.name"))

# there are many analyses we will only do on biallelic pips, so we'll make a separate dataset for them

invmeta_bi <- invmeta[invmeta$PCA.cluster.number==3,]
  
indmeta_bi <- indmeta %>%
  filter(indmeta$PIP.genotype%in%c("AA","AB","BB"))


#####################################################################################################################
# Results 1: overall statistics
#####################################################################################################################

# how many PIPs are there

length(invmeta[,1])

# how many PIPs per species

sptab_all <- data.frame(table(invmeta$Species))
sptab_LD <- data.frame(table(invmeta$Species[invmeta$LD.supported==TRUE]))

colnames(sptab_all) <- c("Species", "Freq")
colnames(sptab_LD) <- c("Species", "Freq")

sptab_all
sptab_LD

# how many species have PIPs
length(sptab_all[,1])

# summary stats across species
mean(sptab_all$Freq)
sd(sptab_all$Freq)
range(sptab_all$Freq)
hist(sptab_all$Freq)

mean(sptab_LD$Freq)
sd(sptab_LD$Freq)
range(sptab_LD$Freq)
hist(sptab_LD$Freq)


# how many PIPs produce consistent patterns on LD and PCA

table(invmeta$LD.supported)


#####################################################################################################################
# Results 2: PIP lengths and proximity to HORs
#####################################################################################################################

# how long are the PIPs 
mean(invmeta$PIP.length[invmeta$LD.supported==T])
sd(invmeta$PIP.length[invmeta$LD.supported==T])
range(invmeta$PIP.length[invmeta$LD.supported==T])

mean(invmeta$PIP.length[invmeta$LD.supported==F])
sd(invmeta$PIP.length[invmeta$LD.supported==F])
range(invmeta$PIP.length[invmeta$LD.supported==F])


range(invmeta$PIP.length)

# distances between PIP breakpoints and HOR annotations

table(inv$min_HOR_dist[!is.na(invmeta$min_HOR_dist)]<50000)

# See 3_CentroAnno.txt for code for information about annotating HORs

#####################################################################################################################
# Results 3: characteristics of PIP genotypes -- frequencies and HWE
#####################################################################################################################

# how many pips are "biallelic"?

length(invmeta_bi[,1]) # 133

# as proportion of total?
length(invmeta_bi[,1]) / length(invmeta[,1])

# overall
table(invmeta$PCA.cluster.number)

# number of biallelic pips with LD support

table(invmeta_bi$LD.supported)

# ALLELE FREQUENCIES

# biallelic only!

# all the little math bits make it so the number reflects the number of B's in the genotype 

indmeta_bi$num_geno <- ((-as.numeric(factor(indmeta_bi$PIP.genotype))+1)*-1)
indmeta_bi$u.record <- paste(indmeta_bi$Species, indmeta_bi$PIP.name)
invmeta_bi$u.record <- paste(invmeta_bi$Species, invmeta_bi$PIP.name)
invmeta_bi$MAF = NA

for(i in unique(indmeta_bi$u.record)){
  
  temp <- indmeta_bi[indmeta_bi$u.record==i,] 
  n <- length(temp$num_geno)*2
  q <- sum(temp$num_geno)/n
  p <- 1-q
  
  invmeta_bi$p[invmeta_bi$u.record==i] <- p
  invmeta_bi$q[invmeta_bi$u.record==i] <- q
  
  if(p>q|p==q) {
    invmeta_bi$MAF[invmeta_bi$u.record==i] = q
  } else{NULL}
  
  if(q>p) {
    invmeta_bi$MAF[invmeta_bi$u.record==i] = p
  } else{NULL}
  
  
  
}

# mean allele frequencies
mean(invmeta_bi$MAF[invmeta_bi$LD.supported==T])
sd(invmeta_bi$MAF[invmeta_bi$LD.supported==T])

mean(invmeta_bi$MAF[invmeta_bi$LD.supported==F])
sd(invmeta_bi$MAF[invmeta_bi$LD.supported==F])



# HWE

invmeta_bi$HWE_pval <- NA

for(i in unique(indmeta_bi$u.record)){
  
  
  temp <- indmeta_bi[indmeta_bi$u.record==i,] 
  
  
  invmeta_bi$HWE_pval[invmeta_bi$u.record==i] <- HWExact(table(temp$PIP.genotype), verbose = TRUE)$pval
  
  
}
table(invmeta_bi$HWE_pval<0.05)
123/133 # proportion that do not deviate from HWE

sigHWE <- invmeta_bi$u.record[invmeta_bi$HWE_pval<0.05]

# to print the D values for the significant ones
# this is silly,  but we have to just look at the output as printed on the console
# because of how this function produces its output
for(i in sigHWE){
  
  temp <- indmeta_bi[indmeta_bi$u.record==i,] 
  
  print(i)
  HWExact(table(temp$PIP.genotype), verbose = TRUE)$pval
  
  
}

#####################################################################################################################
# Results 4: Heterozygosity by genotype class
#####################################################################################################################

# See scripts in folder 1_Bioinformatic_code for information about how individuals are genotyped and how 
# heterozygosity is estimated.


invmeta$u.record <- paste(invmeta$Species, invmeta$PIP.name)
indmeta$u.record <- paste(indmeta$Species, indmeta$PIP.name)

# create a useful column for the loop below
invmeta$ignore.comp <- NA
for(i in invmeta$u.record){
  
  invmeta$ignore.comp[invmeta$u.record==i] <- all(is.na(indmeta$COMP.heterozygosity[indmeta$u.record==i])) 
  
}
# and some other ones where we will skip the "comp" heterozygosity calculations because we don't have good estimates from them
invmeta$ignore.comp[invmeta$u.record%in%c("Catharus ustulatus CM020371.1A")] <- T
invmeta$ignore.comp[invmeta$u.record%in%c("Catharus ustulatus CM020371.1B")] <- T
invmeta$ignore.comp[invmeta$u.record%in%c("Catharus ustulatus CM020371.1C")] <- T


table(invmeta$ignore.comp)


hetmatlist <- matrix(list(), nrow=length(invmeta$u.record), ncol=3)
hetmatlist[,1] <- invmeta$u.record

invmeta$PIPhet_elevated <- NA
invmeta$PIPhom_depressed <- NA
invmeta$PIPCD_diff <- NA
invmeta$COMPhet_elevated <- NA
invmeta$COMPhom_depressed <- NA
invmeta$COMPCD_diff <- NA


for(i in invmeta$u.record){
  
  
  temp <- indmeta[indmeta$u.record==i,]
  
  if(length(unique(temp$PIP.genotype))<4){
    
    
    
    pipem <- emmeans(lm(PIP.heterozygosity ~ PIP.genotype, data=temp), pairwise ~ PIP.genotype)
    
    
    # for the 3-cluster ones
    if(length(unique(temp$PIP.genotype))==3){
      
      
      pipmdf <- as.data.frame(pipem$emmeans)
      pipcdf <- as.data.frame(pipem$contrasts)
      
      AAmean <- pipmdf$emmean[pipmdf$PIP.genotype=="AA"]
      ABmean <- pipmdf$emmean[pipmdf$PIP.genotype=="AB"]
      BBmean <- pipmdf$emmean[pipmdf$PIP.genotype=="BB"]
      
      pAA.AB <- pipcdf$p.value[pipcdf$contrast=="AA - AB"]
      pAA.BB <- pipcdf$p.value[pipcdf$contrast=="AA - BB"]
      pAB.BB <- pipcdf$p.value[pipcdf$contrast=="AB - BB"]
      
      invmeta$PIPhet_elevated[invmeta$u.record==i] <- ABmean>AAmean&ABmean>BBmean&pAA.AB<0.05&pAB.BB<0.05
      
      
      # the "het elevated" condition is true if p value for AA-AB and AB-BB are both less than 0.05 AND
      # AB's emmean is greater than both AA and BB
      
      
      # the "hom depressed" condition is true if one of two conditions is met:
      # AA-AB is significant and AA is low AND AB-BB is not significant
      # or
      # AB-BB is significant and BB is low AND AA-AB is not significant
      
      invmeta$PIPhom_depressed[invmeta$u.record==i] <- pAA.AB<0.05&AAmean<ABmean&pAB.BB>0.05|pAB.BB<0.05&BBmean<ABmean&pAA.AB>0.05
      
      
      if(invmeta$ignore.comp[invmeta$u.record==i]==F){
        compem <- emmeans(lm(COMP.heterozygosity ~ PIP.genotype, data=temp), pairwise ~ PIP.genotype)
        
        
        compmdf <- as.data.frame(compem$emmeans)
        compcdf <- as.data.frame(compem$contrasts)
        
        cAAmean <- compmdf$emmean[compmdf$PIP.genotype=="AA"]
        cABmean <- compmdf$emmean[compmdf$PIP.genotype=="AB"]
        cBBmean <- compmdf$emmean[compmdf$PIP.genotype=="BB"]
        
        cpAA.AB <- compcdf$p.value[compcdf$contrast=="AA - AB"]
        cpAA.BB <- compcdf$p.value[compcdf$contrast=="AA - BB"]
        cpAB.BB <- compcdf$p.value[compcdf$contrast=="AB - BB"]
        
        invmeta$COMPhet_elevated[invmeta$u.record==i] <- cABmean>cAAmean&cABmean>cBBmean&cpAA.AB<0.05&cpAB.BB<0.05
        invmeta$COMPhom_depressed[invmeta$u.record==i] <- cpAA.AB<0.05&cAAmean<cABmean&cpAB.BB>0.05|cpAB.BB<0.05&cBBmean<cABmean&cpAA.AB>0.05
        
      }else{NULL}
      
      
    }else{NULL}
    
    
    # for the 2-cluster ones
    if(length(unique(temp$PIP.genotype))==2){
      
      
      pipcdf <- as.data.frame(pipem$contrasts)
      
      
      invmeta$PIPCD_diff[invmeta$u.record==i] <- pipcdf$p.value<0.05
      
      if(invmeta$ignore.comp[invmeta$u.record==i]==F){
        compem <- emmeans(lm(COMP.heterozygosity ~ PIP.genotype, data=temp), pairwise ~ PIP.genotype)
        compcdf <- as.data.frame(compem$contrasts)
        invmeta$COMPCD_diff[invmeta$u.record==i] <- compcdf$p.value<0.05
      }else{NULL}
      
      
    } else{NULL}
    
    
    # save the objects
    
    hetmatlist[[which(hetmatlist[,1]==i),2]] <- pipem
    hetmatlist[[which(hetmatlist[,1]==i),3]] <- compem
    
    
    
  }else{NULL}
  
  
  
}



table(invmeta$PIPhet_elevated[invmeta$PCA.cluster.number==3])
table(invmeta$PIPhom_depressed[invmeta$PCA.cluster.number==3])

table(invmeta$PIPCD_diff[invmeta$PCA.cluster.number==2])

#####################################################################################################################
# Results 5: Estimation of IBD slopes
#####################################################################################################################

# biallelic only!

# this piece of code estimates IBD. Note that the IBD slope is already in the invmeta data provided.

invIBDlist <- matrix(list(), nrow=length(unique(indmeta_bi$u.record)), ncol=2)
invIBDlist[,1] <- unique(indmeta_bi$u.record)


for(i in unique(indmeta_bi$u.record)){
  
  temp <- indmeta_bi[indmeta_bi$u.record==i,] 
  
  distmat <- as.matrix(dist(temp$num_geno, upper=T, diag=T))
  row.names(distmat) <- temp$Sample
  colnames(distmat) <- temp$Sample
  
  
  
  
  genmatdf <- setNames(reshape2::melt(distmat), c("Ind1", "Ind2", "invhap_dist"))
  
  genmatdf$Ind1 <- as.character(genmatdf$Ind1)
  genmatdf$Ind2 <- as.character(genmatdf$Ind2)
  
  
  genmatdf$Lat1 <- NA
  genmatdf$Lat2 <- NA
  genmatdf$Lon1 <- NA
  genmatdf$Lon2 <- NA
  
  
  
  for(k in unique(c(genmatdf$Ind1, genmatdf$Ind2))){
    
    
    
    genmatdf$Lat1[genmatdf$Ind1==k] <- unique(indmeta_bi$Latitude[indmeta_bi$Sample==k])
    genmatdf$Lat2[genmatdf$Ind2==k] <- unique(indmeta_bi$Latitude[indmeta_bi$Sample==k])
    genmatdf$Lon1[genmatdf$Ind1==k] <- unique(indmeta_bi$Longitude[indmeta_bi$Sample==k])
    genmatdf$Lon2[genmatdf$Ind2==k] <- unique(indmeta_bi$Longitude[indmeta_bi$Sample==k])
    
  }
  
  
  
  genmatdf$geodist <- distGeo(  cbind(as.numeric(genmatdf$Lon1), as.numeric(genmatdf$Lat1)),  
                                cbind(as.numeric(genmatdf$Lon2), as.numeric(genmatdf$Lat2)))
  
  
  
  
  genmatdf <- genmatdf[genmatdf$Ind1!=genmatdf$Ind2,] 
  
  
  genmatdf$distKm <- genmatdf$geodist/1000
  genmatdf$dist1kKm <- genmatdf$distKm/1000
  
  genmatdf$u.record <- i
  genmatdf$IBDslope <- coef(lm(genmatdf$invhap_dist ~ genmatdf$distKm))[2]
  
  
  invIBDlist[[which(invIBDlist[,1]==i),2]] <- genmatdf
  
  
  
}


# how this is added to the main dataset
# allIBD <- do.call(rbind, invIBDlist[,2])
# 
# invmeta$IBD.slope <- NA
# for(i in unique(allIBD$u.record)){
#   invmeta$IBD.slope[invmeta$u.record==i] <- unique(allIBD$IBDslope[allIBD$u.record==i])
# 
# }

range(invmeta$IBD.slope, na.rm=T)





#####################################################################################################################
# Results 6: Trans-specific polymorphisms
#####################################################################################################################

# how many potentially shared pips did we find?
unique(sharemeta$Shared.PIP.tree)

# See scripts in folder 1_Bioinformatic_code for information about how gene trees of possibly shared PIPs are created


