library(data.table)

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
inds <- read.table(args[2])
outcsv <- args[3]
outpdf <- args[4]



C <- as.matrix(read.table(input))
e <- eigen(C)

e$values/sum(e$values)

annotated <- cbind(inds, e$vectors)

 
write.csv(annotated, file=outcsv)

pdf(file=outpdf)
plot(e$vectors[,1], e$vectors[,2],
     xlab="PC1", ylab="PC2")
dev.off()