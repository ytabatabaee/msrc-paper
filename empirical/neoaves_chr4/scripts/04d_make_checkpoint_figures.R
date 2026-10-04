#!/usr/bin/env Rscript
d <- read.delim("empirical/neoaves_chr4/results/stage4d_figure4_reproduction.tsv", check.names=FALSE)
draw <- function() {
  old <- par(mfrow=c(1,2), mar=c(3,4,2,1), oma=c(2,0,0,0), las=1)
  on.exit(par(old))
  for (tr in c("T0","T_PNAS")) {
    x <- d[d$treatment==tr,]
    plot(seq_len(nrow(x)), x$mean, ylim=range(c(x$mean-x$se,x$mean+x$se,x$median,0)),
         pch=16, cex=.55, col="#2040DE", xaxt="n", xlab="", ylab=if(tr=="T0") "Delta quartet score (S2024 - J2014)" else "",
         main=if(tr=="T0") "All loci" else "Published outlier loci removed")
    arrows(seq_len(nrow(x)),x$mean-x$se,seq_len(nrow(x)),x$mean+x$se,
           angle=90,code=3,length=.025,col="#2040DE")
    points(seq_len(nrow(x)),x$median,pch=16,cex=.55,col="#20A060")
    abline(h=0,col="#B02020",lwd=.8)
    if(tr=="T0") legend("bottomleft",c("mean +/- SE","median"),pch=16,col=c("#2040DE","#20A060"),bty="n",cex=.75)
  }
  mtext("Original Figure-4 taxon-removal order",side=1,outer=TRUE,line=.5)
}
pdf("empirical/neoaves_chr4/figures/stage4d_figure_A_original_figure4_reproduction.pdf",width=9,height=3.8,useDingbats=FALSE); draw(); dev.off()
png("empirical/neoaves_chr4/figures/stage4d_figure_A_original_figure4_reproduction.png",width=2700,height=1140,res=300); draw(); dev.off()
