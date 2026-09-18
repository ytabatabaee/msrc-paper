root <- normalizePath(".")
results <- file.path(root, "empirical/neoaves_chr4/results")
figures <- file.path(root, "empirical/neoaves_chr4/figures")
d <- read.delim(file.path(results, "stage4b_inference_treatments.tsv"), check.names=FALSE)
loo <- read.delim(file.path(results, "stage4b_leave_one_event_out.tsv"), check.names=FALSE)
p <- d[d$structural_mask == "primary", ]
clades <- c("Columbea", "N61", "N62")
treatments <- c("T0", "T1", "T2", "T3")
cols <- c(q1="#2563a6", q2="#d16a32", q3="#6b7280")

draw_a <- function() {
  old <- par(mfrow=c(1,3), mar=c(4,4,2,1), oma=c(0,0,0,4), las=1)
  for (clade in clades) {
    z <- p[p$clade == clade, ]; z <- z[match(treatments, z$treatment), ]
    barplot(t(as.matrix(z[,c("q1","q2","q3")])), col=cols, names.arg=treatments,
            ylim=c(0,1), main=clade, xlab="Treatment", ylab=if(clade==clades[1]) "Weighted topology support" else "")
  }
  legend("right", inset=c(-0.18,0), legend=names(cols), fill=cols, bty="n", xpd=NA)
  par(old)
}
pdf(file.path(figures, "stage4b_figure_A_treatment_support.pdf"), width=10.5, height=3.5); draw_a(); dev.off()
png(file.path(figures, "stage4b_figure_A_treatment_support.png"), width=3150, height=1050, res=300); draw_a(); dev.off()

draw_b <- function() {
  plot(NA, xlim=c(.7,4.3), ylim=range(c(p$M_species_ci_low,p$M_species_ci_high)), xaxt="n",
       xlab="Treatment", ylab="M_species = q1 - max(q2, q3)", bty="l")
  abline(h=0, lwd=.8); axis(1, 1:4, treatments)
  offsets <- c(-.16,0,.16); shapes <- c(16,15,17); ccols <- c("#2563a6","#d16a32","#4b5563")
  for (i in seq_along(clades)) {
    z <- p[p$clade == clades[i], ]; z <- z[match(treatments,z$treatment), ]; x <- 1:4+offsets[i]
    arrows(x,z$M_species_ci_low,x,z$M_species_ci_high,angle=90,code=3,length=.035,col=ccols[i])
    lines(x,z$M_species,type="b",pch=shapes[i],col=ccols[i])
  }
  legend("topright",clades,pch=shapes,col=ccols,lty=1,bty="n")
}
pdf(file.path(figures, "stage4b_figure_B_species_margin.pdf"), width=7.2, height=4.2); draw_b(); dev.off()
png(file.path(figures, "stage4b_figure_B_species_margin.png"), width=2160, height=1260, res=300); draw_b(); dev.off()

draw_c <- function() {
  z <- loo[loo$clade == "N61", ]; ids <- unique(z$event_id_removed); short <- sub(".*_", "", ids)
  yr <- range(z$M_q1_q2)
  plot(NA,xlim=c(1,length(ids)),ylim=yr,xaxt="n",xlab="Frozen event removed",ylab="N61 q1 - q2",bty="l")
  abline(h=0,lwd=.8); axis(1,1:length(ids),short,las=2,cex.axis=.75)
  for (i in 1:2) { tr <- c("T2","T3")[i]; zz <- z[z$treatment==tr, ]; zz <- zz[match(ids,zz$event_id_removed), ]; lines(1:length(ids),zz$M_q1_q2,type="b",pch=c(16,15)[i],col=c("#2563a6","#d16a32")[i]) }
  legend("topright",c("T2","T3"),pch=c(16,15),col=c("#2563a6","#d16a32"),lty=1,bty="n")
}
pdf(file.path(figures, "stage4b_figure_C_N61_leave_one_event_out.pdf"), width=8.5, height=4.2); draw_c(); dev.off()
png(file.path(figures, "stage4b_figure_C_N61_leave_one_event_out.png"), width=2550, height=1260, res=300); draw_c(); dev.off()
