root <- normalizePath(".")
results <- file.path(root, "empirical/neoaves_chr4/results")
figures <- file.path(root, "empirical/neoaves_chr4/figures")
d <- read.delim(file.path(results,"stage4c_independent_validation.tsv"),check.names=FALSE)
ctrl <- read.delim(file.path(results,"stage4c_control_chromosomes.tsv"),check.names=FALSE)
clades <- c("Columbea","N61","N62"); treatments <- c("T0","T1","T2","T3")
cols <- c(q1="#2563a6",q2="#d16a32",q3="#6b7280")

draw_a <- function() {
  old <- par(mfrow=c(1,3),mar=c(5,4,2,1),oma=c(0,0,0,4),las=1)
  for (clade in clades) {
    z <- d[d$clade==clade,]; z <- z[match(treatments,z$treatment),]
    values <- rbind(c(z$reference_q1[1],z$treatment_q1),c(z$reference_q2[1],z$treatment_q2),c(z$reference_q3[1],z$treatment_q3))
    barplot(values,col=cols,names.arg=c("non-chr4",treatments),las=2,ylim=c(0,1),main=clade,ylab=if(clade==clades[1]) "Topology support" else "")
  }
  legend("right",inset=c(-.18,0),legend=names(cols),fill=cols,bty="n",xpd=NA); par(old)
}
pdf(file.path(figures,"stage4c_figure_A_reference_comparison.pdf"),10.5,3.8); draw_a(); dev.off()
png(file.path(figures,"stage4c_figure_A_reference_comparison.png"),3150,1140,res=300); draw_a(); dev.off()

draw_b <- function() {
  z <- d[d$clade=="N61",]; z <- z[match(treatments,z$treatment),]
  vals <- c(z$M_reference,z$reference_margin[1]); labs <- c(treatments,"non-chr4\nreference")
  plot(1:5,vals,type="b",pch=16,col="#2563a6",xaxt="n",xlab="Descriptive comparison",ylab="N61 M_reference",bty="l",ylim=range(c(vals,0)))
  axis(1,1:5,labs); abline(h=0,lwd=.8); arrows(1:4,vals[1:4],2:5,vals[2:5],length=.08,col="#6b7280")
  points(1:5,vals,pch=c(16,17,15,18,16),cex=1.2,col=c("#d16a32","#2563a6","#4b849f","#7b6ca8","#111827"))
  text(3,.95*max(vals),"Arrows denote comparisons, not causation",cex=.8)
}
pdf(file.path(figures,"stage4c_figure_B_N61_focal.pdf"),7.5,4.5); draw_b(); dev.off()
png(file.path(figures,"stage4c_figure_B_N61_focal.png"),2250,1350,res=300); draw_b(); dev.off()

draw_c <- function() {
  main <- ctrl[ctrl$focal_clade=="N61",]; main <- main[order(as.numeric(sub("chr","",main$chromosome))),]
  x <- seq_len(nrow(main)); col <- ifelse(main$chromosome=="chr4","#d16a32","#6b7280")
  plot(x,main$delta_M_reference_block_vs_window,pch=16,col=col,xaxt="n",xlab="Main chromosome",ylab="N61 Delta M_reference (block - window)",bty="l")
  axis(1,x,sub("chr","",main$chromosome),cex.axis=.65); abline(h=0,lwd=.8); legend("bottomright",c("chr4","other chromosomes"),pch=16,col=c("#d16a32","#6b7280"),bty="n")
}
pdf(file.path(figures,"stage4c_figure_C_control_context.pdf"),9,4.5); draw_c(); dev.off()
png(file.path(figures,"stage4c_figure_C_control_context.png"),2700,1350,res=300); draw_c(); dev.off()
