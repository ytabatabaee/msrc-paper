root <- normalizePath(".")
dat <- file.path(root, "data/neoaves_chr4/processed")
res <- file.path(root, "empirical/neoaves_chr4/results")
fig <- file.path(root, "empirical/neoaves_chr4/figures")
dir.create(fig, recursive=TRUE, showWarnings=FALSE)
cols <- c(q1="#2673B8", q2="#D65F5F", q3="#4A9B68")
clades <- c("Columbea", "N61", "N62")

drawA <- function() {
  loci <- read.delim(file.path(dat, "locus_table_chr4.tsv"))
  bp <- read.delim(file.path(dat, "denovo_breakpoint_sets.tsv"))
  bp <- bp[tolower(bp$in_primary) %in% c("true", "1", "yes"), ]
  par(mfrow=c(3,1), mar=c(2.4,4.3,1.6,1), oma=c(2.6,0,3.2,0), xaxs="i")
  for (cl in clades) {
    x <- loci[loci$clade == cl, ]
    plot(x$midpoint/1e6, x$dominant_support, type="n", ylim=c(0,1), xlim=c(0,91.31), xaxt=if(cl=="N62") "s" else "n", xlab="", ylab=paste(cl, "support"))
    usr <- par("usr")
    for (p in bp$reference_position/1e6) rect(max(0,p-.25), usr[3], min(91.31,p+.25), usr[4], col=adjustcolor("grey50", alpha.f=.14), border=NA)
    for (q in names(cols)) {
      z <- x$dominant_topology == q
      points(x$midpoint[z]/1e6, x[[q]][z], pch=16, cex=.25, col=adjustcolor(cols[q], alpha.f=.5))
    }
    if (cl == "Columbea") legend("topright", names(cols), col=cols, pch=16, bty="n", horiz=TRUE)
  }
  mtext("Chromosome 4 position (Mb)", side=1, outer=TRUE, line=1)
  mtext("Topology state/support and frozen primary structural neighborhoods", outer=TRUE, side=3, line=1.5, cex=1.15)
  mtext("Shading denotes association; causal identification is not implied", outer=TRUE, side=3, line=.2, cex=.85)
}

drawBars <- function(file, filter_col, levels, title, subtitle, ylabel) {
  x <- read.delim(file)
  par(mfrow=c(1,3), mar=c(4,4,2.2,.6), oma=c(0,0,3.2,0))
  for (cl in clades) {
    z <- x[x$clade == cl, ]
    z <- z[match(levels, z[[filter_col]]), ]
    m <- rbind(c(z$q_1[1],z$q_2[1],z$q_3[1]), c(z$q_1[2],z$q_2[2],z$q_3[2]))
    if (!"q_1" %in% names(z)) m <- rbind(c(z$frequency_q1[1],z$frequency_q2[1],z$frequency_q3[1]), c(z$frequency_q1[2],z$frequency_q2[2],z$frequency_q3[2]))
    barplot(m, beside=TRUE, names.arg=c("q1","q2","q3"), ylim=c(0,max(.6,m)*1.12), col=c("#5B8DB8","#D69A4C"), border=NA, main=cl, ylab=if(cl=="Columbea") ylabel else "")
    if (cl == "N62") legend("topright", gsub("structural-", "", levels), fill=c("#5B8DB8","#D69A4C"), bty="n", cex=.8)
  }
  mtext(title, outer=TRUE, side=3, line=1.5, cex=1.15)
  mtext(subtitle, outer=TRUE, side=3, line=.2, cex=.85)
}

devices <- list(
  list(file="stage4a_figure_A_chr4_structural_topology", fun=drawA),
  list(file="stage4a_figure_B_structural_category_frequencies", fun=function() {
    s <- read.delim(file.path(res,"stage4a_structural_topology_support.tsv")); s <- s[s$structural_set=="primary", ]; tmp <- tempfile(fileext=".tsv"); write.table(s,tmp,sep="\t",row.names=FALSE,quote=FALSE)
    drawBars(tmp,"structural_category",c("structural-associated","structural-background"),"Topology frequencies by primary structural category","Association only; causal identification is not implied","Frequency")
  }),
  list(file="stage4a_figure_C_window_vs_block", fun=function() drawBars(file.path(res,"stage4a_window_vs_block_support.tsv"),"weighting_scheme",c("window-weighted","block-normalized"),"Window-weighted and block-normalized topology support","Descriptive analogue of competing-topology fractions","Evidence fraction"))
)
for (d in devices) {
  pdf(file.path(fig,paste0(d$file,".pdf")), width=11, height=if(grepl("figure_A",d$file)) 7 else 4.2, useDingbats=FALSE); d$fun(); dev.off()
  png(file.path(fig,paste0(d$file,".png")), width=3300, height=if(grepl("figure_A",d$file)) 2100 else 1260, res=300); d$fun(); dev.off()
}
