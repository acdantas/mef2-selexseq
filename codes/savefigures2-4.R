###############################################################################################
#                                   SAVE Figures
###############################################################################################


# MOTIF ANALYSIS - STRIPCHART / PLOT / HEATMAP ------


# Figure 2 ------------------------------


png("Fig2-AB-Stripchart-Plot.png", width = 180, height = 65, units = 'mm', res=300)
layout(matrix(c(1,1,2,2,3), 1, 5, byrow = TRUE),
       widths=c(2,2,1)) #, heights=c(1,2))
par(mgp=c(2,0.75,0), mar=c(8,4,2,1), xpd=TRUE)
plotstripchart(dat=dfByGrp)
plotMotifAff(dat=merged, myxlab="CTAW4TAG")
dev.off()

png("Fig2-PanelC.png", width = 180, height = 65, units = 'mm', res=300 )
get.Kmerby3mer(dat=l.trimer ,
               grp="(3mer)W4TAR",
               col2="kmer", col1="leftKmer")
dev.off()


# KMER ANALYSIS - FIGURE 3 and S4


# Figure 3 ------------------------


# Figure3A - mean relative affinity

png(paste("meanRelAff-heatmap-CTANNNNTAG.png"), width=100, height=100, units = 'mm', res=600)
get.heatmap.myaffMeans(dat=my.affMeans, mylabels=c4.lab)
dev.off()


# Make boxplot --------

# Figure3 B-C

png("Fig3-AffinityBxp.png", width = 120, height = 80, units = 'mm', res=600)
layout(matrix(c(1,1,2,3,4,5,6,7), 2, 4, byrow = TRUE)) 
op <- par(mar=c(3,2.5,2,1)+0.1,  mgp=c(1.5,0.5,0), cex.axis=0.5, cex.lab=0.5)
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')

makebxp.affBy.pos(dat.list = affByMutList, 
                  to.compare=c("A", "G"),
                  mylabels=main.bs.number[-c(1,2,3,8,9,10)])

makebxp.affBy.pos(dat.list = affByMutList, 
                  to.compare=c("A", "C"),
                  mylabels=main.bs.number[-c(1,2,3,8,9,10)])

makebxp.affBy.nucleotide(dat=affByMutList,
                         ylab="Relative affinity")
par(op)
dev.off()

# Supplementary Figure 4 --------------

# Panel A

#PSAM


# Panel B

png(paste("S4-RelAffByPosition.png"), width=183, height=100, units = 'mm', res=600)
par(mfrow=c(2,4),
    pty="s",
    mgp=c(2,1,0),
    cex=0.5,
    cex.lab=1,
    cex.axis=1,
    oma=c(2,2,1,1),
    mar=c(3,3,2,1),
    xpd=TRUE)
make.AffPlot.MutKmer(affByMutList, by.Nuc = TRUE, xlab="position", ylab="")
make.AffPlot.MutKmer(affByMutList, by.Kmer = TRUE, xlab="position", ylab="")
mtext(text="Relative affinity of reference sequence", side=1, line = 0.5, xpd=TRUE, outer=TRUE, cex=0.6)
mtext(text="Relative affinity of alternative sequence", side=2, line = 0, xpd=TRUE, outer=TRUE, cex=0.6)
dev.off()



# SHAPE FIGURES ----------------------


# Make shape heatmap (all bins equal)

plotHeatmapTogether(arr=binned.shape2, shape.par="") -> gl
png("Figure4-Shape-PanelA.png", width = 183, height= 90, res=300, units='mm')
grid.arrange(grobs=gl, ncol=4, clip=TRUE)
dev.off()


# Compare shape across groups

png("Fig4-Shape-PanelsBC.png", width=180, height=150,  res=300, units='mm')
layout(matrix(c(1:8,9,9,9,9), nrow = 3, ncol = 4, byrow = TRUE), heights = c(0.4, 0.4, 0.2))
par(mgp=c(1.45,0.65,0), mar=c(3,3,1.75,0.25)) #,cex.axis=0.65)
getShapeBxpComparison(shape=all.shape[c(1:4)], region=central4.cols)
lapply(names(allshape.df1[1:4]), plot.high.vs.low.aff.sites)
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottomleft", xpd=TRUE, horiz = TRUE,
       legend = c("High affinity (W4)", "Low affinity (W4)",
                  "High affinity (non-W4)", "Low affinity (non-W4)"),
       lty=c(1,3,1,3), lwd=2, col = c( "deeppink3", "deeppink3", "cyan4", "cyan4"))


dev.off()







