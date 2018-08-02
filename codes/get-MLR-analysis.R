# Analyze-mlr-summary


# SET PARAMETERS FOR ANALYSIS ------------
setwd(subDir$ana.path)
mlrdir <-   "mlrresult"
mlr    <- as.matrix(read.table(file.path(
                    subDir$inpath, mlrdir , "summary.txt"), 
                    header = TRUE, row.names = 1))


source(file.path(subDir$code.path, "mlr-functions.R"))


# FILTER DATA TO PLOT

# select data to plot
#nosy <- mlr[- grep(pattern='symm', rownames(mlr)), ] # non-symmetric data
toplot <- mlr[- grep("fixed", rownames(mlr)),]
rownames(toplot) <- sub("^mef2b(.*)(.)aff(0.*)M._(.*)_10_2$", "\\4_\\2", rownames(toplot))
rownames(toplot)


# Compare mlr results for selected columns for specified sample.
op <- par(pty='s', mfrow=c(1,2))
getSeqvsShape(mymlr = toplot, models=c("X1mer", "X1n2shape"))
getSeqvsShape(mymlr = toplot, models=c("X1mer", "X1mer.1n2shape"))
par(op)


# Barplot for multiple motifs.
getMLRbarplot(dat=toplot)





# Final MLR Figures ------------------------------


# Figure 7 ---


png(file.path(subDir$ana.path, "Fig7-MLR-FeatureSelec.png"), width=180, height = 100, res=300, units = 'mm')
layout(matrix(c(1,2,3,2,4,5), nrow = 2, ncol = 3, byrow = FALSE))
par(mar=c(4,4,1,2), mgp=c(2,0.7,0))

getSeqvsShape(mymlr = toplot, models=c("X1mer", "X1mer.1n2shape"))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
courier.legend(dat = toplot, myinset=c(-0.1,0), ncol=3)

op <- par(family='mono')
legend("topright", inset = c(0.4,0.5),
       legend=c("WWWW", "NNNN", "AAAA"), pch=c(0, 1,2), bty="n")
par(op)

getSeqvsShape(mymlr = toplot, models=c("X1mer", "X1n2shape"))

getPlotByPositionMLR(toplot = toplot, model.index=6) # YTAN4TAR_2 = 6

dev.off()




# Figure S7

png(file.path(subDir$ana.path, "S7-MLRmodels.png"), width=180, height = 140, res=300, units = 'mm')
layout(matrix(c(1,2,3,4,4,4), nrow = 2, ncol = 3, byrow = TRUE),
       heights = c(0.4, 0.4), widths = c(0.2, 0.2, 0.1))
par(mar=c(6,5,3,2))

getComparisonMLR(mlr=toplot, model.index = 1)
mtext(text= myFigLab[1], side=3, line=1, at=-0.7)
op <- par(pty='s', mgp=c(4,1,0), mar=c(5,5,3,2))

getPerformanceCountplot(dat=toplot)
mtext(text= myFigLab[3], side=3,  line=1, at=-0.3)
par(op)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
courier.legend(dat=toplot, myinset=c(-0.9,0))

getMLRbarplot(dat=toplot)
mtext(text= myFigLab[2], side=3, line=1, at=-0.9)

dev.off()



# S8 - Feature selection

png(file.path(subDir$ana.path, "S8-featSelection.png"), width=180, height = 180, res=300, units = 'mm')
par(mfcol=c(2,2), oma=c(0,2,1,0))

getPlotByPositionMLR(toplot = toplot, model.index=1)
mtext(text= "A", side=3, outer=TRUE, at=0 )#adj=0)

getPlotByPositionMLR(toplot=toplot, model.index=6)
mtext(text="B", side=3, outer=TRUE, at = 0.5)

dev.off()



