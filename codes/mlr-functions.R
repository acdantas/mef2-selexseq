# MLR functions

mlr.lab <- list( X1mer = "1mer",
                 X2mer = "2mer",
                 X3mer = "3mer",
                 X1n2shape= "shape",
                 X1mer.1n2shape = "1mer + shape",
                 X2mer.1n2shape = "2mer + shape",
                 X3mer.1n2shape = "3mer + shape"
                 )



##
plotFeatureByPosition <- function(
  dat, 
  sample.name, 
  minus1mer=FALSE, 
  remove2Shape=FALSE, 
  minus1mer2shape=FALSE) {
  # remove or add 1st+2nd order shape position by position
  # minus1mer = add shape by position; remove 1mer by position
  # removeShape = from removed shape position by position,
  
  #op <- par( cex.axis=1.2, cex.lab=1.2, cex.main=0.8, mar=c(4.5,6,4,1.5)+0.1)
  op <- par( mar=c(2,4,1,1)+0.1)
  if(minus1mer==TRUE) {
    
    # X1mer.1n2shape.Px - X1mer
    my.positions <- dat[ grep(paste0("X1mer.1n2shape", "P"), colnames(mlr))]
    seq <- dat["X1mer" ]
    print(my.positions)
    print(seq)
    toplot <- my.positions - seq
    bpNames<- c("-3", "-2", "-1", "+1", "+2", "+3")
    print("CEX plotfeatbypostion")
    print(par('cex.axis'))
    barplot(toplot, 
            ylim=c(0,0.07),
            las=1, main=sample.name, names=bpNames, 
            col="honeydew", cex.names=par('cex.axis')-0.1,
            ylab="", cex.main=par('cex.axis'))
    mtext(expression(Delta*italic('R')^"2"~"("*italic('R')^"2"*"" ["seq+shape"]*
                                    "" [""[italic('i')]]~"-"
                                  ~italic('R')^"2"*"" ["seq"]*")"),
          side=2, line=3.5, cex=par('cex'), xpd=TRUE)
   # mtext(c("-3", "-2", "-1", "+1", "+2", "+3"), 
         
  }
  
  if (remove2Shape==TRUE) {
    
    # X1n2shape.noP - X1n2shape
    my.positions <- dat[ grep("X1n2shape.noP", colnames(mlr))]
    seq <- dat["X1n2shape" ]
    toplot <- my.positions - seq
    print(my.positions)
    print(seq)
    bpNames<- sub("X1n2shape.no", "", names(my.positions))
    barplot(toplot, 
            ylim=c(-0.05, 0),
            las=2,  names=bpNames, 
            col="honeydew", xaxt='n', 
            ylab="")
    mtext(expression(Delta*italic('R')^"2"~"("*italic('R')^"2"*"" ["shape-shape"]*
                                     "" [""[italic('i')]]~"-"
                                   ~italic('R')^"2"*"" ["shape"]*")"),
    side=2, line=3.5, cex=par('cex'))
                       

    #ylab=expression(paste(Delta, italic('R')^"2", "", " (", italic('R')^"2", "" ["shape - shape"], " " [italic('i')] , ' - ', italic('R')^2, "" ["shape"], " )")))
    #mtext(expression(paste(Delta, italic('R')^"2", "", " (", italic('R')^"2", "" ["shape - shape"], " " [italic('i')] , ' - ', italic('R')^2, "" ["shape"], " )")), side = 2, line=4)
    
  }
  
  if (minus1mer2shape==TRUE) {
    
    # X1mer.1n2shape.Px - X1mer.1n2shape
    my.positions <- dat[ grep(paste0("X1mer.1n2shape", "P"), colnames(mlr))]
    seq <- dat["X1mer.1n2shape" ]
    toplot <- my.positions - seq
    print("To plot is: ")
    print(toplot)
    bpNames<- sub("X1mer.1n2shape", "", names(my.positions))
    print(bpNames)
    barplot(toplot, 
            las=2,  main=sample.name, names=bpNames, 
            col="honeydew", xaxt='n', ylab='')
    #ylab=expression(paste(Delta, italic('R')^"2", "", " (", italic('R')^"2", "" ["1mer + shape i"] , ' - ', italic('R')^2, "" ["1mer + shape"], " )")))
    #mtext(expression(paste(Delta, italic('R')^"2", "", " (", italic('R')^"2", "" ["1mer + shape i"] , ' - ', italic('R')^2, "" ["1mer + shape"], " )")), side = 2)

  }
  par(op)
}
#..
####


getPlotByPositionMLR <- function (toplot, model.index) {
  
  model       <- toplot[model.index,]
  model.name  <- rownames(toplot)[model.index]
  
  plotFeatureByPosition(dat=model, sample.name=model.name, minus1mer=TRUE)
  plotFeatureByPosition(dat=model, sample.name=model.name, remove2Shape=TRUE)
  model["seqCount"]
  
}
#..


plotComparisonMLR <- function(dat, sample, mylist, singleplot, mytitle) {
  
  if (singleplot==TRUE) {
    op <- par(mar=c(5,5,2,2)+0.1)
    print("PlotComparsion cex")
    print(par('cex'))
    dat <- dat[mylist]
    barplot(t(dat), 
            names.arg = as.vector(unlist(mlr.lab[mylist])),  
            ylab=expression(paste(italic('R')^"2")), 
            las=2, col="honeydew") #,main=mytitle)
    mtext(mytitle, side=2, line=4, xpd=TRUE, cex=par('cex'))
    par(op)
  }
  
  else {
    dat <-  dat[,mylist]
    pdf(paste0(sample, ".pdf"))
    barplot(t(dat), beside=TRUE, ylab=expression(paste(italic('R')^"2")), las=1,
            legend.text=TRUE, angle=120,col="honeydew", main=mytitle,
            args.legend= list(title = "Features", x = "topleft", cex = 1))
    # bp <- barplot(t(dat), beside=TRUE, ylab="R2", las=1, legend.text=TRUE, angle=120,args.legend= list(title = "Features", x = "topleft", cex = 1))
    #text(bp, 0, round(t(dat), 2),cex=1,pos=3)
    dev.off()
  }
}
# ..







##
getComparisonMLR <- function (mlr, model.index, ...) {
  
  model       <- mlr[model.index,] # m1
  model.name  <- rownames(mlr)[model.index]
  model.feat  <- c("X1mer", "X2mer", "X3mer", "X1mer.1n2shape",
                   "X2mer.1n2shape", "X3mer.1n2shape")
  
  model <- model[names(model) %in% model.feat] # removes non-existing features
  
  plotComparisonMLR(model,
                    sample="wt-seqKmer",
                    mylist= names(model),    # model.feat,
                    singleplot=TRUE,
                    mytitle=model.name)
  model["seqCount"]
  
}
# ..



##
getMLRbarplot <- function(dat) {
  
  barplot(t(dat[, c("X1mer", "X1mer.1n2shape", "X1n2shape")]),
          col = rep(rainbow(nrow(dat)), each=3),
          beside=T, ylab=rsq.lab,
          names=rep(c("1mer", "1mer+shape", "shape"),nrow(dat)),
          las=2)#,  cex.names=0.7)
  # legend("topright", legend = rownames(mlr), fill = rainbow(nrow(mlr)), bty = "n",
  #        cex=0.7 , xpd=TRUE, inset=c(-0.25,-0.35))
}
#..


##
getPerformanceCountplot <- function(dat) {
  
  plot(dat[,"X1mer.1n2shape"], dat[,1],
       pch=16,
       las=2,
       yaxs='i', xaxs='i',
       col=rainbow(nrow(dat)),
       xlim=c(0,1), ylim=c(0, 30000),
       xlab="Model performance (1mer+shape)",
       ylab="Sequence count")
  
  #text((mlr[,"X1mer.1shape"]+0.05), mlr[,1], labels = rownames(mlr), cex=0.4, las=3)
}
#..

getSeqvsShape <- function(mymlr, models) {
  
  mydat <- mymlr[, models]
  mycol <- rainbow(nrow(mymlr))
  
  op <- par(pty='s')
  plot(mydat,
       col = mycol , pch=16,
       yaxs='i', xaxs='i',
       xlab=bquote( italic('R')^"2"~"("*.(mlr.lab[[models[1]]])*")" ),
       ylab= bquote( italic('R')^"2"~"("*.(mlr.lab[[models[2]]])*")" ),
       ylim=c(0,1), xlim=c(0,1)
  )
  
  points( mymlr[grepl("NNNN", row.names(mymlr)),][, models]   , pch=1, col="black")
  points( mymlr[grepl("WWWW", row.names(mymlr)),][, models]   , pch=0, col="black")
  points( mymlr[grepl("AAAA", row.names(mymlr)),][, models]   , pch=2, col="black")
  #  legend("bottomright", legend=c("WWWW", "NNNN"), pch=c(0, 1), cex=1, bty="n")
  
  abline(a=0,b=1,col=1)# ,lwd=2)
  par(op)
  
}




courier.legend <- function(dat=mlr, myinset, ncol=1, cex=1) {
  mycol <- rainbow(nrow(dat))
  
  op <- par(family="mono", xpd=NA)
  legend("topleft",
         legend = rownames(dat),
         # inset=c(-0.9,0),
         inset=myinset,
         bty='n',
         ncol=ncol,
         cex=cex,
         border = "white",
         fill = mycol)
  par(op)
  
}


