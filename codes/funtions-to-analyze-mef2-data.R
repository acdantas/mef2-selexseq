library(seqinr)
library(seqLogo)
library(Biostrings)
library(DNAshapeR)
library(lattice)
library(gplots)
library(gridGraphics)
library(grid)
library(gridExtra)
library(seqinr)
library(RColorBrewer)

mergeAllBySample <- function(n, arg1, arg2) {
  # x = one of sample names
  # arg1 = aff data list for each sample/round
  # arg2 = which output is desired (count or afifnity)
  # order by Affinity -R2 -R1
  # get counts/aff for kmer from aff.data
  # merge count.data based on sample "k04e, wt..."
  
  indexes <- grep(pattern=n, names(arg1), ignore.case=TRUE)
  to.plot <- arg1[indexes]
  merged <- Reduce(function(x, y) 
    merge(x, y, by="Kmer", all=TRUE), to.plot) 
  merged <- merged [with (merged, order(-merged[7], -merged[2])), ] # order by affinity
  print(head(merged))
  merged <- merged[, c(1, arg2, (arg2+5) )] 
  colnames(merged) <- as.vector(paste(c("Kmer", rounds)))
  return(merged)
}




filterAffTb <- function(round, name, dat, aff=0, regex, mymismatch, get.tbR=FALSE) {
  
  toMatch.list <- paste0(round, ".*", name)
  lapply(toMatch.list, function (toMatch) {

    my.name <-  filterByName(toMatch=toMatch, dat=dat)
    my.df <- getAfftbByName(dat=dat, df.name=my.name,ColNames=c("Kmer", "Affinity","ObservedCount"), aff=aff)
    ref.df <- filterByAgrep(dat=my.df, regex=regex, sample.name=my.name, aff=aff, mymismatch=mymismatch, get.tbR=get.tbR )
    if (get.tbR==TRUE) {
      return(ref.df)
      
    }
    
  })
  
}


filterByAgrep <- function(dat, regex, sample.name, aff, mymismatch, get.tbR) {
  lapply(seq_along(regex), function(y) { 
    by.kmer <- dat[ agrep(regex[y], dat[,"Kmer"] , 
                          list(cost=mymismatch, substitutions=mymismatch, 
                               deletions=0, insertions=0),  fixed=FALSE) ,]
    print(regex[[y]])
    
    if (grepl( "\\^", regex[[y]]) == TRUE) {
      print("Won't find complement")
    }
    else {
      rev.comp.regex <- as.character(reverseComplement(DNAString(regex[[y]]))) #[1]
     
      if (regex[[y]] != rev.comp.regex) {
        by.kmer.comp <-  dat[ agrep( rev.comp.regex , dat[,"Kmer"] , 
                                     list(cost=mymismatch, substitutions=mymismatch, 
                                          deletions=0, insertions=0),  fixed=FALSE) ,]
        
        by.kmer <- unique(rbind(by.kmer, by.kmer.comp))
        
      }
      
    }
    
    towrite <- by.kmer[order(- by.kmer["Affinity"]) ,]
    file.name <- paste(analysis.name, sample.name, paste0("mismatch",mymismatch , "aff", sub("\\.", "", aff), "M", y),  names(regex)[y], "10", "2.txt" , sep="_")
    if (get.tbR==TRUE) { 
      return(towrite)
    }
    print(paste("Writing file: " , file.path(subDir$ana.path, "output", file.name)))
    write.table(x=towrite, file=file.path("./", file.name), quote=FALSE, sep=" ", row.names=F, col.names=F)
    
  })
}

getAfftbByName <- function (dat, df.name, ColNames= c("Kmer", "Affinity","ObservedCount"), aff=0) {
  df <- dat[[df.name]][, ColNames]
  df <- df[df$Affinity > aff ,]
  return(df)
}

filterByName <- function(toMatch, dat) {
  my.match <- grep(toMatch, names(dat), value=TRUE, ignore.case=T)
  return(my.match)
}



getShapePred <- function(dat) {
  fa.file <- paste0("seq" ,".fa")
  # print(dat$Kmer)
  write.fasta(sequences=as.list(dat$Kmer), 
              names=paste0(dat$Kmer, "_",  dat$Affinity), 
              file.out=fa.file)
  prediction <- getShape(filename=fa.file)
  return(prediction)
}




plotstripchart <- function(dat) {
  
  stripchart(dat$Affinity ~ dat$group,
             method="jitter", ylab="Relative affinity",
             xaxt='n',
             jitter=0.25,  pch=16, cex=0.5, vertical=TRUE, las=2,
             col=rainbow(length(unique(dat$group))),
             ylim=c(0.1,1)) #, cex.lab=0.7)
  
  
  
  my.labels <-  lapply(unique(dat$group), function(x) {
    re.l1.names[[x]]
  })
  str(my.labels)
  my.labels <- do.call(c, my.labels)
  print(my.labels)
  
  
  axis(1, at=seq(1, length(unique(dat$group)), by=1),
       las=2, #cex=0.7,
       family="mono",
       labels= my.labels, xpd=TRUE )

  
}

bin.shape <- function(x, n.bin) {
  f <- rep( c(1:(nrow(x)/n.bin+1)), each=n.bin, length.out=nrow(x))
  binned <- apply(x, 2, function(t) tapply(t, f, mean, na.rm=TRUE))
  return(binned)
}


getShapeHeatMap <- function(dat, name, n.bins=0, return.binned.data=F, per.bin=45) {
  # use n.bin=40 for 40 bins total; (if comparing heatmaps, different bin size could occur)
  # use n.bin=0 and per.bin=100, for 100 sequence in each bin, regardless of  dataset size 
  
  if (nrow(dat$MGW) >= 200 & n.bins!=0) {
    print(paste(name, "is greater than 1000: ", nrow(dat$MGW), ". Will plot per bin"))
    #bins=500 #500 in each bin ; #
    in.bins= (nrow(dat$MGW)/n.bins)
    #in.bins= per.bin
    print(paste("will return", in.bins, "sequences per bin"))
    binned.data <- sapply(dat, function(t) bin.shape(t, n.bin=in.bins))  
    if(return.binned.data==T) { return(binned.data)}
    plotMyHeatmaps(binned.data, sample.name=name)
  }
  
  if (nrow(dat$MGW) >= 200 & n.bins==0) {
    print(paste(name, "is greater than 1000: ", nrow(dat$MGW), ". Will plot per bin"))
    #bins=500 #500 in each bin ; #in.bins= (nrow(dat$MGW)/n.bins)
    in.bins= per.bin
    binned.data <- sapply(dat, function(t) bin.shape(t, n.bin=in.bins))  
    if(return.binned.data==T) { return(binned.data)}
    plotMyHeatmaps(binned.data, sample.name=name)
  }
  if (nrow(dat$MGW) < 200 & n.bins==0) {
    print(paste(name, "is smal: ", nrow(dat$MGW)))
    if(return.binned.data==T) { return(binned.data)}
    plotMyHeatmaps(dat, sample.name=name)
  }
  if (nrow(dat$MGW) < 200 & n.bins!=0) {
    print(paste(name, "is smal: ", nrow(dat$MGW) , " but will BIN"))
    in.bins= (nrow(dat$MGW)/n.bins) 
    binned.data <- sapply(dat, function(t) bin.shape(t, n.bin=in.bins)) 
    if(return.binned.data==T) { return(binned.data)}
    plotMyHeatmaps(binned.data, sample.name=name)
  }
  # else {
  #   print(paste(name, "is smal: ", nrow(dat$MGW)))
  #   plotMyHeatmaps(dat, sample.name=name)
  # }
}
#..


plotMyHeatmaps <- function(each.shape.list, sample.name, save.plot=T) {
  if (save.plot==TRUE){
    
    #mycols <- head(rainbow(1000), n=600)
    mycols <- head(heat.colors(1000), n=950)
    
    lapply(seq_along(each.shape.list), function(x) {
      
      my.par <- names(each.shape.list)[x]
      my.dat <- each.shape.list[[x]]
      print(my.par)
      print(head(my.dat))
      # width = 46, height= 90, res=300, units='mm'
      
      png(paste0(sample.name, "-", my.par, ".png", sep="" ),
          width = 50, height= 90, res=1200, units='mm')
      
      #par( mgp=c(1.5, 0.5, 0), las=0)
      heatmap.2(my.dat, dendrogram ='none',
                Rowv=FALSE, Colv=FALSE, labRow="",
                labCol = x.axisShapeLabels[[my.par]],
                margins=c(3,0), cexCol = 0.7, 
                symm=FALSE, symkey=F, symbreaks=F,
                col=mycols,
                density.info='none', trace='none',
                
                # key.par=list(mgp=c(1, 0.25, 0),
                #               cex.axis=0.6, cex.lab=0.6,
                #               mar=c(3, 5.5, 3,5.5 )), 
                # 
                # key.title= NA,  key.xlab=short.shapelabs[[my.par]],
                # lmat = rbind(c(3,4),c(2,1)) , lhei=c(1.5, 4), lwid=c(0.2,4),
                # 
                key.par=list(mgp=c(0.75, 0.25, 0),
                             cex.axis=0.5, cex.lab=0.6,
                             # mar=c(2.5, 5, 3.5, 5 )),
                             mar=c(2.5, 3, 3.5, 3 )), 
                
                key.title= NA,  key.xlab= short.shapelabs[[my.par]], 
                lmat = rbind(c(3,4),c(2,1)) , lhei=c(1.5, 4), lwid=c(0.2,4),
                
                colsep=1:10)
      dev.off()
      
    }) 
    
  }
  
}



grab_grob <- function(){
  grid.echo()
  grid.grab()
}


plotHeatmapTogether <- function(arr, shape.par) {
  lapply(seq_along(arr), function (i) {
    
    my.par <- names(arr)[[i]]
    mydf <- arr[[i]]
    
    par( mgp=c(1.5, 0.5, 0), las=0)
    #  par(cex.axis=0.5, cex.names=0.5, cex.lab=0.5)
    heatmap.2(mydf, dendrogram ='none',
              Rowv=FALSE, Colv=FALSE, labRow="",
              margins=c(3,0), cexCol = 0.7, 
              symm=FALSE, symkey=F, symbreaks=F,
              #col=head(rainbow(1000), n=600),
              labCol = x.axisShapeLabels[[my.par]],
              density.info='none', trace='none',
              key.par=list(mgp=c(1, 0.5, 0),
                           cex.axis=0.7, cex.lab=0.7,
                           mar=c(6, 15, 1, 15)), #mar=c(2.5, 2.5, 1, 0))
              key.title= NA,  key.xlab=short.shapelabs[[names(arr)[i]]],
              # keysize=1.2,
              lmat = rbind(c(3,4),c(2,1)) , lhei=c(1.5, 4), lwid=c(0.2,4),
              colsep=1:10)
    
    #   title(main=paste(names(arr)[i]), cex.main=0.7)
    grab_grob()
  })
}


getShapeBxpComparison <- function (shape, region) {
  lapply(seq_along(shape), function (x) {
    
    my.par <- names(shape)[[x]]
    my.dat <- shape[[x]]
    plotRowMeansShape(high.aff=head(my.dat, n=1000),
                      low.aff=tail(my.dat, n=1000),
                      n=my.par,
                      region=region, with.pvalue=T, myplot="boxplot")
    
    # mtext(text = shape.labels[[my.par]], side = 2, line=2, cex=1.2)
  })
}



rm(getBoxplotComparison)

# get.Pvalues.U <- function(x,y) {
#   #With mann-whitney test
#   pvals.u <- c()
#   for (i in 1:ncol(x)) {
#     pvals.u <- c(pvals.u, wilcox.test(x[,i], y[,i], paired=FALSE)$p.value)
#   }
#   return(pvals.u)
# }

getShapeBxpComparison <- function (shape, region) {
  lapply(seq_along(shape), function (x) {
    
    my.par <- names(shape)[[x]]
    my.dat <- shape[[x]]
    plotRowMeansShape(high.aff=head(my.dat, n=1000),
                      low.aff=tail(my.dat, n=1000),
                      n=my.par,
                      region=region, with.pvalue=T, myplot="boxplot")
    
    # mtext(text = shape.labels[[my.par]], side = 2, line=2, cex=1.2)
  })  
}

plotRowMeansShape <- function (high.aff, low.aff, with.pvalue=FALSE, region, myplot, n) {
  

  if (n %in% c("HelT", "Roll")) {
    region <- region[- length(region)] 
  }
  
  print(paste(n, " ", region))
  high.aff <- rowMeans(high.aff[, region])
  low.aff <- rowMeans(low.aff[, region])
  
  if(myplot=="boxplot") {
    print(par('cex'))
    boxplot(cbind(high.aff, low.aff ),
            boxwex=0.8, outline=F, #col=c("mistyrose1", "lightcyan1"),
            ylab= short.shapelabs[[n]], #cex.axis=1.2, cex.lab=1.2,
            xlab="",
            names= c("High", "Low")
    )

    
  }

  
  if (with.pvalue==TRUE) {
    p.val <- wilcox.test(high.aff, low.aff, paired=FALSE)$p.value
    print(paste("pval1:", p.val))
    prettypvalue(pval=p.val)
    

  }
  
  
  
}


prettypvalue <- function(pval)  {
  
  print(str(pval))
  mypval1 <- format.pval(as.numeric(format(pval, scientific = TRUE, digits=3)))
  print(paste("mypval1: ", mypval1))
  print(str(mypval1))
  
  if (pval < 2.2e-16) {
    print("Will print smaller than 2.2x10^16")
    mtext( bquote(italic('P')~'<'~'2.22 x '*10^-16),
           side=3, line=0.25, at=1.5, cex=par('cex')*0.725 )
    
  }
  
  if (pval > 2.2e-16 & pval < 0.001) {
    print("still significant and will print")
    mypval2 <- sub(pattern='e(.*)', replacement = ' x ', x=mypval1, perl = TRUE)
    myex <- sub(pattern='(.*)e(.*)', replacement = '\\2', x=mypval1, perl=TRUE)
    mtext( bquote( italic('P')~'='~.(mypval2)*10^.(myex)), side=3, line=0.25, at=1.5, cex=par('cex')*0.725 )
  }
  
  if(pval > 0.001) {
    mtext( "N.S.", side=3, line=0.25, at=1.5, cex=par('cex')*0.725 )
    
  }
  print("DONE")
  
}


ana.shape.for.df <- function(dat, name, get.binned.data=F ) {
  # dat = df with Kmer anf Affinity cols
  # name = name of motif or how to save file
  
  print( head(dat[order(- dat["Affinity"]),], head=30))
  print(tail(dat[order(- dat["Affinity"]),], n=30))
  getShapePred(dat[order(- dat["Affinity"]),]) -> shape.dat
  
  # getShapeHeatMap(dat = shape.dat, name = paste0(name, "-ShapeHeatmap") ,
  #                 return.binned.data = get.binned.data,
  #                 # n.bins=0, per.bin=25 # fixed number of seq per bin
  #                 n.bins=35, per.bin=0 # fixed number of bins
  # )
  if(get.binned.data == FALSE ) {
    return(shape.dat)
    
  }
  
}

make.legend <- function(my.location,
                        my.text,
                        my.col,
                        lwd, lty,
                        ...) {
  legend(my.location,
         legend=my.text,
         fill=my.col,
         inset=0,
         xpd=TRUE,
         lty=lty,
         lwd=lwd,
         horiz = TRUE)
  
  
  
}

addShapeAxis <- function(mypar, mylabels) {
  
  
  if (mypar %in% c("HelT", "Roll")) {
    #mylabels <- x.axis.shortShapeLabels[[3]]
    axis(side = 1, labels = mylabels, at=1.5:(length(mylabels)+0.5), las=2)
  }
  
  if(mypar %in% c("MGW", "ProT")) {
    mylabels <- mylabels[c(-1, -(length(mylabels)))]
    axis(side = 1, labels = mylabels, at=1:(length(mylabels)), las=2)
  }
  
  
}

main.bs.number <- c(paste("\U2212", seq(5,1,-1), sep=""),
                    paste("+", seq(1,5,1), sep=""))


plot.high.vs.low.aff.sites <- function(x,
                                       df1=allshape.df1,
                                       df2 = allshape.df2,
                                       df.names = c("W4 at center", "non-W4 at center"),
                                       n.seq=100,
                                       with.legend=FALSE)  {
  my.ldw <- 2
  mycol <- c( "deeppink3", "cyan4")
  mypar <- x
  
  df1 <- df1[[x]]
  df2 <- df2[[x]]
  
  
  if(mypar %in% c("MGW", "ProT")) {
    df1 <- df1[, colSums(is.na(df1)) != nrow(df1)]
    df2 <- df2[, colSums(is.na(df2)) != nrow(df2)]
    
    h.df1 <- head(df1, n=n.seq)
    l.df1 <- tail(df1, n=n.seq)
    h.df2 <- head(df2, n=n.seq)
    l.df2 <- tail(df2, n=n.seq)
    print(paste0("Ncol of ", mypar, "is", ncol(h.df1)))
    
    
    plot(colMeans(h.df1, na.rm = T),
         type="l", col=mycol[1], lwd=my.ldw,
         ylim = ylim.shape[[x]],
         xaxt="n",
         xlab= "Binding site position",
         ylab = short.shapelabs[[x]])
    
    rm.rows <- c(1,6) # use shape of central 4 "WWWW"
    print(rm.rows)
    pval.df1 <- wilcox.test(rowMeans(h.df1[, - rm.rows], na.rm = T), rowMeans(l.df1[, - rm.rows], na.rm = T), paired=FALSE)$p.value #within W4
    pval.df2 <- wilcox.test(rowMeans(h.df2[, - rm.rows], na.rm = T), rowMeans(l.df2[, - rm.rows], na.rm = T), paired=FALSE)$p.value # within non-W4
    pval.HighAff <- wilcox.test(rowMeans(h.df1[, - rm.rows], na.rm = T), rowMeans(h.df2[, - rm.rows], na.rm = T), paired=FALSE)$p.value #b/e high aff
    pval.LowAff <- wilcox.test(rowMeans(l.df1[, - rm.rows], na.rm = T), rowMeans(l.df2[, - rm.rows], na.rm = T), paired=FALSE)$p.value # b/e low aff
    
  }
  
  if (mypar %in% c("HelT", "Roll")) {
    h.df1 <- head(df1, n=n.seq)
    l.df1 <- tail(df1, n=n.seq)
    h.df2 <- head(df2, n=n.seq)
    l.df2 <- tail(df2, n=n.seq)
    print(paste0("Ncol of ", mypar, "is", ncol(h.df1)))
    
    plot(colMeans(h.df1, na.rm = T),
         type="l", col=mycol[1], lwd=my.ldw,
         ylim = ylim.shape[[x]],
         xaxt="n",
         xlab= "Binding site position",
         xaxs='i',
         ylab = short.shapelabs[[x]])
    
    rm.rows <- c(1,2,3,7,8,9) # use shape of central 4 "WWWW"
    print(rm.rows)
    pval.df1 <- wilcox.test(rowMeans(h.df1[, - rm.rows], na.rm = T), rowMeans(l.df1[, - rm.rows], na.rm = T), paired=FALSE)$p.value #within W4
    pval.df2 <- wilcox.test(rowMeans(h.df2[, - rm.rows], na.rm = T), rowMeans(l.df2[, - rm.rows], na.rm = T), paired=FALSE)$p.value # within non-W4
    pval.HighAff <- wilcox.test(rowMeans(h.df1[, - rm.rows], na.rm = T), rowMeans(h.df2[, - rm.rows], na.rm = T), paired=FALSE)$p.value #b/e high aff
    pval.LowAff <- wilcox.test(rowMeans(l.df1[, - rm.rows], na.rm = T), rowMeans(l.df2[, - rm.rows], na.rm = T), paired=FALSE)$p.value # b/e low aff
    
  }
  
  
  addShapeAxis(mypar = mypar,
               mylabels=main.bs.number[-c(1,10)])
  
  # low.aff W4
  lines(colMeans(l.df1, na.rm = T),
        col=mycol[1], lty=3, lwd=my.ldw)
  
  # Add noW4 daa to plot-----
  
  lines(colMeans(h.df2, na.rm = T),   # high.aff noW4
        col=mycol[2], lwd=my.ldw)
  
  lines(colMeans(l.df2, na.rm = T),   # low.aff noW4
        col=mycol[2], lty=3, lwd=my.ldw)
  
  # pval.df1 <- wilcox.test(rowMeans(h.df1, na.rm = T), rowMeans(l.df1, na.rm = T), paired=FALSE)$p.value #within W4
  # pval.df2 <- wilcox.test(rowMeans(h.df2, na.rm = T), rowMeans(l.df2, na.rm = T), paired=FALSE)$p.value # within non-W4
  # pval.HighAff <- wilcox.test(rowMeans(h.df1, na.rm = T), rowMeans(h.df2, na.rm = T), paired=FALSE)$p.value #b/e high aff
  # pval.LowAff <- wilcox.test(rowMeans(l.df1, na.rm = T), rowMeans(l.df2, na.rm = T), paired=FALSE)$p.value # b/e low aff
  #
  
  pval.list <- list(pval.df1, pval.df2, pval.HighAff, pval.LowAff)
  pval.names <- c("High vs. Low (W4)", "High vs. Low (non-W4)", "High W4 vs. non-W4", "Low W4 vs. non-W4")
  
  pval.lines <- list(one=c(1,3), two=c(1,3), three=c(1,1), four=c(3,3))
  pval.cols <- list(one=c(mycol[1]), two=mycol[2], three=mycol, four=mycol)
  
  print(mypar)
  print(pval.list)
  
  my.count <- 0
  
  legend("bottomleft", bty='n',
         legend = bquote( italic('P')~'< 0.001'),
         cex=par('cex')*0.725, inset=c(0,0.03*5) )
  
  lapply(seq_along(pval.list), function(y) {
    
    my.val <- pval.list[y]
    
    if (my.val < 0.001) {
      
      my.count <<- my.count +1
      #legend("bottomleft", paste(pval.names[y], "*"), bty='n', inset = c(0,(0.05*my.count)))
      legend("bottomleft",
             legend=c("", ""),
             lty = pval.lines[[y]] ,
             col = pval.cols[[y]],
             horiz = TRUE, cex=0.5,
             lwd = my.ldw,  bty='n', inset = c(0,(0.03*my.count)))
    }
  })
  
  #legend("topleft", "A)", bty="n")
  if(with.legend == TRUE) {
    
    legend("bottomright",
           legend=df.names,
           inset = c(0, 0.12),
           bty='n', fill=c(mycol))
    legend("bottomright",
           legend=c("High aff. sites", "Low aff. sites"),
           inset = c(0.03,0), xpd=TRUE,
           bty='n', lty=c(1,3), lwd=my.ldw)
    
  }
  
  
}
