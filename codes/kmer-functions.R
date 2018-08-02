library(RColorBrewer)


groupKmerRegex <- function(dat, regex, title, r=F)  {
  # For a given regex list, find kmers that match regex
  
  each.df <- dat
  print(paste0("Initial dataframe  ", title,  "  nrow is: ", nrow(each.df)))
  toPlot <- data.frame()
  
  lapply(seq_along(regex), function(y) { # look for all different regex
    by.kmer <- subset(each.df, grepl(regex[[y]], Kmer,  perl=TRUE))
    by.kmer$group <- names(regex)[y]
    toPlot <<- rbind(toPlot, by.kmer)
  })
  
  #Final dataset to plot separated by groups
  if (r==T) {
    return(toPlot)
  }
  #plotstripchartAffinitybyGrp(df=toPlot, regex=regex, sample.name=title)
  
}


plotKmerStripChart <- function(df, Round, sample.name) {
  # plot Kmer Stripchart for
  # df with $Kmer, $R2 (affinity), and sample name
  
  df$group<-factor(df$group, levels=names(regex.l)) # set order of levels
  rainbowcols <- rev(rainbow(length(names(regex.l)) , s = 0.7, alpha=0.9)) # rainbow colors to be used
  stripchart(get(Round) ~ group, data=df, xaxt="n",  # stripchart
             method="jitter", #main=sample.name, 
             jitter=0.25,  pch=20, vertical=TRUE, las=2,
             col=rainbowcols, ylim=c(0,1), ylab=toupper(sample.name))
}
#..

##
setPlotKmer <- function(dat="", regex="", nsamples) {
  png(paste("stripChart-Kmer", regex, '.png', sep=""), width=600, height=nsamples*(300))
  par(mfcol=c(nsamples,1), mar = c(2,6,2,2) + 0.1, oma=c(12,5,1,0), #mar = c(1,8,1,2) oma=c(20,12,5,1)
            cex.axis=nsamples-(nsamples*0.3), cex.lab=nsamples-(nsamples*0.3), 
            # mgp = c(6, 1, 0),
             cex.main=nsamples-(nsamples*0.3) 

            )
}
#..

##
finishePlotKmer <- function(regex) {
  mtext(text="Relative Binding Affinity", side=2, line=0, 
        outer=TRUE, cex=1.5, xpd=TRUE) #line=4 changes position
  axis(1, at=seq(1, length(regex), by=1), labels = regex, las=2)
  dev.off()
}
#..

##
plotstripchartAffinitybyGrp <- function(df, sample.name="", regex){
  df$group<-factor(df$group, levels=names(regex)) 
  rainbowcols <- rev(rainbow(length(names(regex)) , s = 0.7, alpha=0.9)) 
  stripchart(Affinity ~ group, data=df, xaxt="n",  
             method="jitter", #main=sample.name, 
             jitter=0.25,  pch=20, vertical=TRUE, las=2,
             col=rainbowcols, ylim=c(0,1), ylab=toupper(sample.name))
}
#.. 

## 
getKmerRegexgroup <- function(dat, regex, title, r=F, makestripchart=F)  {
  
  lapply(seq_along(dat), function(y) { #Apply function to each dataframe in list 
    each.df <- dat[[y]]
    toPlot <- data.frame()
    
    lapply(seq_along(regex), function(x) { # look for all different regex
      by.kmer <- subset(each.df, grepl(regex[[x]], Kmer,  perl=TRUE))
      each.df <<- subset(each.df, !grepl(regex[[x]], Kmer,  perl=TRUE))
      by.kmer$group <- names(regex)[x]
      toPlot <<- rbind(toPlot, by.kmer) 
    })
    
    #Final dataset to plot separated by groups
    if(makestripchart==FALSE) {
      print("wont make plot")
    }

    
    if(makestripchart==TRUE) {
      
      plotstripchartAffinitybyGrp(df=toPlot, regex=regex, sample.name=title)
      
    }
    if (r==T) {
      print("Will return data")
      return(toPlot)
    }    
    
    
  })
} 



analyzeBS.ByKmer <- function(dat) {
  
  df <- dat
  df$leftKmer <- gsub(pattern="^(.{3}).{4}.{3}$", replacement = "\\1", x = df$Kmer)
  df$kmer <- gsub(pattern="^.{3}(.{4}).{3}$", replacement = "\\1", x = df$Kmer) #central kmer 4mer
  df$rightKmer <- gsub(pattern="^.{3}.{4}(.{3})$", replacement = "\\1", x = df$Kmer)
  str(df)
  head(df)
  return(df)
  
}

make.3merStripchart <- function(df, myname) {

  stripchart(df[,2] ~ df[,5], xlab="Kmer 5' to central 4mer",
             ylab="Relative affinity", xaxt='n',
             method="jitter", las=2, #cex.axis=0.7,
             jitter=0.25,  pch=16, vertical=TRUE)
  axis(side=1, cex.axis=0.6, las=2,
       at = seq_along(levels(factor(df[,5]))), 
       labels = levels(factor(df[,5])) )
  


  myplot <- stripchart(df[,2] ~ df[,7], xlab="Kmer 3' to central 4mer",
             ylab="Relative affinity", xaxt='n',
             method="jitter", las=2, #cex.axis=0.7,
             jitter=0.25,  pch=16, vertical=TRUE)
  print(names(df[,7]))
  axis(side=1, cex.axis=0.6, las= 2,
       at = seq_along(levels(factor(df[,7]))), 
       labels = levels(factor(df[,7])) )
}

make.smallKmer.Stripchart <- function(dat, myname, mycol) {

  x <- dat[, mycol]
  y <- dat[,"Affinity"]
  print(str(x))
  print(str(y))
  boxplot(y ~ x, xlab= mycol, las=2,
             ylab="Relative affinity"
            
             )

}


plotMotifAff <- function(dat,myxlab) {
  
  par(pty='s')
  darkcols <- brewer.pal(8, "Dark2")
  
  plot(x=NULL,  
       xlim=c(min(dat[,2]), max(dat[,2])),ylab="Alternative Motif",
       ylim=c(min(dat[,-1]), max(dat[,-1])), xlab=myxlab)
  
  lapply(1:ncol(dat[,-1]), function (x) {
    
    mydat <- dat[, (x+1)] 
    points(dat[, 2], mydat, col=darkcols[x], pch=19)
    
  })
  
  legend("topleft", 
         bty='n',
         legend=c("CTAW4TAG", "CTAW4TAA","TTAW4TAA"),
         fill=darkcols[1:ncol(dat[,-1])],
         inset=c(0,-0.4), xpd=TRUE
  )
}

#####


get.Kmerby3mer <- function(dat, grp, col1, col2, ...) {
  
  mylist <- split(dat[c("Affinity", col1)], dat[, col2])
  print(mylist)
  toplot <- Reduce(function(...) merge(..., by=c(col1), all=T), mylist)
  colnames(toplot)[2:ncol(toplot)]  <- names(mylist)
  rownames(toplot) <- toplot[,1]
  toplot <- as.matrix(toplot[,-1])
  print(toplot)
  
  #op <- par(cex.main=0.8)
  makeAffinityHeatmap3(t(toplot))
  #
  # par(op)
  
}



makeAffinityHeatmap3 <- function(x) {
  
  op <- par(family="mono")
  heatmap.2(x,
            cexCol = 0.7, cexRow= 0.7,
            col =colorRampPalette(brewer.pal(9,"GnBu"))(100),
            breaks=seq(0.1,1, length.out = 101),
            key.par = list(mar=c(3,2,0.5,3),
                           mgp=c(1.2,0.5,0),
                           cex.axis=0.7,
                           cex.lab=0.9,
                           
                           family="sans"),
            key.title= "",
            key.xlab="Relative affinity",
            margins=c(3,9),
            density.info='none', trace='none'
            
  )
  par(op)
}



# FUNCTIONS - Figure 3 & S4

getAffForMutKmer <- function(dat=df2, bs.region=c(4:7)) {
  
  myAT.rich.df <- subset(dat, grepl("^.{3}[AT][AT][AT][AT].{3}$", Kmer, perl=T))
  print(paste0("Total size of df is: ", nrow(dat), ". Head of dat is:"))
  
  mynucleotides <- c("A", "T", "C", "G")
  
  lapply(bs.region, function(j) { # for each position within binding site
    
    myregex <- paste0("^(.{", j-1, "})","(.)", "(.{", n-j , "})$" )
    
    lapply(mynucleotides, function(b) { # for each nucletode at position "j"
      
      sapply(myAT.rich.df$Kmer, function(eachkmer) {
        mut <- sub(pattern=myregex,
                   replacement=paste0("\\1", b, "\\3"),
                   x=eachkmer, perl=TRUE)
      }) -> mut.kmers
      
      
      # Select kmers with alternative nucleotide
      dfByNucPos <- dat[dat$Kmer %in% mut.kmers ,]
      dfByNucPos$Kmer <- sub(pattern=myregex,
                             replacement=paste0("\\1", "-", "\\3"),
                             x=dfByNucPos$Kmer, perl=TRUE)
      return(dfByNucPos)
    }) -> mut.df.list
    
    # process list of delta aff to return merged.data for each position
    merged.data.frame <-  Reduce(function(...) merge(..., by=c("Kmer"), all=F), mut.df.list)
    colnames(merged.data.frame)  <- c("Kmer", mynucleotides)
    return(merged.data.frame)
    
  }) -> my.result
  names(my.result) <- paste("Position", bs.region)
  return(my.result )
  
}


make.AffPlot.MutKmer <- function (dat.list,
                                  by.Nuc=FALSE,
                                  by.Kmer=FALSE,
                                  xlab=bquote("Higher aff. CTAnnnnTAG"),
                                  ylab=bquote("Alternative CTAnnnnTAG"),
                                  ...) {
  
  colnames.toPlot <- c("A", "T", "C", "G")
  mypos <- paste("Position ", main.bs.number[-c(1,2,3,8,9,10)])
 
   # 1. For each position in the list
  lapply(seq_along(dat.list), function(i) {
    
    dat <- dat.list[[i]]
    # start plot area
    
    if ( xlab == "position") {
      xlab <- mypos[[i]]
    }
    plot(x=NULL,
         xlim=c(0.8,1),
         ylim=c(0,1),
         xlab=xlab,
         ylab=ylab
    )

    
    # 2. For each mutation
    # Plot by nucleotide OR ...
    if (by.Nuc==TRUE) {
      sapply(1:nrow(dat), function(X) {
        
        myrow <- dat[X,][, names(dat) %in% colnames.toPlot]
        print(myrow)
        points(x=rep(max(myrow), length(myrow)),
               y=myrow,
               pch= c(15,16,17,18),
               col=c(1:4)
        )
        print("Done with plot")
        
      })
      legend("topleft", legend = colnames.toPlot, col = c(1:4),
             pch=c(15,16,17,18), bty='n')
    }
    
    # 3.
    # ... or plot by kmer
    if(by.Kmer==TRUE) {
      sapply(1:nrow(dat), function(X) {
        
        myrow <- dat[X,][, names(dat) %in% colnames.toPlot]
        print(myrow)
        points(x=rep(max(myrow), length(myrow)),
               y=myrow,
               pch= c(0,1,2,5),
               col=c(4:(4+nrow(dat)))[X]
        )
        print("Done with plot")
        
      })
      legend("topleft", legend = colnames.toPlot, pch=c(0,1,2,5), bty='n', cex=1)
      legend("bottomleft", legend = dat[,1], fill=c(4:(4+nrow(dat))), bty='n', cex=0.6)
      
    }
    
  })
  
  
}


get.heatmap.myaffMeans <- function(dat, mylabels) {
  
  heatmap.2(dat, labCol=mylabels,
            col =colorRampPalette(brewer.pal(9,"GnBu"))(100),
            key.title= "", key.xlab="Mean relative affinity",
            density.info='none', trace='none',
            Colv = FALSE, dendrogram ="none",
            key.par=list(mgp=c(1,0.1,0),
                         cex.axis=0.5,
                         cex.lab=0.7,
                         mar=c(4,1,2,1)+0.1)
  )
}



makebxp.affBy.nucleotide <- function (dat.list,
                                      ylab="") {
  
  
  mypos <- paste("Base at position ", main.bs.number[-c(1,2,3,8,9,10)])
    
    lapply(seq_along(dat.list), function(x) {
      
      dat <- dat.list[[x]]
      boxplot(dat[,-1],
              ylim=c(0.2,1),
              xlab=mypos[[x]],
              ylab=ylab)
      
    })
  
}



makebxp.affBy.pos <- function (dat.list, ylab="Change in affinity", mylabels, to.compare) {
  
    ylab <- bquote(""~"["* .(to.compare[1]) %->% .(to.compare[2])*"]")
    ylim <- c(0,1)
    lapply(seq_along(mylabels), function(x) {
      dat <- dat.list[[x]]
      dat <- dat[,-1]
      dg <- apply(dat[, to.compare ], 1, function(z) {
        abs(diff(z))
        
      })
      return(dg)
      
    }) -> to.plot
  
  boxplot(to.plot,
          xpd=TRUE,
          xlab=names(mylabels),
          names=mylabels,
          ylim=ylim,
          ylab= ylab
  )

     mtext(text="Position", side=1,line=1, xpd=TRUE, cex=0.3)
}
