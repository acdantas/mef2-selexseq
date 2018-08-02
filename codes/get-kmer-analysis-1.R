# MAKE FIGURES FOR OVERVIEW SELEX-SEQ RESULTS -------- 
# Main and Supp Figures -------------

# A) Information Gain 
# B) R2 vs R1 plot
# C) PWM (run MEME)
# D) Shape profile



# Functions for plotting results ---------------

getInfoGainFig <- function(x) {
  
  barplot(height=infoGain.data[[x]]$InformationGain, 
          names.arg =  infoGain.data[[x]]$K,
          ylab="Information Gain", cex.names=1, las=2,
          xlab=expression( italic(k) * "mer length"))
  
}


makeR1vsR2Fig <- function (x){
  
  par(pty='s')
  plot(R1 ~ R2, data=m.aff[[x]], pch=20, cex=0.5,
       xlab="Relative Affinity (R2 vs. R0)",
       ylab="Relative Affinity (R1 vs. R0)")
  
}


################################## MAIN CODE ####################################


# Open output file with analysis details


# parameters
mysample <- "wt"
myaff <- 0.7
mm <- 0 # mismatch
  

# 1. SELEX-seq results -------

# PWM and Shape Means ---------

# Data processing based on parameters

topseqs <- unlist(unlist(filterAffTb(round="r2", 
                                     name=mysample, 
                                     dat=aff.data, 
                                     aff=myaff, 
                                     mymismatch=mm,
                                     get.tbR=TRUE,
                                     regex=list(YTANNNNTAR =".")),
                         recursive=F), 
                  recursive=F)

topseqs <- as.data.frame(topseqs)
head(topseqs, 100)
nrow(topseqs)


# Remove shifted motifs
 topseqs <- subset(topseqs,  
                   !grepl("(.{1,5}CTA[AT][TA])|(.*[TA][TA]TAG.{1,5})", 
                          Kmer,  perl=TRUE))



# write fasta
mytop.file <- file.path(subDir$ana.path, "results", "topkmers.fa")
write.fasta(sequences = as.list(topseqs$Kmer), 
            file.out = mytop.file,
            names=c(paste0(topseqs$Kmer,"_", topseqs$Affinity )))


# OBS: Use topseqs$kmer for making PWM with MEME

#  Make Shape profile
# use same fasta file used for pwm generation
setwd(file.path(subDir$ana.path, "results"))
topseq.pred <- getShape(filename="./topkmers.fa")

#get shape means
ave.pred <- lapply(topseq.pred, colMeans)
ave.pred


# Paper figure : PLOT shape heatmap line
# with levelplot

library(lattice)
lapply(seq_along(ave.pred), function(x) {
  lattice.options(
    layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
    layout.widths=list(left.padding=list(x=5), right.padding=list(x=0))
  )
    myplot <- levelplot(as.matrix(ave.pred[[x]]),
                      col.regions=heat.colors(25), 
                      ylab= list(short.shapelabs[[names(ave.pred)[x]]], cex=1.2, vjust=-4),
                      xlab=list(x.axisShapeLabels[[names(ave.pred)[x]]]), 
                      margins=TRUE,
                     # xlab=list( label= shape.labels[[names(ave.pred)[x]]], cex=1.2), #, hjust=8.6) , 
                      colorkey=list(labels=list(cex=0.7), space="left", hjust=2), #labels=list(cex=0.8),
                      pretty=T, labels=F,  scales=list(draw=F)
  )

  
  return(myplot)
}) -> plot.list



#################################################################################################
# Make Figures

# Figure 1 --------------

# panel B : make with meme suite

# panel C: 
tiff("Fig2-panelC-shapeMeanTopBS.tiff")
print(plot.list[[3]], split = c(1, 1, 1, 4), more = TRUE)
print(plot.list[[1]], split = c(1, 2, 1, 4), more = TRUE)
print(plot.list[[2]], split = c(1, 3, 1, 4), more = TRUE)
print(plot.list[[4]], split = c(1, 4, 1, 4), more = FALSE)
dev.off()


# Figure S3
png("FigS3-BCpanels.png", width=180, height=90, res=600, units='mm')
par(mfrow=c(1,2), cex.axis=1.2, cex.lab=1.2, mar=c(5,5,1,2))
getInfoGainFig(mysample)
makeR1vsR2Fig(mysample)
dev.off()