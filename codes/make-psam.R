
# Required
# install.packages("ggseqlogo")

require(ggplot2)
require(ggseqlogo)

# 

# FUNCTIONS


getPSAM <- function(dat, mystring) {

  con.aff <- unlist(dat[grepl(mystring, dat$Kmer),]["R2"])
  psam <- matrix(nrow=4, ncol=width(mystring))
  for (i in 0:(nchar(mystring)-1)) {

        all.by.pos <- paste0(substr(mystring, 1 , 0+i ) , 
                             c("A", "C", "G", "T") , 
                             substr(mystring, (2+i) ,nchar(mystring)))

    psa <- dat[match(all.by.pos, dat$Kmer),]$R2
    psam[,(i+1)] <- psa
  }
  
  rownames(psam) <- c("A", "C", "G", "T")
  psam.ori <- psam
  psam <- log(psam)



  myfinal <- apply(psam, 2, function(b) {

    mymean <- mean(b) 
    round((b - mymean ), 2) 
  } )
  
   return(myfinal)

}





# Get PSAM based on highest rel. aff. sequence
mypsam2 <- getPSAM(dat= m.aff[[1]], mystring="CTAAAAATAG")

# Make Logo representation
p1 <- ggseqlogo(mypsam2, method='custom', seq_type='dna', facet="grid") + 
  ylab(bquote("-"*Delta*Delta*G ~"/RT")) +
  xlab('Binding site position') + 
  ylim(-1,1)


# Make Figure
png(paste0("psam.png" ), width=120, height=80, units = 'mm', res=600)
gridExtra::grid.arrange(p1)
dev.off()


# Repeat analysis for palindromic sequence
mypsam3 <- getPSAM(dat= m.aff[[1]], mystring="CTATTAATAG")

p2 <- ggseqlogo(mypsam3, method='custom', seq_type='dna', facet="grid") + 
  ylab(bquote("-"*Delta*Delta*G ~"/RT")) +
  xlab('Binding site position') + 
  ylim(-1,1)


png(paste0("psam3.png" ), width=120, height=80, units = 'mm', res=600)
gridExtra::grid.arrange(p2)
dev.off()
