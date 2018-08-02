#  Analyze sequence and shape features of MEF2 sites
#  Kmer stripchart; kmer heatmaps; sequence by position
#  shape heatmap; high vs. low affinity sites

###############################################################################################

# Required libraries and codes:

setwd(subDir$ana.path)

source("../codes/regexList.R")
source("../codes/kmer-functions.R")
source("../codes/refine_selex_data.R")
source( "../codes/make-psam.R")

# refine_selex_data refines data and generates:
# 1) df (all seq based on specified motif) &
# 2) dfByGrp (sequences with grouped motifs)
# df1="W4tract"
# df2="noW4tract"


# Set sample preferences
mysample <- "wt"
myround <- "r2"
regex.l <- regex.l1
aff.lim <- 0.2
seq.mismatch <- 2
setwd(file.path(subDir$ana.path, "results"))



############################# MAIN ANALYSIS ############################


# KMER ANALYSIS (Figure 2)  --------------------------------


# Kmer stripchart analysis (Figure 2A)

plotstripchart(dat=dfByGrp)


# Filter specific kmer groups
t1 <-  df1[df1$group == "CTAW4TAG", ]
t2 <-  df1[df1$group == "CTAW4TAA", ]
t3 <-  df1[df1$group == "TTAW4TAA", ]


# Merge affinity  between kmers.
merged <- merge(t1[,c(2,6,7)], t2[,c(2,6,7)], by=c("kmer", "rightKmer")) 
merged <- merge(merged, t3[,c(2,6)], by="kmer", all=TRUE)
merged <- merged[, -which(names(merged) %in% "rightKmer")]
colnames(merged) <- c("kmer", "CTAW4TAG", "CTAW4TAA","TTAW4TAA") 


# Plot affinity  between kmers YTAW4TAR (Figure 2B)
plotMotifAff(dat=merged, myxlab="CTAW4TAG")


# Edges of the core variations based on AT-rich regions.

# Left trimers (5') data
l.trimer <- df1[df1$rightKmer == "TAG",] 


# Make heatmap (Figure 2C)
get.Kmerby3mer(dat=l.trimer ,
               grp="(3mer)W4TAR",
               col2="kmer", col1="leftKmer")


# At the 3' edge. Not shown on paper
# # Right trimers (3') data
# r.trimer <- df1[df1$leftKmer == "CTA",]
# summary(r.trimer$group)
# get.Kmerby3mer(dat=r.trimer ,
#                grp="(3mer)W4TAR",
#                col2="kmer", col1="rightKmer")



# KMER ANALYSIS (Figure 3 & S4)  -------------------------


# Select data matching CTAWWWWTAG with up to 1 mismatch at center
mydf <- dfByGrp[dfByGrp$group == "01" | dfByGrp$group== "05",]
mydf <- mydf[mydf$rightKmer == "TAG" & mydf$leftKmer== "CTA",]
mydf <-  subset(mydf, !grepl("[GC].*[GC]", kmer, perl=TRUE))



# Refine data to be used
mydf2 <- mydf[,c(1,2)]
n <- 10


# Get a list of dataframes containing aff for each "common" kmer 
# with varying nucleotide at position "N"
affByMutList <- getAffForMutKmer(dat=mydf2, bs.region=c(4:7))


# Get affinity means for variations in nucleotide at each position
my.affMeans <- sapply(affByMutList,
                      function(x) { colMeans(x[,-1]) })


# Labels for heatmap for affMeans 
c4.lab <- main.bs.number[-c(1,2,3,8,9,10)]




# make boxplot aff based on nucleotide mutation
makebxp.affBy.nucleotide(dat=affByMutList,
                         ylab="Relative affinity")


makebxp.affBy.pos(dat.list = affByMutList, 
                  to.compare=c("A", "G"),
                  mylabels=main.bs.number[-c(1,2,3,8,9,10)])

makebxp.affBy.pos(dat.list = affByMutList, 
                  to.compare=c("A", "C"),
                  mylabels=main.bs.number[-c(1,2,3,8,9,10)])




# SHAPE ANALYSIS (Figure 4) ------------------------------


# Run DNAShapeR and proscess shape data.

# Get DNA shape for all 
all.shape <- getShapePred(df[order(-df["Affinity"]),])

# Bin Shape Parameters
binned.shape <- getShapeHeatMap(dat = all.shape,
                                name = "shapebyAff",
                                n.bins=0,
                                return.binned.data=T,
                                per.bin=200)


# Every bin with same number of seqs 
binned.shape2 <- lapply(binned.shape[1:4], function(x) {
  head(x, -1) # Remove last row (not 200 per bin)
})

# Clear directory from shape files
file.remove(list.files(path=".", pattern="seq.*fa*"))



# Set central region for shape means
t <- ncol(shape.dat[["ProT"]])
central4.cols <- c(((t/2)-1):((t/2)+2))

# Get shape means plot comparison
getShapeBxpComparison(shape=all.shape[c(1:4)], region=central4.cols)


# Analyze shape between W4 vs non-W4

# Get shape for individual groups
allshape.df1 <- ana.shape.for.df(dat=df1, name=names.dfs[1])
allshape.df2 <- ana.shape.for.df(dat=df2, name=names.dfs[2])


ylim.shape <- list ( "MGW" = c(4,5.5),
                     "HelT" = c(32,36),
                     "ProT" = c(-16, -6),
                     "Roll" = c(-6, 6)
)


lapply(names(allshape.df1[1:4]), plot.high.vs.low.aff.sites, with.legend=FALSE)

setwd(subDir$ana.path)
source( "../codes/savefigures2-4.R")

setwd(subDir$ana.path)
source("../codes/set-SeqforMLR.R")
