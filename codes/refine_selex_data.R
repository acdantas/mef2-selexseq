# Refine selex data to be used
# Always use this for data filtering with desired arguments

###############################################################################################
# Set sample preferences
mysample <- "wt"
myround <- "r2"
regex.l <- regex.l1
aff.lim <- 0
seq.mismatch <- 2

###############################################################################################
#                                       MAIN                                                  #
###############################################################################################

# REFINE DFs TO BE USED --------------

# Filter data based on regex list and other parameters.  

df <- unlist(filterAffTb(round=myround, 
                         name=mysample, 
                         dat=aff.data,
                         aff=aff.lim, 
                         mymismatch=seq.mismatch ,get.tbR=TRUE,
                         regex=list(YTANNNNTAR = "^[CT]TA[ATCG][ATCG][ATCG][ATCG]TA[AG]$")),
                recursive=F)


dfByGrp <- groupKmerRegex(dat=as.data.frame(df),
                          regex=regex.l, 
                          title=mysample,
                          r=T)

dfByGrp <- analyzeBS.ByKmer(dat=dfByGrp)


# Remove partial matches in the middle of BS.

dfByGrp <- subset(dfByGrp,
                  !grepl("($.*CTA[AT][TA].*)|(.*[TA][TA]TAG.*^)",
                         Kmer,  
                         perl=TRUE))



df1 <- getKmerRegexgroup(dat=df, 
                         regex=regex.l2, 
                         title="YTAWWWWTAR",
                         r=T) 
df1 <- as.data.frame(unlist(df1 , recursive = F))


df2 <- getKmerRegexgroup(dat=df,
                         regex=regex.l3,
                         title="YTA(noW4tract)TAR",
                         r=T)
df2 <- as.data.frame(unlist(df2 , recursive = F))

df3 <- getKmerRegexgroup(dat=df,
                         regex= list('YTAN4TAR' = '^[CT]TA.{4}TA[GA]$'),
                         title="YTA(N4)TAR",
                         r=T)
df3 <- as.data.frame(unlist(df3 , recursive = F))

df <- as.data.frame(unlist(df , recursive = F))

# Creates new cols df$leftKmer df$kmer df$rightKmer 

df1 <- analyzeBS.ByKmer(dat=df1)
df2 <- analyzeBS.ByKmer(dat=df2)
df3 <- analyzeBS.ByKmer(dat=df3)
df4 <- df3[(df3$leftKmer == "CTA" & df3$rightKmer == "TAG") ,]



df1
# shape.dat <- getShapePred(df[order(-df["Affinity"]),]) # all shape prediction
# 
# binned.shape <- getShapeHeatMap(dat = shape.dat, name = "shapebyAff" , 
#                 n.bins=0, return.binned.data=T, per.bin=200)

file.remove(list.files(path=".", pattern="seq.*fa*"))


