
# Scrip filters aff table from selex data based on motif and other parameters:
# 1.) aff  = affnity cutoff to be used
# 2.) mismatch = allowed ANYWHERE in the sequence
# Two main motifs in thi case; can then further filter tables to study variations at central vs edges
# Saves tables and output file with command input 


setwd(subDir$inpath)
dir.create("mlrseq")
mlr.inpath <- file.path(mainDir, "mlr", "input", analysis.name)
setwd("./mlrseq")

aff.data
d.aff <- aff.data[["r2wt"]]
d.aff <- d.aff[, c("Kmer", "Affinity", "ObservedCount")]
head(d.aff)



# MAIN COMMAND TO GET SEQs FIR MRL -------------


# Use Filter Table fun() for filering based on :
# Use Filter Table fun() for filering based on :
aff.cutoff <- 0.0



filterByAgrep(dat=d.aff[,c(1,2,3)],
              regex=list(YTAWWWWTAR = "^[CT]TA[AT][AT][AT][AT]TA[AG]$",
                         YTANNNNTAR = "^[CT]TA[ATCG][ATCG][ATCG][ATCG]TA[AG]$",
                         CTANNNNTAR = "^CTA[ATCG][ATCG][ATCG][ATCG]TAG$"),
              sample.name="wt",
              aff=aff.cutoff,get.tbR=FALSE,
              mymismatch=1)

filterByAgrep(dat=d.aff[,c(1,2,3)],
              regex=list(YTAWWWWTAR = "^[CT]TA[AT][AT][AT][AT]TA[AG]$",
                         CTAWWWWTAG = "^CTA[AT][AT][AT][AT]TAG$",
                         YTANNNNTAR = "^[CT]TA[ATCG][ATCG][ATCG][ATCG]TA[AG]$",
                         CTANNNNTAR = "^CTA[ATCG][ATCG][ATCG][ATCG]TAG$",
                         CTAAAAATAG = "CTAAAAATAG"),
              sample.name="wt",
              aff=aff.cutoff,get.tbR=FALSE,
              mymismatch=2)

filterByAgrep(dat=d.aff[,c(1,2,3)],
              regex=list(YTAWWWWTAR = "^[CT]TA[AT][AT][AT][AT]TA[AG]$",
                         CTAWWWWTAG = "^CTA[AT][AT][AT][AT]TAG$",
                         CTANNNNTAR = "^CTA[ATCG][ATCG][ATCG][ATCG]TAG$",
                         CTAAAAATAG = "CTAAAAATAG"),
              sample.name="wt",
              aff=aff.cutoff,get.tbR=FALSE,
              mymismatch=3)

filterByAgrep(dat=d.aff[,c(1,2,3)],
              regex=list(YTAWWWWTAR = "^[CT]TA[AT][AT][AT][AT]TA[AG]$",
                         CTAWWWWTAG = "^CTA[AT][AT][AT][AT]TAG$",
                         CTAAAAATAG = "CTAAAAATAG"),
              sample.name="wt",
              aff=aff.cutoff,get.tbR=FALSE,
              mymismatch=4)






sapply(list.files(path="./", pattern="mismatch[234]"), function (x) {
  
  print(x)
  df <- read.table(x)
  print(nrow(df))
  
  if (nrow(df) > 1000) {
    
    
    print("Split based on variations at the core or edges")
    
    mydf <- df[grep("^...(AAAA)|(TTTT)...$", df[,1]), ]
    myname <- sub(pattern="mismatch", replacement = paste0("fixedA4", "m"), x=x, perl = TRUE)
    write.table(mydf, file=myname, quote=F, col.names=F, row.names=F)
    
    mydf <- df[grep("^...[AT][AT][AT][AT]...$", df[,1]), ]
    myname <- sub(pattern="mismatch", replacement = paste0("fixedW4", "m"), x=x, perl = TRUE)
    write.table(mydf, file=myname, quote=F, col.names=F, row.names=F)
    
    
  }
  
}
)

# Remove duplicates and datasets with low #seqs
system("rm ./mef2b_*fixed*NNNN*")
system("rm -f ./mef2b_*fixedA*CTAAAAATAG* ")
system("find -name '*.txt' | xargs  wc -l | awk '{if($1 < 500  && index($2, 'txt')>0 ) print $2}' | xargs rm")



