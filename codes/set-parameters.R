
#set main sample names:
analysis.name <- "mef2b"
main.names <- c("wt")
rounds <- c("R1", "R2")
kmer= "10mer"
#..

# set dir names
setwd("C:/Users/dantasma/Dropbox/mef2-data-analysis/selex-seq/mef2-NAR2018")
mainDir <-  getwd()  # add main directory name
subDir <- list(inpath = file.path(mainDir,"data/"),  # results from SELEX package
               selexOutpath = file.path(mainDir, paste0("data/", "selex-out/")),
               ana.path = file.path(mainDir, "analysis/"), # where  all analysis are
               code.path = file.path(mainDir, paste0("codes/")),
               subs = c("figures", "results")
               )
#..

# create dir 
dir.create(file.path(subDir$ana.path))                       
setwd(subDir$ana.path)                                       
sapply(subDir$subs, function(x)  dir.create(x) )     
#..


# read in input files

# Read input aff data
aff.files <- grep(kmer, list.files(path = subDir$selexOutpath,  pattern = "aff_.*\\.txt$") , value=TRUE )
aff.data <- lapply(aff.files, function(x) read.table(file.path(subDir$selexOutpath, x )))
names(aff.data) <- gsub("aff_.{2}(.*)his_.*txt$", "\\1" , aff.files)
str(aff.data)
#..

# Read input infoGain data
infoGain.files <- list.files(path = subDir$selexOutpath,  
                             pattern = "infoGain.*\\.txt$")
infoGain.files
infoGain.data <- lapply(infoGain.files, function(x) 
  read.table(file.path(subDir$selexOutpath, x )))
names(infoGain.data) <- gsub("infoGain_(.*)his.*txt$", "\\1" , infoGain.files)



# source funtions used to analyze data
source(file.path(subDir$code.path, "funtions-to-analyze-mef2-data.R"))


# All affinity data by Sample
m.aff <- lapply(main.names, mergeAllBySample, arg1=aff.data, arg2=5) # arg2 = ("Count"=2; "Affinity"=5)
m.aff <- setNames(m.aff, main.names)
head(m.aff)
#..


# Set labels
x.axisShapeLabels <- list( MGW = c("-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5") ,
                           ProT = c("-5", "-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4", "+5"),
                           HelT = c("-5/-4", "-4/-3", "-3/-2", "-2/-1", "-1/+1", "+1/+2", "+2/+3", "+3/+4", "+4/+5"),
                           Roll = c("-5/-4", "-4/-3", "-3/-2", "-2/-1", "-1/+1", "+1/+2", "+2/+3", "+3/+4", "+4/+5")
)

short.shapelabs <- list( MGW = expression(paste("MGW [" , ring(A) , "]")) ,
                         ProT = expression(paste("ProT [ ",degree," ]")),
                         HelT = expression(paste("HelT [ ",degree," ]")),
                         Roll = expression(paste("Roll [ ",degree," ]")))


rsq.lab <- expression(italic('R')^"2")




shape.labels <- list( MGW = expression(paste("Minor groove width [" , ring(A) , "]")) ,
                      ProT = expression(paste("Propeller twist [" ,degree, "]")),
                      HelT = expression(paste("Helical twist [" ,degree, "]")),
                      Roll = expression(paste("Roll [" ,degree, "]")))





x.axis.shortShapeLabels <- list( MGW = c( "-3", "-2", "-1", "+1", "+2", "+3") ,
                           ProT = c( "-3", "-2", "-1", "+1", "+2", "+3"),
                           HelT = c("-4/-3", "-3/-2", "-2/-1", "-1/+1", "+1/+2", "+2/+3", "+3/+4"),
                           Roll = c("-4/-3", "-3/-2", "-2/-1", "-1/+1", "+1/+2", "+2/+3", "+3/+4")
)

xaxis.Motif <- list( MGW = c( "-3", "-2", "-1", "+1", "+2", "+3") )

myFigLab <- c("A", "B", "C", "D", "E")
