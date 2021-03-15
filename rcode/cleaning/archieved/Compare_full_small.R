if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc/")
} else {
  setwd("~/Documents/github/pheno_bc")
}

require(plyr)
require(dplyr)
require(tidyr)
# require(chillR)

rm(list=ls()) 
options(stringsAsFactors = FALSE)

# read in the cleaning phenology data:
#data<-read.csv("input/bc_phenology.csv", header=T, na.strings=c("","NA"))
data<-read.csv("input/day.of.bb.Jan22021.csv", header=T, na.strings=c("","NA"))
data$dataset<-"full"
small<-read.csv("input/day.of.bb.small.csv", header=T, na.strings=c("","NA"))
small$dataset<-"short"

head(data)
data$treatment<-paste(data$chill, data$photo,data$force, sep = ".")
head(small)

data<-data[,c("lab2","population","flask","treatment","species","tbb","latbb80","latbb50","latbb1")]
small<-small[,c("lab2","population","flask","treatment","species","tbb","latbb80","latbb50","latbb1")]

toget<-rbind(data,small)

tar.var<-c("population","flask","treatment","species")
resp.var<-c("tbb","latbb80","latbb50","latbb1")

## subset data to look for duplicates (resp.var are included or most of the subset is duplicated)
trt.sub<-toget[,c(tar.var,resp.var)]

## remove duplicated rows in a simple way
no.dup<-toget[!duplicated(trt.sub),]

dups<-toget[duplicated(trt.sub),]

test<-subset(toget, lab2 == "mp.HC.HP.HF.10.symalb.1")
