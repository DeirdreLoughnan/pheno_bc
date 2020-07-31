# Cleaning phenology data:
#Begining by combining all the days of observations and fixing any misspelt species names, random ' or typos that are very very obvisouly wrong (eg any scale value that is not on the scale like 22)

if(length(grep("deirdreloughnan", getwd())>0)) { 
  setwd("~/Documents/github/pheno_bc") 
} else
  setwd("~/Documents/github/pheno_bc")


rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(plyr)
require(readr)

files <- list.files(path="~/Documents/github/pheno_bc/pheno_data/typos_cleaned", pattern="*.csv", full.names=TRUE)
files
data<-ldply(files, read_csv)

df<-as.data.frame(data)

#Interesting, there are a number of additional columns being added, some with no data

df<-df[order(df$day),]

#These are all misplaced comments
unique(df$X21) 
unique(df$X20)
unique(df$X22)
unique(df$comments) 

names(df)

