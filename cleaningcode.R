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
data<-ldply(files, read_csv, na = "")

df<-as.data.frame(data)

#Interesting, there are a number of additional columns being added, some with no data

df<-df[order(df$day),]

#These are all misplaced comments
unique(df$X21) 
unique(df$X20)
unique(df$X22)
unique(df$comments) 

#subsetting out the sporatic comments and the final percent columns that are not actually filled out
df<-df[, (1:22)]
names(df)

#Double check that there are no random typos
unique(df$day) #looks good

unique(df$population)
#[1] "mp" "sm" "Sm" NA  
df$population[df$population=="Sm"] <- "sm"

unique(df$treatment)
df$treatment[df$treatment=="LHH"] <- "LC_HP_HF"
df$treatment[df$treatment=="LHL"] <- "LC_HP_LF"
df$treatment[df$treatment=="LLH"] <- "LC_LP_HF"
df$treatment[df$treatment=="LLL"] <- "LC_LP_LF"

unique(df$flask)
#look into 2019, NA, 41
temp<-subset(df, flask =="2019");temp
#it's a shecan in HC LP LF, but there is no data for it, so I am going to remove it
df<-subset(df, flask!= "2019")

unique(df$species)
# rupar apipyr-missing poptre-missing spibe popre bepap shecan8 Shecan 
df$species[df$species=="rupar"] <- "rubpar"
df$species[df$species=="spipyr-missing"] <- "spipyr"
df$species[df$species=="poptre-missing"] <- "poptre"
df$species[df$species=="spibe"] <- "spibet"
df$species[df$species=="popre"] <- "poptre"
df$species[df$species=="bepap"] <- "bepap"
df$species[df$species=="shecan8"] <- "shecan"
df$species[df$species=="Shecan"] <- "shecan"

unique(df$bbch.t)
# 4 32 5 8 157 6 18 Na

#HHLsm36amealn has 4 as its bbc.t for day 5 and 6
df$bbch.t[which(df$treatment == "HC_HP_LF" & df$flask == "36" & df$species =="amealn")] <- NA
df$percent.t[which(df$treatment == "HC_HP_LF" & df$flask == "36" & df$species =="amealn")] <- NA

#temp<-subset(df, bbch.t =="6");temp
#HHHsm34 has 6 as its bbc.t for day 3
df$bbch.t[which(df$treatment == "HC_HP_HF" & df$flask == "34" & df$species =="vibedu")] <- NA
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "34" & df$species =="vibedu")] <- NA

#temp<-subset(df, bbch.t =="32");temp
#LHL12amealn had a value of 32 on day
df$bbch.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn")] <- NA
df$percent.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn")] <- NA


temp<-subset(df, bbch.t =="5");temp
#LLL9vibedu has 5 as its bbc.t for day 15
df$bbch.t[which(df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu")] <- NA
df$percent.t[which(df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu")] <- NA

#temp<-subset(df, bbch.t =="8");temp
#LHH1amealn has 8 as its bbc.t for day 23 and 24
df$bbch.t[which(df$treatment == "LC_HP_HF" & df$flask == "1" & df$species =="amealn")] <- NA
df$percent.t[which(df$treatment == "LC_HP_HF" & df$flask == "1" & df$species =="amealn")] <- NA

#temp<-subset(df, bbch.t =="157");temp #should just be 15 for day 31 and 32
df$bbch.t[which(df$treatment == "LC_HP_HF" & temp$flask == "39" & temp$species =="corsto")] <- 15

temp<-subset(df, bbch.t =="18");temp
#HHHspibet3 on day 31
df$bbch.t[which(df$treatment == "HC_HP_HF" & df$flask == "3" & df$species =="spibet")] <- NA
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "3" & df$species =="spibet")] <- NA

#temp<-subset(df, bbch.t =="NA");temp
df$bbch.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn")] <- NA
df$percent.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn")] <- NA


unique(df$bbch.l)
# 50 60 8 20 45 25 19 73 80 100 2 5 75 70

unique(df$bbch2.l)
# 60 50 40 15 45 2 80 30 20 65 75 85 25 6 5 19 8 70 11' 1- 4 16

unique(df$bbch3.l)
#100 50 20 80 30 45 70 25 35 2 175 40 135 60 14 4 8 75 6 

unique(df$bbch4.l)
# 50 70

unique(df$bbch4.t)
# 2 50 70 40 15 980 19 31 4 

#Look into these more closely
temp<-subset(df, day=="NA");temp

temp<-subset(df, comments=="70");temp
