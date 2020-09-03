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
#data<-ldply(files, read_csv, na = c("","NA"))
data<-ldply(files, read_csv, na = "")

df<-as.data.frame(data)
df[df=="<NA>"] = "NA"   
head(df)
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
df$species[df$species=="bepap"] <- "betpap"
df$species[df$species=="shecan8"] <- "shecan"
df$species[df$species=="Shecan"] <- "shecan"
df$species[df$species=="corso"] <- "corsto"

unique(df$bbch.t)
# 4 32 5 8 157 6 18 Na

temp<-subset(df, bbch.t =="4");temp
#HHLsm36amealn has 4 as its bbc.t for day 5 and 6

df$bbch.t[which(df$treatment == "HC_HP_LF" & df$flask == "36" & df$species =="amealn" & df$day == "5")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_LF" & df$flask == "36" & df$species =="amealn" & df$day == "5")] <- "NA"

df$bbch.t[which(df$treatment == "HC_HP_LF" & df$flask == "36" & df$species =="amealn" & df$day == "6")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_LF" & df$flask == "36" & df$species =="amealn" & df$day == "6")] <- "NA"

temp<-subset(df, bbch.t =="4");temp

###########################################################################################
test<-subset(df, treatment == "HC_HP_LF")

#temp<-subset(df, bbch.t =="6");temp
#HHHsm34 has 6 as its bbc.t for day 3
df$bbch.t[which(df$treatment == "HC_HP_HF" & df$flask == "34" & df$species =="vibedu" & df$day == "3")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "34" & df$species =="vibedu" & df$day == "3")] <- "NA"

df$bbch.t[which(df$treatment == "HC_HP_HF" & df$flask == "34" & df$species =="vibedu" & df$day == "5")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "34" & df$species =="vibedu" & df$day == "5")] <- "NA"

###########################################################################################
temp<-subset(df, bbch.t =="32");temp
#LHL12amealn had a value of 32 on day
df$bbch.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn"  & df$day == "15")] <- "NA"
df$percent.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn"  & df$day == "15")] <- "NA"

df$bbch.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn"  & df$day == "16")] <- "NA"
df$percent.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn"  & df$day == "16")] <- "NA"

df$bbch.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn"  & df$day == "17")] <- "NA"
df$percent.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn"  & df$day == "17")] <- "NA"

df$bbch.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn"  & df$day == "18")] <- "NA"
df$percent.t[which(df$treatment == "LC_HP_LF" & df$flask == "12" & df$species =="amealn"  & df$day == "18")] <- "NA"

###########################################################################################
#temp<-subset(df, bbch.t =="5");temp
#LLL9vibedu has 5 as its bbc.t for day 15
df$bbch.t[which(df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu" & df$day == "15")] <- "NA"
df$percent.t[which(df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu" & df$day == "15")] <- "NA"

###########################################################################################
#temp<-subset(df, bbch.t =="8");temp
#LHH1amealn has 8 as its bbc.t for day 23 and 24
df$bbch.t[which(df$treatment == "LC_HP_HF" & df$flask == "1" & df$species =="amealn"  & df$day == "23")] <- "NA"
df$percent.t[which(df$treatment == "LC_HP_HF" & df$flask == "1" & df$species =="amealn" & df$day == "23")] <- "NA"

df$bbch.t[which(df$treatment == "LC_HP_HF" & df$flask == "1" & df$species =="amealn"  & df$day == "24")] <- "NA"
df$percent.t[which(df$treatment == "LC_HP_HF" & df$flask == "1" & df$species =="amealn" & df$day == "24")] <- "NA"

###########################################################################################
temp<-subset(df, bbch.t =="157");temp #should just be 15 for day 31 and 32
df$bbch.t[which(df$treatment == "LC_HP_HF" & temp$flask == "39" & temp$species =="corsto" & df$day == "31")] <- 15

df$bbch.t[which(df$treatment == "LC_HP_HF" & temp$flask == "39" & temp$species =="corsto" & df$day == "32")] <- 15

df$bbch.t[which(df$treatment == "LC_HP_HF" & temp$flask == "39" & temp$species =="corsto" & df$day == "33")] <- 15

###########################################################################################
#temp<-subset(df, bbch.t =="18");temp
#HHHspibet3 on day 31
df$bbch.t[which(df$treatment == "HC_HP_HF" & df$flask == "3" & df$species =="spibet" & df$day == "31")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "3" & df$species =="spibet" & df$day == "31")] <- "NA"

###########################################################################################

temp<-subset(df, bbch.t =="Na");temp
df$bbch.t[df$bbch.t=="Na"] <- "NA"

temp<-subset(df, bbch.t == "NA");temp
unique(temp$bbch.t)
unique(df$bbch.t)

###########################################################################################
###########################################################################################
###########################################################################################
unique(df$bbch.l)
#should be shifted over, bbch in is the count columns for HHF

#temp<-subset(df, bbch.l =="60");temp
temp<-subset(df, bbch.l =="19");temp
#temp<-subset(df, bbch.l =="25");temp
#temp<-subset(df, bbch.l =="80");temp
#temp<-subset(df, bbch.l =="N");temp
#temp<-subset(df, bbch.l == "100");temp
temp<-subset(df, bbch.l =="5");temp
temp<-subset(df, bbch.l =="2");temp
temp<-subset(df, bbch.l =="\\");temp
temp<-subset(df, bbch.l =="75");temp
temp<-subset(df, bbch.l =="73");temp

# # fixing the incorrect 60 values
# df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "33" & df$species =="corsto" & df$day == "15")] <- NA
# 
# df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "33" & df$species =="corsto" & df$day == "16")] <- NA
# 
# df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "33" & df$species =="corsto" & df$day == "17")] <- NA
# 
# df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "33" & df$species =="corsto" & df$day == "18")] <- NA
# 
# df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "33" & df$species =="corsto" & df$day == "19")] <- NA
# 
# df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "33" & df$species =="corsto" & df$day == "22")] <- NA

# fixing the incorrect 19 values, I think it should be 10, but no way of knowing for certain
df$bbch.l[which(df$treatment == "HC_HP_HF" & df$flask == "11" & df$species =="rhoalb" & df$day == "15")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "11" & df$species =="rhoalb"  & df$day == "15")] <- "NA"

df$bbch.l[which(df$treatment == "HC_HP_HF" & df$flask == "11" & df$species =="rhoalb" & df$day == "16")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "11" & df$species =="rhoalb"  & df$day == "16")] <- "NA"

df$bbch.l[which(df$treatment == "HC_HP_HF" & df$flask == "11" & df$species =="rhoalb" & df$day == "17")] <- "NA"
df$percent.t[which(df$treatment == "LHC_HP_HF" & df$flask == "11" & df$species =="rhoalb"  & df$day == "17")] <- "NA"

df$bbch.l[which(df$treatment == "HC_HP_HF" & df$flask == "11" & df$species =="rhoalb" & df$day == "18")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "11" & df$species =="rhoalb"  & df$day == "18")] <- "NA"

df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "32" & df$species =="shecan" & df$day == "54")] <- "NA"
df$percent.t[which(df$treatment == "LC_LP_LF" & df$flask == "32" & df$species =="shecan"  & df$day == "54")] <- "NA"

df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "32" & df$species =="shecan" & df$day == "55")] <- "NA"
df$percent.t[which(df$treatment == "LC_LP_LF" & df$flask == "32" & df$species =="shecan"  & df$day == "55")] <- "NA"

df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "32" & df$species =="shecan" & df$day == "56")] <- "NA"
df$percent.t[which(df$treatment == "LC_LP_LF" & df$flask == "32" & df$species =="shecan"  & df$day == "56")] <- "NA"

#######################################################################
#Fixing the accidental \\
df$bbch.l[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer" & df$day == "33")] <- "NA"
df$percent.t[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer"  & df$day == "33")] <- "NA"

df$bbch.l[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer" & df$day == "34")] <- "NA"
df$percent.t[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer"  & df$day == "34")] <- "NA"

df$bbch.l[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer" & df$day == "35")] <- "NA"
df$percent.t[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer"  & df$day == "35")] <- "NA"

df$bbch.l[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer" & df$day == "36")] <- "NA"
df$percent.t[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer"  & df$day == "36")] <- "NA"

df$bbch.l[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer" & df$day == "37")] <- "NA"
df$percent.t[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer"  & df$day == "37")] <- "NA"

df$bbch.l[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer" & df$day == "38")] <- "NA"
df$percent.t[which(df$treatment == "HC_LP_LF" & df$flask == "10" & df$species =="menfer"  & df$day == "38")] <- "NA"


# ###########################################################################################
# df$bbch.l[which(df$treatment == "HC_HP_HF" & df$flask == "6" & df$species =="spipyr" & df$day == "10")] <- NA
# df$bbch.l[which(df$treatment == "HC_HP_HF" & df$flask == "6" & df$species =="spipyr" & df$day == "11")] <- NA
# df$bbch.l[which(df$treatment == "HC_HP_HF" & df$flask == "6" & df$species =="spipyr" & df$day == "12")] <- NA

# ###########################################################################################
# # bbch.l is incorrectly 73
# df$bbch.l[which(df$treatment == "LC_HP_LF" & df$flask == "21" & df$species =="rubpar" & df$day == "22")] <- NA
# 
# ###########################################################################################
# # bbch.l is incorrectly N
# df$bbch.l[which(df$treatment == "HC_LP_HF" & df$flask == "15" & df$species =="poptre" & df$day == "23")] <- NA
# df$percent.t[which(df$treatment == "HC_LP_HF" & df$flask == "15" & df$species =="poptre" & df$day == "23")] <- NA

###########################################################################################
# bbch.l is incorrectly 5
df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "27" & df$species =="shecan" & df$day == "38")] <- "NA"
df$percent.t[which(df$treatment == "LC_LP_LF" & df$flask == "27" & df$species =="shecan" & df$day == "38")] <- "NA"

unique(df$bbch.l) # Looks great

###########################################################################################
###########################################################################################
###########################################################################################
unique(df$bbch2.l)
###########################################################################################
# bbch2.l is incorrectly 40
df$bbch2.l[which(df$population == "sm" & df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu" & df$day == "28")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu" & df$day == "28")] <- "NA"

df$bbch2.l[which(df$population == "sm" & df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu" & df$day == "29")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu" & df$day == "29")] <- "NA"

df$bbch2.l[which(df$population == "sm" & df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu" & df$day == "31")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "LC_LP_LF" & df$flask == "9" & df$species =="vibedu" & df$day == "31")] <- "NA"

###########################################################################################
# bbch2.l is incorrectly 1-
df$bbch2.l[which(df$treatment == "LC_LP_HF" & df$flask == "5" & df$species =="popbal" & df$day == "40")] <- "NA"
df$percent2.l[which(df$treatment == "LC_LP_HF" & df$flask == "5" & df$species =="popbal" & df$day == "40")] <- "NA"

df$bbch2.l[which(df$treatment == "LC_LP_HF" & df$flask == "5" & df$species =="popbal" & df$day == "41")] <- "NA"
df$percent2.l[which(df$treatment == "LC_LP_HF" & df$flask == "5" & df$species =="popbal" & df$day == "41")] <- "NA"

###########################################################################################
# bbch2.l is incorrectly 5, supposed to be 3
df$bbch2.l[which(df$treatment == "LC_HP_LF" & df$flask == "6" & df$species =="loninv" & df$day == "24")] <- 3
df$bbch2.l[which(df$treatment == "LC_HP_LF" & df$flask == "6" & df$species =="loninv" & df$day == "25")] <- 3

###########################################################################################
# bbch2.l is incorrectly 6, it should be 7
df$bbch2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "19")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "19")] <- "NA"

df$bbch2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "20")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "20")] <- "NA"

df$bbch2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "21")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "21")] <- "NA"

df$bbch2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "22")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "22")] <- "NA"

df$bbch2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "23")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "HC_LP_LF" & df$flask == "27" & df$species =="riblac" & df$day == "23")] <- "NA"

df$bbch2.l[which(df$population == "sm" & df$treatment == "LC_LP_HF" & df$flask == "19" & df$species =="betpap" & df$day == "27")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "LC_LP_HF" & df$flask == "19" & df$species =="betpap" & df$day == "27")] <- "NA"

df$bbch2.l[which(df$population == "sm" & df$treatment == "LC_LP_HF" & df$flask == "19" & df$species =="betpap" & df$day == "28")] <- "NA"
df$percent2.l[which(df$population == "sm" & df$treatment == "LC_LP_HF" & df$flask == "19" & df$species =="betpap" & df$day == "28")] <- "NA"

#Unlear whether the accidental 19 values should be 9 or 10
df$bbch2.l[df$bbch2.l=="70"] <- "NA"

#Unlear whether the accidental 70 values should be 1 or 0
df$bbch2.l[df$bbch2.l=="19"] <- "NA"

#The 4 should be 3
df$bbch2.l[which(df$population == "mp" & df$treatment == "LC_LP_LF" & df$flask == "28" & df$species =="amealn" & df$day == "44")] <- 3
df$bbch2.l[which(df$population == "mp" & df$treatment == "LC_LP_LF" & df$flask == "28" & df$species =="amealn" & df$day == "47")] <- 3

#The 16 should be 15
df$bbch2.l[which(df$population == "mp" & df$treatment == "HC_LP_LF" & df$flask == "19" & df$species =="popbal" & df$day == "84")] <- 15
df$bbch2.l[which(df$population == "mp" & df$treatment == "HC_LP_LF" & df$flask == "19" & df$species =="popbal" & df$day == "88")] <- 15

#The 8 should be 9
df$bbch2.l[which(df$population == "sm" & df$treatment == "HC_LP_HF" & df$flask == "11" & df$species =="corsto" & df$day == "26")] <- 15
df$bbch2.l[which(df$population == "sm" & df$treatment == "HC_LP_HF" & df$flask == "11" & df$species =="corsto" & df$day == "27")] <- 15
df$bbch2.l[which(df$population == "mp" & df$treatment == "LC_LP_LF" & df$flask == "12" & df$species =="loninv" & df$day == "112")] <- 15
df$bbch2.l[which(df$population == "mp" & df$treatment == "LC_LP_LF" & df$flask == "12" & df$species =="loninv" & df$day == "113")] <- 15

#The 11` should just be 11`
df$bbch2.l[df$bbch2.l=="11`"] <- "11"

#The 2 should be 3
df$bbch2.l[df$bbch2.l=="2"] <- "3"


unique(df$bbch2.l)
##################################################################################
unique(df$bbch3.l)

#The two values should be 1's
df$bbch3.l[df$bbch3.l=="2"] <- "1"

#The 135 values should be 1's
df$bbch3.l[df$bbch3.l=="135"] <- "1"

#The 14 values should be 1's
df$bbch3.l[df$bbch3.l=="14"] <- "1"

#The 4 values should be 3's
df$bbch3.l[df$bbch3.l=="4"] <- "3"

#The 8 values should be 9's
df$bbch3.l[df$bbch3.l=="8"] <- "9"

#The 6 values should be 7's
df$bbch3.l[df$bbch3.l=="6"] <- "7"

#The 70 are unclear so changing to NA
df$bbch3.l[df$bbch3.l=="70"] <- "NA"

#I am fairly certain taht these should be 7, since on day 41 it was 3
df$bbch3.l[df$bbch3.l=="5"] <- "7"

# It is unclear what the 25 values should be 
temp$bbch3.l[temp$bbch3.l=="25"] <- "NA"
df$percent3.l[which(df$population == "sm" & df$treatment == "LC_HP_HF" & df$flask == "34" & df$species =="alninc" & df$day == "46")] <- "NA"
df$percent3.l[which(df$population == "sm" & df$treatment == "LC_HP_HF" & df$flask == "34" & df$species =="alninc" & df$day == "47")] <- "NA"
df$percent3.l[which(df$population == "sm" & df$treatment == "LC_HP_HF" & df$flask == "34" & df$species =="alninc" & df$day == "48")] <- "NA"

unique(df$bbch3.l)

########################################################
########################################################
########################################################

unique(df$bbch4.t)

#The two values should be 1's
df$bbch4.t[df$bbch4.t=="2"] <- "1"

#The I am not sure what the 19 value should be, probably 15?
df$bbch4.t[df$bbch4.t=="19"] <- "NA"

#The I am not sure what the 31 value should be
df$bbch4.t[df$bbch4.t=="31"] <- "NA"

#The 4 values should be 3
df$bbch4.t[df$bbch4.t=="4"] <- "3"

# The 980 should be 9 and 90%
df$bbch4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "20")] <- 9
df$percent4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "20")] <-80

df$bbch4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "21")] <- 9
df$percent4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "21")] <-80

df$bbch4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "22")] <- 9
df$percent4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "22")] <-80

df$bbch4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "23")] <- 9
df$percent4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "23")] <-80

## Double checking
unique(df$bbch.t)
unique(df$bbch4.t)
unique(df$bbch.l)
unique(df$bbch2.l)
unique(df$bbch3.l)

# Fixing typos with the pecent values:
#The -100 values should be 100's
df$percent.l[df$percent.l=="-100"] <- "100"

#The NS values should be NA
df$percent2.l[df$percent2.l=="NS"] <- "NA"

#The -4 values should be 4
df$percent2.l[df$percent2.l=="-4"] <- "4"

#The 751 values should be 75
df$percent2.l[df$percent2.l=="751"] <- "75"

#The i values is unknown
df$percent3.l[df$percent3.l=="i"] <- "NA"

#The i values is unknown
df$percent4.t[df$percent4.t=="8-"] <- "80"


###################################################################################################
# Since symalb and menfer did not have terminal buds, all the bbc.t and percent.t should be NA

symalb<-subset(df, species=="symalb")
head(symalb)
unique(symalb$bbch4.t)

temp<-subset(symalb, bbch.t =="15");temp

df$bbch.t[which( df$species =="symalb")] <- "NA"
df$percent.t[which(df$species =="symalb")] <-"NA"
df$bbch4.t[which( df$species =="symalb")] <- "NA"
df$percent4.t[which(df$species =="symalb")] <-"NA"

df$bbch.t[which( df$species =="menfer")] <- "NA"
df$percent.t[which(df$species =="menfer")] <-"NA"
df$bbch4.t[which( df$species =="menfer")] <- "NA"
df$percent4.t[which(df$species =="menfer")] <-"NA"

###################################################################################################
###################################################################################################

# For the sake of cleaning the data I am going to remove the extra columns I added:
names(df)
df1<-df[,c(1:3,5:15,18,19)]
names(df1)
head(df1)

unique(df1$treatment) # there are two rows of poptre from day 36 and 38, flask 20 that have NA instead of treatment, to be conservative I am going to remove these
df1<-subset(df1, treatment!= "NA")

#write.csv(df1,"bc_phenology.csv")
