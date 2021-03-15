# Cleaning phenology data:
#Begining by combining all the days of observations and fixing any misspelt species names, random ' or typos that are very very obvisouly wrong (eg any scale value that is not on the scale like 22)

if(length(grep("deirdreloughnan", getwd())>0)) { 
  setwd("~/Documents/github/pheno_bc") 
} else {
  setwd("~/Documents/github/pheno_bc")
}

rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(tidyr)
require(plyr)
require(readr)

files <- list.files(path="~/Documents/github/pheno_bc/pheno_data/typos_cleaned", pattern="*.csv", full.names=TRUE)
files

#data <- ldply(files, read_csv, na = c("","NA"))
data <- ldply(files, read_csv, na = "")


df <- as.data.frame(data)
df[df=="<NA>"] = "NA"   
head(df)
#Interesting, there are a number of additional columns being added, some with no data

df <- df[order(df$day),]
unique(df$day)
#These are all misplaced comments

head(df)
unique(df$X22)
unique(df$comments) 
names(df)
#subsetting out the sporatic comments and the final percent columns that are not actually filled out
df <- df[ , c("day", "population","treatment","indiv","flask","species","bbch.t","percent.t","bbch.l"
              ,"percent.l","bbch2.l","percent2.l","bbch3.l","percent3.l","bbch4.t","percent4.t", "comments")]


#Double check that there are no random typos
unique(df$day) #looks good, except some NA

unique(df$population)
#[1] "mp" "sm" "Sm" NA  
df$population[df$population=="Sm"] <- "sm"
df <- subset(df, population !="NA")

unique(df$treatment)
df$treatment[df$treatment=="LHH"] <- "LC_HP_HF"
df$treatment[df$treatment=="LHL"] <- "LC_HP_LF"
df$treatment[df$treatment=="LLH"] <- "LC_LP_HF"
df$treatment[df$treatment=="LLL"] <- "LC_LP_LF"

unique(df$flask)
#look into 2019, NA, 41
#temp <- subset(df, flask =="2019")
# temp <- subset(df, is.na(day))
#temp <- subset(df, flask == 41)
#temp <- subset(df, indiv == "unknown")
#temp <- subset(df, species == "spipyr-missing")
temp <- subset(df, species == "poptre-missing")

# There are a few oddities that I will start out my just removing, like a sample that has 2019 instead of the flask, samples in flask 41 that fell out of other flasks and could not be re-homed etc
df <- subset(df, flask!= "2019")
df <- subset(df, !is.na(day))
df <- subset(df, flask != 41)
df <- subset(df, species != "spipyr-missing")
df <- subset(df, species != "poptre-missing")


unique(df$species)
temp <- subset(df, species == "spipyr-missing")
# rupar apipyr-missing poptre-missing spibe popre bepap shecan8 Shecan 
df$species[df$species=="rupar"] <- "rubpar"
df$species[df$species=="spibe"] <- "spibet"
df$species[df$species=="popre"] <- "poptre"
df$species[df$species=="bepap"] <- "betpap"
df$species[df$species=="shecan8"] <- "shecan"
df$species[df$species=="Shecan"] <- "shecan"
df$species[df$species=="corso"] <- "corsto"

unique(df$bbch.t)
# 4 32 5 8 157 6 18 Na

temp <- subset(df, bbch.t =="4");temp
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


df$bbch.t[which(df$treatment == "LC_HP_HF" & df$flask == "39" & df$species =="corsto" & df$day == "31")] <- 15

df$bbch.t[which(df$treatment == "LC_HP_HF" & df$flask == "39" & df$species =="corsto" & df$day == "32")] <- 15

df$bbch.t[which(df$treatment == "LC_HP_HF" & df$flask == "39" & df$species =="corsto" & df$day == "33")] <- 15

###########################################################################################
#temp<-subset(df, bbch.t =="18");temp
#HHHspibet3 on day 31
df$bbch.t[which(df$treatment == "HC_HP_HF" & df$flask == "3" & df$species =="spibet" & df$day == "31")] <- "NA"
df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "3" & df$species =="spibet" & df$day == "31")] <- "NA"

###########################################################################################

#temp<-subset(df, bbch.t =="Na");temp
df$bbch.t[df$bbch.t=="Na"] <- "NA"

temp <- subset(df, bbch.t == "NA");temp
unique(temp$bbch.t)
unique(df$bbch.t)

###########################################################################################
###########################################################################################
###########################################################################################
unique(df$bbch.l)
#should be shifted over, bbch in is the count columns for HHF

temp <- subset(df, bbch.l =="19");temp
temp <- subset(df, bbch.l =="5");temp
temp <- subset(df, bbch.l =="2");temp
temp <- subset(df, bbch.l =="\\");temp
# temp <- subset(df, bbch.l =="75");temp
# temp <- subset(df, bbch.l =="73");temp

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

###########################################################################################
# bbch.l is incorrectly 5
df$bbch.l[which(df$treatment == "LC_LP_LF" & df$flask == "27" & df$species =="shecan" & df$day == "38")] <- "NA"
df$percent.t[which(df$treatment == "LC_LP_LF" & df$flask == "27" & df$species =="shecan" & df$day == "38")] <- "NA"

sort(unique(df$bbch.l)) # Looks great

###########################################################################################
###########################################################################################
###########################################################################################
sort(unique(df$bbch2.l))
#"16", "2",  "4"   "40"  "5"   "6", "70"  "8","1-"
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

#Somehow on day 15 a samrac pop was changed to sm
df$population[which(df$species =="samrac" & df$treatment == "LC_HP_HF")] <- "mp"

unique(df$bbch2.l)
##################################################################################
unique(df$bbch3.l) # "70"  "135" "2"   "14"  "4"   "8"   "5"   "6"  

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
# temp$bbch3.l[temp$bbch3.l=="25"] <- "NA"
# df$percent3.l[which(df$population == "sm" & df$treatment == "LC_HP_HF" & df$flask == "34" & df$species =="alninc" & df$day == "46")] <- "NA"
# df$percent3.l[which(df$population == "sm" & df$treatment == "LC_HP_HF" & df$flask == "34" & df$species =="alninc" & df$day == "47")] <- "NA"
# df$percent3.l[which(df$population == "sm" & df$treatment == "LC_HP_HF" & df$flask == "34" & df$species =="alninc" & df$day == "48")] <- "NA"

unique(df$bbch3.l)

########################################################
########################################################
########################################################
# head(df)
# unique(df$bbch4.t) # I am ignoring this, since I am not going to include it in the analysis
# 
# #The two values should be 1's
# df$bbch4.t[df$bbch4.t=="2"] <- "1"
# 
# #The I am not sure what the 19 value should be, probably 15?
# df$bbch4.t[df$bbch4.t=="19"] <- "NA"
# 
# #The I am not sure what the 31 value should be
# df$bbch4.t[df$bbch4.t=="31"] <- "NA"
# 
# #The 4 values should be 3
# df$bbch4.t[df$bbch4.t=="4"] <- "3"
# 
# # The 980 should be 9 and 90%
# df$bbch4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "20")] <- 9
# df$percent4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "20")] <-80
# 
# df$bbch4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "21")] <- 9
# df$percent4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "21")] <-80
# 
# df$bbch4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "22")] <- 9
# df$percent4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "22")] <-80
# 
# df$bbch4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "23")] <- 9
# df$percent4.t[which(df$population == "mp" & df$treatment == "HC_LP_HF" & df$flask == "29" & df$species =="amealn" & df$day == "23")] <-80

## Double checking
unique(df$bbch.t)



# df$bbch.t[which(df$treatment == "HC_HP_HF" & df$flask == "3" & df$species =="spibet" & df$day == "31")] <- "NA"
# df$percent.t[which(df$treatment == "HC_HP_HF" & df$flask == "3" & df$species =="spibet" & df$day == "31")] <- "NA"

unique(df$bbch4.t)
unique(df$bbch.l)
unique(df$bbch2.l)
unique(df$bbch3.l)

# Fixing typos with the pecent values:
unique(df$percent.t)
unique(df$percent.l)
unique(df$percent2.l)
unique(df$percent3.l)
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

symalb <- subset(df, species=="symalb")
head(symalb)
unique(symalb$bbch4.t)

temp <- subset(symalb, bbch.t =="15");temp

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
df1 <- df[ , c("day", "population","indiv","treatment","flask","species","bbch.t","percent.t","bbch.l"
            ,"percent.l","bbch2.l","percent2.l","bbch3.l","percent3.l","bbch4.t","percent4.t", "comments")]

names(df1)
head(df1)

df1$lab <- paste(df1$population,df1$treatment,df1$flask, df1$species, sep=".")

unique(df1$treatment)
unique(df1$day) 
unique(df1$population)
unique(df1$flask)# there are two rows of poptre from day 36 and 38, flask 20 that have NA instead of treatment, to be conservative I am going to 



#investigate closer
temp  <-  subset(df1, lab == "sm_LC_LP_HF_2_symalb")
temp  <-  subset(df1, lab == "sm_LC_LP_HF_30_symalb")
temp  <-  subset(df1, lab == "sm_LC_HP_LF_13_corsto")

length(unique(df1$lab))
# Removing the individuals that died
# data  <-  df
# temp  <-  df1[, c("day","indiv","population","species","treatment","lab")]
# temp0  <-  subset(temp, day == "0")
# table(temp0$species)
# sum(table(temp0$species))

##############################################################################################
#begin by dividing the treatment names (C_P_F) into their own columns
require(tidyr)
data <- df1 %>%
  separate(treatment, c("chill","photo","force"), "_")
 head(data)
 
 # here I am overwriting the lab, but with only . instead of _
data$lab <- paste(data$population,data$chill,data$photo,data$force,data$flask, data$species, sep=".")
d <- data
head(d)
#start by identifying the samples that are duplicates, demarcated with T or F

d$dup <- duplicated(d[,c("day","lab")]) # there are 212 samples that are duplicated

#Check that it worked the way I wanted
test <- subset(d, dup == "TRUE") # 14752 duplicates 
# there are two flasks that have 3 of the same species in it!
# test <- subset(d, lab =="sm_HC_LP_HF_37_vacmem")
# test <- subset(d, lab =="mp_LC_LP_LF_4_menfer")
head(d)

d <- d %>% 
  group_by(day, lab) %>% 
  mutate(ref=ifelse(dup, "2", "1"))

d$lab2 <- paste(d$lab, d$ref, sep=".")
d <- as.data.frame(d)
head(d)

d$treatment <- paste(d$chill, d$photo, d$force, sep = ".")
#There are an additional four samples that had three per flask:
unique(d$lab2)
d$dup2 <- duplicated(d[,c("day","lab2")])

d$lab2[which(d$lab == "mp.LC.LP.LF.4.menfer" & d$dup2 == "TRUE")] <- "mp.LC.LP.LF.4.menfer.3"
#d$lab2[which(d$lab == "sm.LC.LP.LF.30.shecan" & d$dup2 == "TRUE")] <- "sm.LC.LP.LF.30.shecan.3"
d$lab2[which(d$lab == "sm.HC.LP.HF.37.vacmem" & d$dup2 == "TRUE")] <- "sm.HC.LP.HF.37.vacmem.3"

temp1 <- subset(d, lab2== "mp.LC.LP.LF.4.menfer.3")
length(unique(d$lab2))
head(d)
# d <- d %>% 
#   group_by(day, lab2) %>% 
#   mutate(ref2=ifelse(dup2, "3", ""))
# head(d)
# d$lab3 <- paste(d$lab2, d$ref2, sep=".")
# d <- as.data.frame(d)
# head(d)

# For curiosity, here I am creating a new datset of jus the samples that have pairs of the same species in a flask
# ddups <- vector()
# for(i in 1:length(d$dup)){
#   if(d$dup[i] == "TRUE"){
#     ddups <- rbind(ddups, d[i,])
#   }
# }

#head(ddups)

length(unique(d$lab))  #2400
length(unique(d$lab2)) # 2619


# How many indiv of each sp are there?
d0 <- subset(d, day == "0")
table(d0$species)

d50 <- subset(d, day == "50")
table(d50$species)
head(d0)

tail(sort(unique(d$lab)))
tail(sort(unique(d$lab2)))

d <- as.data.frame(d);head(d)

# Ones that either died, broke, or were put back wrong according to my notes:
discard <- c("sm.HC.HP.LF.8.sorsco.1","sm.HC.LP.LF.13.sorsco.1","sm.HC.LP.HF.2.sorsco.1", "mp.Lc.LP.LF.1.spibet","sm.Hc.LP.LF.28.shecan","sm.LC.HP.HF.33.vibedu.2", "sm.LC.HP.HF.38.spibet.2", "mp.LC.HP.LF.28.menfer.1", "sm.LC.LP.HF.13.poptre.2", "sm.LC.LP.LF.21.sorsco.1", "sm.HC.HP.HF.25.vibedu.1", "mp.LC.HP.HF.15.acegla.1", "sm.LC.LP.HF.23.amealn.1", "mp.LC.LP.HF.15.alninc.1","sm.LC.HP.HF.2.loninv.1","sm.LC.LP.HF.13.alnvir.1","mp.LC.HP.LF.21.symalb.1","mp.HC.HP.HF.1.acegla.1","sm.HC.LP.HF.5.acegla.1","mp.HC.LP.HF.19.popbal.1","mp.HC.LP.HF.33.vacmem.1","mp.HC.LP.LF.38.corsto.2","sm.HC.HP.HF.10.rubpar.1", "mp.HC.HP.LF.27.vacmem.2", "sm.HC.HP.LF.5.rubpar.1",  "sm.HC.HP.LF.29.poptre.1", "sm.HC.HP.HF.31.alninc.1", "sm.LC.HP.HF.6.alnvir.1","sm.LC.HP.LF.12.sorsco.1", "sm.LC.LP.LF.8.spipyr.1",  "sm.LC.HP.HF.13.rubpar.1", "mp.LC.HP.LF.5.riblac.1",  "mp.HC.LP.HF.9.acegla.1",  "mp.HC.LP.HF.32.shecan.1","mp.HC.LP.HF.35.acegla.1", "mp.LC.HP.HF.3.poptre.1",  "mp.LC.HP.HF.20.acegla.1", "sm.LC.HP.HF.36.acegla.1", "sm.LC.HP.HF.38.vacmem.1","sm.LC.HP.LF.35.acegla.1","mp.LC.LP.HF.39.samrac.1","mp.LC.LP.LF.20.acegla.1","mp.HC.HP.HF.3.spibet.1",  "mp.HC.HP.HF.3.spibet.2", "sm.HC.HP.LF.3.sorsco.1","sm.HC.HP.LF.3.sorsco.2",  "sm.HC.HP.LF.34.sorsco.1", "sm.HC.LP.HF.10.sorsco.1", "sm.HC.LP.LF.3.sorsco.1", "sm.HC.LP.LF.20.sorsco.1","sm.HC.LP.LF.26.sorsco.1", "sm.HC.LP.LF.37.sorsco.1", "sm.HC.LP.LF.38.sorsco.1", "sm.HC.HP.LF.39.sorsco.1", "sm.HC.HP.HF.32.sorsco.1","sm.HC.HP.HF.35.sorsco.1", "sm.HC.LP.LF.19.sorsco.1", "sm.HC.LP.LF.24.sorsco.1", "sm.HC.HP.HF.22.sorsco.1")

deadones <- d[d$lab2 %in% discard,]

d <- d[!d$lab2 %in% discard,]

# Something happened on the 10th and the 15 that some rows got duplicated, not a killer mistake, since so few samples were doing anything yet

oops10 <- c("sm.LC.HP.LF.1.corsto.2","sm.LC.HP.LF.39.corsto.2","sm.LC.HP.LF.28.corsto.2","sm.LC.HP.LF.1.acegla.2","sm.LC.LP.HF.11.acegla.2", "mp.LC.LP.LF.1.acegla.2", "sm.LC.HP.LF.28.alninc.2","sm.LC.LP.HF.22.alninc.2","mp.LC.LP.HF.12.alnvir.2","sm.LC.HP.LF.17.loninv.2","sm.LC.LP.HF.28.loninv.2","sm.LC.LP.HF.27.menfer.2", "sm.LC.LP.HF.27.poptre.2","sm.LC.HP.LF.17.riblac.2", "sm.LC.LP.HF.22.riblac.2", "sm.LC.HP.HF.30.rubpar.2","mp.LC.LP.HF.12.shecan.2", "sm.LC.HP.LF.1.spibet.2", "sm.LC.LP.HF.9.spibet.2","sm.LC.LP.HF.22.spibet.2","sm.LC.HP.LF.28.spibet.2","sm.LC.HP.LF.39.spipyr.2","mp.LC.LP.HF.12.spipyr.2","sm.LC.LP.HF.28.spipyr.2","sm.LC.HP.LF.17.vibedu.2","sm.LC.HP.LF.39.vibedu.2","sm.LC.LP.HF.27.vibedu.2","sm.LC.LP.HF.22.vibedu.2","sm.LC.LP.HF.28.vibedu.2","sm.HC.LP.LF.7.vibedu.1","sm.LC.LP.HF.27.shecan.2","mp.LC.LP.HF.12.amealn.2","sm.LC.LP.HF.28.amealn.2","sm.LC.HP.LF.17.amealn.2","sm.LC.LP.HF.9.alninc.2","sm.LC.HP.HF.30.symalb.2","sm.LC.HP.LF.28.sorsco.2","sm.HC.LP.HF.16.sorsco.1","sm.LC.HP.LF.17.vibedu.2","sm.LC.LP.HF.27.vibedu.2","sm.LC.LP.HF.28.vibedu.2","sm.LC.LP.HF.22.vibedu.2","sm.HC.LP.LF.7.vibedu.1")

oops15 <- c("sm.LC.HP.LF.13.corsto.1", "sm.LC.LP.HF.31.acegla.2","sm.LC.HP.HF.31.alninc.2", "sm.LC.LP.LF.1.alninc.2","mp.LC.LP.HF.2.menfer.2","mp.LC.LP.HF.40.menfer.2","sm.LC.LP.HF.32.menfer.2", "mp.LC.HP.HF.20.poptre.1", "sm.LC_HP_HF.20.poptre.1" , "sm.LC.LP.LF.1.poptre.2","sm.LC.LP.LF.10.shecan.1","sm.LC.LP.LF.38.spibet.2", "sm.LC.HP.LF.39.spipyr.2", "mp.LC.LP.HF.21.spipyr.1","mp.LC.HP.HF.20.symalb.1","sm.LC.HP.LF.39.symalb.2","mp.LC.HP.HF.20.vacmem.1","sm.LC.HP.HF.30.vacmem.2","sm.LC.HP.LF.21.vacmem.1","sm.LC.LP.LF.1.vacmem.2","sm.LC.HP.LF.16.vibedu.2","sm.LC.HP.LF.39.vibedu.2","sm.LC.LP.LF.1.vibedu.2","mp.LC.LP.HF.30.sorsco.2", "mp.LC.HP.HF.20.symalb.1","mp.LC.HP.HF.20.sorsco.1"," mp.LC.LP.HF.21.rhoalb.1", "mp.LC.LP.HF.21.menfer.1", "sm.LC.LP.HF.13.poptre.2","sm.LC.LP.LF.30.shecan.2","sm.HC.HP.HF.8.sorsco.1","sm.HC.HP.HF.29.sorsco.1","sm.HC.LP.HF.25.sorsco.1","sm.HC.LP.HF.32.sorsco.1","sm.HC.LP.HF.3.sorsco.2","sm.LC.LP.LF.12.sorsco.2","sm.LC.HP.HF.19.sorsco.1","sm.LC.LP.HF.16.symalb.1","sm.LC.HP.HF.6.symalb.1","sm.LC.HP.HF.10.symalb.2","mp.LC.LP.HF.21.acegla.2","mp.LC.LP.HF.2.loninv.1","mp.LC.LP.HF.2.poptre.1","mp.LC.LP.HF.21.rhoalb.1","mp.HC.LP.HF.33.rubpar.2","sm.LC.HP.HF.30.symalb.1","sm.LC.HP.LF.17.shecan.2","sm.LC.LP.HF.37.shecan.1", "sm.LC.LP.HF.38.shecan.1","sm.LC.LP.HF.5.shecan.1","sm.LC.LP.LF.7.shecan.1","mp.LC.HP.HF.8.spibet.1", "mp.LC.LP.HF.39.spibet.1", "sm.LC.HP.HF.38.symalb.1", "sm.LC.HP.HF.38.symalb.2", "sm.LC.HP.HF.40.symalb.1","sm.HC.HP.LF.33.symalb.1", "sm.HC.HP.LF.33.symalb.2","mp.LC.HP.HF.20.sorsco.1","sm.LC.HP.LF.13.sorsco.1","sm.LC.HP.LF.21.sorsco.1","sm.LC.HP.LF.28.sorsco.2","sm.LC.HP.LF.7.sorsco.1","mp.LC.LP.HF.30.sorsco.2","sm.LC.LP.HF.9.sorsco.2","sm.HC.HP.HF.1.sorsco.1","sm.HC.HP.HF.9.sorsco.1","sm.HC.LP.HF.16.sorsco.1","sm.HC.LP.HF.3.sorsco.1", "sm.HC.LP.HF.3.sorsco.2", "sm.HC.LP.HF.31.sorsco.1","sm.LC.HP.HF.30.vacmem.2","sm.LC.HP.LF.21.vacmem.1","sm.LC.LP.LF.1.vacmem.2","sm.LC.LP.LF.1.vacmem.2","sm.LC.HP.LF.39.vibedu.2","sm.LC.LP.LF.1.vibedu.2","sm.LC.HP.LF.16.vibedu.2","mp.LC.HP.LF.24.sorsco.2","mp.LC.HP.LF.24.sorsco.1")

d <- d[!d$lab2 %in% oops10,]
d <- d[!d$lab2 %in% oops15,]

head(d)

# dtemp <- d[, 1:19]
# df_dat <- split(dtemp, dtemp$day)
# 
# for (i in 1:length(unique(dtemp$day))){
# write.csv(df_dat[i], paste('rcode/cleaning/cleaned_datafiles/day_', i,'.csv'), row.names = FALSE)
# }
# #lapply(


# Done! Writing the final datafile 
# write.csv(d,"input/bc_phenology_Feb52021.csv", row.names=FALSE)


# double checking i didn't miss any 
# c0 <- subset(d, species == "acegla" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "acegla" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "acegla" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "acegla" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "acegla" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "acegla" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "acegla" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "acegla" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "alninc" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "alninc" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "alninc" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "alninc" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "alninc" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "alninc" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "alninc" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "alninc" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "alnvir" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "alnvir" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "alnvir" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "alnvir" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "alnvir" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "alnvir" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "alnvir" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "alnvir" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "amealn" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "amealn" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "amealn" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "amealn" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "amealn" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "amealn" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "amealn" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "amealn" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "betpap" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "betpap" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "betpap" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "betpap" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "betpap" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "betpap" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "betpap" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "betpap" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "corsto" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "corsto" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "corsto" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "corsto" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "corsto" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "corsto" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "corsto" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "corsto" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "loninv" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "loninv" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "loninv" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "loninv" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "loninv" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "loninv" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "loninv" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "loninv" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "menfer" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "menfer" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "menfer" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "menfer" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "menfer" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "menfer" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "menfer" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "menfer" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "popbal" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "popbal" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "popbal" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "popbal" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "popbal" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "popbal" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "popbal" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "popbal" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "poptre" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "poptre" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "poptre" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "poptre" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "poptre" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "poptre" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "poptre" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "poptre" & treatment == "HC.LP.LF");table(c0$lab2)
# #
# c0 <- subset(d, species == "rhoalb" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "rhoalb" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "rhoalb" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "rhoalb" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "rhoalb" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "rhoalb" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "rhoalb" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "rhoalb" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "riblac" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "riblac" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "riblac" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "riblac" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "riblac" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "riblac" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "riblac" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "riblac" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "rubpar" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "rubpar" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "rubpar" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "rubpar" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "rubpar" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "rubpar" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "rubpar" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "rubpar" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "samrac" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "samrac" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "samrac" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "samrac" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "samrac" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "samrac" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "samrac" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "samrac" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "shecan" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "shecan" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "shecan" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "shecan" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "shecan" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "shecan" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "shecan" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "shecan" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "spibet" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "spibet" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "spibet" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "spibet" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "spibet" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "spibet" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "spibet" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "spibet" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "spipyr" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "spipyr" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "spipyr" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "spipyr" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "spipyr" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "spipyr" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "spipyr" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "spipyr" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "symalb" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "symalb" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "symalb" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "symalb" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "symalb" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "symalb" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "symalb" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "symalb" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "sorsco" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "sorsco" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "sorsco" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "sorsco" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "sorsco" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "sorsco" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "sorsco" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "sorsco" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "vacmem" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "vacmem" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "vacmem" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "vacmem" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "vacmem" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "vacmem" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "vacmem" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "vacmem" & treatment == "HC.LP.LF");table(c0$lab2)
# 
# c0 <- subset(d, species == "vibedu" & treatment == "LC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "vibedu" & treatment == "LC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "vibedu" & treatment == "LC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "vibedu" & treatment == "LC.LP.LF");table(c0$lab2)
# c0 <- subset(d, species == "vibedu" & treatment == "HC.HP.HF");table(c0$lab2)
# c0 <- subset(d, species == "vibedu" & treatment == "HC.HP.LF");table(c0$lab2)
# c0 <- subset(d, species == "vibedu" & treatment == "HC.LP.HF");table(c0$lab2)
# c0 <- subset(d, species == "vibedu" & treatment == "HC.LP.LF");table(c0$lab2)

# sorsco <-subset(d, species =="sorsco")
# sort((unique(sorsco$lab2)))
# length(unique(sorsco$lab2))
# 
# symalb <-subset(d, species =="symalb")
# sort((unique(symalb$lab2)))
# length(unique(sorsco$lab2))
# 
# trt.succ <- symalb %>%
#   group_by(treatment) %>%
#   summarise(no_rows = length(unique(lab2)))
# 
# rubpar <-subset(d, species =="rubpar")
# sort((unique(rub$lab2)))
# length(unique(sorsco$lab2))
# 
# trt.succ <- rubpar %>%
#   group_by(treatment) %>%
#   summarise(no_rows = length(unique(lab2)))
# 
# amealn <-subset(d, species =="amealn")
# 
# trt.succ <- amealn %>%
#   group_by(treatment) %>%
#   summarise(no_rows = length(unique(lab2)))
# 
# menfer <-subset(d, species =="menfer")
# 
# trt.succ <- menfer %>%
#   group_by(treatment) %>%
#   summarise(no_rows = length(unique(lab2)))
# 
# rhoalb <-subset(d, species =="rhoalb")
# 
# trt.succ <- rhoalb %>%
#   group_by(treatment) %>%
#   summarise(no_rows = length(unique(lab2)))
# 
# loninv <-subset(d, species =="loninv")
# 
# trt.succ <- loninv %>%
#   group_by(treatment) %>%
#   summarise(no_rows = length(unique(lab2)))
# 
# shecan <-subset(d, species =="shecan")
# 
# trt.succ <- shecan %>%
#   group_by(treatment) %>%
#   summarise(no_rows = length(unique(lab2)))
# 
# poptre <-subset(d, species =="poptre")
# 
# trt.succ <- poptre %>%
#   group_by(treatment) %>%
#   summarise(no_rows = length(unique(lab2)))
# 
# Done! Writing the final datafile 
# write.csv(d,"input/bc_phenology_Feb52021.csv", row.names=FALSE)

# There are a few unexplained extras:
#corsto
# c1 <- subset(d, species == "corsto" & treatment == "LC.HP.LF")
# sp[16]
# for (i in 1:length(unique(c1$flask))){
#   sp <- sort(unique(c1$lab2))
#   sub <- subset(c1, lab2== sp[16])
#   
#   #spplot <- 
#   ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/corsto",i,".pdf"))
# }
# 
# # Poptre
# c1 <- subset(d, species == "shecan" & treatment == "HC.HP.HF")
# c2 <- subset(d, species == "shecan" & treatment == "LC.LP.HF")
# c3 <- subset(d, species == "shecan" & treatment == "LC.LP.LF")
# 
# sp[16]
# for (i in 1:length(unique(c1$lab2))){
#   sp <- sort(unique(c1$lab2))
#   sub <- subset(c1, lab2== sp[i])
#   
#   spplot <- 
#   ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/shecan",i,".pdf"))
# }
# 
# for (i in 1:length(unique(c2$lab2))){
#   sp <- sort(unique(c2$lab2))
#   sub <- subset(c2, lab2== sp[i])
#   
#   spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/shecan",i,".pdf"))
# }
# sp[15]
# 
# 
# for (i in 1:length(unique(c3$lab2))){
#   sp <- sort(unique(c3$lab2))
#   sub <- subset(c3, lab2== sp[i])
#   
#   spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/shecan",i,".pdf"))
# }
# sp[5]
# temp <- subset(d, lab2 == "sm.LC.LP.LF.30.shecan.3")
# 
# #sorsco
# c1 <- subset(d, species == "sorsco" & treatment == "HC.HP.HF") #19
# c2 <- subset(d, species == "sorsco" & treatment == "HC.LP.HF") #21
# c3 <- subset(d, species == "sorsco" & treatment == "LC.LP.HF") #20
# c4 <- subset(d, species == "sorsco" & treatment == "LC.LP.LF") #17
# c5 <- subset(d, species == "sorsco" & treatment == "LC.HP.HF") #17
# 
# table(c1$population)
# 
# sp[16]
# for (i in 1:length(unique(c1$lab2))){
#   sp <- sort(unique(c1$lab2))
#   sub <- subset(c1, lab2== sp[i])
#   
#   spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/sorsco",i,".pdf"))
# }
# 
# for (i in 1:length(unique(c2$lab2))){
#   sp <- sort(unique(c2$lab2))
#   sub <- subset(c2, lab2== sp[i])
#   
#   spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/sorsco",i,".pdf"))
# }
# sp[18]
# 
# 
# for (i in 1:length(unique(c3$lab2))){
#   sp <- sort(unique(c3$lab2))
#   sub <- subset(c3, lab2== sp[i])
#   
#   spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/sorsco",i,".pdf"))
# }
# sp[3]
# 
# for (i in 1:length(unique(c4$lab2))){
#   sp <- sort(unique(c4$lab2))
#   sub <- subset(c4, lab2== sp[i])
#   
#   spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/sorsco4",i,".pdf"))
# }
# sp[4]
# 
# 
# for (i in 1:length(unique(c5$lab2))){
#   sp <- sort(unique(c5$lab2))
#   sub <- subset(c5, lab2== sp[i])
#   
#   spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.t)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/sorsco5",i,".pdf"))
# }
# sp[13]
# 
# #symalb
# c1 <- subset(d, species == "symalb" & treatment == "LC.HP.HF") #21
# c2 <- subset(d, species == "symalb" & treatment == "LC.LP.HF") #20
# d$bbch.l <- as.numeric(d$bbch.l)
# table(c1$treatment)
# 
# sp[16]
# for (i in 1:length(unique(c1$lab2))){
#   sp <- sort(unique(c1$lab2))
#   sub <- subset(c1, lab2== sp[22])
#   
#   #spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.l)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#  #ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/symalblat",i,".pdf"))
# }
# sp[3]
# sp[12]
# 
# for (i in 1:length(unique(c2$lab2))){
#   sp <- sort(unique(c2$lab2))
#   sub <- subset(c2, lab2== sp[i])
#   
#   #spplot <- 
#     ggplot(sub)+
#     aes(x=day, y=bbch.l)+
#     geom_line() +
#     labs(x="day",y="bbch") +
#     xlim(-0.5,113)+
#     ylim(0,20)
#   #ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/checking/sorsco",i,".pdf"))
# }
# sp[18]
# 
# 
# unique(c2$lab2)
# temp <- subset(d, lab2 == "mp.LC.HP.HF.30.symalb.2")

dead <- c("dead","missing","bud driedout","buds dead","I guess not dead…","term dead","buds died", "died","term dead; dead","extra poptre","extra spibet","buds dried out","bud dead","died -buds dead","died -remaining buds deaed","dead, remaining buds dead","dies;", "died, buds dried", "died,","dies","Dead",  "extra spipyr",  "extra shecan","I guess not dead?","buds are dead","remaining buds are dead","dead -remianing buds dead","dies")
dcom <- d[d$comments %in% dead,]

dcom$bbch.l <- as.numeric(as.character(dcom$bbch.l)) 
dcom.t <- subset(dcom, bbch.t < 7)
dcom.tl <- subset(dcom.t, bbch.l < 7)

unique(dcom.tl$lab2)

ex<-c("XX","x", "xx")
dex <- d[d$indiv %in% ex,]

idk <- c("I guess not dead…","I guess not dead?")
didk <- d[d$comments %in% idk,]
dfeb5 <- read.csv("input/bc_phenology_Feb52021.csv", header=TRUE, na.strings=c("","NA"))