# This code will compile the final data file for use in my phenology analysis

# The aim is to extract the day of study on which the terminal bud reached stage 7 first and the day when the lateral buds reached 80% at stage 7

if(length(grep("deirdreloughnan", getwd())>0)) { 
  setwd("~/Documents/github/pheno_bc") 
} else
  setwd("~/Documents/github/pheno_bc")

rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(plyr)
require(dplyr)
require(tidyr)

# read in the cleaning phenology data:
data<-read.csv("input/bc_phenology.csv")
data<-data[,(2:17)] # getting rid of the X count column
# Starting with the terminal buds:

# Would it be useful to have a unique identifying for every sample?
data$lab<-paste(data$population,data$treatment,data$flask, data$species, sep="_")

pheno<-data %>% 
  separate(treatment, c("chill","photo","force"), "_")

head(pheno)

# I have 21 species, but rhoalb, betpap, samrac were only at one site
18*8*8*2+3*8*8 # there should be 2496 samples

length(unique(pheno$lab)) # there are 2485 unique samples

2496-2485 # only had 11 die! 
# Should look into which ones these are...

sort(unique(data$species))
#############################################################
 #    Starting to work with the terminal buds first 
#############################################################

# I am curious if all samples even reached stage 7 or...

max<-pheno %>% 
  group_by(lab) %>%
  slice(which.max(bbch.t))

low<-subset(max, bbch.t<7)
# there are 440 samples for which the terminal bud did not burst
