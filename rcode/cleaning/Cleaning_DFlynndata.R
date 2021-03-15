# I want to merge the DFlynn data with my data, so first I have to rename to columns etc. 

if(length(grep("deirdreloughnan", getwd()) > 0)) {
  setwd("~/Documents/github/pheno_bc/")
} else {
  setwd("~/Documents/github/pheno_bc")
}

require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
# require(chillR)

rm(list = ls()) 
options(stringsAsFactors = FALSE)

# read in the cleaning phenology data:
df <- read.csv("input/dflynn.data.csv", header=TRUE, na.strings=c("","NA"))
head(df)

dl <- read.csv("input/bc_phenology_Feb52021.csv", header=TRUE, na.strings=c("","NA"))
head(dl)

# In the DF data they used their own scale, where bbch 7 is called 3, 

df <- df[, c("id", "sp", "rep","site","ind","treatcode","warm","photo","chill", "tleaf","lleaf", "day")]

# changing column names 
names(dl)
names(df)
