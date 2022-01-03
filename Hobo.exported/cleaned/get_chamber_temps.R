# January 2, 2022
# aim of this code is to correct for HOBO callibration and get the range of temperatures of the chambers used
# in the experiment
rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/github/pheno_bc/Hobo.exported/cleaned")

require(tidyr)
require(plyr)
require(readr)

files <- list.files(path = "~/Documents/github/pheno_bc/Hobo.exported/cleaned", pattern =".csv" )
files

#data <- ldply(files, read_csv, na = c("","NA"))
data <- ldply(files, read_csv, na = "")

head(data)

data$corr <- data$temp - data$callibration

day <- subset(data, day_night == "day")
day_hf <- subset(day, force == "high")
day_lf <- subset(day, force == "low")

night <- subset(data, day_night == "night")
night_hf <- subset(night, force == "high")
night_lf <- subset(night, force == "low")

##what are the ranges of temperatures
range(day_hf$corr); hist(day_hf$corr); mean(day_hf$corr); sd(day_hf$corr); median(day_hf$corr)
range(day_lf$corr); hist(day_lf$corr); mean(day_lf$corr); sd(day_lf$corr); median(day_lf$corr)

range(night_hf$corr); hist(night_hf$corr); mean(night_hf$corr); sd(night_hf$corr); median(night_hf$corr)
range(night_lf$corr); hist(night_lf$corr); mean(night_lf$corr); sd(night_lf$corr); median(night_lf$corr)
