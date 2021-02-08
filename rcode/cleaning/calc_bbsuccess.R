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
source('rcode/cleaning/pheno_bb_calc.R')

#source("rcode/cleaning/cleaningcode.R")

# what does the data look like generally?
# 20 species from mp,20 species from sm; so in theory there should be 2560 samples, but after chilling we had 

# How many indiv of each sp are there?

d88 <- subset(d, day == "88")
surv <- sum(table(d88$species))
table(d88$species) # acegla had the worst survivorship, followed by sorsco, rubpar, vacmem, spibet and spirpyr

initial <- 18 * 128 + 3 * 64 # 2496 samples went into chilling, 2458 went into forcing and survived the experiment

1 - surv/initial # had 1.52 % mortatility 

###### Excluding the dead, how many samples did not budburst? ###############

smpin <- d88 %>%
  group_by(species) %>%
  summarise(no_rows = length(species))

smpbb <- pheno %>%
  group_by(species) %>%
  summarise(no_rows = length(species))

spwithbb <- merge(smpin, smpbb, by = "species")
spwithbb$prop.bb <- spwithbb$no_rows.y/spwithbb$no_rows.x
spbb <- spwithbb[, c("species","prop.bb")]
names(spbb)[names(spbb) == "species"] <- "Species"
names(spbb)[names(spbb) == "prop.bb"] <- "Proportion Budburst"

write.table(spbb, "output/bb.success.ax.trt.csv", sep = ",", row.names = FALSE)

