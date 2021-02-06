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


##### Adding individual ############
# indiv <- read.csv("input/indiv.no.cleaned.csv", na.strings = "")

# source("rcode/cleaning/cleaning.indivno.R")
# 
# head(indiv)
# #subset to just the colns needed
# indiv <- indiv[ , c("lab3","indiv")]
# # indiv[complete.cases(indiv), ]
# # test <- subset(indiv, indiv!= "NA")
# head(indiv)
# 
# 
# dtemp <- merge(d, indiv, by = "lab2", all.x = TRUE) # this is adding rows! 

# mptemp <- subset(indiv, site == "mp");length(unique(mptemp$labtemp)) # 204
# 
# smtemp <- subset(indiv, site == "sm");length(unique(smtemp$labtemp)) #186


####################################################################

begin <- subset(d, day == 0)
count.begin <- table(begin$species)
sum(count.begin)

end <- subset(d, day == 88)
count.end <- table(end$species)
sum(count.end)

d$day <- as.numeric(d$day)
range(d$day, na.rm = TRUE) # this dataset includes the low chill greenhouse days, 
sort(unique(data$species))

#At the experiment, there were 2499 samples (three extra samples snuck in), by the beginning of forcing, there were still 2495 samples alive. 
#############################################################
#    Starting to work with the terminal buds first 
#############################################################
# Since this dataset includes the LC's days in the greenhouse, I need to subset those out and just have the 12 weeks they spent in the growth chambers:
# in total there were 2366/ 2496
gc <- subset(d, day <= 84)
unique(gc$day)
#head(gc)
gc <- as.data.frame(gc)
# I am curious if all samples even reached stage 7 or...

max <- gc %>% 
  group_by(lab2, ref) %>%
  slice(which.max(bbch.t))

low <- subset(max, bbch.t<7)
# there are 457 out of 2241 samples for which the terminal bud did not burst or approximately 17%

count <- table(low$species)
sum(count)
# For species in both pops there were 128 samples max, so rub par 40% of the time the terminal bud did not bb, for ace gla and sorsco it was 25%
#460 samples did not have terminal bb; 
459/2496 #18.4
####################################################################################################################
# This code is taken from Dan Flynn and cleaning the east coast data (found in the budchill repo)
# 1. for both terminal and lateral buds, what is the max stage within a row? Identify which rows are greater or equal to the specific BBCH stage
# 2. now for that individual, find the earliest day at which that stage was reached.
####################################################################################################################

tbb <- nl <- vector() # This creates empty vectors for the bbday, leaf out, the nl yes no vector of whether this event occured (1 means yes, there was leafout; 0 means no leafout)
gc$lab2 <- as.factor(gc$lab2)
#levels(d$lab2)
for(i in levels(gc$lab2)){ # i = levels(d$lab2)[2505] # for each individual clipping. # DL: why did DF have 602, that seems low

  dx <- gc[gc$lab2 == i,]
  bdax <- which(apply(dx[, c("bbch.t","bbch.t")], MARGIN = 1, max, na.rm = TRUE) >3) # margin = 1 means it is applied over rows, takes the maximum value > 3
  if(length(bdax) < 1) {bdax = NA; nl <- c(nl, 0)} else {bdax = dx[min(bdax), 'day']; nl <- c(nl, 1)}

  tbb <- c(tbb, bdax)

}
dx <- gc[match(levels(gc$lab2), gc$lab2),] # with twig id in same order as the loop above
#dx <- dx[,2:ncol(dx)]
dx <- dx[,c("lab2", "population", "treatment", "flask", "species"#, "indiv"
            )]
terminalbb <- data.frame(dx, tbb, nl)

dfly <- terminalbb %>%
  group_by(species) %>%
  summarise(no_rows = length(species))

smpin <- d88 %>%
  group_by(species) %>%
  summarise(no_rows = length(species))

dfwithbb <- merge(smpin, dfly, by = "species")
dfwithbb$prop.bb <- spwithbb$no_rows.y/spwithbb$no_rows.x

sorsco <- subset(d, species == "sorsco")
sorsco.t <- subset( sorsco, bbch.t >= 7)
length(unique(sorsco.t$lab2))

terminalss <- aggregate(sorsco.t["day"],
                        sorsco.t[c("lab2", "population", "treatment", "flask", "species")], 
                        FUN = min)
names(terminalss)[names(terminalss) == "day"] <- "tbb"
head(terminalss)
unique(terminalss$lab2)


alninc <- subset(d, species == "alninc")
alninc.t <- subset( alninc, bbch.t >= 7)
length(unique(alninc.t$lab2))

terminalai <- aggregate(alninc.t["day"],
                        alninc.t[c("lab2", "population", "treatment", "flask", "species")], 
                        FUN = min)
names(terminalai)[names(terminalai) == "day"] <- "tbb"
head(terminalai)
unique(terminalai$lab2)

alninc <- subset(latdaymin1, species == "alninc") ; length(unique(alninc$lab2))
alninc <- subset(pheno, species == "alninc") ; length(unique(alninc$lab2))
alninc.t <- subset( bursted, bbch.t >= 7 & species == "alninc"); length(unique(alninc.t$lab2))
alninc.t <- subset( terminalbb, species == "sorsco"); length(unique(alninc.t$lab2))
