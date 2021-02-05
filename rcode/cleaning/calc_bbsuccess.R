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
d <- read.csv("input/bc_phenology_Feb42021.csv", header=TRUE, na.strings=c("","NA"))

#source("rcode/cleaning/cleaningcode.R")

# what does the data look like generally?
# 20 species from mp,20 species from sm; so in theory there should be 2560 samples, but after chilling we had 
21*8*8*2
alive <- subset(d, bbch.l > 0)
length(unique(alive$lab2)) #2313

dead <- subset(d, bbch.l < 0)
1-(2316/2560)
# so 9.25 of samples did not budburst, either because they were dead, or becuase of insufficient conditions


# How many indiv of each sp are there?
d0 <- subset(d, day == "0")
table(d0$species)

d88 <- subset(d, day == "88")
table(d88$species)

d <- as.data.frame(d)

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

trt.succ <- terminalbb %>%
  group_by(species) %>%
  summarise(no_rows = length(species))

amealn <- subset(terminalbb, species == "amealn")

######################################################
dlat <- gc
#Task is to select bbch.1 7 and above, sum percentages, then get 1st day where percentage above 80%
#1. reshape data so it is in long format 
dlong <- gather(dlat, key = "bbchL", value = "stage", c(bbch.l, bbch2.l, bbch3.l))
#dlong <- gather(dlong, key = "percentL", value = "l.percent", c(percent.l, percent2.l,percent3.l))

#put relevent percentages with stages - a bit clunky but does the job 
dlong$bbchPercent <- dlong$percent.l
dlong$bbchPercent [dlong$bbchL == "bbch2.l"] <- dlong$percent2.l [dlong$bbchL == "bbch2.l"] 
dlong$bbchPercent [dlong$bbchL == "bbch3.l"] <- dlong$percent3.l [dlong$bbchL == "bbch3.l"] 

#head(dlong)
#str(dlong)

#Remove na rows
dlong <- dlong[!is.na(dlong$stage), ]

#1. Select bbch.1 7 and above
dlong7 <- dlong[dlong$stage >= 7, ]

#sum percentages each sample
#sumPercent <- aggregate(dlong7$bbchPercent, by = list(Category = dlong7$lab2), FUN = sum)
sumPercent <- aggregate(dlong7$bbchPercent, by = list(Category = dlong7$lab2, day = dlong7$day), FUN = sum)

names(sumPercent) <- c("lab2", "day", "sumPercent")
#names(sumPercent) <- c("lab2", "day.l.bb", "sumPercent")

dlong7sum <- merge(dlong7, sumPercent, by = c("lab2", 'day'))

dead <- subset(dlong7sum, sumPercent < 0)
