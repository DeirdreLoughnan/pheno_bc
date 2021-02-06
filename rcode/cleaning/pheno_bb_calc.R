# This code will compile the final data file for use in my phenology analysis

# The aim is to extract the day of study on which the terminal bud reached stage 7 first and the day when the lateral buds reached 80% at stage 7

# D. Flynn calcualted the day that both lat and terminal buds bb, but here I want to treat them seperately and look only at one or the other

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
d <- read.csv("input/bc_phenology_Feb52021.csv", header=TRUE, na.strings=c("","NA"))
head(d)
#source("rcode/cleaning/cleaningcode.R")

# what does the data look like generally?
# 20 species from mp,20 species from mp; so in theory there should be 2560 samples, but after chilling we had 

# alive <- subset(data, bbch.l > 0)
# length(unique(alive$lab2)) #2319
# 
# 1-(2319/2560)
# # so 9.4 of samples did not budburst, either because they were dead, or becuase of insufficient conditions
# 
#
# # How many indiv of each sp are there?
# d0 <- subset(d, day == "0")
# table(d0$species)
# 
# d88 <- subset(d, day == "88")
# table(d88$species)
# 
# d <- as.data.frame(d)

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

begin <- subset(d, day == 50)
count.begin <- table(begin$species)
sum(count.begin)

end <- subset(d, day == 88)
count.end <- table(end$species)
sum(count.end)

d$day <- as.numeric(d$day)
range(d$day, na.rm = TRUE) # this dataset includes the low chill greenhouse days, 
sort(unique(d$species))


#############################################################
 #    Starting to work with the terminal buds first 
#############################################################
# Since this dataset includes the LC's days in the greenhouse, I need to subset those out and just have the 12 weeks they spent in the growth chambers:

gc <- subset(d, day <= 84)
unique(gc$day)
#head(gc)
gc <- as.data.frame(gc)
# I am curious if all samples even reached stage 7 or...

max <- gc %>% 
  group_by(lab2, ref) %>%
  slice(which.max(bbch.t))

low <- subset(max, bbch.t<7)
# there are 405 out of 2406 samples for which the terminal bud did not burst or approximately 17%

count <- table(low$species)

#########################################################################

# starting by subsetting the observations that are above the level of bb, which occurs at stage 7
bursted  <- subset(gc, bbch.t >= 7)

# for every unique, label, pop, trt, flask, sp, it is getting the minimum (ie first) day that an obs of 7 or greater was observed
terminalbb <- aggregate(bursted["day"],
                        bursted[c("lab2", "population", "treatment", "flask", "species")], 
                        FUN = min)
names(terminalbb)[names(terminalbb) == "day"] <- "tbb"
head(terminalbb)

# below is the code used to show that the DFlynn and my data were different 
# head(terminalbb) # here is the data for the terminal bud, there are 2406 individuals, with a value for each

# nobb <- subset(terminalbb, nl == 0) # 209 samples didnt bb
# table(nobb$chill) # mostly low chill
# table(nobb$force) # mostly low force
# table(nobb$photo) # mostly low photoperiod
# 
# 
# comp <- merge(flynn, faith, by = c("lab2","population","treatment","flask","species"))
# comp <- comp[, c("lab2","population","treatment","species","tbb","day")]
# head(comp)
# 
# pdf(file = "figures/compar_bbcalc.pdf")
# ggplot(comp) +
#       aes(x = day, y = tbb, color = species) +
#       geom_point() +
#       labs(x = "Faith's",y = "Flynn's") +
#       facet_wrap("treatment") 
# dev.off()
#############################################################
#    Now moving on to the lateral buds 
#############################################################

#What if I would to use a similar approach as above, but start by calculating when the three bbch levels over phase 7 sum to 80%?
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
#head(dlong7sum)

#Select samples with 80 or more percent --> 1068 rows, if it is lower at 50, then there are 1437 rows
dlong7sum80 <- dlong7sum[dlong7sum$sumPercent >= 80, ]

#Select first day with 80% or more
latdaymin80 <- aggregate(dlong7sum80$day, by = list(Category = dlong7sum80$lab2), FUN = min)
names(latdaymin80) <- c("lab2", "latbb80")

#dlong7sum80Daymin <- merge(dlong7sum80, daymin, by = "lab2")
length(unique(latdaymin80$lab2))

#Select first day with 50%, then there are 1450 rows
dlong7sum50 <- dlong7sum[dlong7sum$sumPercent >= 50, ]
latdaymin50 <- aggregate(dlong7sum50$day, by = list(Category = dlong7sum50$lab2), FUN = min)
names(latdaymin50) <- c("lab2", "latbb50")
nrow(latdaymin50)

#Select first day of lateral bb, then there are  rows 1923
dlong7sum1 <- dlong7sum[dlong7sum$sumPercent > 0, ]
latdaymin1 <- aggregate(dlong7sum1$day, by = list(Category = dlong7sum1$lab2), FUN = min)
names(latdaymin1) <- c("lab2", "latbb1")
nrow(latdaymin1)

# latdaymin80 <- aggregate(dlong7sum80["day"], dlong7sum80[c("lab2", "population", "treatment", "flask", "species")], FUN = min)
# names(latdaymin80)[names(latdaymin80) == "day"] <- "latbb80"
# #names(latdaymin80) <- c("lab2", "latbb80")
# 
# #dlong7sum80Daymin <- merge(dlong7sum80, daymin, by = "lab2")
# length(unique(latdaymin80$lab2))
# 
# #Select first day with 50%, then there are 1450 rows
# dlong7sum50 <- dlong7sum[dlong7sum$sumPercent >= 50, ]
# latdaymin50 <- aggregate(dlong7sum50["day"], dlong7sum50[c("lab2", "population", "treatment", "flask", "species")], FUN = min)
# names(latdaymin50)[names(latdaymin50) == "day"] <- "latbb50"
# 
# #Select first day of lateral bb, then there are  rows 1923
# dlong7sum1 <- dlong7sum[dlong7sum$sumPercent > 0, ]
# latdaymin1 <- aggregate(dlong7sum1["day"], dlong7sum1[c("lab2", "population", "treatment", "flask", "species")], FUN = min)
# names(latdaymin1)[names(latdaymin1) == "day"] <- "latbb1"


######################################################
# Combine the terminal and the lateral bb days
#head(terminalbb)
pheno.temp <- merge(latdaymin1, latdaymin50, by = c("lab2"), all.x = TRUE)
pheno.temp2 <- merge(pheno.temp, latdaymin80, by = "lab2", all.x = TRUE) # the all.x = T is telling it that I want all rows from this dataset, even if there isn't a corresponding row in latdaymin
pheno <- merge(pheno.temp2, terminalbb, by = "lab2", all.x = TRUE, all.y = T) 
pheno$lab3 <- pheno$lab2

# pheno <- pheno[,c("lab2", "population.x", "treatment.x", "flask.x", "species.x","latbb1", "latbb50","latbb80","tbb")]
# names(pheno)[names(pheno) == "population.x"] <- "population"
# names(pheno)[names(pheno) == "treatment.x"] <- "treatment"
# names(pheno)[names(pheno) == "flask.x"] <- "flask"
# names(pheno)[names(pheno) == "species.x"] <- "species"


pheno <- pheno %>% separate(lab3, c("population","chill", "photo","force","flask","species", "rep"), convert = T)

#pheno$treatment <- paste(pheno$chill, pheno$photo, pheno$force, sep = ".")
head(pheno)
write.csv(pheno, "input/day.of.bb.Feb52021.csv", row.names = FALSE)
######################################################
# To make it more comparable to the Flynn dataset, I am adding a treatment column, and then try to calculate chill portions...for the terminal bud? 



# making fequncy tables:
# require(plyr)
# 
# term <- pheno[, c("population","species","tbb","treatment","lab2")]; head(term)
# term.cc <- term[complete.cases(term), ]
# 
# term.summ <- term.cc %>%
#   count(treatment, species, population)
# 
# lat80 <- pheno[, c("population","species","latbb80","treatment","lab2")]; head(lat80)
# lat80.cc <- term[complete.cases(lat80), ]
# 
# lat80.summ <- lat80.cc %>%
#   count(treatment, species, population)
# 
# lat50 <- pheno[, c("population","species","latbb50","treatment","lab2")]; head(lat50)
# lat50.cc <- term[complete.cases(lat50), ]
# 
# lat50.summ <- lat50.cc %>%
#   count(treatment, species, population)



# PLOTS
# start by plotting means by species and by treatments:
# phenoplot_sp <- ddply(pheno, .(treatment, species), summarize, termbb = mean(tbb, na.rm = TRUE), lat50bb = mean(latbb50, na.rm = TRUE))
# 
# ggplot(phenoplot, aes(species, termbb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))
# 
# # plotting means by just treatments:
# phenoplot_trt <- ddply(pheno, .(treatment), summarize, termbb = mean(tbb, na.rm = TRUE), lat50bb = mean(latbb50, na.rm = TRUE))
# 
# #terminal dobb
# ggplot(phenoplot_sp, aes(treatment, termbb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) + 
#   scale_y_continuous(limits = c(0, 88))
# 
# #lateral dobb
# ggplot(phenoplot_sp, aes(treatment, lat50bb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) + 
#   scale_y_continuous(limits = c(0, 88))
# 
# # Just the mean treatment value
# # terminal
# ggplot(phenoplot_trt, aes(treatment, termbb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) + 
#   scale_y_continuous(limits = c(0, 88))
# 
# ggplot(phenoplot_trt, aes(treatment, lat50bb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) + 
#   scale_y_continuous(limits = c(0, 88))
# 
# 
#  hhh <- subset(pheno, treatment == "HC_HP_HF")
# # hll <- subset(pheno, treatment == "HC_LP_LF")
# # hlh <- subset(pheno, treatment == "HC_LP_HF")
# # hhl <- subset(pheno, treatment == "HC_HP_LF")
# # lll <- subset(pheno, treatment == "LC_LP_LF")
# # lhh <- subset(pheno, treatment == "LC_HP_HF")
# # llh <- subset(pheno, treatment == "LC_LP_HF")
# # lhl <- subset(pheno, treatment == "LC_HP_LF")
# 
# 
# 
# # Dan Flynn plots
# 
# colz <- c("darkorchid","blue3", "cadetblue","coral3")
# lcol <- alpha(colz, 0.1)
# names(lcol) = levels(pheno$chill)
# 
# 
# 
# d <- pheno
# 
# pdf( width = 8, height = 4)
# par(mfcol = c(1, 3), mar = c(3, 3, 1, 0.5))
# for(spx in levels(d$species)){ # spx = "BETALL"
#   
#   dxx = d[d$species == spx, ]
#   
#   counter = 1
#   for(i in sort(as.character((unique(d$chill))))){# i = "high chill or low chill"
#     
#     dseq = seq(0, max(dxx$day))
#     plot(dseq, seq(0, 7, length = length(dseq)), type = "n", 
#          ylab2 = "Stage",
#          xlab2 = "")
#     if(counter == 1) mtext(spx, line = -2, adj = 0.5)
#     legend("topleft", bty = "n", i, cex = 0.85, inset = 0)
#     xx <- dxx[dxx$time == i,]
#     # calculate mean response by date and chill
#     xt <- tapply(pmax(xx$tleaf, xx$lleaf, na.rm = T), list(xx$day, xx$chill), mean, na.rm = TRUE)
#     
#     for(j in unique(xx$ind)){ #j = unique(xx$ind)[1]
#       xj <- xx[xx$ind == j, ]
#       pcol = lcol[xj$chill]
#       lines(xj$day, xj$tleaf, col = pcol)
#     }
#     lines(rownames(xt), xt[, 1], col = colz[1], lwd = 2)
#     lines(rownames(xt), xt[, 2], col = colz[2], lwd = 2)
#     lines(rownames(xt), xt[, 3], col = colz[3], lwd = 2)
#     lines(rownames(xt), xt[, 4], col = colz[4], lwd = 2)
#     
#     
#     # add a legend
#     if(counter == 3) {    
#       legend("topright", bty = "n",
#              col = colz,
#              lwd = 2,
#              legend = c(1, 2, 4, 8),
#              title = "°C")
#     }
#     counter = counter + 1
#   }
# }
# dev.off()
# system(paste("open '", paste("figures/Trace Plots ", Sys.Date(), ".pdf", sep=""), "' -a /Applications/Preview.app", sep=""))
# 
#write.csv(pheno, "input/day.of.bb.csv")
#

