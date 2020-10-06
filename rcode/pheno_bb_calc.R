# This code will compile the final data file for use in my phenology analysis

# The aim is to extract the day of study on which the terminal bud reached stage 7 first and the day when the lateral buds reached 80% at stage 7
# building off of Dan Flynn's data

if(length(grep("deirdreloughnan", getwd())>0)) { 
  setwd("~/Documents/github/pheno_bc") 
} else {
  setwd("~/Documents/github/pheno_bc")
}

rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(plyr)
require(dplyr)
require(tidyr)

# read in the cleaning phenology data:
data<-read.csv("input/bc_phenology.csv")

# Starting with the terminal buds:

# Would it be useful to have a unique identifying for every sample?
data$lab<-paste(data$population,data$treatment,data$flask, data$species, sep="_")

d<-data %>% 
  separate(treatment, c("chill","photo","force"), "_")

head(d)

begin<-subset(d, day==0)
count.begin<-table(begin$species)
sum(count.begin)

end<-subset(d, day==88)
count.end<-table(end$species)
sum(count.end)

d$day<-as.numeric(d$day)
range(d$day, na.rm=TRUE) # this dataset includes the low chill greenhouse days, 
sort(unique(data$species))

#The final dataset has 2406 unique individuals, while I started with 2560
(2560-2406)/2560
# There was approximately 6% mortality 
#############################################################
 #    Starting to work with the terminal buds first 
#############################################################
# Since this dataset includes the LC's days in the greenhouse, I need to subset those out and just have the 12 weeks they spent in the growth chambers:
head(gc)
gc<-subset(d, day<=84)
unique(gc$day)

# I am curious if all samples even reached stage 7 or...

max<-gc %>% 
  group_by(lab) %>%
  slice(which.max(bbch.t))

low<-subset(max, bbch.t<7)
# there are 405 out of 2406 samples for which the terminal bud did not burst or approximately 17%

sort(unique(low$species))
# 1+ indiviudal for every species that had a terminal bud

count<-table(low$species)

# For species in both pops there were 128 samples max, so rub par 40% of the time the terminal bud did not bb, for ace gla and sorsco it was 25%

####################################################################################################################
# This code is taken from Dan Flynn and cleaning the east coast data (found in the budchill repo)
# 1. for both terminal and lateral buds, what is the max stage within a row? Identify which rows are greater or equal to the specific BBCH stage
# 2. now for that individual, find the earliest day at which that stage was reached.
####################################################################################################################



tbb <-vector() # This creates empty vectors for the bbday, leaf out, the nl yes no vector of whether this event occured (1 means yes, there was leafout; 0 means no leafout)
gc$lab<-as.factor(gc$lab)
#levels(d$lab)
for(i in levels(gc$lab)){ # i=levels(d$lab)[2496] # for each individual clipping. # DL: why did DF have 602, that seems low
  
  dx <- gc[gc$lab == i,]
  bdax <- which(apply(dx[,c("bbch.t","bbch.l")], MARGIN=1, max, na.rm=T) >3) # margin =1 means it is applied over rows, takes the maximum value > 3
  if(length(bdax) < 1) bdax = NA else bdax = dx[min(bdax),'day']

  tbb<- c(tbb, bdax)

}
dx <- gc[match(levels(gc$lab), gc$lab),] # with twig id in same order as the loop above
#dx <- dx[,2:ncol(dx)] 
dx <- dx[,c(3:8,20)]
terminalbb <- data.frame(dx, tbb)

warnings()

head(terminalbb) # here is the data for the terminal bud, there are 2406 individuals, with a value for each

nobb<-subset(terminalbb, nl==0) # 411 samples didnt bb
table(nobb$chill) # mostly low chill
table(nobb$force) # mostly low force
table(nobb$photo) # mostly low photoperiod

#############################################################
#    Now moving on to the lateral buds 
#############################################################

#What if I would to use a similar approach as above, but start by calculating when the three bbch levels over phase 7 sum to 80%?

d<-gc
#Task is to select bbch.1 7 and above, sum percentages, then get 1st day where percentage above 80%

#1. reshape data so it is in long format 
dlong <- gather(d, key = "bbchL", value = "stage", c(bbch.l, bbch2.l, bbch3.l))
#dlong <- gather(dlong, key = "percentL", value = "l.percent", c(percent.l, percent2.l,percent3.l))

#put relevent percentages with stages - a bit clunky but does the job 
dlong$bbchPercent <- dlong$percent.l
dlong$bbchPercent [dlong$bbchL == "bbch2.l"] <- dlong$percent2.l [dlong$bbchL == "bbch2.l"] 
dlong$bbchPercent [dlong$bbchL == "bbch3.l"] <- dlong$percent3.l [dlong$bbchL == "bbch3.l"] 

head(dlong)
str(dlong)

#Remove na rows
dlong <- dlong[!is.na(dlong$stage),]

#1. Select bbch.1 7 and above
dlong7 <- dlong[dlong$stage >= 7,]

#sum percentages each sample
#sumPercent <- aggregate(dlong7$bbchPercent, by=list(Category=dlong7$lab), FUN=sum)
sumPercent <- aggregate(dlong7$bbchPercent, by=list(Category=dlong7$lab, day=dlong7$day), FUN=sum)

names(sumPercent) <- c("lab","day", "sumPercent")
#names(sumPercent) <- c("lab","day.l.bb", "sumPercent")

dlong7sum <- merge(dlong7, sumPercent, by = c("lab",'day'))
head(dlong7sum)

#Select samples with 80 or more percent --> 1068 rows, if it is lower at 50, then there are 1437 rows
dlong7sum80 <- dlong7sum[dlong7sum$sumPercent >= 80,]

#Select first day with 80% or more
latdaymin80<- aggregate(dlong7sum80$day, by=list(Category=dlong7sum80$lab), FUN=min)
names(latdaymin80) <- c("lab", "latbb80")

#dlong7sum80Daymin <- merge(dlong7sum80, daymin, by = "lab")
length(unique(daymin$lab))

#Select first day with 50%, then there are 1437 rows
dlong7sum50 <- dlong7sum[dlong7sum$sumPercent >= 50,]
latdaymin50<- aggregate(dlong7sum50$day, by=list(Category=dlong7sum50$lab), FUN=min)
names(latdaymin50) <- c("lab", "latbb50")
nrow(latdaymin50)

#Select first day of lateral bb, then there are  rows 19323
dlong7sum1 <- dlong7sum[dlong7sum$sumPercent > 0,]
latdaymin1<- aggregate(dlong7sum1$day, by=list(Category=dlong7sum1$lab), FUN=min)
names(latdaymin1) <- c("lab", "latbb1")
nrow(latdaymin1)


######################################################
# Combine the terminal and the lateral bb days

head(terminalbb)
head(latdaymin80)
head(latdaymin50)
head(latdaymin1)

pheno<-merge(terminalbb, latdaymin80, by= "lab", all.x=TRUE) # the all.x=T is telling it that I want all rows from this dataset, even if there isn't a corresponding row in latdaymin

pheno<-merge(pheno, latdaymin50, by= "lab", all.x=TRUE) 
pheno<-merge(pheno, latdaymin1, by= "lab", all.x=TRUE) 

head(pheno)

######################################################
# To make it more comparable to 
##### GOOO ###########################################

# test<-daymin[order(daymin$firstday),]

#Alternative dplyr solution (sorry Lizzie!)
data.frame(dlong7 %>% 
             filter(stage >=7) %>%
             group_by(lab,day) %>%
             dplyr::mutate(sumPercent = sum(bbchPercent) ) %>%
             filter(sumPercent >= 80) %>%
             filter(day == min(day))
) # you need the dplyr:: before mutate otherwise it doesnt work 

amealnfl10<-subset(data, lab =="mp_HC_HP_HF_10_amealn") # yup first day is day 5

alinc14<-subset(data, lab =="sm_LC_HP_HF_14_alninc") # yup firs day is day 36

betpap17<-subset(data, lab =="sm_HC_LP_LF_17_betpap") # yup firs day is day 36

#Yay this worked!! 