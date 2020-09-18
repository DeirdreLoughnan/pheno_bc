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



tbb <- tlf <- nl <- vector() # This creates empty vectors for the bbday, leaf out, the nl yes no vector of whether this event occured (1 means yes, there was leafout; 0 means no leafout)
gc$lab<-as.factor(gc$lab)
#levels(d$lab)
for(i in levels(gc$lab)){ # i=levels(d$lab)[2496] # for each individual clipping. # DL: why did DF have 602, that seems low
  
  dx <- gc[gc$lab == i,]
  bdax <- which(apply(dx[,c("bbch.t","bbch.l")], MARGIN=1, max, na.rm=T) >= 3) # margin =1 means it is applied over rows, takes the maximum value >= 3
  if(length(bdax) < 1) bdax = NA else bdax = dx[min(bdax),'day']
  
  ldax <- which(apply(dx[,c("bbch.t","bbch.l")], 1, max, na.rm=T) >= 6)
  if(length(ldax) < 1) {ldax = NA; nl <- c(nl, 0)} else {ldax = dx[min(ldax),'day']; nl <- c(nl, 1)}
  
  
  tbb<- c(tbb, bdax)
  tlf <- c(tlf, ldax)
}
dx <- gc[match(levels(gc$lab), gc$lab),] # with twig id in same order as the loop above
#dx <- dx[,2:ncol(dx)] 
dx <- dx[,c(1:6,19)]
terminalbb <- data.frame(dx, tlf, tbb, nl)

warnings()


#############################################################
#    Now moving on to the lateral buds 
#############################################################

#What if I would to use a similar approach as above, but start by calculating when the three bbch levels over phase 7 sum to 80%?

head(gc)

# Idea 1: using some sort of if else, if bbch.l, bbch2.l, bbch3.l are greater than 7 then add them together, else give it a value of 0 or NA
#...not going well
bbch.l.sum<-vector
for ( i in nrow(gc)){
  if(gc$bbch.l[i] > 6) { 
  temp<-gc$bbch.l[i]
  #bbch.l.sum =gc$percent.l
}  else {
  temp<-0
 }
bbch.l.sum<-rbind(temp,bbch.l.sum)  
}

# idea 2, should I be using a while statement instead...
while(gc[,c("bbch.l","bbch2.l")] > 7){
  gc$test<-rowSum(gc[,c("bbch.l","bbch2.l")], na.rm=TRUE)
}
gc$bbch.l.sum<-rowSum(gc[,c("bbch.l","bbch2.l","bbch3.l")], na.rm=TRUE)
gc$percent.l.sum<-rowSums(gc[,c("percent.l","percent2.l","percent3.l")], na.rm=TRUE)

  
 llf <- nll <- vector() # This creates empty vectors for the bbday, leaf out, the nl yes no vector of whether this event occured (1 means yes, there was leafout; 0 means no leafout)
gc$lab<-as.factor(gc$lab)
#levels(d$lab)
for(i in levels(gc$lab)){ # i=levels(d$lab)[2496] # for each individual clipping. # DL: why did DF have 602, that seems low
  
  dx <- gc[gc$lab == i,]
#  bdaxlat <- which(apply(dx[,c("bbch.l","bbch2.l","bbch3.l")], MARGIN=1, max, na.rm=TRUE) >= 7) # for each unique identifier, is the bbch >=3?
#  if(length(bdax) < 1) bdax = NA else bdax = dx[min(bdax),'day']

  ldax <- which(apply(dx[,c("bbch.l","bbch2.l","bbch3.l")], 1, max, na.rm=T) >= 7 & dx$percent.l.sum >=80)
  if(length(ldax) < 1) {ldax = NA; nll <- c(nll, 0)} else {ldax = dx[min(ldax),'day']; nll <- c(nll, 1)}
  
  llf <- c(llf, ldax)

}

dxl <- gc[match(levels(gc$lab), gc$lab),] # with twig id in same order as the loop above
#dx <- dx[,2:ncol(dx)] 
dxl <- dxl[,c(1:6,19)]
lateralbb <- data.frame(dxl, llf, nl)

head(lateralbb)

#######################################################################################
#idea 3: subset the data so you get only rows that have either 2 or 3 of the phases greater than phase 7 and then sum them to see how many samples even reach 80%

# I dont want to use this method, it is so hack and I don't think it gives the right answer

first<-subset(gc, bbch.l>=7)
second<-subset(first, bbch2.l>=7)

second$percent.l.sum<-rowSums(second[,c("percent.l","percent2.l")], na.rm=TRUE)
second80<-subset(second, percent.l.sum>=80)

third<-subset(second, bbch3.l>=7)

third$percent.l.sum<-rowSums(third[,c("percent.l","percent2.l","percent3.l")], na.rm=TRUE)
third80<-subset(third, percent.l.sum>=80)

fin80<-rbind(second80,third80)

############################################################
fin80$lab<-as.factor(fin80$lab)

llf <- nll <- vector()
for(i in levels(fin80$lab)){ # i=levels(d$lab)[2496] # for each individual clipping. # DL: why did DF have 602, that seems low
  
  dx <- fin80[fin80$lab == i,]
 
  ldax <- which( dx$percent.l.sum >=80)
  if(length(ldax) < 1) {ldax = NA; nll <- c(nll, 0)} else {ldax = dx[min(ldax),'day']; nll <- c(nll, 1)}
  
  llf <- c(llf, ldax)
  
}

dxl <- fin80[match(levels(fin80$lab), fin80$lab),] # with twig id in same order as the loop above
#dx <- dx[,2:ncol(dx)] 
dxl <- dxl[,c(1:6,19)]
lateralbb <- data.frame(dxl, llf, nl)

latbb<-subset(lateralbb, day>=0)
unique(latbb$lab)


## <><><><><><><><><><><><><><><><><><><><><>
#Code sent from Faith on September 17, 2020, hugely helpful! I modified it slightly by adding grouping by day but overall it was perfect! 
#-------------------------------------

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
#Select samples with 80 or more percent 
dlong7sum80 <- dlong7sum[dlong7sum$sumPercent >= 80,]

#Select first day with 80% or more
daymin<- aggregate(dlong7sum80$day, by=list(Category=dlong7sum80$lab), FUN=min)
names(daymin) <- c("lab", "firstday")
dlong7sum80Daymin <- merge(dlong7sum80, daymin, by = "lab")

nrow(daymin) # 1101
length(unique(daymin$lab))

test<-daymin[order(daymin$firstday),]

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