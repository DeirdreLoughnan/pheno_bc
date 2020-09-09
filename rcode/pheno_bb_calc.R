# This code will compile the final data file for use in my phenology analysis

# The aim is to extract the day of study on which the terminal bud reached stage 7 first and the day when the lateral buds reached 80% at stage 7
# building off of Dan Flynn's data

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

d<-data %>% 
  separate(treatment, c("chill","photo","force"), "_")

head(d)

# I have 21 species, but rhoalb, betpap, samrac were only at one site
18*8*8*2+3*8*8 # there should be 2496 samples, but I actually have more...2546

begin<-subset(d, day==0)
count.begin<-table(begin$species)

end<-subset(d, day==88)
count.end<-table(end$species)

2496-2406 # only had 90 die! 
count<-table(d$species)
# Should look into which ones these are...

d$day<-as.numeric(d$day)
range(d$day, na.rm=TRUE)
sort(unique(data$species))

end<-subset(d, day==88)
count<-table(end$species)
sum(count)
samrac<-subset(end, species=="samrac")

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
# there are 440 samples for which the terminal bud did not burst
sort(unique(low$species))
# 1+ indiviudal for every species that had a terminal bud

count<-table(low$species)

# For species in both pops there were 128 samples max, so rub par 40% of the time the terminal bud did not bb, for ace gla and sorsco it was 25%

####################################################################################################################
# This code is taken from Dan Flynn and cleaning the east coast data (found in the budchill repo)
# 1. for both terminal and lateral buds, what is the max stage within a row? Identify which rows are greater or equal to the specific BBCH stage
# 2. now for that individual, find the earliest day at which that stage was reached.
####################################################################################################################

bday <- lday <- nl <- vector() # This creates empty vectors for the bbday, leaf out, the nl yes no vector of whether this event occured (1 means yes, there was leafout; 0 means no leafout)
d$lab<-as.factor(gc$lab)
#levels(d$lab)
for(i in levels(gc$lab)){ # i=levels(d$lab)[2496] # for each individual clipping. # DL: why did DF have 602, that seems low
  
  dx <- gc[gc$lab == i,]
  bdax <- which(apply(dx[,c("bbch.t","bbch.l")], 1, max, na.rm=T) >= 3) # for each unique identifier, is the bbch >=3?
  if(length(bdax) < 1) bdax = NA else bdax = dx[min(bdax),'day']
  
  ldax <- which(apply(dx[,c("bbch.t","bbch.l")], 1, max, na.rm=T) >= 6)
  if(length(ldax) < 1) {ldax = NA; nl <- c(nl, 0)} else {ldax = dx[min(ldax),'day']; nl <- c(nl, 1)}
  
  
  bday <- c(bday, bdax)
  lday <- c(lday, ldax)
}
dx <- gc[match(levels(gc$lab), gc$lab),] # with twig id in same order as the loop above
dx <- dx[,2:ncol(dx)]
dx <- data.frame(dx, lday, bday, nl)
head(dx)
warnings()

####################################################################################################################

head(gc)
#What if I would to use a similar approach as above, but selecting for values that are greater than the requried sum?
gc$bbch.l.sum<-rowSums(gc[,c("bbch.l","bbch2.l","bbch3.l")], na.rm=TRUE)
gc$percent.l.sum<-rowSums(gc[,c("percent.l","percent2.l","percent3.l")], na.rm=TRUE)
tail(gc)
  
bday <- lday <- nl <- vector() # This creates empty vectors for the bbday, leaf out, the nl yes no vector of whether this event occured (1 means yes, there was leafout; 0 means no leafout)
d$lab<-as.factor(gc$lab)
#levels(d$lab)
for(i in levels(gc$lab)){ # i=levels(d$lab)[2496] # for each individual clipping. # DL: why did DF have 602, that seems low
  
  dx <- gc[gc$lab == i,]
  bdax <- which(apply(dx[,c("bbch.t","bbch.l")], 1, max, na.rm=T) >= 3) # for each unique identifier, is the bbch >=3?
  if(length(bdax) < 1) bdax = NA else bdax = dx[min(bdax),'day']
  
  ldax <- which(apply(dx[,c("bbch.t","bbch.l")], 1, max, na.rm=T) >= 6)
  if(length(ldax) < 1) {ldax = NA; nl <- c(nl, 0)} else {ldax = dx[min(ldax),'day']; nl <- c(nl, 1)}
  
  
  bday <- c(bday, bdax)
  lday <- c(lday, ldax)
}
dx <- gc[match(levels(gc$lab), gc$lab),] # with twig id in same order as the loop above
dx <- dx[,2:ncol(dx)]
dx <- data.frame(dx, lday, bday, nl)
head(dx)
