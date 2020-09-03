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

# I am curious if all samples even reached stage 7 or...

max<-pheno %>% 
  group_by(lab) %>%
  slice(which.max(bbch.t))

low<-subset(max, bbch.t<7)
# there are 440 samples for which the terminal bud did not burst
sort(unique(low$species))
# 1+ indiviudal for every species that had a terminal bud

count<-table(low$species)

#acegla alninc alnvir amealn  bepap betpap  corso corsto loninv popbal poptre rhoalb riblac rubpar samrac shecan sorsco 
#33      1      7      6     48      1     23      2     37      4     14      1     39     52      8      4     33 
#spibet spipyr vacmem vibedu 
#23     38     59      7 
# For species in both pops there were 128 samples max, so rub par 40% of the time the terminal bud did not bb, for ace gla and sorsco it was 25%

# Now need to subset the day to get doy when bb occured
bday <- lday <- nl <- vector()
d$lab<-as.factor(d$lab)
levels(d$lab)
for(i in levels(d$lab)){ # i=levels(d$lab)[2496] # for each individual clipping. # DL: why did DF have 602, that seems low
  
  dx <- d[d$lab == i,]
  
  # 1. for both terminal and lateral buds, what is the max stage within a row? Identify which rows are greater or equal to the specific BBCH stage
  # 2. now for that individual, find the earliest day at which that stage was reached.
  # Lizzie understand 'nl' to mean a yes/no (1/0) on whether an individual twig made it to leafout (1 means yes, there was leafout; 0 means no leafout)
  bdax <- which(apply(dx[,c("bbch.t","bbch.l")], 1, max, na.rm=T) >= 3)
  if(length(bdax) < 1) bdax = NA else bdax = dx[min(bdax),'day']
  
  ldax <- which(apply(dx[,c("bbch.t","bbch.l")], 1, max, na.rm=T) >= 6)
  if(length(ldax) < 1) {ldax = NA; nl <- c(nl, 0)} else {ldax = dx[min(ldax),'day']; nl <- c(nl, 1)}
  
  
  bday <- c(bday, bdax)
  lday <- c(lday, ldax)
}
dx <- d[match(levels(d$lab), d$lab),] # with twig id in same order as the loop above
dx <- dx[,2:ncol(dx)]
dx <- data.frame(dx, lday, bday, nl)
head(dx)
warnings()
