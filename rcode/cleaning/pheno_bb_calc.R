# This code will compile the final data file for use in my phenology analysis

# The aim is to extract the day of study on which the terminal bud reached stage 7 first and the day when the lateral buds reached 80% at stage 7
# building off of Dan Flynn's data

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc/")
} else {
  setwd("~/Documents/github/pheno_bc")
}

require(plyr)
require(dplyr)
require(tidyr)
# require(chillR)

rm(list=ls()) 
options(stringsAsFactors = FALSE)

# read in the cleaning phenology data:
data<-read.csv("input/bc_phenology.csv", header=T, na.strings=c("","NA"))
#source("rcode/cleaning/cleaningcode.R")




# what does the data look like generally?
# 20 species from mp,20 species from mp; so in theory there should be 2560 samples, but after chilling we had 2539
temp<-data[, c("day","population","species","treatment", "lab")]
temp0<-subset(temp, day == "0")
table(temp0$species)
sum(table(temp0$species))

temp88<-subset(temp, day == "88")
table(temp88$species)


#begin by dividing the treatment names (C_P_F) into their own columns
data<-data %>% 
  separate(treatment, c("chill","photo","force"), "_")
head(data)
data$lab<-paste(data$population,data$chill,data$photo,data$force,data$flask, data$species, sep=".")
d<-data
#start by identifying the samples that are duplicates, demarcated with T or F

d$dup<-duplicated(d[,c("day","lab")])

#Check that it worked the way I wanted
#test<-subset(d, dup == "TRUE") # 17131
# there are two flasks that have 3 of the same species in it!
# test<-subset(d, lab =="sm_HC_LP_HF_37_vacmem")
# test<-subset(d, lab =="mp_LC_LP_LF_4_menfer")
head(d)

d<-d %>% 
  group_by(day, lab) %>% 
  mutate(ref=ifelse(dup, "2", "1"))
head(d)
d$lab2<-paste(d$lab, d$ref, sep=".")
d<-as.data.frame(d)
head(d)

d$dup2<-duplicated(d[,c("day","lab2")])
d<-d %>% 
  group_by(day, lab2) %>% 
  mutate(ref2=ifelse(dup2, "3", ""))
head(d)
d$lab3<-paste(d$lab2, d$ref2, sep=".")
d<-as.data.frame(d)
head(d)

# For curiosity, here I am creating a new datset of jus the samples that have pairs of the same species in a flask
# ddups <- vector()
# for(i in 1:length(d$dup)){
#   if(d$dup[i] == "TRUE"){
#     ddups<-rbind(ddups, d[i,])
#   }
# }

head(ddups)

length(unique(d$lab)) 
length(unique(d$lab2))
length(unique(d$lab3)) 

# How many indiv of each sp are there?
d0<-subset(d, day == "0")
table(d0$species)

d50<-subset(d, day == "50")
table(d50$species)

head(d0)

tail(sort(unique(d$lab)))
tail(sort(unique(d$lab2)))

d<-as.data.frame(d);head(d)

# There are a few species that have some extras
goop<-subset(d, day == 50 & species == "riblac")
sort(unique(goop$lab2))


#there should actually be 2539 samples

# Starting with the terminal buds:


##### Adding individual ############
# indiv<-read.csv("input/indiv.no.cleaned.csv", na.strings = "")

# source("rcode/cleaning/cleaning.indivno.R")
# 
# head(indiv)
# #subset to just the colns needed
# indiv<-indiv[,c("lab3","indiv")]
# # indiv[complete.cases(indiv),]
# # test<-subset(indiv, indiv!= "NA")
# head(indiv)
# 
# 
# dtemp<-merge(d, indiv, by= "lab2", all.x=T) # this is adding rows! 

# mptemp<-subset(indiv, site == "mp");length(unique(mptemp$labtemp)) # 204
# 
# smtemp<-subset(indiv, site == "sm");length(unique(smtemp$labtemp)) #186


####################################################################

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
head(gc)
gc<-as.data.frame(gc)
# I am curious if all samples even reached stage 7 or...

max<-gc %>% 
  group_by(lab,ref) %>%
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



tbb <-nl<-vector() # This creates empty vectors for the bbday, leaf out, the nl yes no vector of whether this event occured (1 means yes, there was leafout; 0 means no leafout)
gc$lab<-as.factor(gc$lab)
#levels(d$lab)
for(i in levels(gc$lab)){ # i=levels(d$lab)[2496] # for each individual clipping. # DL: why did DF have 602, that seems low
  
  dx <- gc[gc$lab == i,]
  bdax <- which(apply(dx[,c("bbch.t","bbch.l")], MARGIN=1, max, na.rm=T) >3) # margin =1 means it is applied over rows, takes the maximum value > 3
  if(length(bdax) < 1) {bdax = NA; nl <- c(nl, 0)} else {bdax = dx[min(bdax),'day']; nl <- c(nl, 1)}
  
  tbb<- c(tbb, bdax)

}
dx <- gc[match(levels(gc$lab), gc$lab),] # with twig id in same order as the loop above
#dx <- dx[,2:ncol(dx)] 
dx <- dx[,c("lab", "population", "chill", "photo", "force", "flask", "species"#, "indiv"
            )]
terminalbb <- data.frame(dx, tbb,nl)


head(terminalbb) # here is the data for the terminal bud, there are 2406 individuals, with a value for each

#
nobb<-subset(terminalbb, nl==0) # 209 samples didnt bb
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
length(unique(latdaymin80$lab))

#Select first day with 50%, then there are 1450 rows
dlong7sum50 <- dlong7sum[dlong7sum$sumPercent >= 50,]
latdaymin50<- aggregate(dlong7sum50$day, by=list(Category=dlong7sum50$lab), FUN=min)
names(latdaymin50) <- c("lab", "latbb50")
nrow(latdaymin50)

#Select first day of lateral bb, then there are  rows 1923
dlong7sum1 <- dlong7sum[dlong7sum$sumPercent > 0,]
latdaymin1<- aggregate(dlong7sum1$day, by=list(Category=dlong7sum1$lab), FUN=min)
names(latdaymin1) <- c("lab", "latbb1")
nrow(latdaymin1)


######################################################
# Combine the terminal and the lateral bb days

pheno<-merge(terminalbb, latdaymin80, by= "lab", all.x=TRUE) # the all.x=T is telling it that I want all rows from this dataset, even if there isn't a corresponding row in latdaymin
pheno<-merge(pheno, latdaymin50, by= "lab", all.x=TRUE) 
pheno<-merge(pheno, latdaymin1, by= "lab", all.x=TRUE) 

head(pheno)


######################################################
# To make it more comparable to the Flynn dataset, I am adding a treatment column, and then try to calculate chill portions...for the terminal bud? 

pheno$treatment<-paste(pheno$chill, pheno$photo, pheno$force, sep = "_")
head(pheno)
#Calculating chill portions

# plots
# start by plotting means by species and by treatments:
# phenoplot_sp<-ddply(pheno, .(treatment, species), summarize, termbb=mean(tbb, na.rm=TRUE), lat50bb=mean(latbb50, na.rm=TRUE))
# 
# head(phenoplot_sp)
# 
# ggplot(phenoplot, aes(species, termbb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
# 
# # plotting means by just treatments:
# phenoplot_trt<-ddply(pheno, .(treatment), summarize, termbb=mean(tbb, na.rm=TRUE), lat50bb=mean(latbb50, na.rm=TRUE))
# 
# #terminal dobb
# ggplot(phenoplot_sp, aes(treatment, termbb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
#   scale_y_continuous(limits = c(0, 88))
# 
# #lateral dobb
# ggplot(phenoplot_sp, aes(treatment, lat50bb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+ 
#   scale_y_continuous(limits = c(0, 88))
# 
# # Just the mean treatment value
# # terminla
# ggplot(phenoplot_trt, aes(treatment, termbb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+ 
#   scale_y_continuous(limits = c(0, 88))
# 
# ggplot(phenoplot_trt, aes(treatment, lat50bb)) +
#   geom_point(aes(color = treatment)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+ 
#   scale_y_continuous(limits = c(0, 88))
# 
# 
#  hhh<-subset(pheno, treatment == "HC_HP_HF")
# # hll<-subset(pheno, treatment == "HC_LP_LF")
# # hlh<-subset(pheno, treatment == "HC_LP_HF")
# # hhl<-subset(pheno, treatment == "HC_HP_LF")
# # lll<-subset(pheno, treatment == "LC_LP_LF")
# # lhh<-subset(pheno, treatment == "LC_HP_HF")
# # llh<-subset(pheno, treatment == "LC_LP_HF")
# # lhl<-subset(pheno, treatment == "LC_HP_LF")
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
# d<-pheno
# 
# pdf( width = 8, height = 4)
# 
# par(mfcol=c(1, 3), mar = c(3,3,1,0.5))
# for(spx in levels(d$species)){ # spx = "BETALL"
#   
#   dxx = d[d$species == spx,]
#   
#   counter = 1
#   for(i in sort(as.character((unique(d$chill))))){# i = "high chill or low chill"
#     
#     dseq = seq(0, max(dxx$day))
#     plot(dseq, seq(0, 7,length=length(dseq)), type = "n", 
#          ylab = "Stage",
#          xlab = "")
#     if(counter == 1) mtext(spx, line = -2, adj = 0.5)
#     legend("topleft",bty="n",i, cex = 0.85, inset = 0)
#     xx <- dxx[dxx$time == i,]
#     # calculate mean response by date and chill
#     xt <- tapply(pmax(xx$tleaf, xx$lleaf,na.rm=T), list(xx$day, xx$chill), mean, na.rm=T)
#     
#     for(j in unique(xx$ind)){ #j=unique(xx$ind)[1]
#       xj <- xx[xx$ind == j,]
#       pcol = lcol[xj$chill]
#       lines(xj$day, xj$tleaf, col = pcol)
#     }
#     lines(rownames(xt), xt[,1], col = colz[1], lwd = 2)
#     lines(rownames(xt), xt[,2], col = colz[2], lwd = 2)
#     lines(rownames(xt), xt[,3], col = colz[3], lwd = 2)
#     lines(rownames(xt), xt[,4], col = colz[4], lwd = 2)
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
#     
#     counter = counter + 1
#   }
#   
# }
# dev.off()
# system(paste("open '", paste("figures/Trace Plots ", Sys.Date(), ".pdf", sep=""), "' -a /Applications/Preview.app", sep=""))
# 

#write.csv(pheno, "input/day.of.bb.csv")
#

##### GOOO ###########################################


##### GOOO ###########################################

# test<-daymin[order(daymin$firstday),]

#Alternative dplyr solution (sorry Lizzie!)
# data.frame(dlong7 %>% 
#              filter(stage >=7) %>%
#              group_by(lab,day) %>%
#              dplyr::mutate(sumPercent = sum(bbchPercent) ) %>%
#              filter(sumPercent >= 80) %>%
#              filter(day == min(day))
# ) # you need the dplyr:: before mutate otherwise it doesnt work 
# 
# amealnfl10<-subset(data, lab =="mp_HC_HP_HF_10_amealn") # yup first day is day 5
# 
# alinc14<-subset(data, lab =="sm_LC_HP_HF_14_alninc") # yup firs day is day 36
# 
# betpap17<-subset(data, lab =="sm_HC_LP_LF_17_betpap") # yup firs day is day 36

#Yay this worked!! 