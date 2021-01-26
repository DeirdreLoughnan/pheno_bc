# checking for any remaining outliers or abnormalities in the data


if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc/")
} else {
  setwd("~/Documents/github/pheno_bc")
}

require(plyr)
require(dplyr)
require(tidyr)
require(ggplot2)
# require(chillR)

rm(list=ls()) 
options(stringsAsFactors = FALSE)

# read in the cleaning phenology data:
data<-read.csv("input/bc_phenology.csv", header=T, na.strings=c("","NA"))
#source("rcode/cleaning/cleaningcode.R")
head(data)

d<-data %>% 
  separate(treatment, c("chill","photo","force"), "_")
d$treatment<-paste(d$chill, d$photo, d$force, sep = ".")
head(d)

#start by identifying the samples that are duplicates, demarcated with T or F

d$dup<-duplicated(d[,c("day","lab")])

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

length(unique(d$lab2))
length(unique(d$lab3))

# now I want to create plots so that I have 168 plots (one for each species) for each treatment and lines within each for each sample

## LC LP LF
lll<-subset(d, treatment =="LC.LP.LF")
for (i in 1:length(unique(lll$species))){
  sp<-sort(unique(lll$species))
sub<-subset(lll, species== sp[i])

spplot<-ggplot(sub)+
  aes(x=day, y=bbch.t, colour = lab3, linetype = population)+
  geom_line() +
  labs(x="day",y="bbch") +
  xlim(-0.5,113)+
  ylim(0,20)
ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/terminal/LLL",i,".pdf"))
}

## LC HP LF
lhl<-subset(d, treatment =="LC.HP.LF")
for (i in 1:length(unique(lhl$species))){
  sp<-sort(unique(lhl$species))
  sub<-subset(lhl, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.t, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/terminal/LHL",i,".pdf"))
}

## LC HP HF
lhh<-subset(d, treatment =="LC.HP.HF")
for (i in 1:length(unique(lhh$species))){
  sp<-sort(unique(lhh$species))
  sub<-subset(lhh, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.t, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/terminal/LHH",i,".pdf"))
}

## LC LP HF
llh<-subset(d, treatment =="LC.LP.HF")
for (i in 1:length(unique(llh$species))){
  sp<-sort(unique(llh$species))
  sub<-subset(llh, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.t, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/terminal/LLH",i,".pdf"))
}

############################################################################################
## HC HP HF
hll<-subset(d, treatment =="HC.LP.LF")
for (i in 1:length(unique(hll$species))){
  sp<-sort(unique(hll$species))
  sub<-subset(hll, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.t, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/terminal/HLL",i,".pdf"))
}

## HC HP LF
hhl<-subset(d, treatment =="HC.HP.LF")
for (i in 1:length(unique(hhl$species))){
  sp<-sort(unique(hhl$species))
  sub<-subset(hhl, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.t, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/terminal/HHL",i,".pdf"))
}

## HC HP HF
hhh<-subset(d, treatment =="HC.HP.HF")
for (i in 1:length(unique(hhh$species))){
  sp<-sort(unique(hhh$species))
  sub<-subset(hhh, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.t, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/terminal/HHH",i,".pdf"))
}

## HLH
hlh<-subset(d, treatment =="HC.LP.HF")
for (i in 1:length(unique(hlh$species))){
  sp<-sort(unique(hlh$species))
  sub<-subset(hlh, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.t, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/terminal/HLH",i,".pdf"))
}

############################################################################################
############################################################################################
# Lateral buds
## LC LP LF
lll<-subset(d, treatment =="LC.LP.LF")
for (i in 1:length(unique(lll$species))){
  sp<-sort(unique(lll$species))
  sub<-subset(lll, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.l, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/LLL",i,".pdf"))
}

## LC HP LF
lhl<-subset(d, treatment =="LC.HP.LF")
for (i in 1:length(unique(lhl$species))){
  sp<-sort(unique(lhl$species))
  sub<-subset(lhl, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.l, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/LHL",i,".pdf"))
}

## LC HP HF
lhh<-subset(d, treatment =="LC.HP.HF")
for (i in 1:length(unique(lhh$species))){
  sp<-sort(unique(lhh$species))
  sub<-subset(lhh, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.l, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/LHH",i,".pdf"))
}

## LC LP HF
llh<-subset(d, treatment =="LC.LP.HF")
for (i in 1:length(unique(llh$species))){
  sp<-sort(unique(llh$species))
  sub<-subset(llh, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.l, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/LLH",i,".pdf"))
}

############################################################################################
## HC HP HF
hll<-subset(d, treatment =="HC.LP.LF")
for (i in 1:length(unique(hll$species))){
  sp<-sort(unique(hll$species))
  sub<-subset(hll, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.l, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/HLL",i,".pdf"))
}

## HC HP LF
hhl<-subset(d, treatment =="HC.HP.LF")
for (i in 1:length(unique(hhl$species))){
  sp<-sort(unique(hhl$species))
  sub<-subset(hhl, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.l, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/HHL",i,".pdf"))
}

## HC HP HF
hhh<-subset(d, treatment =="HC.HP.HF")
for (i in 1:length(unique(hhh$species))){
  sp<-sort(unique(hhh$species))
  sub<-subset(hhh, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.l, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/HHH",i,".pdf"))
}

## HLH
hlh<-subset(d, treatment =="HC.LP.HF")
for (i in 1:length(unique(hlh$species))){
  sp<-sort(unique(hlh$species))
  sub<-subset(hlh, species== sp[i])
  
  spplot<-ggplot(sub)+
    aes(x=day, y=bbch.l, colour = lab3, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,20)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/HLH",i,".pdf"))
}

## To check the lateral buds, I think it would be simplest to plot the sumPercent, ie the total % of lateral bb that is greater than 7, it should show a coninual increase

source('rcode/cleaning/pheno_bb_calc.R')

head(dlong7sum)
dlong7sum$sumPercent<-as.numeric(dlong7sum$sumPercent)
dlong7sum$lab2<-as.character(dlong7sum$lab2)

hhh<-subset(dlong7sum, treatment =="HC.HP.HF")
for (i in 1:length(unique(hhh$species))){
  sp<-sort(unique(hhh$species))
  sub<-subset(hhh, species== sp[i])
  
  spplot<-
    ggplot(sub)+
    aes(x=day, y=sumPercent, colour = lab2, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,100)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/hhhlat",i,".pdf"))
}

hlh<-subset(dlong7sum, treatment =="HC.LP.HF")
for (i in 1:length(unique(hlh$species))){
    sp<-sort(unique(hlh$species))
    sub<-subset(hlh, species== sp[i])

    spplot<-
      ggplot(sub)+
      aes(x=day, y=sumPercent, colour = lab2, linetype = population)+
      geom_line() +
      labs(x="day",y="bbch") +
      xlim(-0.5,113)+
      ylim(0,100)
    ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/hlhlat",i,".pdf"))
}

hll<-subset(dlong7sum, treatment =="HC.LP.LF")
for (i in 1:length(unique(hll$species))){
  sp<-sort(unique(hll$species))
  sub<-subset(hll, species== sp[i])
  
  spplot<-
    ggplot(sub)+
    aes(x=day, y=sumPercent, colour = lab2, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,100)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/hlllat",i,".pdf"))
}

hhl<-subset(dlong7sum, treatment =="HC.HP.LF")
for (i in 1:length(unique(hhl$species))){
  sp<-sort(unique(hhl$species))
  sub<-subset(hhl, species== sp[i])
  
  spplot<-
    ggplot(sub)+
    aes(x=day, y=sumPercent, colour = lab2, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,100)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/hhllat",i,".pdf"))
}

######################################################################################
lhh<-subset(dlong7sum, treatment =="LC.HP.HF")
for (i in 1:length(unique(lhh$species))){
  sp<-sort(unique(lhh$species))
  sub<-subset(lhh, species== sp[i])
  
  spplot<-
    ggplot(sub)+
    aes(x=day, y=sumPercent, colour = lab2, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,100)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/lhhlat",i,".pdf"))
}

llh<-subset(dlong7sum, treatment =="LC.LP.HF")
for (i in 1:length(unique(llh$species))){
  sp<-sort(unique(llh$species))
  sub<-subset(llh, species== sp[i])
  
  spplot<-
    ggplot(sub)+
    aes(x=day, y=sumPercent, colour = lab2, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,100)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/llhlat",i,".pdf"))
}

lll<-subset(dlong7sum, treatment =="LC.LP.LF")
for (i in 1:length(unique(lll$species))){
  sp<-sort(unique(lll$species))
  sub<-subset(lll, species== sp[i])
  
  spplot<-
    ggplot(sub)+
    aes(x=day, y=sumPercent, colour = lab2, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,100)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/llllat",i,".pdf"))
}

lhl<-subset(dlong7sum, treatment =="LC.HP.LF")
for (i in 1:length(unique(lhl$species))){
  sp<-sort(unique(lhl$species))
  sub<-subset(lhl, species== sp[i])
  
  spplot<-
    ggplot(sub)+
    aes(x=day, y=sumPercent, colour = lab2, linetype = population)+
    geom_line() +
    labs(x="day",y="bbch") +
    xlim(-0.5,113)+
    ylim(0,100)
  ggsave(spplot, file=paste("rcode/cleaning/plots_rawdat/lateral/lhllat",i,".pdf"))
}