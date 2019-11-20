#Oct 10, 2019

# Need to randomly assign branches to different flasks, aiming for either 3 or 4 per flask

# 20 species, 2 chilling, 2 forcing, 2 photoperiod; 8 individuals per f/p/c treatment
# We have 64 per species
rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
library(dplyr)
library(tidyr)
library(splitstackshape)
library(reshape)
library(randomizr)

rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("~/Documents/github/pheno_bc")

splist<-read.csv("GC_splist_smithers.csv")

splist<-expandRows(splist, "collected")
by_cyl<-splist %>% group_by(sp.name)
set.seed(343)
exper<-sample_n(by_cyl, 64, replace=FALSE)

#assign to treatments:
tmp<-block_ra(blocks=exper$sp.name, conditions = c("HC_HP_HF","HC_LP_HF","HC_HP_LF","HC_LP_LF","LC_HP_HF","LC_LP_HF","LC_LP_LF","LC_HP_LF"))
exper$assignment<- tmp
exper

#Double checking: do I get 64 rows, 
ag<-filter(exper, sp.name=="acegla")
nrow(ag)
table(ag$assignment)

li<-filter(exper, sp.name=="loninv")
nrow(li)
table(li$assignment)

#exper$sp.trt<-paste(exper$sp.name, exper$assignment, sep ="_")

#HC_HF_HP
unique(exper$assignment)
hchphf<-subset(exper, assignment=="HC_HP_HF")
#Assign to different beakers
X<-complete_ra(N=152, num_arms = 40)
hchphf$flask<-X
hchphf$flask_id<-paste(hchphf$sp.name, hchphf$assignment, hchphf$flask, sep="_")
write.csv(hchphf,"hchphf.csv",row.names = FALSE)

#HC_HF_LP
unique(exper$assignment)
hchplf<-subset(exper, assignment=="HC_HP_LF")
#Assign to different beakers
X<-complete_ra(N=152, num_arms = 40)
hchplf$flask<-X
hchplf$flask_id<-paste(hchplf$sp.name, hchplf$assignment, hchplf$flask, sep="_")
write.csv(hchplf,"hchplf.csv",row.names = FALSE)

#HC_LP_HF
unique(exper$assignment)
hclphf<-subset(exper, assignment=="HC_LP_HF")
#Assign to different beakers
X<-complete_ra(N=152, num_arms = 40)
hclphf$flask<-X
hclphf$flask_id<-paste(hclphf$sp.name, hclphf$assignment, hclphf$flask, sep="_")
write.csv(hclphf,"hclphf.csv",row.names = FALSE)

#HC_LP_LF
unique(exper$assignment)
hclplf<-subset(exper, assignment=="HC_LP_LF")
#Assign to different beakers
X<-complete_ra(N=152, num_arms = 40)
hclplf$flask<-X
hclplf$flask_id<-paste(hclplf$sp.name, hclplf$assignment, hclplf$flask, sep="_")
write.csv(hclplf,"hclplf.csv",row.names = FALSE)


###############
#LC_HF_HP
unique(exper$assignment)
unique(exper$assignment)
LChphf<-subset(exper, assignment=="LC_HP_HF")
#Assign to different beakers
X<-complete_ra(N=152, num_arms = 40)
LChphf$flask<-X
LChphf$flask_id<-paste(LChphf$sp.name, LChphf$assignment, LChphf$flask, sep="_")
write.csv(LChphf,"lchphf.csv",row.names = FALSE)

#LC_HF_LP
unique(exper$assignment)
LChplf<-subset(exper, assignment=="LC_HP_LF")
#Assign to different beakers
X<-complete_ra(N=152, num_arms = 40)
LChplf$flask<-X
LChplf$flask_id<-paste(LChplf$sp.name, LChplf$assignment, LChplf$flask, sep="_")
write.csv(LChplf,"lchplf.csv",row.names = FALSE)

#LC_LP_HF
unique(exper$assignment)
LClphf<-subset(exper, assignment=="LC_LP_HF")
#Assign to different beakers
X<-complete_ra(N=152, num_arms = 40)
LClphf$flask<-X
LClphf$flask_id<-paste(LClphf$sp.name, LClphf$assignment, LClphf$flask, sep="_")
write.csv(LClphf,"lclphf.csv",row.names = FALSE)

#LC_LP_LF
unique(exper$assignment)
LClplf<-subset(exper, assignment=="LC_LP_LF")
#Assign to different beakers
X<-complete_ra(N=152, num_arms = 40)
LClplf$flask<-X
LClplf$flask_id<-paste(LClplf$sp.name, LClplf$assignment, LClplf$flask, sep="_")

write.csv(LClplf,"lclplf.csv",row.names = FALSE)


tags1<-dplyr::select(lclplf, sp.name,assignment, flask,flask_id)
tags2<-dplyr::select(lchplf, sp.name,assignment, flask,flask_id)
tags3<-dplyr::select(lclphf, sp.name,assignment, flask,flask_id)
tags4<-dplyr::select(lchphf, sp.name,assignment, flask,flask_id)
tags5<-dplyr::select(hclplf, sp.name,assignment, flask,flask_id)
tags6<-dplyr::select(hchplf, sp.name,assignment, flask,flask_id)
tags7<-dplyr::select(hclphf, sp.name,assignment, flask,flask_id)
tags8<-dplyr::select(hchphf, sp.name,assignment, flask,flask_id)

tag<-rbind(tags1,tags2,tags3,tags4,tags5,tags6,tags7,tags8)
names(tag)
temp2<-tag[order(flask),]


unique(temp$flask)
write.csv(tag,"final_labels_smithers.csv",row.names = FALSE)
head(tag)
require(reshape2)


#Sorting count for each species so can priorities indiv species and make set up faster?
d<-read.csv("tags_smithers_printing.csv")
head(d)
length(d)


test<-dcast(d, assignment+flask~sp.name, value.var="sp.name")
write.csv(test,"ordering.csv",row.names = FALSE)

###<><><><><><><><><><><><><><><><><><><><>#####

splist<-read.csv("GC_splist_manning.csv")

splist<-expandRows(splist, "collected")
by_cyl<-splist %>% group_by(sp.name)
set.seed(343)
exper<-sample_n(by_cyl, 64, replace=FALSE)

#assign to treatments:
tmp<-block_ra(blocks=exper$sp.name, conditions = c("HC_HP_HF","HC_LP_HF","HC_HP_LF","HC_LP_LF","LC_HP_HF","LC_LP_HF","LC_LP_LF","LC_HP_LF"))
exper$assignment<- tmp
exper

#Double checking: do I get 64 rows, 
ag<-filter(exper, sp.name=="acegla")
nrow(ag)
table(ag$assignment)

li<-filter(exper, sp.name=="loninv")
nrow(li)
table(li$assignment)

#exper$sp.trt<-paste(exper$sp.name, exper$assignment, sep ="_")

#HC_HF_HP
unique(exper$assignment)
hchphf<-subset(exper, assignment=="HC_HP_HF")
#Assign to different beakers
X<-complete_ra(N=160, num_arms = 40)
hchphf$flask<-X
hchphf$flask_id<-paste(hchphf$sp.name, hchphf$assignment, hchphf$flask, sep="_")
write.csv(hchphf,"hchphf_mp.csv",row.names = FALSE)

#HC_HF_LP
unique(exper$assignment)
hchplf<-subset(exper, assignment=="HC_HP_LF")
#Assign to different beakers
X<-complete_ra(N=160, num_arms = 40)
hchplf$flask<-X
hchplf$flask_id<-paste(hchplf$sp.name, hchplf$assignment, hchplf$flask, sep="_")
write.csv(hchplf,"hchplf_mp.csv",row.names = FALSE)

#HC_LP_HF
unique(exper$assignment)
hclphf<-subset(exper, assignment=="HC_LP_HF")
#Assign to different beakers
X<-complete_ra(N=160, num_arms = 40)
hclphf$flask<-X
hclphf$flask_id<-paste(hclphf$sp.name, hclphf$assignment, hclphf$flask, sep="_")
write.csv(hclphf,"hclphf_mp.csv",row.names = FALSE)

#HC_LP_LF
unique(exper$assignment)
hclplf<-subset(exper, assignment=="HC_LP_LF")
#Assign to different beakers
X<-complete_ra(N=160, num_arms = 40)
hclplf$flask<-X
hclplf$flask_id<-paste(hclplf$sp.name, hclplf$assignment, hclplf$flask, sep="_")
write.csv(hclplf,"hclplf_mp.csv",row.names = FALSE)


###############
#LC_HF_HP
unique(exper$assignment)
unique(exper$assignment)
LChphf<-subset(exper, assignment=="LC_HP_HF")
#Assign to different beakers
X<-complete_ra(N=160, num_arms = 40)
LChphf$flask<-X
LChphf$flask_id<-paste(LChphf$sp.name, LChphf$assignment, LChphf$flask, sep="_")
write.csv(LChphf,"lchphf_mp.csv",row.names = FALSE)

#LC_HF_LP
unique(exper$assignment)
LChplf<-subset(exper, assignment=="LC_HP_LF")
#Assign to different beakers
X<-complete_ra(N=160, num_arms = 40)
LChplf$flask<-X
LChplf$flask_id<-paste(LChplf$sp.name, LChplf$assignment, LChplf$flask, sep="_")
write.csv(LChplf,"lchplf_mp.csv",row.names = FALSE)

#LC_LP_HF
unique(exper$assignment)
LClphf<-subset(exper, assignment=="LC_LP_HF")
#Assign to different beakers
X<-complete_ra(N=160, num_arms = 40)
LClphf$flask<-X
LClphf$flask_id<-paste(LClphf$sp.name, LClphf$assignment, LClphf$flask, sep="_")
write.csv(LClphf,"lclphf_mp.csv",row.names = FALSE)

#LC_LP_LF
unique(exper$assignment)
LClplf<-subset(exper, assignment=="LC_LP_LF")
#Assign to different beakers
X<-complete_ra(N=160, num_arms = 40)
LClplf$flask<-X
LClplf$flask_id<-paste(LClplf$sp.name, LClplf$assignment, LClplf$flask, sep="_")

write.csv(LClplf,"lclplf_mp.csv",row.names = FALSE)


#Sorting count for each species so can priorities indiv species and make set up faster?
d<-read.csv("tags_manning_printing.csv")
head(d)
length(d)


test<-dcast(d, assignment+flask~sp.name, value.var="sp.name")
write.csv(test,"ordering.csv",row.names = FALSE)
