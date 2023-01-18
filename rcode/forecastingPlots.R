#started January 11, 2022 by Deirdre

#Aim of this code is to make prjection plots for individual species for the growth chamber study
#based off of forecast_changebb.R from ospree


rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(rstan)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  

load("output/bb_4sites_phylo_contin.Rda")

fit <- mdl.4phylo
fit.sum <- summary(fit)$summary

sp<-c("alninc","betpap")
sp.num<-c(5,11) # sp fact no.
tempforecast<-c(1,2,3,4,5,6,7)

###################################################################
# get the climate data:
chillsm<-read.csv("input/smithersChilling.csv")
chillmp<-read.csv("input/manningparkChilling.csv")

# just do both from 1975
tempsm<-read.csv("input/smithersDaily_1943_2018.csv"); tempsm <- tempsm[tempsm$year>1974,]
tempmp<-read.csv("input/hopeDaily_1975_2018.csv")


tempsm$tempMean <- as.numeric(tempsm$tempMean)
tempmp$tempMean <- as.numeric(tempmp$tempMean)

sprtempSM <- mean(tempsm$tempMean[tempsm$month>2 & tempsm$month<5], na.rm = T)#March 1-April 30 (4 degrees C)
latSM <- 54.8907; longSM <- 126.8918
sprtempMP <- mean(tempmp$tempMean[tempmp$month>2 & tempmp$month<5], na.rm = T)#March 1-April 30 (4 degrees C)
latMP <- 49.0646; longMP <- 120.7816

# Now get phenology data:
pheno <- read.csv("input/dl_allbb_mini.csv")
phenoAlninc <- pheno[pheno$species %in% sp[1], ]
phenoBetpap <- pheno[pheno$species %in% sp[2], ]

bbdoyAlninc <- as.integer(mean(phenoAlninc$bb))
bbdoyBetpap <- as.integer(mean(phenoBetpap$bb))

# Start with a plot for Smithers
daylengthbbdoySMAI <- daylength(latSM, (bbdoyAlninc+59))#$Daylength + julian march 
daylengthbbdoySMBP <- daylength(latSM, (bbdoyBetpap+59))
chillportSM <- mean(chillsm$Chill_portions)
chillSM<-mean(chillsm$Utah_Model)/240

daylengthbbdoyMPAI <- daylength(latMP, (bbdoyAlninc+59))#$Daylength + julian march 
daylengthbbdoyMPBP <- daylength(latMP, (bbdoyBetpap+59))
chillportMP <- mean(chillmp$Chill_portions)
chillMP<-mean(chillmp$Utah_Model)/240

#########################################################################
# Let's start with a simple simulation of Alninc in Smithers:

predicts <- as.data.frame(matrix(NA,ncol=5,nrow=7))
predicts.25per <- as.data.frame(matrix(NA,ncol=5,nrow=7))
predicts.75per <- as.data.frame(matrix(NA,ncol=5,nrow=7))

#ad hoc adj daylength -leaving this for now
# predicts.wdl <- as.data.frame(matrix(NA,ncol=5,nrow=7))
# predicts.25per.wdl <- as.data.frame(matrix(NA,ncol=5,nrow=7))
# predicts.75per.wdl <- as.data.frame(matrix(NA,ncol=5,nrow=7))

colnames(predicts)<-colnames(predicts.25per) <-colnames(predicts.75per) <-
  #colnames(predicts.wdl)<-colnames(predicts.25per.wdl) <-colnames(predicts.75per.wdl) <- 
  c("warming","nowarm","sprwarm","winwarm","bothwarm")

  listofdraws <- rstan::extract(fit)
  
  s<-1
  avgbb <- listofdraws$a_sp[,sp.num[s]] + listofdraws$b_warm[,sp.num[s]]*sprtempSM +
    listofdraws$b_photo[,sp.num[s]]*daylengthbbdoySMAI + listofdraws$b_chill[,sp.num[s]]*chillportSM
  
  warmsprbb <- listofdraws$a_sp[,sp.num[s]] + listofdraws$b_warm[,sp.num[s]]*(sprtempSM + tempforecast[j]) +
    listofdraws$b_photo[,sp.num[s]]*(daylength + daylengthbbdoySMAI) + listofdraws$b_chill[,sp.num[s]]*chillportSM
  
  warmwinbb <- listofdraws$a_sp[,sp.num[s]] + listofdraws$b_warm[,sp.num[s]]*sprtempSM +
    listofdraws$b_photo[,sp.num[s]]*(daylength + daylengthbbdoySMAI) + listofdraws$b_chill[,sp.num[s]]*(chillportSM - (chillportSM/2))
  
  warmsprwinbb <- listofdraws$a_sp[,sp.num[s]] + listofdraws$b_warm[,sp.num[s]]*(sprtempSM + tempforecast[j]) +
    listofdraws$b_photo[,sp.num[s]]*(daylength + daylengthbbdoySMAI) + listofdraws$b_chill[,sp.num[s]]*(chillportSM - (chillportSM/2))
  
  yebbest <- list(avgbb, warmsprbb, warmwinbb, warmsprwinbb)
  return(yebbest)
}


photo.forplotAI <- daylengthbbdoySMAI
warmspring <-tempforecast[j]
warmwinterAI <- mean(chillsm$Chill_portions)-chillportSM

# zeros are bc not altering daylength
bbposteriors <- getspest.bb(fit, sprtempSM, daylengthbbdoySMAI, chillportSM, warmspring, warmwinterAI, 0, 0, 0)

meanz <- unlist(lapply(bbposteriors, mean))

quantz <- lapply(bbposteriors, function(x) quantile(x,  c(0.25, 0.5, 0.75)))

quant25per <- unlist(lapply(bbposteriors, function(x) quantile(x,  c(0.25))))
quant75per <- unlist(lapply(bbposteriors, function(x) quantile(x,  c(0.75))))
daychange.springwarm<-meanz[2]-meanz[1]
daychange.wintwarm<-meanz[3]-meanz[1]
daychange.bothwarm<-meanz[4]-meanz[1]
daylengthchange.springwarm<-daylength(latSM,bbdoyAlninc+daychange.springwarm)-daylengthbbdoySMAI
daylengthchange.wintwarm<- daylength(latSM,bbdoyAlninc+daychange.wintwarm)-daylengthbbdoySMAI
daylengthchange.bothwarm<-daylength(latSM,bbdoyAlninc+daychange.bothwarm)-daylengthbbdoySMAI

bbposteriors.wdaylength <- getspest.bb(fit, sprtempSM, daylengthbbdoySMAI, chillportSM, warmspring, warmwinterAI, daylengthchange.springwarm, daylengthchange.wintwarm, daylengthchange.bothwarm)

predicts[j,]<-c(warmspring,meanz,chillportSM,warmwinterAI)
predicts.25per[j,]<-c(warmspring,quant25per)
predicts.75per[j,]<-c(warmspring,quant75per)

# predicts<-rbind(c(0,predicts$nowarm[1:4],chillportSM,0),predicts)
# predicts<-predicts[,-2]
# predicts.25per<-rbind(c(0,predicts.25per$nowarm[1:4]),predicts.25per)
# predicts.25per<-predicts.25per[,-2]
# predicts.75per<-rbind(c(0,predicts.75per$nowarm[1:4]),predicts.75per)
# predicts.75per<-predicts.75per[,-2]

}
predicts$lat<-latSM
predicts$lon<-longSM
predicts.25per$lat<-latSM
predicts.25per$lon<-longSM
predicts.75per$lat<-latSM
predicts.75per$lon<-longSM
spests<-c()
spests<-rbind(spests,predicts)

ymin = 13#min(predicts[,-1],predicts.25per[,-1],predicts.75per[,-1])
ymax = 35#max(predicts[,-1],predicts.25per[,-1],predicts.75per[,-1])
xlim = c(0, 7)
ylim = c(ymin,ymax)

plot(x=NULL,y=NULL, xlim=xlim, xlab="Amount of warming (°C)", ylim=ylim,
     ylab="Days to budburst", main=maintext, bty="l", cex.axis = 1.2, cex.lab = 1.2)
pos.x <- 0
pos.y <- predicts[1,2]
points(pos.x, pos.y, cex=1.2, pch=19, bg="gray")

for(t in 3:5){
  polygon(c(rev(predicts$warming), predicts$warming), c(rev(predicts.75per[,t-1]), predicts.25per[,t-1]), col = alpha(cols[t-2], 0.2), border = NA)
}

for(t in 3:5){
  lines(predicts$warming, predicts[,t-1],
        col=cols[t-2], lwd=2)}