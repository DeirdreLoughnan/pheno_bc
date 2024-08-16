# Started May 3, 2022 by deirde

# the aim of this code is to generate the model output for my phenology ms
# rm(list=ls())
# options(stringsAsFactors = FALSE)

library(rstan)
#library(shinystan)
#library(reshape2)
#library(bayesplot)
library(ggplot2)
# library(dplyr)
library(plyr)
library(stringr)
library(phytools)



# if(length(grep("deirdreloughnan", getwd()) > 0)) { 
#   setwd("~/Documents/github/pheno_bc") 
# }  

#load("output/final/ew_phylo_output_newpriors.Rda")
dl <- read.csv("..//input/dl_allbb_mini.csv")

temp <- str_split_fixed(dl$trt, "_", 3); head(temp)
dl$chill<- temp[,1]
dl$photo <- temp[,2]
dl$force <- temp[,3]

dl.chill <- read.csv("..//input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")
dl.wchill$transect <- "west"

df <- read.csv("..//input/df_dxb_prepped_data.csv")
df.chill <- read.csv("..//input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill <- df.wchill[, c("population", "chill","force","photo","lab2", "bb","species", "treatment","Chill_portions","Utah_Model")]
df.wchill$transect <- "east"

# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)

head(pheno)
# combined the data has 3197 unique samples
############################################################
# Preping the data for the model
pheno <- subset(pheno, chill != "chill2")
pheno$force.n <- pheno$force
# pheno$force.n[pheno$force.n == "HF"] <- "1"
# pheno$force.n[pheno$force.n == "LF"] <- "0"
# pheno$force.n <- as.numeric(pheno$force.n)
pheno$force.n[pheno$force.n == "HF" & pheno$population == "mp"] <- "15"
pheno$force.n[pheno$force.n == "HF" & pheno$population == "sm"] <- "15"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "mp"] <- "10"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "sm"] <- "10"

pheno$force.n[pheno$force.n == "HF" & pheno$population == "HF"] <- "13.33"
pheno$force.n[pheno$force.n == "HF" & pheno$population == "SH"] <- "13.33"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "HF"] <- "8.33"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "SH"] <- "8.33"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$site.n <- as.character(pheno$population)
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.numeric(pheno$site.n)

pheno$transect.n <- pheno$transect
pheno$transect.n[pheno$transect.n == "east"] <- "1"
pheno$transect.n[pheno$transect.n == "west"] <- "0"
pheno$transect.n <- as.numeric(pheno$transect.n)

head(pheno)
#add dummy/ site level effects:
# pheno <- pheno %>%
#   mutate ( site2 = if_else(site.n == 2, 1, 0),
#            site3 = if_else(site.n == 3, 1, 0),
#            site4 = if_else(site.n == 4, 1, 0))

pheno$site2 <- ifelse(pheno$site.n == "2", "1", pheno$site.n)
pheno$site2 <- ifelse(pheno$site.n == c("1"), "0", pheno$site2)
pheno$site2 <- ifelse(pheno$site.n == c("3"), "0", pheno$site2)
pheno$site2 <- ifelse(pheno$site.n == c("4"), "0", pheno$site2)

pheno$site3 <- ifelse(pheno$site.n == "3", "1", pheno$site.n)
pheno$site3 <- ifelse(pheno$site.n == c("1"), "0", pheno$site3)
pheno$site3 <- ifelse(pheno$site.n == c("2"), "0", pheno$site3)
pheno$site3 <- ifelse(pheno$site.n == c("4"), "0", pheno$site3)

pheno$site4 <- ifelse(pheno$site.n == "4", "1", pheno$site.n)
pheno$site4 <- ifelse(pheno$site.n == c("1"), "0", pheno$site4)
pheno$site4 <- ifelse(pheno$site.n == c("2"), "0", pheno$site4)
pheno$site4 <- ifelse(pheno$site.n == c("3"), "0", pheno$site4)

pheno$site2 <- as.numeric(pheno$site2)
pheno$site3 <- as.numeric(pheno$site3)
pheno$site4 <- as.numeric(pheno$site4)

# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)

pheno$transect.n <- pheno$transect
pheno$transect.n[pheno$transect.n == "east"] <- "1"
pheno$transect.n[pheno$transect.n == "west"] <- "0"
pheno$transect.n <- as.numeric(pheno$transect.n)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
#pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term), ] # 3609

pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
#sort(unique(pheno.t$species.fact)) # 49 

# now get the phylogeny and pair it with species names:
spInfo <- read.csv("..//input/species_list.csv")
# head(spInfo)
# head(pheno.t)
# change the df to lowercase:

pheno.t <- merge(pheno.t, spInfo, by = "species")

head(pheno.t)

westSite <- c("sm", "mp")
westPheno <- pheno.t[pheno.t$population %in% westSite, ]

meanBBWest <- round(mean(westPheno$bb), 1)
meanBBWestU <- format(round(quantile(westPheno$bb, c(0.95)),1), nsmall =1)
meanBBWestL <- format(round(quantile(westPheno$bb, c(0.05)),1), nsmall =1)

eastSite <- c("HF", "SH")
eastPheno <- pheno.t[pheno.t$population %in% eastSite, ]

meanBBEast <- round(mean(eastPheno$bb), 1)
meanBBEastU <- format(round(quantile(eastPheno$bb, c(0.95)),1), nsmall =1)
meanBBEastL <-format(round(quantile(eastPheno$bb, c(0.05)),1), nsmall =1)

#################################################################
# load("output/final/ew_phylo_output_newpriors_allncp.Rda")
# sumew <- summary(mdl.ewphylo)$summary

#load("..//output/bb_4sites_phylo_contphotothermo_zscored_Apr19.Rda")
load("..//output/bb_phylo_contphotothermo_2zscoredMay13.Rda")

sum <- summary(mdl.2z)$summary
fit <- rstan::extract(mdl.2z)

# sum <- summary(mdl.4phyloContWP)$summary 
# fit <- rstan::extract(mdl.4phyloContWP)

# load("..//output/final/dl_phylo_lat1.Rda")
# load("..//output/final/dl_phylo_lat50.Rda")
#############################################
col4fig <- c("mean","sd","25%","50%","75%","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")

mu_params_4 <- c( 
               #"a_z",
               # "lam_interceptsa",
                "mu_b_warm",
                "mu_b_photo",
                "mu_b_chill1",
                "b_site2",
                "b_site3",
                "b_site4",
                "mu_b_inter_wp",
                "mu_b_inter_wc1",
                "mu_b_inter_pc1",
                "mu_b_inter_ws2",
                "mu_b_inter_ps2",
                "mu_b_inter_s2c1",
                "mu_b_inter_ws3",
                "mu_b_inter_ps3",
                "mu_b_inter_s3c1",
                "mu_b_inter_ws4",
                "mu_b_inter_ps4",
                "mu_b_inter_s4c1")

meanz4 <- sum[mu_params_4, col4table]

rownames(meanz4) = c( 
                      #"Root trait intercept", "Lambda",
                      "Forcing",
                      "Photoperiod",
                      "Chilling",
                      "Manning Park",
                       "Harvard Forest",
                       "St. Hippolyte",
                      "Forcing x photoperiod",
                      "Forcing x chilling",
                      "Photoperiod x chilling",
                      "Forcing x Manning Park",
                      "Photoperiod x Manning Park",
                      "Chilling x Manning Park",
                      "Forcing x Harvard Forest",
                      "Photoperiod x Harvard Forest",
                      "Chilling x Harvard Forest",
                      "Forcing x St. Hippolyte",
                      "Photoperiod x St. Hippolyte",
                      "Chilling x St. Hippolyte"            
)

meanz4.table <- sum[mu_params_4, col4table]
row.names(meanz4.table) <- row.names(meanz4)
head(meanz4.table)
#write.table(meanzew.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)

# call table values:

chillCue <- round(meanz4["Chilling",1],1)
chillCueU <- round(quantile(fit$mu_b_chill1, c(0.95)),1)
chillCueL <- round(quantile(fit$mu_b_chill1, c(0.05)),1)

chillCue2.5 <- round(meanz4["Chilling",3],1)
chillCue97.5 <- round(meanz4["Chilling",5],1)

photoCue <- round(meanz4["Photoperiod",1],1)
photoCueU <- format(round(quantile(fit$mu_b_photo, c(0.95)),1), nsmall=1)
photoCueL <- round(quantile(fit$mu_b_photo, c(0.05)),1)
photoCue2.5 <- round(meanz4["Photoperiod",3],1)
photoCue97.5 <- round(meanz4["Photoperiod",5],1)

forceCue <- round(meanz4["Forcing",1],1)
forceCueU <- round(quantile(fit$mu_b_warm, c(0.95)),1)
forceCueL <- round(quantile(fit$mu_b_warm, c(0.05)),1)

intrxnCF <- round(meanz4["Forcing x chilling",1],1)
intrxnCFU <- round(quantile(fit$mu_b_inter_wc1, c(0.95)),1)
intrxnCFL <- round(quantile(fit$mu_b_inter_wc1, c(0.05)),1)
intrxnCF2.5 <- round(meanz4["Forcing x chilling",3],1)
intrxnCF97.5 <- round(meanz4["Forcing x chilling",5],1)

siteHF <- round(meanz4["Harvard Forest",1],1)
siteHFU <- round(quantile(fit$b_site3, c(0.95)),1)
siteHFL <- round(quantile(fit$b_site3, c(0.05)),1)
siteHF2.5 <- round(meanz4["Harvard Forest",3],1)
siteHF97.5 <- round(meanz4["Harvard Forest",5],1)

siteSH <- round(meanz4["St. Hippolyte",1],1)
siteSHU <- round(quantile(fit$b_site4, c(0.95)),1)
siteSHL <- round(quantile(fit$b_site4, c(0.05)),1)
siteSH2.5 <- round(meanz4["St. Hippolyte",3],1)
siteSH97.5 <- round(meanz4["St. Hippolyte",5],1)

siteMP <- round(meanz4["Manning Park",1],1)
siteMPU <- round(quantile(fit$b_site2, c(0.95)),1)
siteMPL <- round(quantile(fit$b_site2, c(0.05)),1)
siteMP2.5 <- round(meanz4["Manning Park",3],1)
siteMP97.5 <- round(meanz4["Manning Park",5],1)

forceHF <- round(meanz4["Forcing x Harvard Forest",1],1)
forceHFU <- round(quantile(fit$mu_b_inter_ws3, c(0.95)),1)
forceHFL <- round(quantile(fit$mu_b_inter_ws3, c(0.05)),1)
forceHF2.5 <- round(meanz4["Forcing x Harvard Forest",3],1)
forceHF97.5 <- round(meanz4["Forcing x Harvard Forest",5],1)

forceSH <- round(meanz4["Forcing x St. Hippolyte",1],1)
forceSHU <- round(quantile(fit$mu_b_inter_ws4, c(0.95)),1)
forceSHL <- round(quantile(fit$mu_b_inter_ws4, c(0.05)),1)
forceSH2.5 <- round(meanz4["Forcing x St. Hippolyte",3],1)
forceSH97.5 <- round(meanz4["Forcing x St. Hippolyte",5],1)

forceMP <- round(meanz4["Forcing x Manning Park",1],1)
forceMPU <- round(quantile(fit$mu_b_inter_ws2, c(0.95)),1)
forceMPL <- round(quantile(fit$mu_b_inter_ws2, c(0.05)),1)
forceMP2.5 <- round(meanz4["Forcing x Manning Park",3],1)
forceMP97.5 <- round(meanz4["Forcing x Manning Park",5],1)

# forceSH <- round(meanz4["Forcing x St. Hippolyte",1],1)
# forceSH2.5 <- round(meanz4["Forcing x St. Hippolyte",3],1)
# forceSH97.5 <- round(meanz4["Forcing x St. Hippolyte",5],1)

chillHF <- round(meanz4["Chilling x Harvard Forest",1],1)
chillHFU <- round(quantile(fit$mu_b_inter_s3c1, c(0.95)),1)
chillHFL <- round(quantile(fit$mu_b_inter_s3c1, c(0.05)),1)
chillHF2.5 <- round(meanz4["Chilling x Harvard Forest",3],1)
chillHF97.5 <- round(meanz4["Chilling x Harvard Forest",5],1)

chillSH <- round(meanz4["Chilling x St. Hippolyte",1],1)
chillSHU <- round(quantile(fit$mu_b_inter_s4c1, c(0.95)),1)
chillSHL <- round(quantile(fit$mu_b_inter_s4c1, c(0.05)),1)
chillSH2.5 <- round(meanz4["Chilling x St. Hippolyte",3],1)
chillSH97.5 <- round(meanz4["Chilling x St. Hippolyte",5],1)

chillMP <- round(meanz4["Chilling x Manning Park",1],1)
chillMPU <- round(quantile(fit$mu_b_inter_s2c1, c(0.95)),1)
chillMPL <- round(quantile(fit$mu_b_inter_s2c1, c(0.05)),1)
chillMP2.5 <- round(meanz4["Chilling x Manning Park",3],1)
chillMP97.5 <- round(meanz4["Chilling x Manning Park",5],1)

photoHF <- round(meanz4["Photoperiod x Harvard Forest",1],1)
photoHFU <- round(quantile(fit$mu_b_inter_ps3, c(0.95)),1)
photoHFL <- round(quantile(fit$mu_b_inter_ps3, c(0.05)),1)
photoHF2.5 <- round(meanz4["Photoperiod x Harvard Forest",3],1)
photoHF97.5 <- round(meanz4["Photoperiod x Harvard Forest",5],1)

photoSH <- round(meanz4["Photoperiod x St. Hippolyte",1],1)
photoSHU <- round(quantile(fit$mu_b_inter_ps4, c(0.95)),1)
photoSHL <- round(quantile(fit$mu_b_inter_ps4, c(0.05)),1)
photoSH2.5 <- round(meanz4["Photoperiod x St. Hippolyte",3],1)
photoSH97.5 <- round(meanz4["Photoperiod x St. Hippolyte",5],1)

photoMP <- round(meanz4["Photoperiod x Manning Park",1],1)
photoMPU <- round(quantile(fit$mu_b_inter_ps2, c(0.95)),1)
photoMPL <- round(quantile(fit$mu_b_inter_ps2, c(0.05)),1)
photoMP2.5 <- round(meanz4["Photoperiod x Manning Park",3],1)
photoMP97.5 <- round(meanz4["Photoperiod x Manning Park",5],1)

# General survival values:

# Total number samples that went into exp:

totChill <- 8 * 20 * 2 * 8 # 8 trt X 2 sites X 20 sp (19 both 2 at just one) X 8 reps
  #length(unique(begin$lab2))

d <- read.csv("..//input/bc_phenology_Feb52021.csv", header=TRUE, na.strings=c("","NA"))
head(d)

end <- subset(d, day == 88)
totExpSurv <- length(unique(end$lab2))

perMort <- round((totChill- totExpSurv)/(totChill)*100, 1)

nobb <- length(unique(pheno$lab2))

perNoBB <- round(((totExpSurv-nobb)/totExpSurv)*100,1)

noTermEnd <- subset(end, bbch.t <7)
noTerm <- length(unique(noTermEnd$lab2))

perNoTerm <- round(((noTerm/totExpSurv)*100),1)

noTermEnd$count <- 1
temp <- aggregate(noTermEnd["count"], noTermEnd[c("species")], FUN = sum)

# what is the min and max average bb?
meanBBData <- aggregate(pheno.t["bb"], pheno.t[c("species")], FUN = mean)

minBB <- round(min(meanBBData$bb), 1)
maxBB <- round(max(meanBBData$bb), 1)

meanBBPop <- aggregate(pheno.t["bb"], pheno.t[c("species", "transect")], FUN = mean)
meanBBEast <- subset(meanBBPop, transect != "west")
meanBBWest <- subset(meanBBPop, transect != "east")

minBBE <- round(min(meanBBEast$bb), 1)
maxBBE <- round(max(meanBBEast$bb), 1)
diffBBEast <- maxBBE-minBBE

minBBW <- round(min(meanBBWest$bb), 1)
maxBBW <- round(max(meanBBWest$bb), 1)
diffBBWest <- maxBBW-minBBW

meanBB <- round(mean(meanBBData$bb),1)
meanBBLowerer <- round(quantile(meanBBData$bb, c(0.05)),1)
meanBBUpper <- format(round(quantile((meanBBData$bb), c(0.95)),1),nsmall =1)
# What are the lambda values?
lam_params <- c( 
  "a_z",
  "lam_interceptsa"
  )
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")

phylo <- sum[lam_params, col4table]

rootT <- round(phylo[1,1],1)

rootTLower  <- round(quantile(fit$a_z, c(0.95)),1)
rootTUpper <- round(quantile(fit$a_z, c(0.05)) ,1)

# rootTLower <- as.numeric(round(HPDI(data.frame(fit$a_z), prob = 0.90)[1],1))
# rootTUpper <- as.numeric(round(HPDI(data.frame(fit$a_z), prob = 0.90)[2],1))

lamT <- round(phylo[2,1],1)

lamTLower <- round(quantile(fit$lam_interceptsa, c(0.95)),1)
lamTUpper <- round(quantile(fit$lam_interceptsa, c(0.05)),1) 

# lamTLower <- as.numeric(round(HPDI(data.frame(fit$lam_interceptsa), prob = 0.90)[1],1))
# lamTUpper <- as.numeric(round(HPDI(data.frame(fit$lam_interceptsa), prob = 0.90)[2],1))

# First to budburst
phenoBB <- subset(pheno.t, bb >0)
first <- subset(phenoBB, bb == min(phenoBB$bb))

last <- subset(phenoBB, bb == max(phenoBB$bb))
#### How different are the lateral vs 50 lat estimates? 
# sum1 <- summary(mdlLat1)$summary
# 
# col4fig <- c("mean","sd","25%","50%","75%","n_eff","Rhat")
# col4table <- c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")
# 
# mu_params_1l <- c(
#   "mu_b_warm", "mu_b_photo", "mu_b_chill1", "b_site","mu_b_inter_wp",
#   "mu_b_inter_wc1","mu_b_inter_pc1", "mu_b_inter_ws",
#   "mu_b_inter_ps","mu_b_inter_sc1")
# 
# meanz1l <- sum1[mu_params_1l, col4fig]
# 
# rownames(meanz1l) = c( #"Root trait intercept","Lambda",
#   "Forcing",
#   "Photoperiod",
#   "Chilling",
#   "Manning Park",
#   "Forcing x photoperiod",
#   "Forcing x chilling",
#   "Photoperiod x chilling",
#   "Forcing x Manning Park",
#   "Photoperiod x Manning Park",
#   "Chilling x Manning Park"
#   
#   
# )
# 
# 
# sum50 <- summary(mdlLat50)$summary
# 
# col4fig <- c("mean","sd","25%","50%","75%","n_eff","Rhat")
# col4table <- c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")
# 
# mu_params_50l <- c(
#   #"a_z","lam_interceptsa",
#   "mu_b_warm", "mu_b_photo", "mu_b_chill1", "b_site", "mu_b_inter_wp",
#   "mu_b_inter_wc1","mu_b_inter_pc1", "mu_b_inter_ws",
#   "mu_b_inter_ps","mu_b_inter_sc1")
# 
# meanz50l <- sum50[mu_params_50l, col4fig]
# 
# rownames(meanz50l) = c( #"Root trait intercept","Lambda",
#   "Forcing",
#   "Photoperiod",
#   "Chilling",
#   "Manning Park",
#   "Forcing x photoperiod",
#   "Forcing x chilling",
#   "Photoperiod x chilling",
#   "Forcing x Manning Park",
#   "Photoperiod x Manning Park",
#   "Chilling x Manning Park"
#   
# )
# 
# (meanz50l[,1] - meanz1l[,1])/meanz50l[,1]

###################################

a_sp = (sum[grep("a_sp", rownames(sum)), 1])
a_spU = round(quantile(fit$a_sp, c(0.95)),1)
a_spL = round(quantile(fit$a_sp, c(0.05)),1)

a_sp97.5 = (sum[grep("a_sp", rownames(sum)), "97.5%"])
a_sp2.5 = (sum[grep("a_sp", rownames(sum)), "2.5%"])

b_photo = sum[grep("b_photo\\[", rownames(sum)), 1]
b_chill = sum[grep("b_chill1\\[", rownames(sum)), 1]
b_force = sum[grep("b_warm\\[", rownames(sum)), 1]

b_photo25 = sum[grep("b_photo\\[", rownames(sum)), "25%"]
b_chill25 = sum[grep("b_chill1\\[", rownames(sum)), "25%"]
b_force25 = sum[grep("b_warm\\[", rownames(sum)), "25%"]

b_photo75 = sum[grep("b_photo\\[", rownames(sum)), "75%"]
b_chill75 = sum[grep("b_chill1\\[", rownames(sum)), "75%"]
b_force75 = sum[grep("b_warm\\[", rownames(sum)), "75%"]

b_photo2.5 = sum[grep("b_photo\\[", rownames(sum)), "2.5%"]
b_chill2.5 = sum[grep("b_chill1\\[", rownames(sum)), "2.5%"]
b_force2.5 = sum[grep("b_warm\\[", rownames(sum)), "2.5%"]

b_photo97.5 = sum[grep("b_photo\\[", rownames(sum)), "97.5%"]
b_chill97.5 = sum[grep("b_chill1\\[", rownames(sum)), "97.5%"]
b_force97.5 = sum[grep("b_warm\\[", rownames(sum)), "97.5%"]

a_sp5 <- vector()
for(i in 1:ncol(fit$a_sp)){
  quantU <- round(quantile(fit$a_sp[,i], c(0.05)),1)
  a_sp5 <- rbind(a_sp5, quantU)
}
colnames(a_sp5) <- c("Int5")

a_sp95 <- vector()
for(i in 1:ncol(fit$a_sp)){
  quantU <- round(quantile(fit$a_sp[,i], c(0.95)),1)
  a_sp95 <- rbind(a_sp95, quantU)
}
colnames(a_sp95) <- c("Int95")

b_chill5 <- vector()
for(i in 1:ncol(fit$b_chill1)){
  quantU <- round(quantile(fit$b_chill1[,i], c(0.05)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5")

b_chill95 <- vector()
for(i in 1:ncol(fit$b_chill1)){
  quantU <- round(quantile(fit$b_chill1[,i], c(0.95)),1)
  b_chill95 <- rbind(b_chill95, quantU)
}
colnames(b_chill95) <- c("chill95")

b_force5 <- vector()
for(i in 1:ncol(fit$b_force)){
  quantU <- round(quantile(fit$b_force[,i], c(0.05)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5")

b_force95 <- vector()
for(i in 1:ncol(fit$b_force)){
  quantU <- round(quantile(fit$b_force[,i], c(0.95)),1)
  b_force95 <- rbind(b_force95, quantU)
}
colnames(b_force95) <- c("force95")

b_photo5 <- vector()
for(i in 1:ncol(fit$b_photo)){
  quantU <- round(quantile(fit$b_photo[,i], c(0.05)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5")

b_photo95 <- vector()
for(i in 1:ncol(fit$b_photo)){
  quantU <- round(quantile(fit$b_photo[,i], c(0.95)),1)
  b_photo95 <- rbind(b_photo95, quantU)
}
colnames(b_photo95) <- c("photo95")

# b_photo5 <- quantile(fit$mu_b_photo, c(0.05))
# b_chill5 <- quantile(fit$mu_b_chill1, c(0.05))
# b_force5 <- quantile(fit$mu_b_warm, c(0.05))
# 
# b_photo95 <- quantile(fit$mu_b_photo, c(0.95))
# b_chill95 <- quantile(fit$mu_b_chill1, c(0.95))
# b_force95 <- quantile(fit$mu_b_warm, c(0.95))


photo <- -0.5041133 #8 h photo
siteSM <- 0
force <- -0.5077191 #15 C trt
chill <- -0.4023109 # 5/15

m <- matrix(nrow = 1000, ncol = 47)

post <- fit

for(sp in 1:47){
  for (it in 1:nrow(m)){
    m[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM +
      post$b_warm[it,sp] * force + post$b_photo[it, sp] * photo + post$b_chill[it,sp] * chill +
      post$b_inter_wp[it,sp] * (force*photo) + post$b_inter_wc1[it,sp] * (force*chill) + post$b_inter_pc1[it,sp] * (photo*chill) +
      post$b_inter_s2c1[it,sp] * (chill*siteSM) + post$b_inter_ws2[it,sp] * (force*siteSM) + post$b_inter_ps2[it,sp] * (photo*siteSM) +
      post$b_inter_s3c1[it,sp] * (chill*siteSM) + post$b_inter_ws3[it,sp] * (force*siteSM) + post$b_inter_ps3[it,sp] * (photo*siteSM) +
      post$b_inter_s4c1[it,sp] * (chill*siteSM) + post$b_inter_ws4[it,sp] * (force*siteSM) + post$b_inter_ps4[it,sp] * (photo*siteSM)
  }
}

photoHigh <- 0.4965051 #8 h photo
siteSM <- 0
forceHigh <- 0.5877121 
chillHigh <- 0.3660412 # high chill for Smithers

mHigh <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(mHigh)){
    mHigh[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM + 
      post$b_warm[it,sp] * forceHigh + post$b_photo[it, sp] * photoHigh + post$b_chill[it,sp] * chillHigh +
      post$b_inter_wp[it,sp] * (forceHigh*photoHigh) + post$b_inter_wc1[it,sp] * (forceHigh*chillHigh) + post$b_inter_pc1[it,sp] * (photoHigh*chillHigh) +
      post$b_inter_s2c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws2[it,sp] * (forceHigh*siteSM) + post$b_inter_ps2[it,sp] * (photoHigh*siteSM) +
      post$b_inter_s3c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws3[it,sp] * (forceHigh*siteSM) + post$b_inter_ps3[it,sp] * (photoHigh*siteSM) +
      post$b_inter_s4c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws4[it,sp] * (forceHigh*siteSM) + post$b_inter_ps4[it,sp] * (photoHigh*siteSM)
  }
}

meanLowBB <- mean(m)
meanHighBB <- mean(mHigh)


#manning park
# mp <- matrix(nrow = 1000, ncol = 47)
# for(sp in 1:47){
#   for (it in 1:nrow(mp)){
#     mp[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * 1 + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM +
#       post$b_warm[it,sp] * force + post$b_photo[it, sp] * photo + post$b_chill[it,sp] * chill +
#       post$b_inter_wp[it,sp] * (force*photo) + post$b_inter_wc1[it,sp] * (force*chill) + post$b_inter_pc1[it,sp] * (photo*chill) +
#       post$b_inter_s2c1[it,sp] * (chill*1) + post$b_inter_ws2[it,sp] * (force*1) + post$b_inter_ps2[it,sp] * (photo*1) +
#       post$b_inter_s3c1[it,sp] * (chill*siteSM) + post$b_inter_ws3[it,sp] * (force*siteSM) + post$b_inter_ps3[it,sp] * (photo*siteSM) +
#       post$b_inter_s4c1[it,sp] * (chill*siteSM) + post$b_inter_ws4[it,sp] * (force*siteSM) + post$b_inter_ps4[it,sp] * (photo*siteSM)
#   }
# }
# 
# #manning park
# hf <- matrix(nrow = 1000, ncol = 47)
# for(sp in 1:47){
#   for (it in 1:nrow(mp)){
#     mp[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * 1 + post$b_site4[it] * siteSM +
#       post$b_warm[it,sp] * force + post$b_photo[it, sp] * photo + post$b_chill[it,sp] * chill +
#       post$b_inter_wp[it,sp] * (force*photo) + post$b_inter_wc1[it,sp] * (force*chill) + post$b_inter_pc1[it,sp] * (photo*chill) +
#       post$b_inter_s2c1[it,sp] * (chill*siteSM) + post$b_inter_ws2[it,sp] * (force*siteSM) + post$b_inter_ps2[it,sp] * (photo*siteSM) +
#       post$b_inter_s3c1[it,sp] * (chill*siteSM) + post$b_inter_ws3[it,sp] * (force*siteSM) + post$b_inter_ps3[it,sp] * (photo*siteSM) +
#       post$b_inter_s4c1[it,sp] * (chill*siteSM) + post$b_inter_ws4[it,sp] * (force*siteSM) + post$b_inter_ps4[it,sp] * (photo*siteSM)
#   }
# }
# 

# now get the order of diff spp bb that I can use to order the figure
#spInfo <- read.csv("..//input/species_list.csv")

spNo <- nrow(spInfo)

eastSpp <- subset(spInfo, transect != "west")
eastSpNo <- nrow(eastSpp)
eastShrub <- subset(eastSpp, type == "shrub")
eastTree <- subset(eastSpp, type == "tree")

eastShrubNo <- nrow(eastShrub)
eastTreeNo <- nrow(eastTree)

westSpp <- subset(spInfo, transect != "east")
westSpNo <- nrow(westSpp)
westShrub <- subset(westSpp, type == "shrub")
westTree <- subset(westSpp, type == "tree")

westShrubNo <- nrow(westShrub)
westTreeNo <- nrow(westTree)

bothSpp <- subset(spInfo, transect == "east/west")
bothSpNo <- nrow(bothSpp)

spInfo <- spInfo[order(spInfo$species),]
head(spInfo)
spInfo$meanBB <- colMeans(m)
colnames(m) <- spInfo$species.name

spInfo$meanBBHigh <- colMeans(mHigh)
colnames(mHigh) <- spInfo$species.name

spInfo$Int <- a_sp
spInfo$Int2.5 <- a_sp2.5
spInfo$Int97.5 <- a_sp97.5
spInfo$Int5 <- a_sp5
spInfo$Int95 <- a_sp95


spInfo$force <- b_force
spInfo$chill <- b_chill
spInfo$photo <- b_photo

spInfo$force2.5 <- b_force2.5
spInfo$chill2.5 <- b_chill2.5
spInfo$photo2.5 <- b_photo2.5

spInfo$force97.5 <- b_force97.5
spInfo$chill97.5 <- b_chill97.5
spInfo$photo97.5 <- b_photo97.5

spInfo$force25 <- b_force25
spInfo$chill25 <- b_chill25
spInfo$photo25 <- b_photo25

spInfo$force75 <- b_force75
spInfo$chill75 <- b_chill75
spInfo$photo75 <- b_photo75

spInfo$force5 <- b_force5
spInfo$chill5 <- b_chill5
spInfo$photo5 <- b_photo5

spInfo$force95 <- b_force95
spInfo$chill95 <- b_chill95
spInfo$photo95 <- b_photo95

QaLlDiff <- round(spInfo[28,"meanBB"] - spInfo[21,"meanBB"] ,1)
east <- subset(spInfo, transect != "west")
eastSp <- unique(east$species.name)

qaChill <- round(spInfo[28,"chill"] ,1)
qaForcing <- round(spInfo[28,"force"] ,1)

llChill <- round(spInfo[21,"chill"] ,1)
llForcing <- round(spInfo[21,"force"] ,1)

dataEast <- spInfo[spInfo$species.name %in% eastSp, ]

overlappingE <- c("aromel","vibcas","betlen","lyolig","rhopri")
spMiniE <- east[!east$species %in% overlappingE,]
spTopE <- east[east$species %in% overlappingE,]

west <- subset(spInfo, transect != "east")
westSp <- unique(west$species.name)

dataWest <- spInfo[spInfo$species.name %in% westSp, ]

overlappingW <- c("spialb","betpap","popbal")
spMiniW <- west[!west$species %in% overlappingW,]
spTopW <- west[west$species %in% overlappingW,]

meanPtW <- aggregate(dataWest[c("meanBB","meanBBHigh", "Int")], dataWest[c("species.name","type","transect")], FUN = mean)

meanPtE <- aggregate(dataEast[c("meanBB","meanBBHigh", "Int")], dataEast[c("species.name","type","transect")], FUN = mean)

diffBAWest <- round(mean((meanPtW$Int/meanPtW$meanBB)*100),1)
diffBAEast <- format(round(mean((meanPtE$Int/meanPtE$meanBB)*100),1),nsmall =1)

diffBAWestHigh <- 100-round(mean((meanPtW$meanBBHigh/meanPtW$Int)*100),1)
diffBAEastHigh <- format(100-round(mean((meanPtE$meanBBHigh/meanPtE$Int)*100),1),nsmall =1)

# Ranking info:
rank <- spInfo[,c("species.name","species","type","transect","meanBB","meanBBHigh","Int")]
rank <- rank[order(rank$Int),]
rank$rankInt <- seq(1:nrow(rank))

rank <- rank[order(rank$meanBB),]
rank$rankLowC <- seq(1:nrow(rank))

rank <- rank[order(rank$meanBBHigh),]
rank$rankHighC <- seq(1:nrow(rank))

# How much to ranks change?
rank$absDiff <- abs(rank$rankHighC-rank$rankLowC)
# eastern ranks
rankE <- subset(rank, transect != "west")

rankE <- rankE[order(rankE$Int),]
rankE$rankInt <- seq(1:nrow(rankE))

rankE <- rankE[order(rankE$meanBB),]
rankE$rankLowC <- seq(1:nrow(rankE))

rankE <- rankE[order(rankE$meanBBHigh),]
rankE$rankHighC <- seq(1:nrow(rankE))

minRankDiffE <- min(rankE$absDiff)
maxRankDiffE <- max(rankE$absDiff)
meanRankDiffE <- round(mean(rankE$absDiff),0)
# Western ranks
rankW <- subset(rank, transect != "east")

rankW <- rankW[order(rankW$Int),]
rankW$rankInt <- seq(1:nrow(rankW))

rankW <- rankW[order(rankW$meanBB),]
rankW$rankLowC <- seq(1:nrow(rankW))

rankW <- rankW[order(rankW$meanBBHigh),]
rankW$rankHighC <- seq(1:nrow(rankW))

minRankDiffW <- min(rankW$absDiff)
maxRankDiffW <- max(rankW$absDiff)
meanRankDiffW <- round(mean(rankW$absDiff),0)
# Subset just the spp that are in both transect:
both <- c("Alnus_incana","Betula_papyrifera","Populus_tremuloides")

rankESub <- rankE[rankE$species.name %in% both,]
rankWSub <- rankW[rankW$species.name %in% both,]

# How many early vs late spp were trees vs shrubs
earlyShrub <- nrow(subset(rank, rankLowC <24 & type == "shrub"))
earlyTree <- nrow(subset(rank, rankLowC <24 & type == "tree"))

lateShrub <- nrow(subset(rank, rankLowC >23 & type == "shrub"))
lateTree <- nrow(subset(rank, rankLowC >23 & type == "tree"))

nShrub <- nrow(subset(rank, type == "shrub"))
nTree <- nrow(subset(rank, type == "tree"))

perEarlyTree <- round((earlyTree/nTree)*100,1)
perLateShrub <- round((lateShrub/nShrub)*100,1)

perLateTree <- round((lateTree/nTree)*100,1)
perEarlyShrub <- round((earlyShrub/nShrub)*100,1)
## Interaction values:


#1 E vs w being overall earlier:
E <- subset(spInfo, transect != "west")
eastMean <- mean(E$meanBB)

W <- subset(spInfo, transect != "east")
westMean <- mean(W$meanBB)

# siteChillNoHigh
# eC <- round(mean(c(siteChillNoHigh[1,1],siteChillNoHigh[2,1],siteChillNoHigh[7,1],siteChillNoHigh[8,1])),1)
# wC <- round(mean(c(siteChillNoHigh[3,1],siteChillNoHigh[4,1],siteChillNoHigh[5,1],siteChillNoHigh[6,1])),1)
# 
# eC2.5 <-round( mean(c(siteChillNoHigh[1,3],siteChillNoHigh[2,3],siteChillNoHigh[7,3],siteChillNoHigh[8,3])),1)
# wC2.5 <- round(mean(c(siteChillNoHigh[3,3],siteChillNoHigh[4,3],siteChillNoHigh[5,3],siteChillNoHigh[6,3])),1)
# 
# eC97.5 <- round(mean(c(siteChillNoHigh[1,4],siteChillNoHigh[2,4],siteChillNoHigh[7,4],siteChillNoHigh[8,4])),1)
# wC97.5 <- round(mean(c(siteChillNoHigh[3,4],siteChillNoHigh[4,4],siteChillNoHigh[5,4],siteChillNoHigh[6,4])),1)
# 
# siteForce
# eF <- round(mean(c(siteForce[1,1],siteForce[2,1],siteForce[7,1],siteForce[8,1])),1)
# wF <- round(mean(c(siteForce[3,1],siteForce[4,1],siteForce[5,1],siteForce[6,1])),1)
# 
# eF2.5 <- round(mean(c(siteForce[1,3],siteForce[2,3],siteForce[7,3],siteForce[8,3])),1)
# wF2.5 <- round(mean(c(siteForce[3,3],siteForce[4,3],siteForce[5,3],siteForce[6,3])),1)
# 
# eF97.5 <- round(mean(c(siteForce[1,4],siteForce[2,4],siteForce[7,4],siteForce[8,4])),1)
# wF97.5 <- round(mean(c(siteForce[3,4],siteForce[4,4],siteForce[5,4],siteForce[6,4])),1)
# 
# sitePhoto <- sitePhoto[order(sitePhoto$site),]
# eP <- round(mean(c(sitePhoto[1,1],sitePhoto[2,1],sitePhoto[7,1],sitePhoto[8,1])),1)
# wP <- round(mean(c(sitePhoto[3,1],sitePhoto[4,1],sitePhoto[5,1],sitePhoto[6,1])),1)
# 
# eP2.5 <- round(mean(c(sitePhoto[1,3],sitePhoto[2,3],sitePhoto[7,3],sitePhoto[8,3])),1)
# wP2.5 <- round(mean(c(sitePhoto[3,3],sitePhoto[4,3],sitePhoto[5,3],sitePhoto[6,3])),1)
# 
# eP97.5 <- round(mean(c(sitePhoto[1,4],sitePhoto[2,4],sitePhoto[7,4],sitePhoto[8,4])),1)
# wP97.5 <- round(mean(c(sitePhoto[3,4],sitePhoto[4,4],sitePhoto[5,4],sitePhoto[6,4])),1)
# 
# # 1. lower lat in eastern diff from western in their forcing
# # eastern
# hfF <- round(mean(c(siteForce[1,1],siteForce[2,1])),1)
# shF <- round(mean(c(siteForce[7,1],siteForce[8,1])),1)
# 
# hfF2.5 <- round(mean(c(siteForce[1,3],siteForce[2,3])),1)
# shF2.5 <- round(mean(c(siteForce[7,3],siteForce[8,3])),1)
# 
# hfF97.5 <- round(mean(c(siteForce[1,4],siteForce[2,4])),1)
# shF97.5 <- round(mean(c(siteForce[7,4],siteForce[8,4])),1)
# 
# # western:
# mpF <- round(mean(c(siteForce[3,1],siteForce[4,1])),1)
# smF <- round(mean(c(siteForce[5,1],siteForce[6,1])),1)
# 
# mpF2.5 <- round(mean(c(siteForce[3,3],siteForce[4,3])),1)
# smF2.5 <- round(mean(c(siteForce[5,3],siteForce[6,3])),1)
# 
# mpF97.5 <- round(mean(c(siteForce[3,4],siteForce[4,4])),1)
# smF97.5 <- round(mean(c(siteForce[5,4],siteForce[5,4])),1)
# 
# 
# 
# a_sp95 <- round(quantile(fit$a_sp, c(0.05)) ,2)
# a_sp95 <- round(quantile(fit$a_sp, c(0.95)) ,2)
# 


# Shrub vs tree values in the main text:

typeValues <- aggregate(spInfo[c("force","chill","photo","force2.5","chill2.5","photo2.5", "force97.5","chill97.5","photo97.5","force5","chill5","photo5", "force95","chill95","photo95")], spInfo[c("type")], FUN = mean)

shrubMeanForce <- typeValues[1,"force"]
shrubMeanForce2.5 <- typeValues[1,"force2.5"]
shrubMeanForce97.5 <- typeValues[1,"force97.5"]
shrubMeanForce5 <- typeValues[1,"force5"]
shrubMeanForce95 <- typeValues[1,"force95"]

shrubMeanChill <- typeValues[1,"chill"]
shrubMeanChill2.5 <- typeValues[1,6]
shrubMeanChill97.5 <- typeValues[1,9]
shrubMeanChill5 <- typeValues[1,"chill5"]
shrubMeanChill95 <- typeValues[1,"chill95"]

shrubMeanPhoto <- typeValues[1,"photo"]
shrubMeanPhoto2.5 <- typeValues[1,7]
shrubMeanPhoto97.5 <- typeValues[1,10]
shrubMeanPhoto5 <- typeValues[1,"photo5"]
shrubMeanPhoto95 <- typeValues[1,"photo95"]

# trees
treeMeanForce <- typeValues[2,2]
treeMeanForce2.5 <- typeValues[2,5]
treeMeanForce97.5 <- typeValues[2,8]
treeMeanForce5 <- typeValues[2,"force5"]
treeMeanForce95 <- typeValues[2,"force95"]

treeMeanChill <- typeValues[2,3]
treeMeanChill2.5 <- typeValues[2,6]
treeMeanChill97.5 <- typeValues[2,9]
treeMeanChill5 <- typeValues[2,"chill5"]
treeMeanChill95 <- typeValues[2,"chill95"]

treeMeanPhoto <- typeValues[2,4]
treeMeanPhoto2.5 <- typeValues[2,7]
treeMeanPhoto97.5 <- typeValues[2,10]
treeMeanPhoto5 <- typeValues[2,"photo5"]
treeMeanPhoto95 <- typeValues[2,"photo95"]

# look at the tree shrub divide within transects:
# westData <- subset(spInfo, transect != "east")
# eastData <- subset(spInfo, transect != "west")
# bothData <-  subset(spInfo, transect == "east/west")
# 
# westTranVal <- aggregate(westData[c("force","chill","photo","force2.5","chill2.5","photo2.5", "force97.5","chill97.5","photo97.5")], westData[c("type")], FUN = mean)
# 
# eastTranVal <- aggregate(eastData[c("force","chill","photo","force2.5","chill2.5","photo2.5", "force97.5","chill97.5","photo97.5")], eastData[c("type")], FUN = mean)
# 
# bothTranVal <- aggregate(bothData[c("force","chill","photo","force2.5","chill2.5","photo2.5", "force97.5","chill97.5","photo97.5")], bothData[c("type")], FUN = mean)

tempC <- aggregate(pheno["bb"], pheno[c("transect","chill")], FUN = mean);tempC

tempF <- aggregate(pheno["bb"], pheno[c("transect","force")], FUN = mean);tempF

tempP <- aggregate(pheno["bb"], pheno[c("transect","photo")], FUN = mean);tempP

tempCS <- aggregate(pheno["bb"], pheno[c("population","chill")], FUN = mean);tempCS

tempFS <- aggregate(pheno["bb"], pheno[c("population","force")], FUN = mean);tempFS

tempPS <- aggregate(pheno["bb"], pheno[c("population","photo")], FUN = mean);tempPS

west <- subset(spInfo, transect != "east")
westAgg <- aggregate(west[c("force","chill","photo","force2.5","chill2.5","photo2.5", "force97.5","chill97.5","photo97.5","force5","chill5","photo5", "force95","chill95","photo95")], west[c("type")], FUN = mean)


# Adding posteriors to get the means!
a_sp = mean(sum[grep("a_sp", rownames(sum)), c("mean")])
a_sp5 <- round(quantile(post$a_sp, c(0.05)),1)
a_sp95 <- round(quantile(post$a_sp, c(0.95)),1)
a_sp25 <- round(quantile(post$a_sp, c(0.25)),1)
a_sp75 <- round(quantile(post$a_sp, c(0.75)),1)
a_sp <- cbind(a_sp, a_sp5,a_sp95, a_sp25,a_sp75)

#a_z = (sum[grep("a_z", rownames(sum)), c("mean")])
mu_b_warm = sum[grep("mu_b_warm", rownames(sum)), c("mean")]
mu_b_warm5 <- round(quantile(post$mu_b_warm, c(0.05)),1)
mu_b_warm95 <- round(quantile(post$mu_b_warm, c(0.95)),1)
mu_b_warm25 <- round(quantile(post$mu_b_warm, c(0.25)),1)
mu_b_warm75 <- round(quantile(post$mu_b_warm, c(0.75)),1)
mu_b_warm <- (cbind(mu_b_warm, mu_b_warm5,mu_b_warm95, mu_b_warm25,mu_b_warm75))

mu_b_photo = sum[grep("mu_b_photo", rownames(sum)), c("mean")]
mu_b_photo5 <- round(quantile(post$mu_b_photo, c(0.05)),1)
mu_b_photo95 <- round(quantile(post$mu_b_photo, c(0.95)),1)
mu_b_photo25 <- round(quantile(post$mu_b_photo, c(0.25)),1)
mu_b_photo75 <- round(quantile(post$mu_b_photo, c(0.75)),1)
mu_b_photo <- (cbind(mu_b_photo, mu_b_photo5,mu_b_photo95, mu_b_photo25,mu_b_photo75))

mu_b_chill1 = sum[grep("mu_b_chill1", rownames(sum)), c("mean")]
mu_b_chill5 <- round(quantile(post$mu_b_chill, c(0.05)),1)
mu_b_chill95 <- round(quantile(post$mu_b_chill, c(0.95)),1)
mu_b_chill25 <- round(quantile(post$mu_b_chill, c(0.25)),1)
mu_b_chill75 <- round(quantile(post$mu_b_chill, c(0.75)),1)
mu_b_chill1 <- (cbind(mu_b_chill1, mu_b_chill5,mu_b_chill95, mu_b_chill25,mu_b_chill75))

mu_b_inter_pc1 = sum[grep("mu_b_inter_pc1", rownames(sum)), c("mean")]
mu_b_inter_pc15 <- round(quantile(post$mu_b_inter_pc1, c(0.05)),1)
mu_b_inter_pc195 <- round(quantile(post$mu_b_inter_pc1, c(0.95)),1)
mu_b_inter_pc125 <- round(quantile(post$mu_b_inter_pc1, c(0.25)),1)
mu_b_inter_pc175 <- round(quantile(post$mu_b_inter_pc1, c(0.75)),1)
mu_b_inter_pc1 <- (cbind(mu_b_inter_pc1, mu_b_inter_pc15,mu_b_inter_pc195, mu_b_inter_pc125,mu_b_inter_pc175))

mu_b_inter_wp = sum[grep("mu_b_inter_wp", rownames(sum)), c("mean")]
mu_b_inter_wp5 <- round(quantile(post$mu_b_inter_wp, c(0.05)),1)
mu_b_inter_wp95 <- round(quantile(post$mu_b_inter_wp, c(0.95)),1)
mu_b_inter_wp25 <- round(quantile(post$mu_b_inter_wp, c(0.25)),1)
mu_b_inter_wp75 <- round(quantile(post$mu_b_inter_wp, c(0.75)),1)
mu_b_inter_wp <- (cbind(mu_b_inter_wp, mu_b_inter_wp5,mu_b_inter_wp95, mu_b_inter_wp25,mu_b_inter_wp75))

mu_b_inter_wc1 = sum[grep("mu_b_inter_wc1", rownames(sum)), c("mean")]
mu_b_inter_wc15 <- round(quantile(post$mu_b_inter_wc1, c(0.05)),1)
mu_b_inter_wc195 <- round(quantile(post$mu_b_inter_wc1, c(0.95)),1)
mu_b_inter_wc125 <- round(quantile(post$mu_b_inter_wc1, c(0.25)),1)
mu_b_inter_wc175 <- round(quantile(post$mu_b_inter_wc1, c(0.75)),1)
mu_b_inter_wc1 <- (cbind(mu_b_inter_wc1, mu_b_inter_wc15,mu_b_inter_wc195, mu_b_inter_wc125,mu_b_inter_wc175))

mu_b_inter_ws2 = sum[grep("mu_b_inter_ws2", rownames(sum)), c("mean")]
mu_b_inter_ws25 <- round(quantile(post$mu_b_inter_ws2, c(0.05)),1)
mu_b_inter_ws295 <- round(quantile(post$mu_b_inter_ws2, c(0.95)),1)
mu_b_inter_ws225 <- round(quantile(post$mu_b_inter_ws2, c(0.25)),1)
mu_b_inter_ws275 <- round(quantile(post$mu_b_inter_ws2, c(0.75)),1)
mu_b_inter_ws2 <- (cbind(mu_b_inter_ws2, mu_b_inter_ws25,mu_b_inter_ws295, mu_b_inter_ws225,mu_b_inter_ws275))

mu_b_inter_s2c1 = sum[grep("mu_b_inter_s2c1", rownames(sum)), c("mean")]
mu_b_inter_s2c15 <- round(quantile(post$mu_b_inter_s2c1, c(0.05)),1)
mu_b_inter_s2c195 <- round(quantile(post$mu_b_inter_s2c1, c(0.95)),1)
mu_b_inter_s2c125 <- round(quantile(post$mu_b_inter_s2c1, c(0.25)),1)
mu_b_inter_s2c175 <- round(quantile(post$mu_b_inter_s2c1, c(0.75)),1)
mu_b_inter_s2c1 <- (cbind(mu_b_inter_s2c1, mu_b_inter_s2c15,mu_b_inter_s2c195, mu_b_inter_s2c125,mu_b_inter_s2c175))

mu_b_inter_ps2 = sum[grep("mu_b_inter_ps2", rownames(sum)), c("mean")]
mu_b_inter_ps25 <- round(quantile(post$mu_b_inter_ps2, c(0.05)),1)
mu_b_inter_ps295 <- round(quantile(post$mu_b_inter_ps2, c(0.95)),1)
mu_b_inter_ps225 <- round(quantile(post$mu_b_inter_ps2, c(0.25)),1)
mu_b_inter_ps275 <- round(quantile(post$mu_b_inter_ps2, c(0.75)),1)
mu_b_inter_ps2 <- (cbind(mu_b_inter_ps2, mu_b_inter_ps25,mu_b_inter_ps295, mu_b_inter_ps225,mu_b_inter_ps275))

mu_b_inter_ws3 = sum[grep("mu_b_inter_ws3", rownames(sum)), c("mean")]
mu_b_inter_ws35 <- round(quantile(post$mu_b_inter_ws3, c(0.05)),1)
mu_b_inter_ws395 <- round(quantile(post$mu_b_inter_ws3, c(0.95)),1)
mu_b_inter_ws325 <- round(quantile(post$mu_b_inter_ws3, c(0.25)),1)
mu_b_inter_ws375 <- round(quantile(post$mu_b_inter_ws3, c(0.75)),1)
mu_b_inter_ws3 <- (cbind(mu_b_inter_ws3, mu_b_inter_ws35,mu_b_inter_ws395, mu_b_inter_ws325,mu_b_inter_ws375))

mu_b_inter_s3c1 = sum[grep("mu_b_inter_s3c1", rownames(sum)), c("mean")]
mu_b_inter_s3c15 <- round(quantile(post$mu_b_inter_s3c1, c(0.05)),1)
mu_b_inter_s3c195 <- round(quantile(post$mu_b_inter_s3c1, c(0.95)),1)
mu_b_inter_s3c125 <- round(quantile(post$mu_b_inter_s3c1, c(0.25)),1)
mu_b_inter_s3c175 <- round(quantile(post$mu_b_inter_s3c1, c(0.75)),1)
mu_b_inter_s3c1 <- (cbind(mu_b_inter_s3c1, mu_b_inter_s3c15,mu_b_inter_s3c195, mu_b_inter_s3c125,mu_b_inter_s3c175))

mu_b_inter_ps3 = sum[grep("mu_b_inter_ps3", rownames(sum)), c("mean")]
mu_b_inter_ps35 <- round(quantile(post$mu_b_inter_ps3, c(0.05)),1)
mu_b_inter_ps395 <- round(quantile(post$mu_b_inter_ps3, c(0.95)),1)
mu_b_inter_ps325 <- round(quantile(post$mu_b_inter_ps3, c(0.25)),1)
mu_b_inter_ps375 <- round(quantile(post$mu_b_inter_ps3, c(0.75)),1)
mu_b_inter_ps3 <- (cbind(mu_b_inter_ps3, mu_b_inter_ps35,mu_b_inter_ps395, mu_b_inter_ps325,mu_b_inter_ps375))

mu_b_inter_ws4 = sum[grep("mu_b_inter_ws4", rownames(sum)), c("mean")]
mu_b_inter_ws45 <- round(quantile(post$mu_b_inter_ws4, c(0.05)),1)
mu_b_inter_ws495 <- round(quantile(post$mu_b_inter_ws4, c(0.95)),1)
mu_b_inter_ws425 <- round(quantile(post$mu_b_inter_ws4, c(0.25)),1)
mu_b_inter_ws475 <- round(quantile(post$mu_b_inter_ws4, c(0.75)),1)
mu_b_inter_ws4 <- (cbind(mu_b_inter_ws4, mu_b_inter_ws45,mu_b_inter_ws495, mu_b_inter_ws425,mu_b_inter_ws475))

mu_b_inter_s4c1 = sum[grep("mu_b_inter_s4c1", rownames(sum)), c("mean")]
mu_b_inter_s4c15 <- round(quantile(post$mu_b_inter_s4c1, c(0.05)),1)
mu_b_inter_s4c195 <- round(quantile(post$mu_b_inter_s4c1, c(0.95)),1)
mu_b_inter_s4c125 <- round(quantile(post$mu_b_inter_s4c1, c(0.25)),1)
mu_b_inter_s4c175 <- round(quantile(post$mu_b_inter_s4c1, c(0.75)),1)
mu_b_inter_s4c1 <- (cbind(mu_b_inter_s4c1, mu_b_inter_s4c15,mu_b_inter_s4c195, mu_b_inter_s4c125,mu_b_inter_s4c175))

mu_b_inter_ps4 = sum[grep("mu_b_inter_ps4", rownames(sum)), c("mean")]
mu_b_inter_ps45 <- round(quantile(post$mu_b_inter_ps4, c(0.05)),1)
mu_b_inter_ps495 <- round(quantile(post$mu_b_inter_ps4, c(0.95)),1)
mu_b_inter_ps425 <- round(quantile(post$mu_b_inter_ps4, c(0.25)),1)
mu_b_inter_ps475 <- round(quantile(post$mu_b_inter_ps4, c(0.75)),1)
mu_b_inter_ps4 <- (cbind(mu_b_inter_ps4, mu_b_inter_ps45,mu_b_inter_ps495, mu_b_inter_ps425,mu_b_inter_ps475))

b_site2 = sum[grep("b_site2", rownames(sum)), c("mean")]
b_site25 <- round(quantile(post$b_site2, c(0.05)),1)
b_site295 <- round(quantile(post$b_site2, c(0.95)),1)
b_site225 <- round(quantile(post$b_site2, c(0.25)),1)
b_site275 <- round(quantile(post$b_site2, c(0.75)),1)
b_site2 <- (cbind(b_site2, b_site25,b_site295, b_site225,b_site275))

b_site3 = sum[grep("b_site3", rownames(sum)), c("mean")]
b_site35 <- round(quantile(post$b_site3, c(0.05)),1)
b_site395 <- round(quantile(post$b_site3, c(0.95)),1)
b_site325 <- round(quantile(post$b_site3, c(0.25)),1)
b_site375 <- round(quantile(post$b_site3, c(0.75)),1)
b_site3 <- (cbind(b_site3, b_site35,b_site395, b_site325,b_site375))

b_site4 = sum[grep("b_site4", rownames(sum)), c("mean")]
b_site45 <- round(quantile(post$b_site4, c(0.05)),1)
b_site495 <- round(quantile(post$b_site4, c(0.95)),1)
b_site425 <- round(quantile(post$b_site4, c(0.25)),1)
b_site475 <- round(quantile(post$b_site4, c(0.75)),1)
b_site4 <- (cbind(b_site4, b_site45,b_site495, b_site425,b_site475))

b_warm = sum[grep("b_warm\\[", rownames(sum)), c("mean")]
mu_b_warm5 <- round(quantile(post$mu_b_warm, c(0.05)),1)
mu_b_warm95 <- round(quantile(post$mu_b_warm, c(0.95)),1)
mu_b_warm25 <- round(quantile(post$mu_b_warm, c(0.25)),1)
mu_b_warm75 <- round(quantile(post$mu_b_warm, c(0.75)),1)
mu_b_warm <- (cbind(mu_b_warm, mu_b_warm5,mu_b_warm95, mu_b_warm25,mu_b_warm75))

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# Get the means for the lowest F, C, P and each site
#

lowForce <- c(-0.3568628,-0.6723508)
lowPhoto <- -0.5033863
lowChill <- c(-0.3546922,-0.3494652,-0.2935465,-0.731798)


site2 <- unique(pheno$site2.z2)
site3 <- unique(pheno$site3.z2)
site4 <- unique(pheno$site4.z2)


# St. Hippolyte plot first for the high forcing - site = 1 is the second value (site4[2])
bb_site4 = a_sp[1:5]   + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[2]  + mu_b_warm[1:5] * lowForce[2] + mu_b_photo[1:5] * lowPhoto + mu_b_chill1[1:5] * lowChill[4] +
  mu_b_inter_wp[1:5] * (lowForce[2]*lowPhoto) +
  mu_b_inter_wc1[1:5] * (lowForce[2]*lowChill[4]) + mu_b_inter_pc1[1:5] * (lowPhoto*lowChill[4]) +
  mu_b_inter_s2c1[1:5] * (lowChill[4]*site2[2]) + mu_b_inter_ws2[1:5] * (lowForce[2]*site2[2]) +mu_b_inter_ps2[1:5] * (lowPhoto*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lowChill[4]*site3[1]) + mu_b_inter_ws3[1:5] * (lowForce[2]*site3[1]) +mu_b_inter_ps3[1:5] * (lowPhoto*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lowChill[4]*site4[2]) + mu_b_inter_ws4[1:5] * (lowForce[2]*site4[2]) +mu_b_inter_ps4[1:5] * (lowPhoto*site4[2])

# Manning park trends - site2[1]
bb_site2 = a_sp[1:5] + b_site2[1:5] * site2[1] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * lowForce[1] + mu_b_photo[1:5] * lowPhoto + mu_b_chill1[1:5] * lowChill[2] +
  mu_b_inter_wp[1:5] * (lowForce[1]*lowPhoto) +
  mu_b_inter_wc1[1:5] * (lowForce[1]*lowChill[2]) + mu_b_inter_pc1[1:5] * (lowPhoto*lowChill[2]) +
  mu_b_inter_s2c1[1:5] * (lowChill[2]*site2[1]) + mu_b_inter_ws2[1:5] * (lowForce[1]*site2[1]) +mu_b_inter_ps2[1:5] * (lowPhoto*site2[1]) +
  mu_b_inter_s3c1[1:5] * (lowChill[2]*site3[1]) + mu_b_inter_ws3[1:5] * (lowForce[1]*site3[1]) +mu_b_inter_ps3[1:5] * (lowPhoto*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lowChill[2]*site4[1]) + mu_b_inter_ws4[1:5] * (lowForce[1]*site4[1]) +mu_b_inter_ps4[1:5] * (lowPhoto*site4[1])

# Smithers trends - site2[2]
bb_site1 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * lowForce[1] + mu_b_photo[1:5] * lowPhoto + mu_b_chill1[1:5] * lowChill[1] +
  mu_b_inter_wp[1:5] * (lowForce[1]*lowPhoto) +
  mu_b_inter_wc1[1:5] * (lowForce[1]*lowChill[1]) + mu_b_inter_pc1[1:5] * (lowPhoto*lowChill[1]) +
  mu_b_inter_s2c1[1:5] * (lowChill[1]*site2[2]) + mu_b_inter_ws2[1:5] * (lowForce[1]*site2[2]) +mu_b_inter_ps2[1:5] * (lowPhoto*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lowChill[1]*site3[1]) + mu_b_inter_ws3[1:5] * (lowForce[1]*site3[1]) +mu_b_inter_ps3[1:5] * (lowPhoto*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lowChill[1]*site4[1]) + mu_b_inter_ws4[1:5] * (lowForce[1]*site4[1]) +mu_b_inter_ps4[1:5] * (lowPhoto*site4[1])

# HarvardForest trends
bb_site3 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[2] + b_site4[1:5] * site4[1] + mu_b_warm[1:5]* lowForce[2] + mu_b_photo[1:5] * lowPhoto + mu_b_chill1[1:5] * lowChill[3] +
  mu_b_inter_wp[1:5] * (lowForce[2]*lowPhoto) +
  mu_b_inter_wc1[1:5] * (lowForce[2]*lowChill[3]) + mu_b_inter_pc1[1:5] * (lowPhoto*lowChill[3]) +
  mu_b_inter_s2c1[1:5] * (lowChill[3]*site2[2]) + mu_b_inter_ws2[1:5] * (lowForce[2]*site2[2]) +mu_b_inter_ps2[1:5] * (lowPhoto*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lowChill[3]*site3[2]) + mu_b_inter_ws3[1:5] * (lowForce[2]*site3[2]) +mu_b_inter_ps3[1:5] * (lowPhoto*site3[2]) +
  mu_b_inter_s4c1[1:5] * (lowChill[3]*site4[1]) + mu_b_inter_ws4[1:5] * (lowForce[2]*site4[1]) +mu_b_inter_ps4[1:5] * (lowPhoto*site4[1])



temp <- rbind(bb_site1,bb_site2,bb_site3,bb_site4)
siteForce <- data.frame(temp, site = c("Smithers","Manning park", "Harvard forest","St.Hippolyte"))
names(siteForce) <- c("mean", "lower1", "upper1","lower2","upper2","site")
siteForce <- siteForce[order(siteForce$site),]
siteForce$temp <- rownames(siteForce)

meanWestEsti <- format(round(mean(siteForce[2,1], siteForce[3,1]), 1), nsmall =1)
meanWestEstiU <- format(round(mean(siteForce[2,5], siteForce[3,5]), 1), nsmall =1)
meanWestEstiL <- format(round(mean(siteForce[2,4], siteForce[3,4]), 1), nsmall =1)

meanEastEsti <- format(round(mean(siteForce[1,1], siteForce[4,1]), 1), nsmall =1)
meanEastEstiU <- format(round(mean(siteForce[1,5], siteForce[4,5]), 1), nsmall =1)
meanEastEstiL <- format(round(mean(siteForce[1,4], siteForce[4,4]), 1), nsmall =1)

expLength <- 86

tot <- read.csv("..//input/phenoMini.csv")

#Only want the total number that reached 7:
below <- subset(tot, bbch.t < 8)
totalObs <- nrow(below)

totalDays <- max(tot$day)
# what spp are different bb but similar 