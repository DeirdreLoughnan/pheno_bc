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
library(rethinking)


# if(length(grep("deirdreloughnan", getwd()) > 0)) { 
#   setwd("~/Documents/github/pheno_bc") 
# }  

#load("output/final/ew_phylo_output_newpriors.Rda")
dl <- read.csv("..//input/dl_allbb.csv")

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

pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$site.n <- pheno$population
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
sort(unique(pheno.t$species.fact)) # 49 

# now get the phylogeny and pair it with species names:
spInfo <- read.csv("..//input/species_list.csv")
# head(spInfo)
# head(pheno.t)
# change the df to lowercase:

pheno.t <- merge(pheno.t, spInfo, by = "species")


#################################################################
# load("output/final/ew_phylo_output_newpriors_allncp.Rda")
# sumew <- summary(mdl.ewphylo)$summary

load("..//output/final/bb_4sites_phylo_mini.Rda")
sum <- summary(mdl.4phyloMini)$summary 

fit <- rstan::extract(mdl.4phyloMini)

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
chillCue2.5 <- round(meanz4["Chilling",3],1)
chillCue97.5 <- round(meanz4["Chilling",5],1)

photoCue <- round(meanz4["Photoperiod",1],1)
photoCue2.5 <- round(meanz4["Photoperiod",3],1)
photoCue97.5 <- round(meanz4["Photoperiod",5],1)

intrxnCF <- round(meanz4["Forcing x chilling",1],1)
intrxnCF2.5 <- round(meanz4["Forcing x chilling",3],1)
intrxnCF97.5 <- round(meanz4["Forcing x chilling",5],1)

siteHF <- round(meanz4["Harvard Forest",1],1)
siteHF2.5 <- round(meanz4["Harvard Forest",3],1)
siteHF97.5 <- round(meanz4["Harvard Forest",5],1)

siteSH <- round(meanz4["St. Hippolyte",1],1)
siteSH2.5 <- round(meanz4["St. Hippolyte",3],1)
siteSH97.5 <- round(meanz4["St. Hippolyte",5],1)

siteMP <- round(meanz4["Manning Park",1],1)
siteMP2.5 <- round(meanz4["Manning Park",3],1)
siteMP97.5 <- round(meanz4["Manning Park",5],1)

forceHF <- round(meanz4["Forcing x Harvard Forest",1],1)
forceHF2.5 <- round(meanz4["Forcing x Harvard Forest",3],1)
forceHF97.5 <- round(meanz4["Forcing x Harvard Forest",5],1)

forceSH <- round(meanz4["Forcing x St. Hippolyte",1],1)
forceSH2.5 <- round(meanz4["Forcing x St. Hippolyte",3],1)
forceSH97.5 <- round(meanz4["Forcing x St. Hippolyte",5],1)
forceHF <- round(meanz4["Forcing x Harvard Forest",1],1)
forceHF2.5 <- round(meanz4["Forcing x Harvard Forest",3],1)
forceHF97.5 <- round(meanz4["Forcing x Harvard Forest",5],1)

forceSH <- round(meanz4["Forcing x St. Hippolyte",1],1)
forceSH2.5 <- round(meanz4["Forcing x St. Hippolyte",3],1)
forceSH97.5 <- round(meanz4["Forcing x St. Hippolyte",5],1)

chillHF <- round(meanz4["Chilling x Harvard Forest",1],1)
chillHF2.5 <- round(meanz4["Chilling x Harvard Forest",3],1)
chillHF97.5 <- round(meanz4["Chilling x Harvard Forest",5],1)

chillSH <- round(meanz4["Chilling x St. Hippolyte",1],1)
chillSH2.5 <- round(meanz4["Chilling x St. Hippolyte",3],1)
chillSH97.5 <- round(meanz4["Chilling x St. Hippolyte",5],1)

chillMP <- round(meanz4["Chilling x Manning Park",1],1)
chillMP2.5 <- round(meanz4["Chilling x Manning Park",3],1)
chillMP97.5 <- round(meanz4["Chilling x Manning Park",5],1)

photoHF <- round(meanz4["Photoperiod x Harvard Forest",1],1)
photoHF2.5 <- round(meanz4["Photoperiod x Harvard Forest",3],1)
photoHF97.5 <- round(meanz4["Photoperiod x Harvard Forest",5],1)

photoSH <- round(meanz4["Photoperiod x St. Hippolyte",1],1)
photoSH2.5 <- round(meanz4["Photoperiod x St. Hippolyte",3],1)
photoSH97.5 <- round(meanz4["Photoperiod x St. Hippolyte",5],1)

photoMP <- round(meanz4["Photoperiod x Manning Park",1],1)
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

meanBB <- round(mean(meanBBData$bb),1)
# What are teh lambda values?
lam_params <- c( 
  "a_z",
  "lam_interceptsa"
  )
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")

phylo <- sum[lam_params, col4table]

rootT <- round(phylo[1,1],1)

rootTLower <- as.numeric(round(HPDI(data.frame(fit$a_z), prob = 0.90)[1],1))
rootTUpper <- as.numeric(round(HPDI(data.frame(fit$a_z), prob = 0.90)[2],1))

lamT <- round(phylo[2,1],1)
lamTLower <- as.numeric(round(HPDI(data.frame(fit$lam_interceptsa), prob = 0.90)[1],1))
lamTUpper <- as.numeric(round(HPDI(data.frame(fit$lam_interceptsa), prob = 0.90)[2],1))

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

spInfo <- spInfo[order(spInfo$species),]
head(spInfo)
spInfo$meanBB <- colMeans(m)
colnames(m) <- spInfo$species.name

spInfo$Int <- a_sp
spInfo$Int2.5 <- a_sp2.5
spInfo$Int97.5 <- a_sp97.5

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



east <- subset(spInfo, transect != "west")
eastSp <- unique(east$species.name)

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

meanPtW <- aggregate(dataWest[c("meanBB", "Int")], dataWest[c("species.name","type","transect")], FUN = mean)

meanPtE <- aggregate(dataEast[c("meanBB", "Int")], dataEast[c("species.name","type","transect")], FUN = mean)

diffBAWest <- round(mean((meanPtW$Int/meanPtW$meanBB)*100),1)
diffBAEast <- round(mean((meanPtE$Int/meanPtE$meanBB)*100),1)

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

typeValues <- aggregate(spInfo[c("force","chill","photo","force2.5","chill2.5","photo2.5", "force97.5","chill97.5","photo97.5")], spInfo[c("type")], FUN = mean)

shrubMeanForce <- typeValues[1,2]
shrubMeanForce2.5 <- typeValues[1,5]
shrubMeanForce97.5 <- typeValues[1,8]

shrubMeanChill <- typeValues[1,3]
shrubMeanChill2.5 <- typeValues[1,6]
shrubMeanChill97.5 <- typeValues[1,9]

shrubMeanPhoto <- typeValues[1,4]
shrubMeanPhoto2.5 <- typeValues[1,7]
shrubMeanPhoto97.5 <- typeValues[1,10]

# trees
treeMeanForce <- typeValues[2,2]
treeMeanForce2.5 <- typeValues[2,5]
treeMeanForce97.5 <- typeValues[2,8]

treeMeanChill <- typeValues[2,3]
treeMeanChill2.5 <- typeValues[2,6]
treeMeanChill97.5 <- typeValues[2,9]

treeMeanPhoto <- typeValues[2,4]
treeMeanPhoto2.5 <- typeValues[2,7]
treeMeanPhoto97.5 <- typeValues[2,10]

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
