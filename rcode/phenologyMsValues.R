# Started May 3, 2022 by deirde

# the aim of this code is to generate the model output for my phenology ms
# rm(list=ls()) 
# options(stringsAsFactors = FALSE)

library(rstan)
library(shinystan)
#library(reshape2)
#library(bayesplot)
library(ggplot2)
library(dplyr)
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
pheno <- pheno %>%
  mutate ( site2 = if_else(site.n == 2, 1, 0),
           site3 = if_else(site.n == 3, 1, 0),
           site4 = if_else(site.n == 4, 1, 0))

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
head(spInfo)
head(pheno.t)
# change the df to lowercase:

pheno.t <- merge(pheno.t, spInfo, by = "species")

#################################################################
# load("output/final/ew_phylo_output_newpriors_allncp.Rda")
# sumew <- summary(mdl.ewphylo)$summary

load("..//output/final/bb_4sites_phylo.Rda")
sum <- summary(mdl.4phylo)$summary 

fit <- rstan::extract(mdl.4phylo)

load("..//output/final/dl_phylo_lat1.Rda")
load("..//output/final/dl_phylo_lat50.Rda")
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

meanz4 <- sum[mu_params_4, col4fig]

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

# General survival values:

# Total number samples that went into exp:
source('..//rcode/cleaning/pheno_bb_calc.R')

totChill <- length(unique(allbb$lab2))

end <- subset(d, day == 88)
totExpSurv <- length(unique(end$lab2))

perMort <- round((totExpSurv-totChill)/(totExpSurv)*100, 1)

nobb <- length(unique(pheno$lab2))

perNoBB <- round(((totExpSurv-nobb)/totExpSurv)*100,1)

noTermEnd <- subset(end, bbch.t <7)
noTerm <- length(unique(noTermEnd$lab2))

perNoTerm <- round(((noTerm/totExpSurv)*100),1)

noTermEnd$count <- 1
temp <- aggregate(noTermEnd["count"], noTermEnd[c("species")], FUN = sum)

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

#### How different are the lateral vs 50 lat estimates? 
sum1 <- summary(mdlLat1)$summary

col4fig <- c("mean","sd","25%","50%","75%","n_eff","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")

mu_params_1l <- c(
  "mu_b_warm", "mu_b_photo", "mu_b_chill1", "b_site","mu_b_inter_wp",
  "mu_b_inter_wc1","mu_b_inter_pc1", "mu_b_inter_ws",
  "mu_b_inter_ps","mu_b_inter_sc1")

meanz1l <- sum1[mu_params_1l, col4fig]

rownames(meanz1l) = c( #"Root trait intercept","Lambda",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Manning Park",
  "Forcing x photoperiod",
  "Forcing x chilling",
  "Photoperiod x chilling",
  "Forcing x Manning Park",
  "Photoperiod x Manning Park",
  "Chilling x Manning Park"
  
  
)


sum50 <- summary(mdlLat50)$summary

col4fig <- c("mean","sd","25%","50%","75%","n_eff","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")

mu_params_50l <- c(
  #"a_z","lam_interceptsa",
  "mu_b_warm", "mu_b_photo", "mu_b_chill1", "b_site", "mu_b_inter_wp",
  "mu_b_inter_wc1","mu_b_inter_pc1", "mu_b_inter_ws",
  "mu_b_inter_ps","mu_b_inter_sc1")

meanz50l <- sum50[mu_params_50l, col4fig]

rownames(meanz50l) = c( #"Root trait intercept","Lambda",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Manning Park",
  "Forcing x photoperiod",
  "Forcing x chilling",
  "Photoperiod x chilling",
  "Forcing x Manning Park",
  "Photoperiod x Manning Park",
  "Chilling x Manning Park"
  
)

(meanz50l[,1] - meanz1l[,1])/meanz50l[,1]


