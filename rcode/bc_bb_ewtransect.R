## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

# Adding standardization: following the methods outlined by Andrew Gelman in his paper

# To properly add site to the model I am going to use the indexing approach used on pg 153 of statistical rethinking, modeling code was drafted by Lizzie

#library(scales)
#library(arm)
library(rstan)
library(shinystan)
#library(reshape2)
library(bayesplot)
library(ggplot2)
#library(RColorBrewer)
library(dplyr)
library(plyr)
library(tidybayes)

options(mc.cores = parallel::detectCores())

rm(list=ls()) 
options(stringsAsFactors = FALSE)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

#source('rcode/cleaning/pheno_bb_calc.R')
# head(pheno)
# length(unique(pheno$lab2))

df <- read.csv("input/day.of.bb.DFlynn.chill0.csv", header=TRUE, na.strings=c("","NA"))
head(df)
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill$transect <- "east"

dl <- read.csv("input/day.of.bb.DL.csv", header=TRUE, na.strings=c("","NA"))
head(dl)
dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")
dl.wchill$transect <- "west"

# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)
#pheno <- dl

pheno$first <- ifelse(pheno$tbb < pheno$latbb1,"t", ifelse (pheno$tbb == pheno$latbb1,"tl", "l"))

head(pheno)
table(pheno$species, pheno$first)
#write.csv(pheno, "input/pheno.w5chill.csv")
#pheno <- read.csv("pheno.w5chill.csv")
############################################################
#convert forcing and photoperiod treatments into binary
pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

# convert the transect into a dummy variable:
pheno$transect.n <- pheno$transect
pheno$transect.n[pheno$transect.n == "east"] <- "1"
pheno$transect.n[pheno$transect.n == "west"] <- "0"
pheno$transect.n <- as.numeric(pheno$transect.n)


pheno$Chill_portions <- as.numeric(pheno$Chill_portions)

# pheno$force.z <- (pheno$force.t-mean(pheno$force.t,na.rm=TRUE))/sd(pheno$force.t,na.rm=TRUE)
# pheno$photo.z <- (pheno$photo.t-mean(pheno$photo.t,na.rm=TRUE))/sd(pheno$photo.t,na.rm=TRUE)
# pheno$chillport.z <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/sd(pheno$Chill_portions,na.rm=TRUE)
# z scoring values, but since some are binary I will use 2SD as per the Gelmen et al paper
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)
pheno$utah.z2 <- (pheno$Utah_Model-mean(pheno$Utah_Model,na.rm=TRUE))/(sd(pheno$Utah_Model,na.rm=TRUE)*2)

#going to split it into analyses of terminal bb and lateral bb
# Starting with the terminal buds:
pheno.term <- pheno[,c("tbb", "force.n", "photo.n", "transect.n", "species", "lab2","Utah_Model","Chill_portions","force.z2", "photo.z2", "chillport.z2", "utah.z2")]

pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 

nrow(pheno.term) - nrow(pheno.t) # That had no terminal bb, 609


# head(pheno.t)
datalist_z <- list( N=nrow(pheno.t),
                    Nsite = nrow(pheno.t),
                    n_sp = length(unique(pheno.t$species.fact)),
                    n_site = length(unique(pheno.t$site.n)),
                    bb = pheno.t$tbb,
                    sp = pheno.t$species.fact,
                    chill = pheno.t$chillport.z2,
                    photo = pheno.t$photo.z2,
                    force = pheno.t$force.z2,
                    site = pheno.t$transect.n)

mdl <- stan("stan/bc.bb.ncpphoto.ncpinter.standardize.old.stan",
              data = datalist_z,
            include = FALSE, pars = c("ypred_new","y_hat"),
              iter = 4000, warmup = 3000, chains=4, control = list(adapt_delta = 0.99))
save(mdl, file="output/tbb_cport_stnd_transect.Rda")

# mdl.i <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_index.stan",
#               data = datalist_z,
#               iter = 4000, chains=4)
# save(mdl.i, file="output/tbb_cport_stnd_index.Rds")
# 
# 
# mdl.d <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_dummy.stan",
#               data = datalist_z,
#               iter = 4000, chains=4)
# save(mdl.d, file="output/tbb_cport_stnd_dummy.Rds")
#######################################################################

load("output/tbb_cport_stnd_transect.Rda")
# #  load("output/tbb_ncp_chillportions_zsc_dl.Rda")
# # #
ssm <-  as.shinystan(mdl)
launch_shinystan(ssm)
# # # 
sum <- summary(mdl)$summary 
summary(mdl)$summary[c("mu_a", "mu_force","mu_chill","mu_photo","mu_site","mu_inter_fp", "mu_inter_fs","mu_inter_ps","mu_inter_fc","mu_inter_pc","mu_inter_sc","sigma_a", "sigma_force","sigma_chill","sigma_photo","sigma_site","sigma_b_inter_fp","sigma_b_inter_fs","sigma_b_inter_ps","sigma_b_inter_fc","sigma_b_inter_pc","sigma_b_inter_sc", "sigma_y"),"mean"]

