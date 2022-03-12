# Fresh start!! Blank slate!! 

#started March 11, 2022

# the aim of this code is to generate test data for my growth chamber study

# First step: 21 species, 2 sites, 2 photoperiod, 2 forcing, 2 chill portions 
# 
# Next steps: 47 species, 4 sites, 2 photoperiod, 2 forcing, 5 chill portions

# Here we go! Note the comments are more a stream of conciousness so I can possibly figure out where my thought process is going wrong. 

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

set.seed(19902)

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc")
} else if(length(grep("Lizzie", getwd())>0)) {
  setwd("~/Documents/git/teaching/stan")
} else if(length(grep("faith", getwd()) > 0)){
  setwd("/home/faith/Documents/mnt/UBC/otherPeople/Deirdre")
}  else {
  setwd("/home/deirdre/pheno_bc") # for midge
}

library(rstan)
require(shinystan)
require(bayesplot)
require(tidybayes)
require(truncnorm)
library(ggplot2)
library(dplyr)
library(plyr)
require(lme4)

# 1. Create a list of the speices and sites, I am starting really small and then will increase the number to be sure the structure is what I think it is:

#Set number of data groups
nsite = 2 #number of sites with trait data
nsp = 3 #number of species - 
nind = 2
nChill <- 2 #This many chilling treatment levels 
nPhoto <- 2
nForce <- 2
nPhenCombo <- nChill*nPhoto*nForce

sites <- 1:nsite
site <- rep(rep(sites, times = nPhenCombo, each = nPhenCombo*nind))

#Makle a vector of species names 
species <- 1:nsp 
sp <- rep(species, times = nPhenCombo, each = nsite*nPhenCombo*nind)
fake <- data.frame(cbind(site, sp))

fake$site <- as.numeric(fake$site)
fake$site[fake$site==1] <- 0
fake$site[fake$site==2] <- 1

fake$photo <- rep(c(0, 1))
fake$warm <- rep(c(0,1), each = 2)
chill.port <- c(10,20)
#fake$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)
fake$chill <-  rep(rep(chill.port, each = nForce*nPhoto), times = nsp*nsite*nind)

#Set Parameters 
#--------------------------
mu_grand = 30
# b_site = 0 
sitediff2 = 1
# sitediff3 = 2
# sitediff4 = 3
mu_force = -10 
mu_photo = -15
mu_chill = -14

sigma_a = 1
sigma_site = 1
sigma_force = 1
sigma_chill = 1
sigma_photo = 1
sigma_y = 1

# sigma_a = 5
# sigma_site = 5
# sigma_force = 5
# sigma_chill = 5
# sigma_photo = 5
# sigma_y = 5

#interactions:
mu_fp <- 3
mu_cp <- 2
mu_fc <- 1

sigma_fp <- 5
sigma_cp <- 5
sigma_fc <- 5

alpha.pheno.sp <- rnorm(nsp, 0, sigma_a) 
fake$alpha.pheno.sp <- rep(alpha.pheno.sp, each =  nsite*nPhenCombo*nind)

alpha.force.sp <- rnorm(nsp, mu_force, sigma_force)
fake$alpha.force.sp <- rep(alpha.force.sp, each =  nsite*nPhenCombo*nind)

alpha.photo.sp <- rnorm(nsp, mu_photo, sigma_photo)
fake$alpha.photo.sp <- rep(alpha.photo.sp, each = nsite*nPhenCombo*nind)

alpha.chill.sp <- rnorm(nsp, mu_chill, sigma_chill)
fake$alpha.chill.sp <- rep(alpha.chill.sp, each =  nsite*nPhenCombo*nind)

alpha.site <- c(0,1)
fake$alpha.site <- rep(alpha.site,  each = nsite*nind*nPhenCombo)

# adding the interaction effects:
alpha.fc <- rnorm(nsp, mu_fc, sigma_fc)
fake$alpha.fc <- rep(alpha.fc, each =  nsite*nPhenCombo*nind)

alpha.cp <- rnorm(nsp, mu_cp, sigma_cp)
fake$alpha.cp <- rep(alpha.cp, each =  nsite*nPhenCombo*nind)

alpha.fp <- rnorm(nsp, mu_fp, sigma_fp)
fake$alpha.fp <- rep(alpha.fp, each =  nsite*nPhenCombo*nind)

fake$gen.var <- rnorm(nrow(fake), 0, sigma_y) 
#hist(rnorm(nrow(fake), 0, sigma_y))
#Simulate continuous cue values

# Quick sanity check:
fake$count <- 1
faket <- aggregate(fake ["count"],
                   fake[c("sp","site","alpha.pheno.sp", "alpha.force.sp","alpha.chill.sp", "alpha.photo.sp","alpha.site","alpha.fc", "alpha.cp", "alpha.fp","warm","chill","photo")],
                   FUN = sum)
# Ok great, every combination has 2...whic makes sense

# Ok I included the interactions above, but will not add them until I am confident that the test data is working 

fake$bb <-  mu_grand  + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$site  +
  #fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
  fake$gen.var

summary(lm(bb ~  warm + photo + chill + site, data = fake)) # 

datalist <- list( N=nrow(fake),
                  n_sp = length(unique(fake$sp)),
                  n_site = length(unique(fake$site)),
                  bb = fake$bb_int,
                  sp = fake$sp,
                  chill = fake$chill,
                  photo = fake$photo,
                  force = fake$warm,
                  site = fake$site)

mdl.full <- stan("stan/test_model_interactions_truepriors.stan",
                 data = datalist,
                 include = FALSE, pars = c("ypred_new","y_hat"),
                 iter = 4000, chains= 4, warmup = 2000)

