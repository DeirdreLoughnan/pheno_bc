# Started November 2021 by DL

# working on test data for phenology model 
# testing for differences between dummy variables and indexing code

# the test data and model seem to being doing well as of Jan3
# Now I will try to add back in more complexity:
#   1. 
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

set.seed(194842)

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

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
#building the required dataframe 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

#Set number of data groups
nsite = 4 #number of sites with trait data
nsp = 8 #number of species - 
nind = 4 
nChill <- 6 #This many chilling treatment levels 
nPhoto <- 2
nForce <- 2
nPhenCombo <- nChill*nPhoto*nForce

#Set Parameters 
#--------------------------
mu_grand = 50
b_site = 0 
sitediff2 = 1
sitediff3 = 2
sitediff4 = 3
mu_force = -10 
mu_photo = -15
mu_chill = -14

sigma_a = 3
sigma_site = 10
sigma_force = 1
sigma_chill = 1
sigma_photo = 1
sigma_y = 3

#interactions:
mu_fp <- 3
mu_cp <- 2
mu_fc <- 1

sigma_fp <- 0.5
sigma_cp <- 0.5
sigma_fc <- 0.5

mu_fsite <- 1
mu_csite <- 5
mu_psite  <- 3

sigma_fsite <- 0.5
sigma_csite <- 0.5
sigma_psite <- 0.5

#Make a vector of site names
sites <- 1:nsite
site <- rep(rep(sites, times = nPhenCombo, each = nPhenCombo*nind))

#Makle a vector of species names 
species <- 1:nsp 
sp <- rep(species, each = nsite*nPhenCombo*nind)

#Make a data frame with ass sites and species laid out 
#---------------------

fake <- data.frame(cbind(site, sp))
head(fake)

# replicating the method for calculating the response variable used in traitors (which I understand the best)
alpha.pheno.sp <- rnorm(nsp, 0, sigma_a) 
fake$alpha.pheno.sp <- rep(alpha.pheno.sp, each =  nsite*nPhenCombo*nind)

alpha.force.sp <- rnorm(nsp, mu_force, sigma_force)
fake$alpha.force.sp <- rep(alpha.force.sp, each =  nsite*nPhenCombo*nind)

alpha.photo.sp <- rnorm(nsp, mu_photo, sigma_photo)
fake$alpha.photo.sp <- rep(alpha.photo.sp, each = nsite*nPhenCombo*nind)

alpha.chill.sp <- rnorm(nsp, mu_chill, sigma_chill)
fake$alpha.chill.sp <- rep(alpha.chill.sp, each =  nsite*nPhenCombo*nind)

alpha.site <- rnorm(nsite, b_site, sigma_site)
fake$alpha.site <- rep(rep(alpha.site, times = nPhenCombo, each = nPhenCombo*nind))

# adding the interaction effects:
alpha.fc <- rnorm(nsp, mu_fc, sigma_fc)
fake$alpha.fc <- rep(alpha.fc, each =  nsite*nForce * nPhoto * nind)

alpha.cp <- rnorm(nsp, mu_cp, sigma_cp)
fake$alpha.cp <- rep(alpha.cp, each =  nsite*nChill * nPhoto *nind)

alpha.fp <- rnorm(nsp, mu_fp, sigma_fp)
fake$alpha.fp <- rep(alpha.fp, each =  nsite*nForce * nPhoto *nind)

#adding interactions between the cues and the 
alpha.fsite <- rnorm(nsite, mu_fsite, sigma_fsite)
fake$alpha.fsite <- rep(alpha.fsite, each =  nPhenCombo*nind)

alpha.csite <- rnorm(nsite, mu_csite, sigma_csite)
fake$alpha.csite <- rep(alpha.csite, each =  nPhenCombo*nind)

alpha.psite <- rnorm(nsite, mu_psite, sigma_psite)
fake$alpha.psite <- rep(alpha.psite, each =  nPhenCombo*nind)

# add dummy/ site level effects:
fake <- fake %>%
  mutate ( d2 = if_else(site == 2, 1, 0),
           d3 = if_else(site == 3, 1, 0),
           d4 = if_else(site == 4, 1, 0))

#general variance

fake$gen.var <- rnorm(nrow(fake), 0, sigma_y) 

#Simulate continuous cue values

fake$photo <- rep(c(0, 1))
fake$warm <- rep(c(0,1), each = 2)
fake$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)


# for dummy variable test data
fake$bb <-  mu_grand + alpha.site[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
  fake$gen.var

fake$bb_int <-  mu_grand + alpha.site[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 + fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) +
  fake$gen.var

# check if works with lmer or brms

summary(lm(bb ~  warm + photo + chill + d2 + d3 + d4, data = fake)) # 

summary(lm(bb_int ~  warm + photo + chill + d2 + d3 + d4 + warm * chill + chill * photo + photo * warm, data = fake)) # 

mu_grand + alpha.site[1]
alpha.site
# The values are very close between trail data and lmer output

  # tbb, f/c/p.n, site.n, species, f/c/p.i, species.fact, d2, d3, d4
datalist <- list( N=nrow(fake),
                    n_sp = length(unique(fake$sp)),
                    n_site = length(unique(fake$site)),
                    bb = fake$bb,
                    sp = fake$sp,
                    chill = fake$chill,
                    photo = fake$photo,
                    force = fake$warm,
                    site = fake$site, 
                    site2 = fake$d2,
                    site3 = fake$d3,
                    site4 = fake$d4)
  
mdl.dum <- stan("stan/test_model_interactions.stan",
                   data = datalist,
                   include = FALSE, pars = c("ypred_new","y_hat"),
                   iter = 4000, chains= 4, warmup = 2000)
save(mdl.dum, file = "output/BBDummy_int.Rda")
  
# mdl.dumint <- stan("stan/test_model_interaction_ncp.stan",
#                    data = datalist,
#                    include = FALSE, pars = c("ypred_new","y_hat"),
#                    iter = 4000, chains= 4, warmup = 2000)
# save(mdl.dumint, file = "output/BBDummy_ncpint_test.Rda")
###################################################################
load("output/BBDummy_int_test.Rda")
# # 
# ssm <-  as.shinystan(mdl.dum)
# launch_shinystan(ssm)
# #
# get_variables(mdl.dum)
# summary(mdl.dum)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","b_site2","b_site3","b_site4","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y"),"mean"]
#    
# load("output/BBDummy_int_test.Rda")
#    
summary(mdl.dumint)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","b_site2","b_site3","b_site4","mu_fc" ,"mu_fp" ,"mu_cp","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y", "sigma_fp" , "sigma_cp" , "sigma_fc", "sigma_fp", "sigma_cp"),"mean"]

 
# summary(mdl.dum)$summary[c("mu_grand","mu_a","mu_force", "mu_chill","mu_photo","b_site2","b_site3","b_site4","mu_b_inter_fp" ,"mu_b_inter_pc", "mu_b_inter_fc","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y","sigma_b_inter_fc","sigma_b_inter_pc","sigma_b_inter_fp"),"mean"]

# mu_grand + alpha.site[1]
# alpha.site
  
  # pairs(mdl.simpdum, pars = c("mu_grand", "mu_force", "mu_chill", "mu_photo","b_site2","b_site3","b_site4","sigma_a", "sigma_force", "sigma_chill", "sigma_force","sigma_y", "lp__")) 
  #dev.off()
  
