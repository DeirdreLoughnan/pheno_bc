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
sigma_y = 4

#interactions:
mu_fp <- 3
mu_cp <- 2
mu_fc <- 1

sigma_fp <- 0.2
sigma_cp <- 0.2
sigma_fc <- 0.2

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
#hist(rnorm(nrow(fake), 0, sigma_y))
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

fake$bb_int_site <-  mu_grand + alpha.site[1] + alpha.psite[1] + alpha.csite[1] + alpha.fsite[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 + fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) + fake$alpha.fsite * (fake$d2*fake$warm) + fake$alpha.fsite * (fake$d3*fake$warm) + fake$alpha.fsite * (fake$d4*fake$warm) + fake$alpha.csite * (fake$d2*fake$chill) + fake$alpha.csite * (fake$d3*fake$chill) + fake$alpha.csite * (fake$d4*fake$chill) + fake$alpha.psite * (fake$d2*fake$photo) + fake$alpha.psite * (fake$d3*fake$photo) + fake$alpha.psite * (fake$d4*fake$photo) + fake$gen.var
# estimates don't look too bad, worse for photo but might need to be ncp?

# check if works with lmer or brms

summary(lm(bb ~  warm + photo + chill + d2 + d3 + d4, data = fake)) # 

summary(lm(bb_int ~  warm + photo + chill + d2 + d3 + d4 + warm * chill + chill * photo + photo * warm, data = fake)) # 

summary(lm(bb_int_site ~  warm + photo + chill + d2 + d3 + d4 + warm * chill + chill * photo + photo * warm + d2 * warm + d2 * chill + d2 * photo + d3 * warm + d3 * chill + d3 * photo + d4 * warm + d4 * chill + d4 * photo , data = fake))

mu_grand + alpha.site[1]
alpha.site
# The values are very close between trail data and lmer output

  # tbb, f/c/p.n, site.n, species, f/c/p.i, species.fact, d2, d3, d4
datalist <- list( N=nrow(fake),
                    n_sp = length(unique(fake$sp)),
                    n_site = length(unique(fake$site)),
                    bb = fake$bb_int,
                    sp = fake$sp,
                    chill = fake$chill,
                    photo = fake$photo,
                    force = fake$warm,
                    site = fake$site, 
                    site2 = fake$d2,
                    site3 = fake$d3,
                    site4 = fake$d4)
  
mdl.full <- stan("stan/test_model_interactions_trueprior_tightprior.stan",
                   data = datalist,
                   include = FALSE, pars = c("ypred_new","y_hat"),
                   iter = 4000, chains= 4, warmup = 2000)
save(mdl.full, file = "output/BBDummy_int_tightprior.Rda")

# fix: issues with sigma_y, sigma_fc, sigma_cp

# mdl.dumint <- stan("stan/test_model_interaction_ncp.stan",
#                    data = datalist,
#                    include = FALSE, pars = c("ypred_new","y_hat"),
#                    iter = 4000, chains= 4, warmup = 2000)
# save(mdl.dumint, file = "output/BBDummy_ncpint_test.Rda")
###################################################################
#  load("output/BBDummy_int.Rda")
# # # # #
# ssm <-  as.shinystan(mdl.full)
# launch_shinystan(ssm)
# # # sigma_fp looks really bad, of all the sigma it is the one that is most squished up to one side
# 
# get_variables(mdl.full)
summary(mdl.full)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","b_site2","b_site3","b_site4","mu_fc" ,"mu_fp" ,"mu_cp","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y", "sigma_fp" , "sigma_cp" , "sigma_fc"),"mean"]

# summary(mdl.full)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","b_site2","b_site3","b_site4","mu_inter_fc" ,"mu_inter_fp" ,"mu_inter_pc","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y", "sigma_b_inter_fp" , "sigma_b_inter_pc" , "sigma_b_inter_fc"),"mean"]
#
# pairs(mdl.full, pars = c("sigma_a", "sigma_force", "sigma_chill", "sigma_photo","sigma_y","sigma_fp" , "sigma_cp" , "sigma_fc", "sigma_fp", "sigma_cp", "lp__"))
# #dev.off()
# # 
# #ext<-rstan::extract(mdl.full)
# # 
# # 
# # 
# # #Plot model posteriors and simulated values
# # #-------------------------------------------------
# # 
# par(mfrow=c(1,1))
# 
# plot(hist(ext$mu_grand))
# abline(v=mu_grand + alpha.site[1],col="red")
# 
# plot(hist(ext$mu_force))
# abline(v=mu_force,col="red")
# 
# plot(hist(ext$mu_chill))
# abline(v=mu_chill,col="red")
# 
# plot(hist(ext$mu_photo))
# abline(v=mu_photo,col="red")
# 
# plot(hist(ext$mu_fc))
# abline(v = mu_fc, col = "red")
# 
# plot(hist(ext$mu_cp))
# abline(v = mu_cp, col = "red")
# 
# plot(hist(ext$mu_fp))
# abline(v = mu_fp, col = "red")
# 
# plot(hist(ext$sigma_y), xlim = c(0,6))
# abline(v = sigma_y,col="red")
# 
# plot(hist(ext$sigma_force))
# abline(v = sigma_force,col="red")
# 
# plot(hist(ext$sigma_chill))
# abline(v = sigma_chill,col="red")
# 
# plot(hist(ext$sigma_photo))
# abline(v = sigma_photo,col="red")
# 
# plot(hist(ext$b_site2))
# abline(v = alpha.site[2],col="red")
# 
# plot(hist(ext$b_site3))
# abline(v = alpha.site[3],col="red")
# 
# plot(hist(ext$b_site4))
# abline(v = alpha.site[4],col="red")
# 
# plot(hist(ext$sigma_fp))
# abline(v = sigma_fp,col="red")
# 
# plot(hist(ext$sigma_fc))
# abline(v = sigma_fc,col="red")
# 
# plot(hist(ext$sigma_cp))
# abline(v = sigma_cp,col="red")

#plot posterior expected values againt real values

# ext<-rstan::extract(mdl.full)
# y <- fake$bb
# y.ext <- ext$ypred_new
# pdf("output/yvsypred.pdf", width =3, height = 3)
# ppc_dens_overlay(y, y.ext[1:50, ])
# dev.off()
# #------------------------------------------------------
# 
# 
# 
