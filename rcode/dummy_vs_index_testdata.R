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
alpha.fsite2 <- rnorm(nsite, mu_fsite, sigma_fsite)
fake$alpha.fsite2 <- rep(alpha.fsite2, each =  nPhenCombo*nind)

alpha.fsite3 <- rnorm(nsite, mu_fsite, sigma_fsite)
fake$alpha.fsite3 <- rep(alpha.fsite3, each =  nPhenCombo*nind)

alpha.fsite4 <- rnorm(nsite, mu_fsite, sigma_fsite)
fake$alpha.fsite4 <- rep(alpha.fsite4, each =  nPhenCombo*nind)

alpha.csite2 <- rnorm(nsite, mu_csite, sigma_csite)
fake$alpha.csite2 <- rep(alpha.csite2, each =  nPhenCombo*nind)

alpha.csite3 <- rnorm(nsite, mu_csite, sigma_csite)
fake$alpha.csite3 <- rep(alpha.csite3, each = nPhenCombo*nind)

alpha.csite4 <- rnorm(nsite, mu_csite, sigma_csite)
fake$alpha.csite4 <- rep(alpha.csite4, each =  nPhenCombo*nind)

alpha.psite2 <- rnorm(nsite, mu_psite, sigma_psite)
fake$alpha.psite2 <- rep(alpha.psite2, each =  nPhenCombo*nind)

alpha.psite3 <- rnorm(nsite, mu_psite, sigma_psite)
fake$alpha.psite3 <- rep(alpha.psite3, each =  nPhenCombo*nind)

alpha.psite4 <- rnorm(nsite, mu_psite, sigma_psite)
fake$alpha.psite4 <- rep(alpha.psite4, each =  nPhenCombo*nind)

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

fake$bb_int_site <-  mu_grand + alpha.site[1] + alpha.csite2[1]+ fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
  fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) +
  fake$alpha.fsite2 * (fake$d2 * fake$warm) + fake$alpha.fsite2 * (fake$d3 * fake$warm) + fake$alpha.fsite2 * (fake$d4 * fake$warm) +
  fake$alpha.csite2 * (fake$d2 * fake$chill) + fake$alpha.csite2 * (fake$d3 * fake$chill) + fake$alpha.csite2 * (fake$d4 * fake$chill) +
  fake$alpha.psite2 * (fake$d2 * fake$photo) + fake$alpha.psite2 * (fake$d3 * fake$photo) + fake$alpha.psite2 * (fake$d4 * fake$photo) +
  fake$gen.var

fake$bb_int_site <-  mu_grand + alpha.site[1] + alpha.csite2[1]+ fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
  fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) +
  fake$alpha.fsite2 * (fake$d2 * fake$warm) + fake$alpha.fsite2 * (fake$d3 * fake$warm) + fake$alpha.fsite2 * (fake$d4 * fake$warm) +
  fake$alpha.csite2 * (fake$d2 * fake$chill) + fake$alpha.csite2 * (fake$d3 * fake$chill) + fake$alpha.csite2 * (fake$d4 * fake$chill) +
  fake$alpha.psite2 * (fake$d2 * fake$photo) + fake$alpha.psite2 * (fake$d3 * fake$photo) + fake$alpha.psite2 * (fake$d4 * fake$photo) +
  fake$gen.var

fake$bb_int_site2 <-  mu_grand + alpha.site[1] + alpha.csite3[1] + alpha.csite4[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm +
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo +
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
  fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) +
  fake$alpha.fsite2 * (fake$d2 * fake$warm) + fake$alpha.fsite3 * (fake$d3 * fake$warm) + fake$alpha.fsite4 * (fake$d4 * fake$warm) +
  fake$alpha.csite2 * (fake$d2 * fake$chill) + fake$alpha.csite3 * (fake$d3 * fake$chill) + fake$alpha.csite4 * (fake$d4 * fake$chill) +
  fake$alpha.psite2 * (fake$d2 * fake$photo) + fake$alpha.psite3 * (fake$d3 * fake$photo) + fake$alpha.psite4 * (fake$d4 * fake$photo) +
  fake$gen.var

# check if works with lmer or brms

summary(lm(bb ~  warm + photo + chill + d2 + d3 + d4, data = fake)) # 

summary(lm(bb_int ~  warm + photo + chill + d2 + d3 + d4 + warm * chill + chill * photo + photo * warm, data = fake)) # 

summary(lm(bb_int_site ~  warm + photo + chill + d2 + d3 + d4 +
             warm * chill + chill * photo + photo*warm +
             d2 * warm + d3 * warm + d4 * warm + d2 * chill + d3 * chill + d4 * chill + d2 * photo + d3 * photo + d4 * photo
             , data = fake)) # 

summary(lm(bb_int_site2 ~  warm + photo + chill + d2 + d3 + d4 +
             warm * chill + chill * photo + photo*warm +
             d2 * warm + d3 * warm + d4 * warm + d2 * chill + d3 * chill + d4 * chill + d2 * photo + d3 * photo + d4 * photo
           , data = fake)) # 



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
  
  mdl.simpdum <- stan("stan/test_model_interactions_simple.stan",
                      data = datalist,
                      # include = FALSE, pars = c("ypred_new","y_hat"),
                      iter = 4000, chains= 4, warmup = 2000)
  
  #save(mdl.simpdum, file = "output/simBBPosterior.Ddata")
  save(mdl.simpdum, file = "output/simBBDummy_interactions.Rda")



#You should make a list using your original values in lines 35-83 rather than redefine values here. It saves confusion. 
#param <- list(mu_grand = 50, mu_force = -20,
#              mu_chill = -20,  mu_photo = -14, b_site = 10, sigma_a = 5,
#              sigma_force = 1,sigma_chill = 1, sigma_photo =1,  sigma_y = 5, site2 = alpha.site[2], site3 = alpha.site[3], site4 = alpha.site[4])
  
  # load("output/simBBDummy_interactions.Rda")
  # 
  # 
  
  summary(mdl.simpdum)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","b_site2","b_site3","b_site4","mu_inter_fc" ,"mu_inter_fp" ,"mu_inter_pc","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y", "sigma_b_inter_fp" , "sigma_b_inter_pc" , "sigma_b_inter_fc", "sigma_b_inter_fp", "sigma_b_inter_pc"),"mean"]
  
  mu_grand + alpha.site[1]
  alpha.site
  
  # pairs(mdl.simpdum, pars = c("mu_grand", "mu_force", "mu_chill", "mu_photo","b_site2","b_site3","b_site4","sigma_a", "sigma_force", "sigma_chill", "sigma_force","sigma_y", "lp__")) 
  #dev.off()
  
  # ext<-rstan::extract(mdl.simpdum)
  # 
  # 
  # 
  # #Plot model posteriors and simulated values
  # #-------------------------------------------------
  # 
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
  # plot(hist(ext$sigma_y))
  # abline(v=sigma_y,col="red")
  # 
  # plot(hist(ext$sigma_force))
  # abline(v=sigma_force,col="red")
  # 
  # plot(hist(ext$sigma_chill))
  # abline(v=sigma_chill,col="red")
  # 
  # plot(hist(ext$sigma_photo))
  # abline(v=sigma_photo,col="red")
  # 
  # plot(hist(ext$b_site2))
  # abline(v=alpha.site[2],col="red")
  # 
  # plot(hist(ext$b_site3))
  # abline(v=alpha.site[3],col="red")
  # 
  # plot(hist(ext$b_site4))
  # abline(v=alpha.site[4],col="red")
  # 
  # #plot posterior expected values againt real values
  # load("output/simBBDummy.Rda")
  # 
  # ext<-rstan::extract(mdl.simpdum)
  # y <- fake$bb
  # y.ext <- ext$ypred_new
  # pdf("output/yvsypred.pdf", width =3, height = 3)
  # ppc_dens_overlay(y, y.ext[1:50, ])
  # dev.off()
#   #------------------------------------------------------
#   
#   
#   #Add in modelled Values
#   posteriorMean <- fake #make a new df with the same size as fake
#   alpha.pheno.sp.post <- colMeans(data.frame(ext$a_sp)) # get mean species intercepts 
#   posteriorMean$alpha.pheno.sp <- rep(alpha.pheno.sp.post, each = ntot) # put mean species values in the posterior df 
#   
#   alpha.force.sp.post <- colMeans(data.frame(ext$b_force)) # get mean species intercepts 
#   posteriorMean$alpha.force.sp <- rep(alpha.force.sp.post, each = ntot) # put mean species values in the posterior df 
#   
#   pairs(mdl.simpdum, pars = c("mu_grand", "mu_force", "mu_chill", "mu_photo","sigma_a", "sigma_force", "sigma_chill", "sigma_force","sigma_y", "lp__")) 
#   
# 
# # prior predictive check added by Deirdre
# #Make a vector of site names
# sites <- 1:nsite
# site <- rep(rep(sites, times = nPhenCombo, each = nPhenCombo*nind))
# 
# #Makle a vector of species names 
# species <- 1:nsp 
# sp <- rep(species, each = nsite*nPhenCombo*nind)
# 
# #Make a data frame with ass sites and species laid out 
# #---------------------
# 
# prior <- data.frame(cbind(site, sp))
# head(prior)
# 
# prior$photo <- rep(c(0, 1))
# prior$warm <- rep(c(0,1), each = 2)
# prior$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)
# prior$simRep <- rep(1:nrow(prior),times =1)
# 
# 
# prior <- prior %>%
#   mutate ( d2 = if_else(site == 2, 1, 0),
#            d3 = if_else(site == 3, 1, 0),
#            d4 = if_else(site == 4, 1, 0))
# 
# prior$mu_grand <- rep(rnorm(n= nsp, mean =50, sd =5), each = nind*nPhenCombo*nsite )
# prior$mu_force <- rep(rnorm(nsp, 0, 1), each = nind*nPhenCombo*nsite*nsp )
# prior$mu_photo <- rnorm(nrep,0, 1)
# prior$mu_chill <- rnorm(nrep, 0, 1)
# prior$sigma_force <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
# prior$sigma_chill <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
# prior$sigma_photo <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
# 
# prior$sigma_y <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
# prior$sigma_a <- rtruncnorm(nrep, a=0, mean = 3, sd = 1)
# 
# prior$sigma_site <- rnorm(nrep, 3,1)
# prior$b_site <- 0
# 
# # define priors:
# nrep <- 10
# output <- data.frame(rep(1:nrep,times = nPhenCombo*nind*nsite*nsp), rep(site, times = nrep), rep(species, times = nrep)) # df that will help save data
# names(output) <- c("rep","site","species")
# output$photo <- rep(c(0, 1))
# output$warm <- rep(c(0,1), each = 2)
# output$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)
# 
# output <- output %>%
#   mutate ( d2 = if_else(site == 2, 1, 0),
#            d3 = if_else(site == 3, 1, 0),
#            d4 = if_else(site == 4, 1, 0))
# 
# output$yhat.prior <- NA
# 
# # now get the species specific values:
# 
# output$mu_force <- rep(rnorm(nsp, 0, 1), times = nsite*nPhenCombo*nind)
# output$mu_photo <- rep(rnorm(nsp,0, 1), times = nsite*nPhenCombo*nind)
# output$mu_chill <- rep(rnorm(nsp, 0, 1), times = nsite*nPhenCombo*nind)
# output$sigma_force <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)
# output$sigma_chill <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)
# output$sigma_photo <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)
# 
# output$sigma_y <- rtruncnorm(nrow(output), a=0, mean = 0, sd = 5)
# output$sigma_a <- rep(rtruncnorm(nsp, a=0, mean = 3, sd = 1), times = nsite*nPhenCombo*nind)
# sigma_site <- rnorm(nsite, 3,1)
# output$sigma_site <- rep(sigma_site, each = nPhenCombo*nind)
# 
# mu_force_sim <- rep(rnorm(nsp, 0, 1), times = nsite*nPhenCombo*nind)
# mu_photo_sim <- rep(rnorm(nsp,0, 1), times = nsite*nPhenCombo*nind)
# mu_chill_sim <- rep(rnorm(nsp, 0, 1), times = nsite*nPhenCombo*nind)
# sigma_force_sim <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)
# sigma_chill_sim <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)
# sigma_photo_sim <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)
# 
# sigma_y_sim <- rtruncnorm(nrow(output), a=0, mean = 0, sd = 5)
# sigma_a_sim <- rep(rtruncnorm(nsp, a=0, mean = 3, sd = 1), times = nsite*nPhenCombo*nind)
# sigma_site <- rnorm(nsite, 3,1)
# sigma_site_sim <- rep(sigma_site, each = nPhenCombo*nind)
# 
# d1_sim <- rep(c(0,0,0,0), each = nrep)
# d2_sim <- rep(c(0,1,0,0), each = nrep)
# d3_sim <- rep(c(0,0,1,0), each = nrep)
# d4_sim <- rep(c(0,0,0,1), each = nrep)
# 
# force.i <- rep(1:2, times = nrep)
# chill.i <- rep(rnorm(n = nChill, mean = 0, sd = 10), times = nrep)
# photo.i <- rep(1:2, times = nrep)
# 
# 
# temp <- vector()
# for(r in 1:nrep){
#   
#   alpha.pheno.sp <- rnorm(nrep, 0, sigma_a_sim[r])
#   
#   alpha.force.sp <- rnorm(nrep, mu_force_sim[r], sigma_force_sim[r])
#   
#   alpha.photo.sp <- rnorm(nrep, mu_photo_sim[r], sigma_photo_sim[r])
#   
#   alpha.chill.sp <- rnorm(nrep, mu_chill_sim[r], sigma_chill_sim[r])
#   
#   alpha.site <- rnorm(nrep, 0, sigma_site_sim[r])
#   
#   gen.var <- rnorm(nrep, 0, sigma_y_sim[1])
#   
#   doy.prior <- mu_grand +  alpha.pheno.sp + alpha.pheno.sp + alpha.force.sp * force.i +
#     alpha.chill.sp * chill.i + alpha.photo.sp * photo.i + alpha.site*d1_rep +
#     alpha.site*d2_rep  + alpha.site*d3_rep  + alpha.site*d4_rep +
#     gen.var
#   
#   temp <- rbind(temp, doy.prior)
# }
# 
# hist(temp)

##########################################
# let's try again
# define priors:
nrep <- 10
output <- data.frame(rep(1:nrep,times = nPhenCombo*nind*nsite*nsp), rep(site, times = nrep), rep(species, times = nrep)) # df that will help save data
names(output) <- c("rep","site","species")
output$iter <- c(1:nrow(output))
output$photo <- rep(c(0, 1))
output$warm <- rep(c(0,1), each = 2)
output$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)

output <- output %>%
  mutate ( d1 = if_else(site == 1, 0, 0),
           d2 = if_else(site == 2, 1, 0),
           d3 = if_else(site == 3, 1, 0),
           d4 = if_else(site == 4, 1, 0))

output$yhat.prior <- NA

# now generate the mu and sigmas
for( ir in 1:nrep){
mu_grand <- rnorm(n= nrep, mean =50, sd =5)
output$mu_grand[output$rep == ir] <- mu_grand[ir]

mu_force <- rnorm(n= nrep, mean =0, sd =35)
output$mu_force[output$rep == ir] <- mu_force[ir]

mu_chill <- rnorm(n= nrep, mean =0, sd =35)
output$mu_chill[output$rep == ir] <- mu_chill[ir]

mu_photo <- rnorm(n= nrep, mean =0, sd =35)
output$mu_photo[output$rep == ir] <- mu_photo[ir]

b_site <- rnorm(n= nrep, mean =0, sd =5)
output$b_site[output$rep == ir] <- b_site[ir]

sigma_force <- rtruncnorm(n = nrep, a = 0, mean = 1, sd = 5)
output$sigma_force [output$rep == ir] <- sigma_force[ir]

sigma_chill <-rtruncnorm(n = nrep, a = 0, mean = 1, sd = 5)
output$sigma_chill [output$rep == ir] <- sigma_chill[ir]

sigma_photo <- rtruncnorm(n = nrep, a = 0, mean = 1, sd = 5)
output$sigma_photo [output$rep == ir] <- sigma_photo[ir]

sigma_y <- rtruncnorm(n = nrep, a = 0, mean = 0, sd = 5)
output$sigma_y [output$rep == ir] <- sigma_y[ir]

sigma_a <- rtruncnorm(n = nrep, a = 0, mean = 0, sd = 5)
output$sigma_a [output$rep == ir] <- sigma_a[ir]
}


for(r in 1:nrow(output)){
  alpha.pheno.sp <- rnorm(1, 0, output$sigma_a[r])
  output$alpha.pheno.sp[output$iter == r] <- alpha.pheno.sp
  
  alpha.force.sp <- rnorm(1, output$mu_force[r], output$sigma_force[r])
  output$alpha.force.sp[output$iter == r] <- alpha.force.sp
  
  alpha.photo.sp <- rnorm(1, output$mu_photo[r],output$sigma_photo[r])
  output$alpha.photo.sp[output$iter == r] <- alpha.photo.sp
  
  alpha.chill.sp <- rnorm(1, output$mu_chill[r], output$sigma_chill[r])
  output$alpha.chill.sp[output$iter == r] <- alpha.chill.sp
  
  alpha.site <- rnorm(1, 0, 10)
  output$alpha.site[output$iter == r] <- alpha.site
  
  gen.var <- rnorm(1, 0, output$sigma_y[r])
  output$gen.var[output$iter == r] <- gen.var
  
  
}

output$yhat.prior <- output$mu_grand +  output$alpha.pheno.sp + output$alpha.pheno.sp + output$alpha.force.sp * output$warm +
  output$alpha.chill.sp * output$chill +output$alpha.photo.sp * output$photo + output$alpha.site* output$d1 +
  output$alpha.site* output$d2  + output$alpha.site* output$d3  + output$alpha.site* output$d4 +
  output$gen.var
pdf("output/priorpc_histyhat.pdf")
hist(output$yhat.prior)
dev.off()

write.csv(output, "output/priorpredcheck.csv", row.names = F)