# Started November 2021 by DL

# working on test data for phenology model 
# testing for differences between dummy variables and indexing code

# The test data is still having issues, I am going to try simplying it to two sites, two chilling levels
# sigmas should probably be bigger than I had
#The way I was adding site was also redundant - should fix this too

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

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
#building the required dataframe 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

#Set number of data groups
nsite = 2 #number of sites with trait data
nsp = 8 #number of species - 
nind = 8
nChill <- 2 #This many chilling treatment levels 
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

sigma_a = 10
sigma_site = 5
sigma_force = 10
sigma_chill = 10
sigma_photo = 10
sigma_y = 20

#interactions:
mu_fp <- 3
mu_cp <- 2
mu_fc <- 1

sigma_fp <- 5
sigma_cp <- 5
sigma_fc <- 5

# mu_fsite <- 1
# mu_csite <- 5
# mu_psite  <- 3
# 
# sigma_fsite <- 0.5
# sigma_csite <- 0.5
# sigma_psite <- 0.5

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

#alpha.site <- rnorm(nsite, b_site, sigma_site)
alpha.site <- c(0,1)
fake$alpha.site <- rep(rep(alpha.site, times = nPhenCombo, each = nPhenCombo*nind))

# adding the interaction effects:
alpha.fc <- rnorm(nsp, mu_fc, sigma_fc)
fake$alpha.fc <- rep(alpha.fc, each =  nsite*nPhenCombo*nind)

alpha.cp <- rnorm(nsp, mu_cp, sigma_cp)
fake$alpha.cp <- rep(alpha.cp, each =  nsite*nPhenCombo*nind)

alpha.fp <- rnorm(nsp, mu_fp, sigma_fp)
fake$alpha.fp <- rep(alpha.fp, each =  nsite*nPhenCombo*nind)

#adding interactions between the cues and the 
# alpha.fsite <- rnorm(nsite, mu_fsite, sigma_fsite)
# fake$alpha.fsite <- rep(alpha.fsite, each =  nPhenCombo*nind)
# 
# alpha.csite <- rnorm(nsite, mu_csite, sigma_csite)
# fake$alpha.csite <- rep(alpha.csite, each =  nPhenCombo*nind)
# 
# alpha.psite <- rnorm(nsite, mu_psite, sigma_psite)
# fake$alpha.psite <- rep(alpha.psite, each =  nPhenCombo*nind)

# add dummy/ site level effects:
# fake <- fake %>%
#   mutate ( d2 = if_else(site == 2, 1, 0))
           #, d3 = if_else(site == 3, 1, 0),
           # d4 = if_else(site == 4, 1, 0))



#general variance

fake$gen.var <- rnorm(nrow(fake), 0, sigma_y) 
#hist(rnorm(nrow(fake), 0, sigma_y))
#Simulate continuous cue values

fake$photo <- rep(c(0, 1))
fake$warm <- rep(c(0,1), each = 2)
#chill.port <- c(20, 40, 60, 80, 100, 120)
fake$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)
#fake$chill<- rep(rep(chill.port, each = nForce*nPhoto), times = nsp*nsite*nind)

#check to make things are being replicated the same number of times:
fake$count <- 1
faket <- aggregate(fake ["count"],
                           fake[c("sp","site","warm", "photo","alpha.photo.sp")],
                           FUN = sum)
#View(faket)

# for dummy variable test data
fake$bb <-  mu_grand  + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$site  +
  #fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
  fake$gen.var

fake$bb_int <-  mu_grand + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$site  + 
  #fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 + 
  fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) +
  fake$gen.var

# fake$bb_int_site <-  mu_grand + alpha.site[1] + alpha.psite[1] + alpha.csite[1] + alpha.fsite[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
#   fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + fake$alpha.site*fake$site  + 
#   #fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 + 
#   fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) + fake$alpha.fsite * (fake$site*fake$warm) + 
#   #fake$alpha.fsite * (fake$d3*fake$warm) + fake$alpha.fsite * (fake$d4*fake$warm) + 
#   fake$alpha.csite * (fake$site*fake$chill) + 
#   #fake$alpha.csite * (fake$d3*fake$chill) + fake$alpha.csite * (fake$d4*fake$chill) + 
#   fake$alpha.psite * (fake$site*fake$photo) + 
#   #fake$alpha.psite * (fake$d3*fake$photo) + fake$alpha.psite * (fake$d4*fake$photo) + 
#   fake$gen.var
# estimates don't look too bad, worse for photo but might need to be ncp?

# check if works with lmer or brms

summary(lm(bb ~  warm + photo + chill + site, data = fake)) # 
# 
# summary(lm(bb_int ~  warm + photo + chill + site + warm * chill + chill * photo + photo * warm, data = fake)) # 

mu_grand 
alpha.site
alpha.fc
alpha.cp
alpha.fp
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
                    site = fake$site)
  
mdl.full <- stan("stan/test_model_interactions_truepriors.stan",
                   data = datalist,
                   include = FALSE, pars = c("ypred_new","y_hat"),
                   iter = 4000, chains= 4, warmup = 2000)
save(mdl.full, file = "output/BBDummy_2chill2site_trueprior.Rda")


# fix: issues with sigma_y, sigma_fc, sigma_cp

# mdl.dumint <- stan("stan/test_model_interaction_ncp.stan",
#                    data = datalist,
#                    include = FALSE, pars = c("ypred_new","y_hat"),
#                    iter = 4000, chains= 4, warmup = 2000)
# save(mdl.dumint, file = "output/BBDummy_ncpint_test.Rda")
###################################################################
  load("output/BBDummy_2chill2site_trueprior.Rda")
# # # # #
 ssm <-  as.shinystan(mdl.df)
 launch_shinystan(ssm)
# # # sigma_fp looks really bad, of all the sigma it is the one that is most squished up to one side
#
 sum <- summary(mdl.full)$summary
# get_variables(mdl.full)
 summary(mdl.full)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","b_site","mu_fc" ,"mu_fp" ,"mu_cp","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y", "sigma_fp" , "sigma_cp" , "sigma_fc"),"mean"]

# summary(mdl.df)$summary[c("mu_a","mu_b_warm", "mu_b_chill","mu_b_photo","mu_b_site","mu_b_inter_wc" ,"mu_b_inter_wp" ,"mu_b_inter_pc","mu_b_inter_wc" ,"mu_b_inter_wp" ,"mu_b_inter_ws","sigma_b_warm","sigma_b_chill","sigma_b_photo","sigma_b_site","sigma_y"),"mean"]
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
 ext<-rstan::extract(mdl.full)
par(mfrow=c(1,1))

plot(hist(ext$mu_grand))
abline(v=mu_grand + alpha.site[1],col="red")

plot(hist(ext$mu_force))
abline(v=mu_force,col="red")

plot(hist(ext$mu_chill))
abline(v=mu_chill,col="red")

plot(hist(ext$mu_photo))
abline(v=mu_photo,col="red")

plot(hist(ext$mu_fc))
abline(v = mu_fc, col = "red")

plot(hist(ext$mu_cp))
abline(v = mu_cp, col = "red")

plot(hist(ext$mu_fp))
abline(v = mu_fp, col = "red")

plot(hist(ext$sigma_y))
abline(v = sigma_y,col="red")

plot(hist(ext$sigma_force))
abline(v = sigma_force,col="red")

plot(hist(ext$sigma_chill))
abline(v = sigma_chill,col="red")

plot(hist(ext$sigma_photo))
abline(v = sigma_photo,col="red")

plot(hist(ext$b_site))
abline(v = alpha.site[2],col="red")

plot(hist(ext$sigma_fp))
abline(v = sigma_fp,col="red")

plot(hist(ext$sigma_fc))
abline(v = sigma_fc,col="red")

plot(hist(ext$sigma_cp))
abline(v = sigma_cp,col="red")

#plot posterior expected values againt real values

# ext<-rstan::extract(mdl.full)
# y <- fake$bb
# y.ext <- ext$ypred_new
# pdf("output/yvsypred.pdf", width =3, height = 3)
# ppc_dens_overlay(y, y.ext[1:50, ])
# dev.off()
# #------------------------------------------------------

# compare to fake data:
intSp <- as.data.frame(sum[grep("a_sp", rownames(sum)), ]) 
slopeF <- as.data.frame(sum[grep("b_force", rownames(sum)), ]) 
slopeC <- as.data.frame(sum[grep("b_chill", rownames(sum)), ]) 
slopeP <- as.data.frame(sum[grep("b_photo", rownames(sum)), ]) 
slopeFC <- as.data.frame(sum[grep("b_fc", rownames(sum)), ]) 
slopeCP <- as.data.frame(sum[grep("b_cp", rownames(sum)), ]) 
slopeFP <- as.data.frame(sum[grep("b_fp", rownames(sum)), ]) ; slopeFP <- slopeFP[c(9:16),]


par(mfrow = c(4,2))
plot(alpha.pheno.sp ~ intSp$mean, pch = 19)
plot(alpha.force.sp ~ slopeF$mean, pch = 19)
plot(alpha.chill.sp ~ slopeC$mean, pch = 19)
plot(alpha.photo.sp ~ slopeP$mean, pch = 19)
plot(alpha.fc ~ slopeFC$mean, pch = 19)
plot(alpha.cp ~ slopeCP$mean, pch = 19)
plot(alpha.fp ~ slopeFP$mean, pch = 19)
