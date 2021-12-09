# Started November 2021 by DL

# working on test data for phenology model 
# testing for differences between dummy variables and indexing code

if(length(grep("deirdreloughnan", getwd()) > 0)) {
  setwd("~/Documents/github/pheno_bc")
} else{
  setwd("/home/deirdre/pheno_bc")
}

library(rstan)
require(shinystan)
require(bayesplot)
require(tidybayes)
require(truncnorm)
library(ggplot2)
library(dplyr)
library(plyr)



rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

nsite = 4 #number of sites with trait data
nsp = 5 #number of species
nind = 3

rep = 1 # making it greater for the test data

#number of treatment levels
nwarm = 2
nphoto = 2
nchill = 2

ntot<-nsite*nwarm*nphoto*nchill*nind #48 
ntotsp<-nsite*nwarm*nphoto*nchill*nsp*nind
# code below to loop through species and indiviudals 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
#building the required dataframe 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

site = gl(nsite, rep, length = ntot)
warm = gl(nwarm, rep*nsite, length = ntot)
photo = gl(nphoto, rep*nsite*nwarm, length = ntot)
chill = gl(nchill, rep*nsite*nwarm*nphoto, length = ntot)

mugrand = 50  
sitediff = 2
sitediff2 = 0.5
sitediff3 = 0.75
sitediff4 = 1
warmdiff = -20 # days earlier from 1 to 2
photodiff = -14
chilldiff = -20 
warmphoto = 3.5
warmchill = -3
chillphoto = -2

mugrand.sd = 5
sitediff.sd = 0.5
warmdiff.sd = 1 
photodiff.sd = 1
chilldiff.sd =1
#sitewarm.sd = 1
# site2photo.sd = 1
# site3photo.sd = 1
# site4photo.sd = 1
warmphoto.sd = 1
warmchill.sd = 1
chillphoto.sd = 1

fake <- vector()

for(i in 1:(nsp)){ 
  fakex <- data.frame(sp = i, site, warm,chill, photo)
  fake <- rbind(fake, fakex)  
}
head(fake)

# replicating the method for calculating the response variable used in traitors (which I understand the best)
alpha.pheno.sp <- rnorm(nsp, mugrand, mugrand.sd) 
fake$alpha.pheno.sp <- rep(alpha.pheno.sp, each = ntot)

alpha.force.sp <- rnorm(nsp, warmdiff, warmdiff.sd)
fake$alpha.force.sp <- rep(alpha.force.sp, each = ntot)

alpha.photo.sp <- rnorm(nsp, photodiff, photodiff.sd)
fake$alpha.photo.sp <- rep(alpha.photo.sp, each = ntot)

alpha.chill.sp <- rnorm(nsp, chilldiff, chilldiff.sd)
fake$alpha.chill.sp <- rep(alpha.chill.sp, each = ntot)

# add dummy/ site level effects:
fake <- fake %>%
  mutate ( d2 = if_else(site == 2, 1, 0),
           d3 = if_else(site == 3, 1, 0),
           d4 = if_else(site == 4, 1, 0)) 

fake$site.effect <- fake$site
fake$site.effect[fake$site.effect==1] <- 0.1
fake$site.effect[fake$site.effect==2] <- 0.5
fake$site.effect[fake$site.effect==3] <- 0.75
fake$site.effect[fake$site.effect==4] <- 1


#general variance
sigma.gen <- 5
gen.var <- rnorm(ntotsp, 0, sigma.gen) 
fake$gen.er <- gen.var

fake$warm <- as.numeric(fake$warm)
fake$warm[fake$warm == 2] <- 0
fake$photo <- as.numeric(fake$photo)
fake$photo[fake$photo == 2] <- 0
fake$chill <- as.numeric(fake$chill)
fake$site <- as.numeric(fake$site)

# now set up the data to be z-scored:
fake$force.z2 <- (fake$warm-mean(fake$warm,na.rm=TRUE))/(sd(fake$warm,na.rm=TRUE)*2)
fake$photo.z2 <- (fake$photo-mean(fake$photo,na.rm=TRUE))/(sd(fake$photo,na.rm=TRUE)*2)
fake$chill.z2 <- (fake$chill-mean(fake$chill,na.rm=TRUE))/(sd(fake$chill,na.rm=TRUE)*2) 


#"run" the full model to simulate data 
# fake$bb <-  fake$alpha.pheno.sp + fake$alpha.force.sp * fake$force.z2 + fake$alpha.chill.sp * fake$chill.z2 + fake$alpha.photo.sp * fake$photo.z2 + fake$gen.er + fake$site.effect2 + fake$site.effect3 + fake$site.effect4 

fake$bb <-  fake$alpha.pheno.sp + fake$alpha.force.sp * fake$warm + fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + fake$gen.er + fake$site.effect

# check if works with lmer or brms
require(lme4)
summary(lmer(bb ~ site + warm + photo + chill + (1|sp), data = fake)) # sanity check 


head(fake)

# Still having issues with the test data: doing a more indepth ppc:

# have site as a single column, with the site effects substituted for the site number
# try running the model in brms/lmer

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
datalist$n_site
datalist$force
head(fake)


mdl.simp2 <- stan("stan/bc.bb.stan",
              data = datalist,
              iter = 4000, chains=1)

mdl.i <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_index.stan",
              data = datalist,
              iter = 4000, chains=4, control = list(adapt_delta = 0.99))
save(mdl.i, file = "output/test_index_01.Rds")

mdl.d <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_dummy.stan",
              data = datalist,
              iter = 4000, chains=4,
              control = list(adapt_delta = 0.99))

save(mdl.simp2, file = "output/test_dummy.Rds")
# Look at model output:
sm.i <- summary(mdl.simp2)$summary
# ext<-rstan::extract(mdl.i)
# 
ssm<- as.shinystan(mdl.simp2)
launch_shinystan(ssm)

pairs(mdl.simp2, pars = c("mu_grand","mu_a", "mu_force", "mu_chill", "mu_photo","mu_site", "lp__")) 


load("output/test_index_01.Rds")
ssm<- as.shinystan(mdl.i)
launch_shinystan(ssm)

sm.i <- summary(mdl.i)$summary

pairs(mdl.i, pars = c("mu_grand","mu_a", "mu_force", "mu_chill", "mu_photo","mu_inter_fp","mu_inter_fc","mu_inter_pc", "lp__")) 

pairs(mdl.i, pars = c("sigma_force", "sigma_chill","sigma_photo","sigma_a","sigma_y", "lp__")) 

ext<-rstan::extract(mdl.i)
# get_variables(mdl.i)
# 
h1 <- hist(rnorm(1000, 50,10))
h2 <- hist(ext$mu_grand)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(0, 400))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,20))
h2 <- hist(ext$mu_force)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50,0))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,35))
h2 <- hist(ext$mu_chill)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-500, 500))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,35))
h2 <- hist(ext$mu_photo)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-500, 500))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,35))
h2 <- hist(ext$mu_inter_fc)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-500, 500))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,35))
h2 <- hist(ext$mu_inter_pc)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-500, 500))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,35))
h2 <- hist(ext$mu_inter_fc)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-500, 500))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,10))
h2 <- hist(ext$sigma_y)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,10))
h2 <- hist(ext$sigma_force)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,30))
h2 <- hist(ext$sigma_chill)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,10))
h2 <- hist(ext$sigma_photo)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,10))
h2 <- hist(ext$sigma_photo)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,10))
h2 <- hist(ext$sigma_b_inter_fp)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,10))
h2 <- hist(ext$sigma_b_inter_pc)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,10))
h2 <- hist(ext$sigma_b_inter_fc)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

y<-pheno.dat$bb

y.ext<-ext$y_hat # I want this to be a matrix, which it is, with one element for each data point in y

ppc_dens_overlay(y, y.ext[1:50, ])


