# Started November 2021 by DL

# working on test data for phenology model 
# testing for differences between dummy variables and indexing code

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc")
} else if(length(grep("Lizzie", getwd())>0)) {
  setwd("~/Documents/git/teaching/stan")
} else{
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

mu_grand = 50
b_site = 10
sitediff2 = 1
sitediff3 = 2
sitediff4 = 3
mu_force = -20 
mu_photo = -14
mu_chill = -20 
# warmphoto = 3.5
# warmchill = -3
# chillphoto = -2

sigma_a = 5
sigma_site = 0.5
sigma_force = 20
sigma_chill = 50
sigma_photo = 20
sigma_y = 5
#sitewarm.sd = 1
# site2photo.sd = 1
# site3photo.sd = 1
# site4photo.sd = 1
# warmphoto.sd = 1
# warmchill.sd = 1
# chillphoto.sd = 1

fake <- vector()

for(i in 1:(nsp)){ 
  fakex <- data.frame(sp = i, site, warm,chill, photo)
  fake <- rbind(fake, fakex)  
}
head(fake)

# replicating the method for calculating the response variable used in traitors (which I understand the best)
alpha.pheno.sp <- rnorm(nsp, 0, sigma_a) 
fake$alpha.pheno.sp <- rep(alpha.pheno.sp, each = ntot)

alpha.force.sp <- rnorm(nsp, mu_force, sigma_force)
fake$alpha.force.sp <- rep(alpha.force.sp, each = ntot)

alpha.photo.sp <- rnorm(nsp, mu_photo, sigma_photo)
fake$alpha.photo.sp <- rep(alpha.photo.sp, each = ntot)

alpha.chill.sp <- rnorm(nsp, mu_chill, sigma_chill)
fake$alpha.chill.sp <- rep(alpha.chill.sp, each = ntot)

fake$site <- as.numeric(fake$site)
alpha.site <- rnorm(nsite, b_site, sigma_site)
fake$alpha.site <- fake$site
fake$alpha.site[fake$alpha.site==1] <- alpha.site[1]
fake$alpha.site[fake$alpha.site==2] <- alpha.site[2]
fake$alpha.site[fake$alpha.site==3] <- alpha.site[3]
fake$alpha.site[fake$alpha.site==4] <- alpha.site[4]

# add dummy/ site level effects:
fake <- fake %>%
  mutate ( d1 = if_else(site == 1, 1, 0),
           d2 = if_else(site == 2, 1, 0),
           d3 = if_else(site == 3, 1, 0),
           d4 = if_else(site == 4, 1, 0))

#general variance

gen.var <- rnorm(ntotsp, 0, sigma_y) 
fake$gen.er <- gen.var

fake$warm <- as.numeric(fake$warm)
fake$warm[fake$warm == 2] <- 0
fake$photo <- as.numeric(fake$photo)
fake$photo[fake$photo == 2] <- 0
fake$chill <- as.numeric(fake$chill)

# now set up the data to be z-scored:
fake$force.z2 <- (fake$warm-mean(fake$warm,na.rm=TRUE))/(sd(fake$warm,na.rm=TRUE)*2)
fake$photo.z2 <- (fake$photo-mean(fake$photo,na.rm=TRUE))/(sd(fake$photo,na.rm=TRUE)*2)
fake$chill.z2 <- (fake$chill-mean(fake$chill,na.rm=TRUE))/(sd(fake$chill,na.rm=TRUE)*2) 


#"run" the full model to simulate data
# calcualte test data for two sites:
# fake$site[fake$site == 2] <- 0
#fake$site <- as.numeric(fake$site)
# fake$bb <-  mu_grand + fake$alpha.pheno.sp + fake$alpha.force.sp * fake$warm + fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + fake$gen.er + fake$alpha.site*fake$site

# for dummy variable test data
# Faith suggested the issue might be that the effect of site 1 is not being included; bb without site 1 bb2 is with site 1 effects
fake$bb <-  mu_grand + fake$alpha.pheno.sp + fake$alpha.force.sp * fake$warm + fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + fake$gen.er  + fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4

fake$bb2 <-  mu_grand + fake$alpha.pheno.sp + fake$alpha.force.sp * fake$warm + fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + fake$gen.er  + alpha.site*fake$d1 + alpha.site*fake$d2  + alpha.site*fake$d3  + alpha.site*fake$d4
 
fake$site2 <- fake$alpha.site*fake$d2
fake$site3 <- fake$alpha.site*fake$d3
fake$site4 <- fake$alpha.site*fake$d4
# check if works with lmer or brms

summary(lmer(bb ~  warm + photo + chill + alpha.site + (1|sp), data = fake)) 
summary(lmer(bb2 ~  warm + photo + chill + alpha.site + (1|sp), data = fake)) 

summary(lm(bb ~  warm + photo + chill + d2 + d3 + d4 , data = fake)) # 
summary(lmer(bb ~  warm + photo + chill + d2 + d3 + d4 + (1|sp), data = fake)) # 

# sanity check that it is not just that all the values are the same:
fake$bb3 <-  mu_grand + fake$alpha.pheno.sp + fake$alpha.force.sp * fake$warm + fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + fake$gen.er  + sitediff2*fake$d2  + sitediff3*fake$d3  + sitediff4*fake$d4

summary(lmer(bb3 ~  warm + photo + chill + d2 + d3 + d4 + (1|sp), data = fake)) # 
# The values are ok, but not stellar

# tbb, f/c/p.n, site.n, species, f/c/p.i, species.fact, d2, d3, d4
datalist <- list( N=nrow(fake),
                    n_sp = length(unique(fake$sp)),
                    n_site = length(unique(fake$site)),
                    bb = fake$bb,
                    sp = fake$sp,
                    chill = fake$chill,
                    photo = fake$photo,
                    force = fake$warm,
                    site = fake$site
                    , site2 = fake$d2,
                    site3 = fake$d3,
                    site4 = fake$d4)

# mdl.simp <- stan("stan/bc.bb.stan",
#                     data = datalist,
#                     include = FALSE, pars = c("ypred_new","y_hat"),
#                     iter = 4000, chains= 4)

mdl.simpdum <- stan("stan/test_model.stan",
                  data = datalist,
                  include = FALSE, pars = c("ypred_new","y_hat"),
                  iter = 4000, chains= 4)

sm <- summary(mdl.simpdum)$summary

param <- list(mu_grand = 50, mu_force = -20,
              mu_chill = -20,  mu_photo = -14, b_site = 10, sigma_a = 5,
              sigma_force = 20,sigma_chill = 50, sigma_photo =20,  sigma_y = 5, site2 = alpha.site[2], site3 = alpha.site[3], site4 = alpha.site[4])

# summary(mdl.simp)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","mu_site","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_site","sigma_y"),"mean"]; t(param)

summary(mdl.simpdum)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","b_site2","b_site3","b_site4","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y"),"mean"]; t(param)

ssm<- as.shinystan(mdl.simpdum)
launch_shinystan(ssm)

# mdl.i <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_index.stan",
#               data = datalist,
#               iter = 4000, chains=4, control = list(adapt_delta = 0.99))
# save(mdl.simp, file = "output/test_simple2sites.Rds")
# 
# mdl.d <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_dummy.stan",
#               data = datalist,
#               iter = 4000, chains=4,
#               control = list(adapt_delta = 0.99))
# 
# save(mdl.simp, file = "output/test_2sites.Rds")
# Look at model output:
# sm.i <- summary(mdl.simp2)$summary
# ext<-rstan::extract(mdl.i)
load("output/test_simple2sites.Rds")
ssm<- as.shinystan(mdl.simpdum)
launch_shinystan(ssm)

pairs(mdl.simp, pars = c("mu_grand", "mu_force", "mu_chill", "mu_photo","mu_site","sigma_a", "sigma_force", "sigma_chill", "sigma_force","sigma_y", "lp__")) 
#pdf("simp_dum_pairs.pdf", height = 5, width = 5)

pairs(mdl.simpdum, pars = c("mu_grand", "mu_force", "mu_chill", "mu_photo","b_site2","b_site3","b_site4","sigma_a", "sigma_force", "sigma_chill", "sigma_force","sigma_y", "lp__")) 
#dev.off()

ext<-rstan::extract(mdl.simpdum)
# get_variables(mdl.i)
# 
h1 <- hist(rnorm(1000, 50,10))
h2 <- hist(ext$mu_grand)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(0, 400))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,20))
h2 <- hist(ext$mu_force)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-500,0))
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

h1 <- hist(rnorm(1000, 1,5))
h2 <- hist(ext$b_site2)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

y<-pheno.dat$bb

y.ext<-ext$y_hat # I want this to be a matrix, which it is, with one element for each data point in y

ppc_dens_overlay(y, y.ext[1:50, ])

### # Still having issues with the test data: starting to do a more indepth ppc:


ntot_prior<-nsite*nwarm*nphoto*nchill*nind #48 
ntotsp<-nsite*nwarm*nphoto*nchill*nsp*nind #480
# code below to loop through species and indiviudals 

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
#building the required dataframe 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

site = gl(nsite, rep, length = ntot)
warm = gl(nwarm, rep*nsite, length = ntot)
photo = gl(nphoto, rep*nsite*nwarm, length = ntot)
chill = gl(nchill, rep*nsite*nwarm*nphoto, length = ntot)

#priors:
mu_grand <- rnorm(1000, 50, 5)
mu_force <- rnorm(1000, 0, 50)
mu_photo <- rnorm(1000,0, 35)
mu_chill <- rnorm(1000,0, 35)
b_site <- rnorm(1000,0,0.5)
sigma_force <- rnorm(1000,20, 1)
sigma_photo <- rnorm(1000,0, 10)
sigma_chill <- rnorm(1000,50, 30)
sigma_y <- rnorm(1000,0,10)
sigma_a <- rnorm(1000, 30, 20)

b_force <- rnorm(1000, mu_force[1], sigma_force[1]);
b_chill <- rnorm(1000, mu_chill[1], sigma_chill[1]);
b_photo <- rnorm(1000, mu_photo[1], sigma_photo[1]);

a_sp <- rnorm(1000, 0,0.1);

mu_bb <- mu_grand[1] + a_sp[1] + b_force[1] * force + b_chill[1] * chill + b_photo * photo + b_site * site2 + b_site * site3 + b_site * site4
yhat <- rnorm(mu_bb, sigma_y[1])

post <- rstan::extract(mdl.simpdum)
plot(density(sigma_a)); lines(density(post$sigma_a), col = "red", lwd = 3)
plot(density(sigma_y)); lines(density(post$sigma_y), col = "red", lwd = 3)
plot(density(sigma_force)); lines(density(post$sigma_force), col = "red", lwd = 3)
plot(density(sigma_chill)); lines(density(post$sigma_chill), col = "red", lwd = 3)
plot(density(sigma_photo)); lines(density(post$sigma_photo), col = "red", lwd = 3)


prior.
