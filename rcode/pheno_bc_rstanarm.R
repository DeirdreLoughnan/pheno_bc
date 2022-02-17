#started Feb 1, 2022 by deirdre

# trying to get the full dummy model with interactions to work in stan lmer:

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
require(rstanarm)
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
nind = 10 
nChill <- 5 #This many chilling treatment levels 
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
site <- rep(rep(sites, times = nsite, each = nPhenCombo*nind))

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
fake$alpha.site <- rep(rep(alpha.site, times = nsite, each = nPhenCombo*nind))

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

fake$bb_int_nochill <-  mu_grand + alpha.site[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 + fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) +
  fake$gen.var


# check if works with lmer or brms

 summary(lm(bb ~  warm + photo + chill + d2 + d3 + d4, data = fake)) # 
# 
 summary(lm(bb_int ~  warm + photo + chill + d2 + d3 + d4 + warm * chill + chill * photo + photo * warm, data = fake)) # 

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

# stan lmer

# mdl 1: site not included, just species as a random intercepts and slopes and interactions
# mdl1 <-stan_lmer(bb ~ warm + photo + chill +  warm * photo + chill * photo + warm * chill +#main effects
#                               (1| sp) + # random intercepts
#                    (sp|warm) + (sp|photo) + (sp|chill) , #random slopes
#                             data = fake,
#                             chains = 2, cores = 2,iter = 2500, warmup=1500)
# save(mdl1, file="rstanarm_nosite.Rda")

# mdl 2: site included as random slope, species as random slope and intercept, interactions for cues
# error message that can't use site as grouping factor more than once, suggests using || to fix this issue
# mdl2 <-stan_lmer(bb ~ warm + photo + chill +  warm * photo + chill * photo + warm * chill +#main effects
#                    (1| sp) + # random intercepts
#                    (sp||warm) + (sp||photo) + (sp||chill) + (site||warm) + (site||photo) + (site||chill), #random slopes
#                  data = fake,
#                  chains = 2, cores = 2,iter = 2500, warmup=1500)
# save(mdl2, file="rstanarm_siteRandInt.Rda")

# 59 divergent transitions after warmup
# mdl 3: site included as a dummy variable,  species as a random intercepts and slopes and interactions for cues
#In this case sp is used as a grouping factor multiple times and (Intercept) is included multiple times.
#Consider using || or -1 in your formulas to prevent this from happening.
# mdl3 <- stan_lmer(bb ~ warm + photo + chill + d2 + d3 + d4 +  
#                          # warm*photo + warm*chill +
#                          # photo*chill +
#                         # + (1|sp) + # random intercepts
#                           #(warm|sp) + (photo|sp) + (chill|sp),   # random slopes
#                     ( warm+ photo + chill|sp),
#                         data=fake,
#                 chains = 2, cores = 2,iter = 2500, warmup=1500)
# 
# mdl3 <- stan_lmer(bb ~ ((warm + photo + chill)|sp) + d2 + d3 + d4+
#                   warm*photo + warm*chill + photo*chill, # + as.factor(site)
#                   data=fake,
#                   chains = 2, cores = 2,iter = 2500, warmup=1500, adapt_delta = 0.99)
# save(mdl3, file="rstanarm_siteDummInt.Rda")


# mdl3b <- stan_lmer(bb ~ ((warm* photo * chill)^2|sp) + d2 + d3 + d4,
#                   data=fake,
#                   chains = 2, cores = 2,iter = 2500, warmup=1500)

mdl6 <- stan_lmer(bb ~ ((warm + photo + chill)|sp) + as.factor(site) 
                    # + warm*photo + warm*chill + photo*chill, 
                  ,data=fake,
                  chains = 4, cores = 4,iter = 4000, warmup=3000)
save(mdl6, file="rstanarm_siteFactNoInt.Rda")



# mdl 4: site included as a dummy variable,  species as a random intercepts and slopes, and interactions for cues and sites
# mdl4 <-stan_lmer(bb ~ warm + photo + chill + d2 + d3 + d4 +
#                  warm*photo + warm*chill + photo*chill +  
#                    d2* warm + d3* warm + d4*warm +
#                    d2* chill + d3* chill + d4*chill +
#                    d2* photo + d3* photo + d4*photo
#                    + (1|sp) + # random intercepts
#                    (warm|sp) + (photo|sp) + (chill|sp),   # random slopes
#                  data=fake,
#                  chains = 2, cores = 2,iter = 2500, warmup=1500)
# 
# save(mdl4, file="rstanarm_siteDummIntCueSite.Rda")

#(photo*chill||genus)
#load("output/rstanarm_siteDummInt.Rda")
# load("output/rstanarm_siteFact.Rda")
load("output/rstanarm_siteFactIterAdapt.Rda")
ssm <-  as.shinystan(mdl4)
launch_shinystan(ssm)

sumt <- summary(mdl4)
#
range(sumt[, "n_eff"]) #
range(sumt[, "Rhat"])
#
# pdf("pairs_dummAll.pdf")
pairs(mdl5, pars = c("(Intercept)","warm","chill","photo", "as.factor(site)2","as.factor(site)3", "as.factor(site)4", "sigma","log-posterior"))
# dev.off()

# get the values to plot against the fake data:
temp <- as_tibble(sumt, row.names(sumt))


mdl.out <- mdl4 %>%
  spread_draws(b[(Intercept),sp], sigma) %>%
  median_qi()
mdl.out  <- as.data.frame(mdl.out )
colnames(mdl.out)[colnames(mdl.out) == "(Intercept)"] <- "param"

intSp <- subset(mdl.out, param == "(Intercept)")
intF <- subset(mdl.out, param == "warm")
intC <- subset(mdl.out, param == "chill")
intP <- subset(mdl.out, param == "photo")

par(mfrow = c(2,2))
plot(alpha.pheno.sp ~ intSp$b, pch =19)
plot(alpha.force.sp ~ intF$b, pch =19)
plot(alpha.chill.sp ~ intC$b, pch =19)
plot(alpha.photo.sp ~ intP$b, pch =19)


# load("output/rstanarm_siteFactInt_B.Rda")
# ssm <-  as.shinystan(mdl3b)
# launch_shinystan(ssm)
# 
# sumb <- summary(mdl3b)
# #
# range(sumb[, "n_eff"])
# range(sumb[, "Rhat"])
# # #
# # pdf("pairs_dummAll.pdf")
# pairs(mdl3, pars = c("(Intercept)","warm","chill","photo", "d2", "d3","d4", "sigma","log-posterior"))


