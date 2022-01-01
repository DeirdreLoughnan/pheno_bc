 # Started November 2021 by DL

# working on test data for phenology model 
# testing for differences between dummy variables and indexing code


    #Notes for Deirdre

    # You need rm(list=ls()) to be the first line of code you run, else you just wipe away all your libraries and wd

    # Your dummy variables in STAN should have upper and lower bounds set at 0 and 1 

    #WHy are the sigma_cue slopes' priors so far from sinulated values? They are centred around 20, with a sd of only 1. 
    #SO the model is forced to chose close to 20 as the slope values. But your simulated sigma_force is set at 1. 
     # sigma_force ~ normal(20, 1); 
     # sigma_photo ~ normal(20, 1); 
     # sigma_chill ~ normal(50, 1);

    #What is "rep", and how is it different from "ind" when simulating data?

    #I'm not convinced teh way you were simulating cue values was right. 

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

set.seed(194842)

#Flags

Midge <-  TRUE

if(Midge == TRUE){
  if(length(grep("faith", getwd())>0)) {
    setwd("/home/faith/CatPhen")
  } else {setwd("/home/deirdre/pheno_bc")}


} else if (Midge == FALSE){

    if(length(grep("deirdreloughnan", getwd())>0)) {
    setwd("~/Documents/github/pheno_bc")
  } else if(length(grep("Lizzie", getwd())>0)) {
    setwd("~/Documents/git/teaching/stan")
  } else if(length(grep("faith", getwd()) > 0)){
    setwd("/home/faith/Documents/mnt/UBC/otherPeople/Deirdre")
  }  else {
    setwd("/home/deirdre/pheno_bc") # for midge
  }
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
nsp = 8 #number of species - Its really hard to estimate a random effect with only 5 groups 
nind = 4 
nChill <- 6 #This many chilling treatment levels 
nPhoto <- 2
nForce <- 2
nPhenCombo <- nChill*nPhoto*nForce

#rep = 5 # making it greater for the test data. WHat's the difference between ind and rep? You don't ahve a hierarchical effect of individual in the model 
#ntot<-nsite*nwarm*nphoto*nchill*nind #48 #Why do you use thsi value? I don't understand why you multiply by nwarm*nphoto*nchill
#ntot <- nsite*nind*rep
#ntotsp<-ntot*nsp

#site = gl(nsite, rep, length = nind*rep*nsp)
#warm = gl(nwarm, rep*nsite, length = ntot)
#photo = gl(nphoto, rep*nsite*nwarm, length = ntot)
#chill = gl(nchill, rep*nsite*nwarm*nphoto, length = ntot)

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
fake$chill<- rep(rep(rnorm(n = nChill, mean = 50, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)


# for dummy variable test data

fake$bb <-  mu_grand + alpha.site[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
    fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
    fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
    fake$gen.var

hist(fake$bb )

# check if works with lmer or brms

summary(lmer(bb ~  warm + photo + chill + alpha.site + (1|sp), data = fake)) 
summary(lm(bb ~  warm + photo + chill + d2 + d3 + d4 , data = fake)) # 



# The values are very close between trail data and lmer output

if(Midge == TRUE){ #Run stan if I'm on Midge 
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


  mdl.simpdum <- stan("stan/test_model.stan",
                  data = datalist,
                 # include = FALSE, pars = c("ypred_new","y_hat"),
                  iter = 4000, chains= 4, warmup = 2000)

  #save(mdl.simpdum, file = "output/simBBPosterior.Ddata")
  save(mdl.simpdum, file = "output/simBBDummy.Rda")


}


#You should make a list using your original values in lines 35-83 rather than redefine values here. It saves confusion. 
#param <- list(mu_grand = 50, mu_force = -20,
#              mu_chill = -20,  mu_photo = -14, b_site = 10, sigma_a = 5,
#              sigma_force = 1,sigma_chill = 1, sigma_photo =1,  sigma_y = 5, site2 = alpha.site[2], site3 = alpha.site[3], site4 = alpha.site[4])


if(Midge == FALSE){

    load("output/simBBDummy.Rda")



    summary(mdl.simpdum)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","b_site2","b_site3","b_site4","sigma_a","sigma_force","sigma_chill","sigma_photo","sigma_y"),"mean"]


    pairs(mdl.simpdum, pars = c("mu_grand", "mu_force", "mu_chill", "mu_photo","b_site2","b_site3","b_site4","sigma_a", "sigma_force", "sigma_chill", "sigma_force","sigma_y", "lp__")) 
    #dev.off()

    ext<-rstan::extract(mdl.simpdum)
  


    #Plot model posteriors and simulated values
    #-------------------------------------------------

    par(mfrow=c(1,1))

    plot(hist(ext$mu_grand))
    abline(v=mu_grand + alpha.site[1],col="red")

    plot(hist(ext$mu_force))
    abline(v=mu_force,col="red")

    plot(hist(ext$mu_chill))
    abline(v=mu_chill,col="red")

    plot(hist(ext$mu_photo))
    abline(v=mu_photo,col="red")

    plot(hist(ext$sigma_y))
    abline(v=sigma_y,col="red")

    plot(hist(ext$sigma_force))
    abline(v=sigma_force,col="red")

    plot(hist(ext$sigma_chill))
    abline(v=sigma_chill,col="red")

    plot(hist(ext$sigma_photo))
    abline(v=sigma_photo,col="red")

    plot(hist(ext$b_site2))
    abline(v=alpha.site[2],col="red")

    plot(hist(ext$b_site3))
    abline(v=alpha.site[3],col="red")

    plot(hist(ext$b_site4))
    abline(v=alpha.site[4],col="red")

    #plot posterior expected values againt real values
    load("output/simBBDummy.Rda")
    
    ext<-rstan::extract(mdl.simpdum)
    y <- fake$bb
    y.ext <- ext$ypred_new
    pdf("output/yvsypred.pdf", width =3, height = 3)
    ppc_dens_overlay(y, y.ext[1:50, ])
    dev.off()
    #------------------------------------------------------


    #Add in modelled Values
    posteriorMean <- fake #make a new df with the same size as fake
    alpha.pheno.sp.post <- colMeans(data.frame(ext$a_sp)) # get mean species intercepts 
    posteriorMean$alpha.pheno.sp <- rep(alpha.pheno.sp.post, each = ntot) # put mean species values in the posterior df 

    alpha.force.sp.post <- colMeans(data.frame(ext$b_force)) # get mean species intercepts 
    posteriorMean$alpha.force.sp <- rep(alpha.force.sp.post, each = ntot) # put mean species values in the posterior df 
    
    pairs(mdl.simpdum, pars = c("mu_grand", "mu_force", "mu_chill", "mu_photo","sigma_a", "sigma_force", "sigma_chill", "sigma_force","sigma_y", "lp__")) 


    


}

# prior predictive check added by Deirdre
#Make a vector of site names
sites <- 1:nsite
site <- rep(rep(sites, times = nPhenCombo, each = nPhenCombo*nind))

#Makle a vector of species names 
species <- 1:nsp 
sp <- rep(species, each = nsite*nPhenCombo*nind)

#Make a data frame with ass sites and species laid out 
#---------------------

prior <- data.frame(cbind(site, sp))
head(prior)

prior$photo <- rep(c(0, 1))
prior$warm <- rep(c(0,1), each = 2)
prior$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)
prior$simRep <- rep(1:nrow(prior),times =1)


prior <- prior %>%
  mutate ( d2 = if_else(site == 2, 1, 0),
           d3 = if_else(site == 3, 1, 0),
           d4 = if_else(site == 4, 1, 0))

prior$mu_grand <- rep(rnorm(n= nsp, mean =50, sd =5), each = nind*nPhenCombo*nsite )
prior$mu_force <- rep(rnorm(nrep, 0, 1), each = nind*nPhenCombo*nsite*nsp )
prior$mu_photo <- rnorm(nrep,0, 1)
prior$mu_chill <- rnorm(nrep, 0, 1)
prior$sigma_force <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
prior$sigma_chill <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
prior$sigma_photo <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)

prior$sigma_y <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
prior$sigma_a <- rtruncnorm(nrep, a=0, mean = 3, sd = 1)

prior$sigma_site <- rnorm(nrep, 3,1)
prior$b_site <- 0

# define priors:
#priors:
# each time goes over the loop overwrite, repeat the alpha sp values, alpah force
pheno.ppc <- data.frame(cbind(rep(1:1000, times = length(nForce)), rep(nForce, times = 1000), rep(nChill, times = 1000), rep(nChill, times = 1000),rep(nPhoto, times = 1000), rep(nPhoto, times = 1000)), rep(nsite, times = 1000), rep(nsite, times = 1000))

nrep <- 10
output <- data.frame(rep(1:nrep,times = nPhenCombo*nind*nsite*nsp), rep(site, times = nrep), rep(species, times = nrep)) # df that will help save data
names(output) <- c("rep","site","species")
output$photo <- rep(c(0, 1))
output$warm <- rep(c(0,1), each = 2)
output$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)

output <- output %>%
  mutate ( d2 = if_else(site == 2, 1, 0),
           d3 = if_else(site == 3, 1, 0),
           d4 = if_else(site == 4, 1, 0))

output$yhat.prior <- NA

for (ir in 1:nrep){
  mu_grand <- rnorm(nsp, 70, 5)
  output$mu_grand[output$ir == ir] <- (rep(mu_grand, each = nPhenCombo*nind*nsite ))
  
  mu_force <- rnorm(nrep, 0, 1)
  mu_photo <- rnorm(nrep,0, 1)
  mu_chill <- rnorm(nrep, 0, 1)
  sigma_force <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
  sigma_chill <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
  sigma_photo <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)

  sigma_y <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
  sigma_a <- rtruncnorm(nrep, a=0, mean = 3, sd = 1)
  sigma_site <- rnorm(nrep, 3,1)
  b_site <- 0
  
  b_force <- rnorm(nsp, mu_force[ir], sigma_force[ir])
  output$b_force[output$rep == ir] <- rep(b_force, each = nsp)
  
  b_chill <- rnorm(nsp, mu_chill[ir], sigma_chill[ir])
  output$b_chill[output$rep == ir] <- rep(b_chill, each = nsp)
  
   b_photo <- rnorm(nsp, mu_photo[ir], sigma_photo[ir])
   output$b_photo[output$rep == ir] <- rep(b_photo, each = nsp)
   

    b_site <- rnorm(nsite, b_site, sigma_site[ir])
    output$b_site[output$rep == ir] <- rep(b_site, each = nsite)
    

    alpha.pheno.sp <- rnorm(nsp, 0, sigma_a[ir])
    output$alpha.pheno.sp[output$rep == ir] <- rep(alpha.pheno.sp, each = nsp)
    

    alpha.force.sp <- rnorm(nsp, mu_force[ir], sigma_force[ir])
    output$alpha.force.sp[output$rep == ir] <- rep(alpha.force.sp, each = nsp)
    

    alpha.photo.sp <- rnorm(nsp, mu_photo[ir], sigma_photo[ir])
    output$alpha.photo.sp[output$rep == ir] <- rep(alpha.photo.sp, each = nsp)
    

    alpha.chill.sp <- rnorm(nsp, mu_chill[ir], sigma_chill[ir])
    output$alpha.chill.sp[output$rep == ir] <- rep(alpha.chill.sp, each = nind)
    

    alpha.site <- rnorm(nsite, b_site, sigma_site[ir])
    output$alpha.site[output$rep == ir] <- rep(alpha.site, each = nsite)
    

    output$gen.var[output$rep == ir] <- rnorm(nrow(prior), 0, sigma_y)
}


    output$mu_bb <- output$mu_grand[1] +  output$alpha.site[1] + output$alpha.pheno.sp[1] + output$alpha.force.sp[1] * output$warm +
      output$alpha.chill.sp[1] * output$chill + output$alpha.photo.sp[1] * output$photo +
      output$alpha.site[1]*output$d2  + output$alpha.site[1]*output$d3  + output$alpha.site[1]*output$d4 +
      output$gen.var[1]
    output$yhat <- outout$mu_bb + output$sigma_y


nrep <- 10

# mu_grand <- rnorm(nrep, 70, 5)
# mu_force <- rnorm(nrep, 0, 1)
# mu_photo <- rnorm(nrep,0, 1)
# mu_chill <- rnorm(nrep, 0, 1)
# sigma_force <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
# sigma_chill <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
# sigma_photo <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
# 
# sigma_y <- rtruncnorm(nrep, a=0, mean = 0, sd = 5)
# sigma_a <- rtruncnorm(nrep, a=0, mean = 3, sd = 1)

# sigma_force <- rnorm(nrep,5,1)
# sigma_photo <- rnorm(nrep,5,1)
# sigma_chill <- rnorm(nrep,5,1)
# sigma_y <- rnorm(nrep,5,1)
# sigma_a <- rnorm(nrep, 3, 1)
# sigma_site <- rnorm(nrep, 3,1)
# b_site <- 0
# names(pheno.ppc) <- c("iteration","forcing", "chilling","photo")
# pheno.ppc$doy.prior <- NA

nrep <- 10
output <- data.frame(rep(1:nrep,times = nPhenCombo*nind*nsite*nsp), rep(site, times = nrep), rep(species, times = nrep)) # df that will help save data
names(output) <- c("rep","site","species")
output$photo <- rep(c(0, 1))
output$warm <- rep(c(0,1), each = 2)
output$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)

output <- output %>%
  mutate ( d2 = if_else(site == 2, 1, 0),
           d3 = if_else(site == 3, 1, 0),
           d4 = if_else(site == 4, 1, 0))

output$yhat.prior <- NA

# now get the species specific values:

output$mu_force <- rep(rnorm(nsp, 0, 1), times = nsite*nPhenCombo*nind)
output$mu_photo <- rep(rnorm(nsp,0, 1), times = nsite*nPhenCombo*nind)
output$mu_chill <- rep(rnorm(nsp, 0, 1), times = nsite*nPhenCombo*nind)
output$sigma_force <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)
output$sigma_chill <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)
output$sigma_photo <- rep(rtruncnorm(nsp, a=0, mean = 0, sd = 5), times = nsite*nPhenCombo*nind)

output$sigma_y <- rtruncnorm(nrow(output), a=0, mean = 0, sd = 5)
output$sigma_a <- rep(rtruncnorm(nsp, a=0, mean = 3, sd = 1), times = nsite*nPhenCombo*nind)
output$sigma_site <- rep(rnorm(nsite, 3,1), times = nPhenCombo*nind*nsp)

############################################################################
  b_force <- rnorm(nrep, mu_force, sigma_force)
  b_chill <- rnorm(nrep, mu_force, sigma_force)
  b_photo <- rnorm(nrep, mu_force, sigma_force)

  b_site <- rnorm(nsite, b_site, sigma_site)

  alpha.pheno.sp <- rnorm(nrep, 0, sigma_a)

  alpha.force.sp <- rnorm(nrep, mu_force, sigma_force)

  alpha.photo.sp <- rnorm(nrep, mu_photo, sigma_photo)

  alpha.chill.sp <- rnorm(nrep, mu_chill, sigma_chill)

  alpha.site <- rnorm(nrep, b_site, sigma_site)

  gen.var <- rnorm(nrep, 0, sigma_y)

  for(i in 1:nrep){
    doy.prior <- mu_grand +  alpha.site + alpha.pheno.sp + alpha.force.sp * output$force +
      alpha.chill.sp * output$chill + alpha.photo.sp * output$photo +
      alpha.site*output$d2  + alpha.site*output$d3  + alpha.site*output$d4 +
      gen.var
  
    output$doy.prior[output$rep == i] <- doy.prior[i]
  }
  
  hist(output$doy.prior)
