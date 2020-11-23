## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

library(scales)
library(arm)
library(rstan)
library(shinystan)
# library(sjPlot)
library(rstanarm)
library(RColorBrewer)

options(mc.cores = parallel::detectCores())

rm(list=ls()) 
options(stringsAsFactors = FALSE)

if(length(grep("deirdreloughnan", getwd())>0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("~/deirdre/") # for midge
}

source('rcode/cleaning/pheno_bb_calc.R')
head(pheno)

############################################################
# Preping the data for the model
pheno$species.fact<-as.numeric(as.factor(pheno$species))

pheno<-pheno[complete.cases(pheno),]

datalist<-with(pheno,
                    list( N=nrow(pheno),
                          n_sp =length(unique(pheno$species)),
                          bb = tbb,
                          sp =species.fact,
                          chill = chill,
                          photo = photo,
                          force = force,
                          site=population
                    ))

mdl<-stan("stan/bc_bb.stan",
            data= datalist
            ,iter=2000, chains=1, seed=1235)

