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
#1. converting species to a factor
pheno$species.fact<-as.numeric(as.factor(pheno$species))

#2. Adding columns of treatments as numeric values
pheno$chill.n<-pheno$chill
pheno$chill.n[pheno$chill.n=="HC"] <- "1"
pheno$chill.n[pheno$chill.n=="LC"] <- "0"
pheno$chill.n<-as.numeric(pheno$chill.n)

pheno$force.n<-pheno$force
pheno$force.n[pheno$force.n=="HF"] <- "1"
pheno$force.n[pheno$force.n=="LF"] <- "0"
pheno$force.n<-as.numeric(pheno$force.n)

pheno$photo.n<-pheno$photo
pheno$photo.n[pheno$photo.n=="HP"] <- "1"
pheno$photo.n[pheno$photo.n=="LP"] <- "0"
pheno$photo.n<-as.numeric(pheno$photo.n)

pheno$site.n<-pheno$population
pheno$site.n[pheno$site.n=="sm"] <- "1"
pheno$site.n[pheno$site.n=="mp"] <- "0"
pheno$site.n<-as.numeric(pheno$site.n)

head(pheno)

#going to split it into analysis of terminal bb and lateral bb
pheno.term<-pheno[,c("tbb","species.fact","chill.n","force.n","photo.n","site.n","species")]
pheno.t<-pheno.term[complete.cases(pheno.term),]


nrow(pheno.term)-nrow(pheno.t)
temp<-subset(pheno.term, is.na(tbb));head(temp); unique(temp$species)
# there were 204 samples that did not have terminal bb

datalist<-with(pheno.t,
                    list( N=nrow(pheno.t),
                          n_sp =length(unique(pheno.t$species.fact)),
                          n_site =length(unique(pheno.t$site.n)),
                          bb = tbb,
                          sp =species.fact,
                          chill = chill.n,
                          photo = photo.n,
                          force = force.n,
                          site=site.n
                    ))

mdl<-stan("stan/bc_bb.stan",
            data= datalist
            ,iter=2000, chains=4, seed=1235)

sm.sum <- summary(mdl)$summary
sm.sum
ssm<- as.shinystan(mdl)
launch_shinystan(ssm)

