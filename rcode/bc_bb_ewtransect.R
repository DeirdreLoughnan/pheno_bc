## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

# Adding standardization: following the methods outlined by Andrew Gelman in his paper

# changing from site to transect, with transect added as a 0/1 dummy variable, west = 0, east = 1

#Feb 28: model consistently has 7 divergent transitions, these persisted when I increased the iterations to 3000 warmup, 4000 iterations, and increased the adapt delta to 0.99 - from the log posterior vs parameter plots, the worst are sigma_chill with them all butting up against one side, sigma_b_inter_pc -- I am going to try ncp for chilling as well as photo and see if that helps

# March 1: the ncp for chilling produced 79 divergent transitions. The fact that this made it worse got me thinking that maybe the issue is the model has more complexity then it needs, so I removed the ncp on photo

# the model with no ncp on photo or chilling produced dt, running it for 5000 warmup and 6000 iterations produced only 1 dt but the rhat were bad at 1.13

# next I tried running no ncp on photo and chill, but with a bigger adapt delta, 4000:3000 iterations, there were no dt, but the rhat was really bad at 1.39 and the max tree depth exceeded
# 
#library(scales)
#library(arm)
library(rstan)
library(shinystan)
#library(reshape2)
library(bayesplot)
library(ggplot2)
#library(RColorBrewer)
library(dplyr)
library(plyr)
#library(tidybayes)

options(mc.cores = parallel::detectCores())

rm(list=ls()) 
options(stringsAsFactors = FALSE)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

#source('rcode/cleaning/pheno_bb_calc.R')
# head(pheno)
# length(unique(pheno$lab2))

df <- read.csv("input/day.of.bb.DFlynn.chill0.csv", header=TRUE, na.strings=c("","NA"))
head(df)
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill$transect <- "east"

dl <- read.csv("input/day.of.bb.DL.csv", header=TRUE, na.strings=c("","NA"))
head(dl)
dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")
dl.wchill$transect <- "west"

# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)
#pheno <- dl

pheno$first <- ifelse(pheno$tbb < pheno$latbb1,"t", ifelse (pheno$tbb == pheno$latbb1,"tl", "l"))

head(pheno)
table(pheno$species, pheno$first)
#write.csv(pheno, "input/pheno.w5chill.csv")
#pheno <- read.csv("pheno.w5chill.csv")
############################################################
#convert forcing and photoperiod treatments into binary
pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

# convert the transect into a dummy variable:
pheno$transect.n <- pheno$transect
pheno$transect.n[pheno$transect.n == "east"] <- "1"
pheno$transect.n[pheno$transect.n == "west"] <- "0"
pheno$transect.n <- as.numeric(pheno$transect.n)


pheno$Chill_portions <- as.numeric(pheno$Chill_portions)

# pheno$force.z <- (pheno$force.t-mean(pheno$force.t,na.rm=TRUE))/sd(pheno$force.t,na.rm=TRUE)
# pheno$photo.z <- (pheno$photo.t-mean(pheno$photo.t,na.rm=TRUE))/sd(pheno$photo.t,na.rm=TRUE)
# pheno$chillport.z <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/sd(pheno$Chill_portions,na.rm=TRUE)
# z scoring values, but since some are binary I will use 2SD as per the Gelmen et al paper
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)
pheno$utah.z2 <- (pheno$Utah_Model-mean(pheno$Utah_Model,na.rm=TRUE))/(sd(pheno$Utah_Model,na.rm=TRUE)*2)

#going to split it into analyses of terminal bb and lateral bb
# Starting with the terminal buds:
pheno.term <- pheno[,c("tbb", "force.n", "photo.n", "transect.n", "species", "lab2","Utah_Model","Chill_portions","force.z2", "photo.z2", "chillport.z2", "utah.z2")]

pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 

nrow(pheno.term) - nrow(pheno.t) # That had no terminal bb, 609


# head(pheno.t)
datalist_z <- list( N=nrow(pheno.t),
                    Nsite = nrow(pheno.t),
                    n_sp = length(unique(pheno.t$species.fact)),
                    n_site = length(unique(pheno.t$site.n)),
                    bb = pheno.t$tbb,
                    sp = pheno.t$species.fact,
                    chill = pheno.t$chillport.z2,
                    photo = pheno.t$photo.z2,
                    force = pheno.t$force.z2,
                    site = pheno.t$transect.n)

mdl <- stan("stan/bc.bb.noncp.standardize.old.stan",
            data = datalist_z,
            include = FALSE, pars = c("ypred_new","y_hat"),
            iter = 6000, warmup = 5000, chains=4,
            control = list(adapt_delta= 0.99))
save(mdl, file="output/tbb_cport_stnd_transect_noncp_moreit.Rda")
# 
mdl <- stan("stan/bc.bb.noncp.standardize.old.stan",
            data = datalist_z,
            include = FALSE, pars = c("ypred_new","y_hat"),
            iter = 4000, warmup = 3000, chains=4,
            control = list(adapt_delta= 0.99, max_treedepth = 15))
save(mdl, file="output/tbb_cport_stnd_transect_noncpadapt.Rda")

# mdl <- stan("stan/bc.bb.ncpphotochill.ncpinter.standardize.old.stan",
#               data = datalist_z,
#             include = FALSE, pars = c("ypred_new","y_hat"),
#               iter = 4000, warmup = 3000, chains=4)
# save(mdl, file="output/tbb_cport_stnd_transect_ncpchillphoto.Rda")

# mdl.i <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_index.stan",
#               data = datalist_z,
#               iter = 4000, chains=4)
# save(mdl.i, file="output/tbb_cport_stnd_index.Rds")
# 
# 
# mdl.d <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_dummy.stan",
#               data = datalist_z,
#               iter = 4000, chains=4)
# save(mdl.d, file="output/tbb_cport_stnd_dummy.Rds")
#######################################################################

# load("output/tbb_cport_stnd_transect_ncpchillphoto.Rda")
# # #  load("output/tbb_ncp_chillportions_zsc_dl.Rda")
# # # #
# ssm <-  as.shinystan(mdl)
# launch_shinystan(ssm)
# # #
sum <- summary(mdl)$summary 
summary(mdl)$summary[c("mu_a", "mu_force","mu_chill","mu_photo","mu_site","mu_inter_fp", "mu_inter_fs","mu_inter_ps","mu_inter_fc","mu_inter_pc","mu_inter_sc","sigma_a", "sigma_force","sigma_chill","sigma_photo","sigma_site","sigma_b_inter_fp","sigma_b_inter_fs","sigma_b_inter_ps","sigma_b_inter_fc","sigma_b_inter_pc","sigma_b_inter_sc", "sigma_y"),"mean"]

# pairs(mdl, pars = c("mu_a", "mu_force","mu_chill","mu_photo","mu_site","mu_inter_fp", "mu_inter_fs","mu_inter_ps","mu_inter_fc","mu_inter_pc","mu_inter_sc", "lp__"))
# 
# pairs(mdl, pars = c("sigma_a", "sigma_force","sigma_chill","sigma_photo","sigma_site","sigma_b_inter_fp","sigma_b_inter_fs","sigma_b_inter_ps","sigma_b_inter_fc","sigma_b_inter_pc","sigma_b_inter_sc", "sigma_y", "lp__"))



# post <- rstan::extract(mdl)
# # histograms
# par(mfrow = c(1,1))
# h1 <- hist(rnorm(1000, 1, 35))
# h2 <- hist(post$mu_a)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_force
# h1 <- hist(rnorm(1000, 0, 50))
# h2 <- hist(post$mu_force)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-150,150))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_chill
# h1 <- hist(rnorm(1000, 0, 35))
# h2 <- hist(post$mu_chill)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_photo
# h1 <- hist(rnorm(1000, 0, 35))
# h2 <- hist(post$mu_photo)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_site
# h1 <- hist(rnorm(1000, 1, 35))
# h2 <- hist(post$mu_site)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_inter_fp
# h1 <- hist(rnorm(1000, 0, 35))
# h2 <- hist(post$mu_inter_fp)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_inter_fs
# h1 <- hist(rnorm(1000, 0, 35))
# h2 <- hist(post$mu_inter_fs)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_inter_fc
# h1 <- hist(rnorm(1000, 0, 35))
# h2 <- hist(post$mu_inter_fc)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_inter_ps
# h1 <- hist(rnorm(1000, 0, 35))
# h2 <- hist(post$mu_inter_ps)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_inter_pc
# h1 <- hist(rnorm(1000, 0, 35))
# h2 <- hist(post$mu_inter_pc)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # mu_inter_sc
# h1 <- hist(rnorm(1000, 0, 35))
# h2 <- hist(post$mu_inter_sc)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # sigma_a
# h1 <- hist(rnorm(1000, 0, 10))
# h2 <- hist(post$sigma_a)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-50,50))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # sigma_force
# h1 <- hist(rnorm(1000, 0, 10))
# h2 <- hist(post$sigma_force)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-50,50))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # sigma_chill
# h1 <- hist(rnorm(1000, 0, 40))
# h2 <- hist(post$sigma_chill)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # sigma_photo
# h1 <- hist(rnorm(1000, 0, 10))
# h2 <- hist(post$sigma_photo)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-50,50))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # sigma_site
# h1 <- hist(rnorm(1000, 0, 40))
# h2 <- hist(post$sigma_site)
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# # sigma_interaction
# h1 <- hist(rnorm(1000, 1, 10))
# h2 <- hist(post$sigma_b_inter_fp, xlim =c(-10,10))
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# h1 <- hist(rnorm(1000, 1, 10))
# h2 <- hist(post$sigma_b_inter_fs, xlim =c(-10,10))
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# h1 <- hist(rnorm(1000, 1, 10))
# h2 <- hist(post$sigma_b_inter_ps, xlim =c(-10,10))
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# h1 <- hist(rnorm(1000, 1, 10))
# h2 <- hist(post$sigma_b_inter_fc, xlim =c(-10,10))
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# h1 <- hist(rnorm(1000, 1, 10))
# h2 <- hist(post$sigma_b_inter_pc, xlim =c(-10,10))
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# h1 <- hist(rnorm(1000, 1, 20))
# h2 <- hist(post$sigma_b_inter_sc, xlim =c(-10,10))
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
# h1 <- hist(rnorm(1000, 1, 10))
# h2 <- hist(post$sigma_y, xlim =c(-10,10))
# plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
# plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
# 
