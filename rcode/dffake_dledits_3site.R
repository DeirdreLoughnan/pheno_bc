# Started in April 2016 #
# Fake data creation by Dan Flynn, updates by Lizzie #

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
# Set up: same as experiment, with two sites, 28 species, two levels each of warming and photoperiod, and three levels of chilling. 2016-04-01 adding interactions. This ends up generating expected differences, but variation in effect sizes across species is minimal currently.
# modifying it for western analysis:
#1. change chilling to chill portions
#2. change site number, initially try having 3 sites
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>

nsite = 2
nsp = 35

nwarm = 2
nphoto = 2
nchill = 5

rep = 12 # within each combination of treatments. 

(ntot = nsite*nwarm*nphoto*nchill*rep) # 792 rows; 22k rows across species

# Build up the data frame
site = gl(nsite, rep, length = ntot)

warm = gl(nwarm, rep*nsite, length = ntot)
photo = gl(nphoto, rep*nsite*nwarm, length = ntot)
chill = gl(nchill, rep*nsite*nwarm*nphoto, length = ntot)
#chillP <- rep(c(10,20), each = rep*nsite*nwarm*nphoto)
chillP <- rep(c(10,15,20,25,30), each = rep*nsite*nwarm*nphoto)

#chill1 = ifelse(chill == 2, 1, 0) 
#chill2 = ifelse(chill == 3, 1, 0) 


treatcombo = paste(warm, photo, chill, sep = "_")

# setting up the dummy variables for site:
# site2 = ifelse(site == 2, 1, 0)
# site3 = ifelse(site == 3, 1, 0)

#d <- data.frame(site, warm, photo, chill, chill1, chill2, site2, site3, treatcombo) # critical coding error here!
d <- data.frame(site, warm, photo, chill,  treatcombo) # critical coding error here!

head(d)


###### Set up differences for each level
sitediff = 2 
warmdiff = -20 # days earlier from 1 to 2
photodiff = -14
chill1diff = -5
# chill2diff = -19

# interactions. 9 two-way interactions
sitewarm = 0
sitephoto = 0
sitechill1 = -1 # similar to stan results
# sitechill2 = -2
warmphoto = 3.5 # positive 3.5. So at the warm level, the effect of longer days is muted by 3.5 days.
warmchill1 = 11 # both positive ~ 10. 
# warmchill2 = 9
photochill1 = 0.1 # from stan results
# photochill2 = 1


######## SD for each treatment
sitediff.sd = 1.5 
warmdiff.sd = 1 
photodiff.sd = 1
chill1diff.sd = 1.5
# chill2diff.sd = 2

# interactions. 9 two-way interactions
sitewarm.sd = 1
sitephoto.sd = 1
sitechill1.sd = 2 
# sitechill2.sd = 2
warmphoto.sd = 1
warmchill1.sd = 1.5
# warmchill2.sd = 1.5
photochill1.sd = 1
# photochill2.sd = 1


##### Again, now with species now.

baseinter = 35 # baseline intercept across all species 
spint <- baseinter + c(1:nsp)-mean(1:nsp) # different intercepts by species

mm <- model.matrix(~(site+ warm + photo + chillP)^2, data.frame(site, warm, photo, chillP))
# remove last column, chill1 x chill2, empty; remove site2:site3
#mm <- mm[,-grep("chill1:chill2", colnames(mm))]
#mm <- mm[,-grep("site2:site3", colnames(mm))]
colnames(mm)
fake <- vector()
values <- vector()
head(mm)

for(i in 1:nsp){ # loop over species, as these are the random effect modeled
  
  # Give species different difference values, drawn from normal. Could make dataframe of diffs and diff.sds, and use apply..
  
  coeff <- c(spint[i], 
             rnorm(1, sitediff, sitediff.sd),
             #rnorm(1, sitediff, sitediff.sd),
             rnorm(1, warmdiff, warmdiff.sd),
             rnorm(1, photodiff, photodiff.sd), 
             rnorm(1, chill1diff, chill1diff.sd),
             # rnorm(1, chill2diff, chill2diff.sd), 
             rnorm(1, sitewarm, sitewarm.sd), 
            # rnorm(1, sitewarm, sitewarm.sd), 
             rnorm(1, sitephoto, sitephoto.sd),
            # rnorm(1, sitephoto, sitephoto.sd),
             rnorm(1, sitechill1, sitechill1.sd),
             # rnorm(1, sitechill2, sitechill2.sd),
             rnorm(1, warmphoto, warmphoto.sd),
             rnorm(1, warmchill1, warmchill1.sd),
             # rnorm(1, warmchill2, warmchill2.sd),
             rnorm(1, photochill1, photochill1.sd)
             # rnorm(1, photochill2, photochill2.sd)
  )
  values <- rbind(values, coeff)
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 0.1)
  
  fakex <- data.frame(bb, sp = i, site, warm, photo, chillP)
  #fakex <- data.frame(bb, sp = i, site2, site3, warm, photo, chill)
  
  fake <- rbind(fake, fakex)  
}
head(values)
head(fake)

values <- as.data.frame(values)
names(values) <- c("int","b.site","b.force", "b.photo","b.chill1","b.fs","b.ps","b.cs","b.fp","b.fc", "b.pc")
#summary(lm(bb ~ (site + warm + photo + chillP)^2, data = fake)) # sanity check 
# 
# plot(fake$bb ~ fake$chillP)
# abline(35,-5)

# fake$count <- 1
# faket <- aggregate(fake ["count"],
#                    fake[c("sp","photo","chillP","warm")],
#                    FUN = sum)

# now fix the levels to 0/1 (not 1/2) as R does
fake$site <- as.numeric(fake$site)
fake$site[fake$site==1] <- 0
fake$site[fake$site==2] <- 1

fake$warm <- as.numeric(fake$warm)
fake$warm[fake$warm==1] <- 0
fake$warm[fake$warm==2] <- 1

fake$photo <- as.numeric(fake$photo)
fake$photo[fake$photo==1] <- 0
fake$photo[fake$photo==2] <- 1

#summary(lm(bb ~ (site+warm+photo+chill)^2, data = fake)) # double sanity check 

#summary(lmer(bb ~ (site|sp) + (warm|sp) + (photo|sp) + (chill1|sp) + (chill2|sp), data = fake)) # too hard for lmer.

#save(list=c("fake"), file = "Fake Budburst.RData")


datalist <- list( N=nrow(fake),
                  n_sp = length(unique(fake$sp)),
                  n_site = length(unique(fake$site)),
                  lday = fake$bb,
                  sp = as.numeric(fake$sp),
                  chill1 = as.numeric(fake$chill),
                  photo = as.numeric(fake$photo),
                  warm = as.numeric(fake$warm),
                   site = as.numeric(fake$site))
                  #,site2 = as.numeric(fake$site2),
                  # site3 = as.numeric(fake$site3))


#datalist$site

mdl.full <- stan("stan/df_mdl.stan",
                 data = datalist,
                 include = FALSE, pars = c("ypred_new","y_hat"),
                 iter = 4000, chains= 4, warmup = 2000)
save(mdl.full, file = "output/df_mdl_morespp.Rda")

mdl.ind <- stan("stan/df_interdirect.stan",
                 data = datalist,
                 include = FALSE, pars = c("ypred_new","y_hat"),
                 iter = 4000, chains= 4, warmup = 2000)
save(mdl.ind, file = "output/df_interdirect_morespp.Rda")

#load("output/df_interdirect.Rda")
# #  load("output/tbb_ncp_chillportions_zsc_dl.Rda")
# # #
# ssm <-  as.shinystan(mdl.full)
# launch_shinystan(ssm)
# # # 
sum <- summary(mdl.full)$summary 
summary(mdl.full)$summary[c("mu_a", "mu_b_warm","mu_b_chill1","mu_b_photo","mu_b_site","mu_b_inter_wp","mu_b_inter_ps","mu_b_inter_wc1","mu_b_inter_pc1","mu_b_inter_sc1","sigma_b_warm","sigma_b_chill1","sigma_b_photo","sigma_b_site","sigma_a","sigma_b_inter_wp","sigma_b_inter_ps","sigma_b_inter_wc1","sigma_b_inter_pc1","sigma_b_inter_sc1", "sigma_y"),"mean"]

# compare to fake data:

intSp <- as.data.frame(sum[grep("a_sp", rownames(sum)), ]) 
slopeF <- as.data.frame(sum[grep("b_warm", rownames(sum)), ]); slopeF <- slopeF[c(1:28),]
slopeC <- as.data.frame(sum[grep("b_chill1", rownames(sum)), ]) ; slopeC <- slopeC[c(1:28),]
slopeP <- as.data.frame(sum[grep("b_photo", rownames(sum)), ]) ; slopeP <- slopeP[c(1:28),]
slopeSite <- as.data.frame(sum[grep("b_site", rownames(sum)), ]) ; slopeSite <- slopeSite[c(1:28),]
slopeFP <- as.data.frame(sum[grep("b_inter_wp", rownames(sum)), ]) ; slopeFP <- slopeFP[c(1:28),]
slopeFC <- as.data.frame(sum[grep("b_inter_wc1", rownames(sum)), ]) ; slopeFC <- slopeFC[c(1:28),]
slopeCP <- as.data.frame(sum[grep("b_inter_pc1", rownames(sum)), ]) ; slopeCP <- slopeCP[c(1:28),]
slopeFS <- as.data.frame(sum[grep("b_inter_ws", rownames(sum)), ]) ; slopeFS <- slopeFS[c(1:28),]
slopePS <- as.data.frame(sum[grep("b_inter_ps", rownames(sum)), ]) ; slopePS <- slopePS[c(1:28),]
slopeCS <- as.data.frame(sum[grep("b_inter_sc1", rownames(sum)), ]) ; slopeCS <- slopeCS[c(1:28),]

pdf("fake_esti_correlations_norspp.pdf", height = 10, width = 5)
par(mfrow = c(6,2))
plot(values$int ~ intSp$mean, pch =19)
plot( values$b.force ~ slopeF$mean, pch = 19)
plot(values$b.chill1 ~ slopeC$mean, pch = 19)
plot(values$b.photo ~ slopeP$mean, pch = 19)
plot(values$b.fc ~ slopeFC$mean, pch = 19)
plot(values$b.pc ~ slopeCP$mean, pch = 19)
plot(values$b.fp ~ slopeFP$mean, pch = 19)
plot(values$b.fs ~ slopeFS$mean, pch = 19)
plot(values$b.cs ~ slopeCS$mean, pch = 19)
plot(values$b.ps ~ slopePS$mean, pch = 19)
plot(values$b.site ~ slopeSite$mean, pch = 19)
dev.off()

#check that values fall in the intervals:
b_site_int <- sum[grep("mu_b_site", rownames(sum)), ]
sigma_b_site <- sum[grep("sigma_b_site", rownames(sum)), ]
b_intr_sc <- sum[grep("mu_b_inter_sc1", rownames(sum)), ]

sum[c(309:331),]
