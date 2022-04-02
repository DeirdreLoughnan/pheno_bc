# Started in April 2016 #
# Fake data creation by Dan Flynn, updates by Lizzie #

# March 28: 5 chill, 4 sites - no site interactions

# March 31 still not working, going back to basics - 4 sites, 5 chill no interactions at all or ncp
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

nsite = 4
nsp = 15

nwarm = 2
nphoto = 2
nchill = 5

rep = 8 # within each combination of treatments. 

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
site2 = ifelse(site == 2, 1, 0)
site3 = ifelse(site == 3, 1, 0)
site4 = ifelse(site == 4, 1, 0)

#d <- data.frame(site, warm, photo, chill, chill1, chill2, site2, site3, treatcombo) # critical coding error here!
d <- data.frame(site, warm, photo, chill, site2, site3, site4, treatcombo) # critical coding error here!

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

mm <- model.matrix(~site2 + site3 + site4 +(site2 + site3 + site4 +chillP +warm + photo)^2 + (warm + photo + chillP)^2, data.frame(site2, site3, site4, warm, photo, chillP))
# remove last column, chill1 x chill2, empty; remove site2:site3
#mm <- mm[,-grep("chill1:chill2", colnames(mm))]
mm <- mm[,-grep("site2:site3", colnames(mm))]
mm <- mm[,-grep("site2:site4", colnames(mm))]
mm <- mm[,-grep("site3:site4", colnames(mm))]
colnames(mm)
fake <- vector()
values <- vector()
head(mm)

for(i in 1:nsp){ # loop over species, as these are the random effect modeled
  
  # Give species different difference values, drawn from normal. Could make dataframe of diffs and diff.sds, and use apply..
  
  coeff <- c(spint[i], 
             rnorm(1, sitediff, sitediff.sd),
             rnorm(1, sitediff, sitediff.sd),
             rnorm(1, sitediff, sitediff.sd),
             rnorm(1, chill1diff, chill1diff.sd),
             rnorm(1, warmdiff, warmdiff.sd),
             rnorm(1, photodiff, photodiff.sd), 
             rnorm(1, sitechill1, sitechill1.sd),
             rnorm(1, sitewarm, sitewarm.sd),
             rnorm(1, sitephoto, sitephoto.sd),
             rnorm(1, sitechill1, sitechill1.sd),
             rnorm(1, sitewarm, sitewarm.sd),
             rnorm(1, sitephoto, sitephoto.sd),
             rnorm(1, sitechill1, sitechill1.sd),
             rnorm(1, sitewarm, sitewarm.sd),
             rnorm(1, sitephoto, sitephoto.sd),
             rnorm(1, warmchill1, warmchill1.sd),
             rnorm(1, photochill1, photochill1.sd),
             rnorm(1, warmphoto, warmphoto.sd)
             
  )
  values <- rbind(values, coeff)
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 0.1)
  
  #fakex <- data.frame(bb, sp = i, site, warm, photo, chillP)
  fakex <- data.frame(bb, sp = i, site2, site3, site4, chillP, warm, photo)
  
  fake <- rbind(fake, fakex)  
}
head(values)
head(fake)

values <- as.data.frame(values)
names(values) <- c("int","b.site2","b.site3","b.site4","b.chill1","b.force", "b.photo",
                  "b.cs2","b.fs2","b.ps2","b.cs3","b.fs3","b.ps3","b.cs4","b.fs4","b.ps4",
                  #"b.fs2","b.ps2","b.fs3","b.ps3","b.fs4","b.ps4",
                  "b.fc", "b.pc","b.fp"); head(values)
#summary(lm(bb ~ (site + warm + photo + chillP)^2, data = fake)) # sanity check 
write.csv(values, "values_4sites_simple_allint.csv", row.names = F) 
#values <-read.csv( "values_4sites_simple_allint.csv")
# plot(fake$bb ~ fake$chillP)
# abline(35,-5)

# fake$count <- 1
# faket <- aggregate(fake ["count"],
#                    fake[c("sp","photo","chillP","warm")],
#                    FUN = sum)

# now fix the levels to 0/1 (not 1/2) as R does
# fake$site <- as.numeric(fake$site)
# fake$site[fake$site==1] <- 0
# fake$site[fake$site==2] <- 1

fake$warm <- as.numeric(fake$warm)
fake$warm[fake$warm==1] <- 0
fake$warm[fake$warm==2] <- 1

fake$photo <- as.numeric(fake$photo)
fake$photo[fake$photo==1] <- 0
fake$photo[fake$photo==2] <- 1

#summary(lm(bb ~ (site2 + site3 + site4 + warm+photo+chillP)^2, data = fake)) # double sanity check 

# summary(lm(bb ~ warm+photo+chillP+site2 + site3, data = fake)) # double sanity check 

#summary(lmer(bb ~ (site|sp) + (warm|sp) + (photo|sp) + (chill1|sp) + (chill2|sp), data = fake)) # too hard for lmer.

#save(list=c("fake"), file = "Fake Budburst.RData")


# datalist <- list( N=nrow(fake),
#                   n_sp = length(unique(fake$sp)),
#                   n_site = length(unique(fake$site)),
#                   bb = fake$bb,
#                   sp = as.numeric(fake$sp),
#                   chill = as.numeric(fake$chill),
#                   photo = as.numeric(fake$photo),
#                   force = as.numeric(fake$warm),
#                    # site = as.numeric(fake$site))
#                   site2 = as.numeric(fake$site2),
#                   site3 = as.numeric(fake$site3))
datalist <- list( N=nrow(fake),
                  n_sp = length(unique(fake$sp)),
                  n_site = length(unique(fake$site)),
                  lday = fake$bb,
                  sp = as.numeric(fake$sp),
                  chill1 = as.numeric(fake$chillP),
                  photo = as.numeric(fake$photo),
                  warm = as.numeric(fake$warm),
                  site = as.numeric(fake$site),
                  site2 = as.numeric(fake$site2),
                  site3 = as.numeric(fake$site3),
                  site4 = as.numeric(fake$site4))


mdl.full <- stan("stan/df_mdl_4sites_again_allint.stan",
                 data = datalist,
                 include = FALSE, pars = c("ypred_new","y_hat"),
                 iter = 2000, chains= 4, warmup = 1000)
save(mdl.full, file = "output/df_noncp_4sites_simple_allint.Rda")


# load("output/df_noncp_4sites_simple_chillforceint.Rda")
# # 
# ssm <-  as.shinystan(mdl.full)
# launch_shinystan(mdl.full)


summary(mdl.full)$summary[c("mu_a", "mu_b_warm","mu_b_chill1","mu_b_photo","b_site2","b_site3","b_site4",
                            "mu_b_inter_wp","mu_b_inter_wc1","mu_b_inter_pc1","mu_b_inter_s2c1","mu_b_inter_s3c1","mu_b_inter_s4c1","mu_b_inter_ws2","mu_b_inter_ws3","mu_b_inter_ws4","mu_b_inter_ps2","mu_b_inter_ps3","mu_b_inter_ps4",
                            "sigma_b_warm","sigma_b_chill1","sigma_b_photo","sigma_a",
                         "sigma_b_inter_wp","sigma_b_inter_wc1","sigma_b_inter_pc1",
                         "sigma_b_inter_ws2","sigma_b_inter_ws3","sigma_b_inter_ws4","sigma_b_inter_ps2","sigma_b_inter_ps3","sigma_b_inter_ps4","sigma_b_inter_s2c1","sigma_b_inter_s3c1","sigma_b_inter_s3c1",
                            "sigma_y"),c("mean","25%","75%")]
#  
sum <- summary(mdl.full)$summary 
intSp <- as.data.frame(sum[grep("a_sp", rownames(sum)), ])
slopeF <- as.data.frame(sum[grep("b_warm", rownames(sum)), ]); slopeF <- slopeF[c(2:16),]
slopeC <- as.data.frame(sum[grep("b_chill1", rownames(sum)), ]) ; slopeC <- slopeC[c(2:16),]
slopeP <- as.data.frame(sum[grep("b_photo", rownames(sum)), ]) ; slopeP <- slopeP[c(2:16),]

slopeFP <- as.data.frame(sum[grep("b_inter_wp", rownames(sum)), ]) ; slopeFP <- slopeFP[c(2:16),]
slopeFC <- as.data.frame(sum[grep("b_inter_wc1", rownames(sum)), ]) ; slopeFC <- slopeFC[c(2:16),]
slopeCP <- as.data.frame(sum[grep("b_inter_pc1", rownames(sum)), ]) ; slopeCP <- slopeCP[c(2:16),]
slopeFS2 <- as.data.frame(sum[grep("b_inter_ws2", rownames(sum)), ]) ; slopeFS2 <- slopeFS2[c(2:16),]
slopePS2 <- as.data.frame(sum[grep("b_inter_ps2", rownames(sum)), ]) ; slopePS2 <- slopePS2[c(2:16),]
 slopeCS2 <- as.data.frame(sum[grep("b_inter_s2c1", rownames(sum)), ]) ; slopeCS2 <- slopeCS2[c(2:16),]
slopeFS3 <- as.data.frame(sum[grep("b_inter_ws3", rownames(sum)), ]) ; slopeFS3 <- slopeFS3[c(2:16),]
slopePS3 <- as.data.frame(sum[grep("b_inter_ps3", rownames(sum)), ]) ; slopePS3 <- slopePS3[c(2:16),]
 slopeCS3 <- as.data.frame(sum[grep("b_inter_s3c1", rownames(sum)), ]) ; slopeCS3 <- slopeCS3[c(2:16),]
 slopeCS4 <- as.data.frame(sum[grep("b_inter_s4c1", rownames(sum)), ]) ; slopeCS4 <- slopeCS4[c(2:16),]
 slopeFS4 <- as.data.frame(sum[grep("b_inter_ws4", rownames(sum)), ]) ; slopeFS4 <- slopeFS4[c(2:16),]
 slopePS4 <- as.data.frame(sum[grep("b_inter_ps4", rownames(sum)), ]) ; slopePS4 <- slopePS4[c(2:16),]
 
# 
# 
slopeSite2 <- as.data.frame(sum[grep("b_site2", rownames(sum)), "mean"]) ; mean(values$b.site2); slopeSite2
slopeSite3 <- as.data.frame(sum[grep("b_site3", rownames(sum)), "mean"]) ; mean(values$b.site3); slopeSite3
slopeSite4 <- as.data.frame(sum[grep("b_site4", rownames(sum)), "mean"]) ; mean(values$b.site4); slopeSite4

pdf("fake_esti_correlations_4site_someint_allint.pdf", height = 10, width = 15)
par(mfrow = c(4,5))
plot(values$int ~ intSp$mean, pch =19);abline(0,1)
plot( values$b.force ~ slopeF$mean, pch = 19);abline(0,1)
plot(values$b.chill1 ~ slopeC$mean, pch = 19);abline(0,1)
plot(values$b.photo ~ slopeP$mean, pch = 19);abline(0,1)
plot(values$b.fc ~ slopeFC$mean, pch = 19);abline(0,1)
plot(values$b.pc ~ slopeCP$mean, pch = 19);abline(0,1)
plot(values$b.fp ~ slopeFP$mean, pch = 19);abline(0,1)
plot(values$b.fs2 ~ slopeFS2$mean, pch = 19);abline(0,1)
 plot(values$b.cs2 ~ slopeCS2$mean, pch = 19);abline(0,1)
plot(values$b.ps2 ~ slopePS2$mean, pch = 19);abline(0,1)
plot(values$b.fs3 ~ slopeFS3$mean, pch = 19);abline(0,1)
 plot(values$b.cs3 ~ slopeCS3$mean, pch = 19);abline(0,1)
plot(values$b.ps3 ~ slopePS3$mean, pch = 19);abline(0,1)
plot(values$b.fs4 ~ slopeFS4$mean, pch = 19);abline(0,1)
 plot(values$b.cs4 ~ slopeCS4$mean, pch = 19);abline(0,1)
plot(values$b.ps4 ~ slopePS4$mean, pch = 19);abline(0,1)
dev.off()

# #check that values fall in the intervals:
# b_site_int <- sum[grep("mu_b_site", rownames(sum)), ]
# sigma_b_site <- sum[grep("sigma_b_site", rownames(sum)), ]
# b_intr_sc <- sum[grep("mu_b_inter_sc1", rownames(sum)), ]
# 
# sum[c("mu_a","mu_b_warm","mu_b_chill1","mu_b_photo","mu_b_site","mu_b_inter_wp","mu_b_inter_wc1","mu_b_inter_pc1","mu_b_inter_sc1","mu_b_inter_ws","mu_b_inter_ps","sigma_a","sigma_b_warm","sigma_b_chill1","sigma_b_photo","sigma_b_site","sigma_b_inter_wp","sigma_b_inter_wc1","sigma_b_inter_pc1","sigma_b_inter_sc1","sigma_b_inter_ws","sigma_b_inter_ps","sigma_y"),c(1,5,7)]
# 
# sum[c("mu_grand","mu_force","mu_chill1","mu_photo","b_site","mu_fp","mu_fc","mu_cp","sigma_a","sigma_force","sigma_chill1","sigma_photo","sigma_fp","sigma_fc","sigma_cp","sigma_y"),c(1,5,7)]
# 
# # Sigma's corrected for ncp?:
# 
# slopeFP <- as.data.frame(sum[grep("b_inter_wp", rownames(sum)), ]) ; slopeFP <- slopeFP[c(38:72),]
# slopeFC <- as.data.frame(sum[grep("b_inter_wc1", rownames(sum)), ]) ; slopeFC <- slopeFC[c(38:72),]
# slopeCP <- as.data.frame(sum[grep("b_inter_pc1", rownames(sum)), ]) ; slopeCP <- slopeCP[c(38:72),]
# slopeFS <- as.data.frame(sum[grep("b_inter_ws", rownames(sum)), ]) ; slopeFS <- slopeFS[c(38:72),]
# slopePS <- as.data.frame(sum[grep("b_inter_ps", rownames(sum)), ]) ; slopePS <- slopePS[c(38:72),]
# slopeCS <- as.data.frame(sum[grep("b_inter_sc1", rownames(sum)), ]) ; slopeCS <- slopeCS[c(38:72),]
# 
# ncpslopeFP <- as.data.frame(sum[grep("b_inter_wp_ncp", rownames(sum)), ]) 
# ncpslopeFC <- as.data.frame(sum[grep("b_inter_wc1_ncp", rownames(sum)), ])
# ncpslopeCP <- as.data.frame(sum[grep("b_inter_pc1_ncp", rownames(sum)), ]) 
# ncpslopeFS <- as.data.frame(sum[grep("b_inter_ws_ncp", rownames(sum)), ]) 
# ncpslopePS <- as.data.frame(sum[grep("b_inter_ps_ncp", rownames(sum)), ]) 
# ncpslopeCS <- as.data.frame(sum[grep("b_inter_sc1_ncp", rownames(sum)), ]) 
# 
# sigmaFP <- as.data.frame(sum[grep("sigma_b_inter_wp", rownames(sum)),]) 
# sigmaFC <- as.data.frame(sum[grep("sigma_b_inter_wc1", rownames(sum)), ]) 
# sigmaCP <- as.data.frame(sum[grep("sigma_b_inter_pc1", rownames(sum)), ]) 
# sigmaFS <- as.data.frame(sum[grep("sigma_b_inter_ws", rownames(sum)), ]) 
# sigmaPS <- as.data.frame(sum[grep("sigma_b_inter_ps", rownames(sum)), ]) 
# sigmaCS <- as.data.frame(sum[grep("sigma_b_inter_sc1", rownames(sum)), ]) 
# 
# muFP <- as.data.frame(sum[grep("mu_b_inter_wp", rownames(sum)),]) 
# sigmaFC <- as.data.frame(sum[grep("sigma_b_inter_wc1", rownames(sum)), ]) 
# sigmaCP <- as.data.frame(sum[grep("sigma_b_inter_pc1", rownames(sum)), ]) 
# sigmaFS <- as.data.frame(sum[grep("sigma_b_inter_ws", rownames(sum)), ]) 
# sigmaPS <- as.data.frame(sum[grep("sigma_b_inter_ps", rownames(sum)), ]) 
# sigmaCS <- as.data.frame(sum[grep("sigma_b_inter_sc1", rownames(sum)), ]) 
# 
# 
# full_sigma_fp <- ncpslopeFP$mean * sigmaFP["mean",]; mean(full_sigma_fp)
# full_sigma_fc <- ncpslopeFC$mean * sigmaFC["mean",]; mean(full_sigma_fc)
# 
# posterior <- rstan::extract(mdl.full)

# pdf("figures/lnc_prior_post_dist.pdf", width = 15, height = 25)
# par(mfrow = c(4,4))
# #plot priors against posteriors
# h1 <- hist(rnorm(1000, 0,35), col = rgb(1,0,1,1/4))
# hist(posterior$mu_b_warm,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,35), col = rgb(1,0,1,1/4))
# hist(posterior$mu_b_chill1,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,35), col = rgb(1,0,1,1/4))
# hist(posterior$mu_b_photo,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,35), col = rgb(1,0,1,1/4))
# hist(posterior$mu_b_site,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_warm,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_chill1,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_photo,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_site,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_inter_wp,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_inter_wc1,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_inter_pc1,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_inter_sc1,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_inter_ws,add=T,col=rgb(0,0,1,1/4))
# 
# h1 <- hist(rnorm(1000, 0,10), col = rgb(1,0,1,1/4))
# hist(posterior$sigma_b_inter_ps,add=T,col=rgb(0,0,1,1/4))
# # dev.off()
# 

## Posterior predictive checks: http://avehtari.github.io/BDA_R_demos/demos_rstan/ppc/poisson-ppc.html
fit <- mdl.full
y_rep <- as.matrix(fit, pars = "y_rep")

y <- fake$bb
#ppc_hist(y, y_rep[1:8, ], binwidth = 1)

pdf("yvsypred_4site_5chill_siteint.pdf")
ppc_dens_overlay(y, y_rep[1:100, ])
dev.off()

posterior <- as.matrix(fit)
mcmc_areas(posterior, 
           pars = c("mu_a","mu_b_warm","mu_b_chill1","mu_b_photo","mu_b_site","mu_b_inter_wp","mu_b_inter_wc1","mu_b_inter_pc1","mu_b_inter_sc1","mu_b_inter_ws","mu_b_inter_ps"),
           prob = 0.8) 

mcmc_areas(posterior, 
           pars = c("mu_grand","mu_force","mu_chill1","mu_photo","b_site","mu_fp","mu_fc","mu_cp"),
           prob = 0.8) 

mcmc_areas(posterior, 
           pars = c("sigma_a","sigma_force","sigma_chill1","sigma_photo","sigma_fp","sigma_fc","sigma_cp","sigma_y"),
           prob = 0.8) 


#,