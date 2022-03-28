## Started March 2016 ##
## By Dan Flynn ##
## Updates by Lizzie in 2017 ##

# Fake data testing of pheno budburst experiment 2015

library(rstan)
library(xtable)
library(ggplot2)
library(shinystan)
library(tibble)
#library(tidybayes)

setwd("/home/deirdre/pheno_bc")
# setwd("~/Documents/git/buds/analyses")
# setwd("~/Documents/git/projects/treegarden/budexperiments/analyses")
# 
# source('stan/savestan.R')
# source('source/plotlet.R')
# get latest .Rdata file

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

nsite = 2
nsp = 28

nwarm = 2
nphoto = 2
nchill = 3

rep = 10 # within each combination of treatments. 

(ntot = nsite*nwarm*nphoto*nchill*rep) # 792 rows; 22k rows across species

# Build up the data frame
site = gl(nsite, rep, length = ntot)

warm = gl(nwarm, rep*nsite, length = ntot)
photo = gl(nphoto, rep*nsite*nwarm, length = ntot)
chill = gl(nchill, rep*nsite*nwarm*nphoto, length = ntot)

chill1 = ifelse(chill == 2, 1, 0) 
#chill2 = ifelse(chill == 3, 1, 0) 

treatcombo = paste(warm, photo, chill1, sep = "_")

(d <- data.frame(site, warm, photo, chill1, treatcombo)) # critical coding error here!

###### Set up differences for each level
sitediff = 2 
warmdiff = -20 # days earlier from 1 to 2
photodiff = -14
chill1diff = -20


# interactions. 9 two-way interactions
sitewarm = 0
sitephoto = 0
sitechill1 = -1 # similar to stan results
warmphoto = 3.5 # positive 3.5. So at the warm level, the effect of longer days is muted by 3.5 days.
warmchill1 = 11 # both positive ~ 10. 
photochill1 = 0.1 # from stan results

######## SD for each treatment
sitediff.sd = 1.5 
warmdiff.sd = 1 
photodiff.sd = 1
chill1diff.sd = 1.5

# interactions. 9 two-way interactions
sitewarm.sd = 1
sitephoto.sd = 1
sitechill1.sd = 2 
warmphoto.sd = 1
warmchill1.sd = 1.5
photochill1.sd = 1


mm <- model.matrix(~(site+warm+photo+chill1)^2, data.frame(site, warm, photo))
# remove last column, chill1 x chill2, empty
# mm <- mm[,-grep("chill1:chill2", colnames(mm))]
colnames(mm)

coeff <- c(1, sitediff, warmdiff, photodiff, chill1diff, 
           sitewarm, sitephoto, sitechill1, 
           warmphoto, warmchill1, 
           photochill1
)


bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 1) # should be able to do sd = mm %*% sd.coeff as well, with a different sd for each parameter.

(fake <- data_frame(bb, site, warm, photo, chill1))

summary(lm(bb ~ (site+warm+photo+chill1)^2, data = fake)) # sanity check 

##### Again, now with species now.

baseinter = 35 # baseline intercept across all species 
spint <- baseinter + c(1:nsp)-mean(1:nsp) # different intercepts by species

fake <- vector()
values <- vector()
head(mm)
for(i in 1:nsp){ # loop over species, as these are the random effect modeled
  
  # Give species different difference values, drawn from normal. Could make dataframe of diffs and diff.sds, and use apply..
  
  coeff <- c(spint[i], 
             rnorm(1, sitediff, sitediff.sd),
             rnorm(1, warmdiff, warmdiff.sd),
             rnorm(1, photodiff, photodiff.sd), 
             rnorm(1, chill1diff, chill1diff.sd),
             rnorm(1, sitewarm, sitewarm.sd), 
             rnorm(1, sitephoto, sitephoto.sd),
             rnorm(1, sitechill1, sitechill1.sd),
             rnorm(1, warmphoto, warmphoto.sd),
             rnorm(1, warmchill1, warmchill1.sd),
             rnorm(1, photochill1, photochill1.sd)
  )
  values <- rbind(values, coeff)
  
  bb <- rnorm(n = length(warm), mean = mm %*% coeff, sd = 0.1)
  
  fakex <- data.frame(bb, sp = i, site, warm, photo, chill1)
  
  fake <- rbind(fake, fakex)  
}
head(fake)

values <- as.data.frame(values)
names(values) <- c("int","b.site","b.force", "b.photo","b.chill1","b.fs","b.ps","b.c1s","b.fp","b.fc1", "b.pc1")

#save(list=c("fake"), file = "Fake_Budburst.RData")
write.csv(values, "df_values_2chill_sitefixed_noncp.csv", row.names = F)
# summary(lm(bb ~ (site+warm+photo+chi,"b.fc1"ll1+chill2)^2, data = fake)) # sanity check 

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

summary(lm(bb ~ (site+warm+photo+chill1)^2, data = fake))
# load("Stan Output 2016-04-04 Fake Interax.RData")
# fakeout <- dir()[grep("Stan Output", dir())[is.na(match(grep("Stan Output", dir()), grep("Fake", dir())))]]
# load(sort(realout, T)[1])
# ls() 

# <><><><><> Unpooled intercepts (Dan's original model) <><><><> #

# To Stan!
datalist.f <- list(lday = fake$bb, # budburst as respose 
                   warm = as.numeric(fake$warm), 
                   site = as.numeric(fake$site), 
                   sp = as.numeric(fake$sp), 
                   photo = as.numeric(fake$photo), 
                   chill1 = as.numeric(fake$chill1), 
                   N = nrow(fake), 
                   n_site = length(unique(fake$site)), 
                   n_sp = length(unique(fake$sp))
                  )


# To Stan! Pooled intercepts!
doym <- stan('stan/lday_site_sp_chill_inter_poola_NOncp_2chill_sitefixed.stan', data = datalist.f,
               iter = 4000,
               include = FALSE, pars = c("ypred_new","y_hat"))
#               # control = list(adapt_delta = 0.9,
#               #                max_treedepth = 15))
save(doym, file="output/df_2site_2chill_sitefixed_noncp.Rda")

# doym <- stan('stan/lday_site_sp_chill_inter_poola_ncp_2chill.stan', data = datalist.f,
#              iter = 4000,
#              include = FALSE, pars = c("ypred_new","y_hat"))
# #               # control = list(adapt_delta = 0.9,
# #               #                max_treedepth = 15))
# save(doym, file="output/df_2site_2chill.Rda")
# 
# values <- read.csv("df_values_2site.csv")
# load("output/df_2site_2chill.Rda")
# 

summary(doym)$summary[c( "mu_b_warm","mu_b_chill1", "mu_b_photo","b_site","mu_b_inter_wp","mu_b_inter_ps","mu_b_inter_wc1","mu_b_inter_pc1","mu_b_inter_sc1","sigma_b_warm","sigma_b_chill1","sigma_b_photo","sigma_b_inter_wp","sigma_b_inter_ps","sigma_b_inter_wc1","sigma_b_inter_pc1","sigma_b_inter_sc1", "sigma_y"),"mean"]

pdf("pairs_mu_2chill_sitefixed_noncp.pdf", width = 5, height =5)
pairs(doym, pars = c("mu_b_warm","mu_b_chill1", "mu_b_photo","b_site","mu_b_inter_wp","mu_b_inter_ps","mu_b_inter_wc1","mu_b_inter_pc1","mu_b_inter_sc1", "lp__")) 
dev.off()

pdf("pairs_sigma_2chill_sitefixed_noncp.pdf", width = 10, height =10)
pairs(doym, pars = c("sigma_b_warm","sigma_b_chill1","sigma_b_photo","sigma_b_inter_wp","sigma_b_inter_ps","sigma_b_inter_wc1","sigma_b_inter_pc1","sigma_b_inter_sc1", "sigma_y", "lp__")) 
dev.off()

sum <- summary(doym)$summary 
sum[c(526:554),c(1,5,7)]

intSp <- as.data.frame(sum[grep("a_sp", rownames(sum)), ])
slopeF <- as.data.frame(sum[grep("b_warm", rownames(sum)), ]); slopeF <- slopeF[c(1:28),]
slopeC <- as.data.frame(sum[grep("b_chill1", rownames(sum)), ]) ; slopeC <- slopeC[c(1:28),]
slopeP <- as.data.frame(sum[grep("b_photo", rownames(sum)), ]) ; slopeP <- slopeP[c(1:28),]
#slopeSite <- as.data.frame(sum[grep("b_site", rownames(sum)), ]) ; slopeSite <- slopeSite[c(1:28),]
slopeFP <- as.data.frame(sum[grep("b_inter_wp", rownames(sum)), ]) ; slopeFP <- slopeFP[c(1:28),]
slopeFC <- as.data.frame(sum[grep("b_inter_wc1", rownames(sum)), ]) ; slopeFC <- slopeFC[c(1:28),]
slopeCP <- as.data.frame(sum[grep("b_inter_pc1", rownames(sum)), ]) ; slopeCP <- slopeCP[c(1:28),]
slopeFS <- as.data.frame(sum[grep("b_inter_ws", rownames(sum)), ]) ; slopeFS <- slopeFS[c(1:28),]
slopePS <- as.data.frame(sum[grep("b_inter_ps", rownames(sum)), ]) ; slopePS <- slopePS[c(1:28),]
slopeCS <- as.data.frame(sum[grep("b_inter_sc1", rownames(sum)), ]) ; slopeCS <- slopeCS[c(1:28),]


# intSp <- as.data.frame(sum[grep("a_sp", rownames(sum)), ])
# slopeF <- as.data.frame(sum[grep("b_warm", rownames(sum)), ]); slopeF <- slopeF[c(1:28),]
# slopeC <- as.data.frame(sum[grep("b_chill1", rownames(sum)), ]) ; slopeC <- slopeC[c(1:28),]
# slopeP <- as.data.frame(sum[grep("b_photo", rownames(sum)), ]) ; slopeP <- slopeP[c(1:28),]
# #slopeSite <- as.data.frame(sum[grep("b_site", rownames(sum)), ]) ; slopeSite <- slopeSite[c(1:28),]
# slopeFP <- as.data.frame(sum[grep("b_inter_wp", rownames(sum)), ]) ; slopeFP <- slopeFP[c(31:58),]
# slopeFC <- as.data.frame(sum[grep("b_inter_wc1", rownames(sum)), ]) ; slopeFC <- slopeFC[c(31:58),]
# slopeCP <- as.data.frame(sum[grep("b_inter_pc1", rownames(sum)), ]) ; slopeCP <- slopeCP[c(31:58),]
# slopeFS <- as.data.frame(sum[grep("b_inter_ws", rownames(sum)), ]) ; slopeFS <- slopeFS[c(31:58),]
# slopePS <- as.data.frame(sum[grep("b_inter_ps", rownames(sum)), ]) ; slopePS <- slopePS[c(31:58),]
# slopeCS <- as.data.frame(sum[grep("b_inter_sc1", rownames(sum)), ]) ; slopeCS <- slopeCS[c(31:58),]

pdf("fake_esti_correlations_dflynn2chillsitefixedNOncp.pdf", height = 10, width = 5)
par(mfrow = c(6,2))
plot(values$int ~ intSp$mean, pch =19); abline(0,1)
plot( values$b.force ~ slopeF$mean, pch = 19); abline(0,1)
plot(values$b.chill1 ~ slopeC$mean, pch = 19); abline(0,1)
plot(values$b.photo ~ slopeP$mean, pch = 19);abline(0,1)
plot(values$b.fc1 ~ slopeFC$mean, pch = 19);abline(0,1)
plot(values$b.pc1 ~ slopeCP$mean, pch = 19);abline(0,1)
plot(values$b.fp ~ slopeFP$mean, pch = 19);abline(0,1)
plot(values$b.fs ~ slopeFS$mean, pch = 19);abline(0,1)
plot(values$b.c1s ~ slopeCS$mean, pch = 19);abline(0,1)
plot(values$b.ps ~ slopePS$mean, pch = 19);abline(0,1)
plot(values$b.site ~ slopeSite$mean, pch = 19);abline(0,1)
dev.off()

range(sum[,"Rhat" ])
# Below runs!
# doym.fpoola.ncp <- stan('stan/lday_site_sp_chill_inter_poola_ncp.stan', data = datalist.f, 
#                iter = 3000)
# save(doym.fpoola.ncp, file="output/lday_site_sp_chill_inter_poola_ncpfull.Rda")

# The below runs deadly slowly and never converges
# doym.fpoola.ncpfull  <- stan('stan/lday_site_sp_chill_inter_poola_ncpfull.stan', data = datalist.f, 
#               iter = 3000)
# Ditto for lday_site_sp_chill_inter_poola_ncpmore.stan, which I waited on for hours and gave up at 300 iter. 

# sf.poola <- summary(doym.fpoola)$summary
# sf.poola[grep("mu_", rownames(sf.poola)),]
# 
# # sanity checks
# summary(lm(bb ~ (site+warm+photo+chill1+chill2)^2, data = fake))
# # library(lme4)
# # summary(lmer(bb ~ (site+warm+photo+chill1+chill2)^2 + (1|sp), data = fake))
# 
# 
# 
# save(sf.poola, file="stan/lday_site_sp_chill_inter_poola.Rda")
# 
# save(doym.fpoola.ncp, file="stan/lday_site_sp_chill_inter_poola_ncpfull.Rda")
# # savestan("Fake Interax poola") # Not working! Saves a corrupted file.
# 
# # <> R stanarm <> #
# library(rstanarm)
# 
# arm.doym.f <- stan_lmer(lday ~ warm + photo + chill1 + chill2 + site +
#        warm*photo + warm*chill1 + warm*chill2 + warm*site + 
#        photo*chill1 + photo*chill2 + photo*site + site*chill1 + site*chill2
#        + (1|sp) +
#        (sp|warm) + (sp|photo) + (sp|chill1) + (sp|chill2) + (sp|site), 
#        data=datalist.f, 
#        prior=normal(), prior_intercept=normal(0,30), prior_aux=cauchy(0,10))
# 
# # <><><><><> End pooled intercepts <><><><> #
# 
# 
# ##
# library(rethinking)
# goober <- map2stan(
#      alist(
#          lday ~ dnorm(yhat, sigma_y),
#          yhat ~ a[sp] + b_site*sitis b_warm*warm + b_photo*photo + b_chill1*chill1 + b_chill2*chill2 + b_interwc*warm*photo,
#          a[sp] ~ dnorm(mu_a_sp, sigma_a_sp),
#          sigma_y ~ dnorm(0, 30),
#          mu_a_sp ~ dnorm(0, 50),
#          sigma_a_sp ~ dnorm(0, 30),
#          b_site~ dnorm(0, 30),
#          b_warm~ dnorm(0, 30),
#          b_photo~ dnorm(0, 30),
#          b_chill1~ dnorm(0, 30),
#          b_chill2~ dnorm(0, 30),
#          b_interwc ~ dnorm(0, 30)
#      ) ,
#      data=datalist.f, iter=2000)
# ##
# stop("Lizzie has not worked with the below code .... ")
# 
# #setwd("~/Documents/git/buds/analyses")
# 
# # Load lastest fake data output. Grep for both Fake and Stan Output.
# if(!exists('doym.f')){
#   
#   fakes <- na.omit(grep("Stan Output", dir())[match(grep("Fake", dir()), grep("Stan Output", dir()))])
# 
#   load(sort(dir()[fakes], T)[1])
# }
# 
# sf <- summary(doym.f)$summary
# 
# plotlet("b_warm", "b_photo", 
#         #xlab = "Advance due to 30d 4° chilling", 
#         #ylab = "Advance due to 30d 1.5° chilling", 
#         dat = sf)
# 
# plotlet("b_inter_wc2", "b_inter_wc1", 
#         #xlab = "Advance due to 30d 4° chilling", 
#         #ylab = "Advance due to 30d 1.5° chilling", 
#         dat = sf)
# 
# bees <- sf[grep("mu_b", rownames(sf)),]
# di <- sf[grep("mu_b_inter", rownames(sf)),]
# 
# plot(seq(min(di[,"mean"]-di[,"sd"]*1.5), max(di[,"mean"]+di[,"sd"]*1.5), length.out = nrow(di)),
#      1:nrow(di), type ="n")


# # now with fixed site and fixed sp
# 
# datalist.f <- list(lday = fake$bb, # budburst as respose 
#                    warm = as.numeric(fake$warm), 
#                    site = as.numeric(fake$site), 
#                    sp = as.numeric(fake$sp), 
#                    photo = as.numeric(fake$photo), 
#                    chill = as.numeric(fake$chill), 
#                    N = nrow(fake), 
#                    n_site = length(unique(fake$site)), 
#                    n_sp = length(unique(fake$sp))
# )
# 
# doym.f2 <- stan('stan/lday0_fixedsite_fixedsp.stan', data = datalist.f, iter = 4000, chains = 4) 
# 
# ssm.f <- as.shinystan(doym.f2)
# #launch_shinystan(ssm.f2) 
# 
# (sumer <- summary(doym.f)$summary)
# 
# setwd("~/Dropbox")
# 
# savestan()
# 
# setwd("~/Documents/git/buds/analyses")
# 
# # Tips for speeding up, from google groups
# set_cppo(mode = "fast")
# # For finding part of code that is slow
# dir(tempdir())

# posterior <- extract(doym.f)
# 
# sitediff = 2 
# warmdiff = -20 # days earlier from 1 to 2
# photodiff = -14
# chill1diff = -20
# chill2diff = -19
# 
# # interactions. 9 two-way interactions
# sitewarm = 0
# sitephoto = 0
# sitechill1 = -1 # similar to stan results
# sitechill2 = -2
# warmphoto = 3.5 # positive 3.5. So at the warm level, the effect of longer days is muted by 3.5 days.
# warmchill1 = 11 # both positive ~ 10. 
# warmchill2 = 9
# photochill1 = 0.1 # from stan results
# photochill2 = 1
# 
# ######## SD for each treatment
# sitediff.sd = 1.5 
# warmdiff.sd = 1 
# photodiff.sd = 1
# chill1diff.sd = 1.5
# chill2diff.sd = 2
# 
# # interactions. 9 two-way interactions
# sitewarm.sd = 1
# sitephoto.sd = 1
# sitechill1.sd = 2 
# sitechill2.sd = 2
# warmphoto.sd = 1
# warmchill1.sd = 1.5
# warmchill2.sd = 1.5
# photochill1.sd = 1
# photochill2.sd = 1
# 
# #pdf("figures/lnc_prior_post_dist.pdf", width = 15, height = 25)
# par(mfrow = c(4,4))
# hist(posterior$mu_b_warm,col=rgb(0,0,1,1/4))
# abline(v = warmdiff, col = "red", lty = 1, lwd = 5)
# 
# hist(posterior$mu_b_chill1,col=rgb(0,0,1,1/4))
# abline(v = chill1diff, col = "red", lty = 1, lwd = 5)
