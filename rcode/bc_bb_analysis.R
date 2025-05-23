## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

# modified by DL on Feb 8, 2022, summarizing the results from various ran models
# 1. Running my data only, site dummy, chill dummy, with interactions
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(scales)
library(arm)
library(rstan)
library(shinystan)
library(reshape2)
library(bayesplot)
library(ggplot2)
require(tidybayes)
library(dplyr)
library(plyr)

options(mc.cores = parallel::detectCores())

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  
#else{
#  setwd("~/deirdre/") # for midge
#}

#source('rcode/cleaning/pheno_bb_calc.R')
head(pheno)
length(unique(pheno$lab2))

df <- read.csv("input/day.of.bb.DFlynn.csv", header=TRUE, na.strings=c("","NA"))
head(df)

dl <- read.csv("input/day.of.bb.DL.csv", header=TRUE, na.strings=c("","NA"))
head(dl)
dl$lab3 <- dl$lab2
dl$lab2 <- paste(dl$species, dl$population, dl$rep, sep = "_")

# mergeing the my data with DF
#pheno <- rbind.fill(dl, df)
pheno <- dl
head(pheno)

# combined the data has 3197 unique samples
############################################################
# Preping the data for the model
#1. converting species to a factor
# colnames(pheno)[colnames(pheno) == "day"] <- "tbb"
# pheno <- pheno %>% separate(treatment, c("chill", "photo","force")); pheno <- as.data.frame(pheno)
#2. Adding columns of treatments as numeric values
pheno$chill.n <- pheno$chill
pheno$chill.n[pheno$chill.n == "HC"] <- "1"
pheno$chill.n[pheno$chill.n == "LC"] <- "0"
pheno$chill.n <- as.numeric(pheno$chill.n)

pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "0"
pheno$site.n[pheno$site.n == "mp"] <- "1"
# pheno$site.n[pheno$site.n == "HF"] <- "2"
# pheno$site.n[pheno$site.n == "SH"] <- "3"
pheno$site.n <- as.numeric(pheno$site.n)

head(pheno)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 

pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 30 species 

nrow(pheno.term) - nrow(pheno.t) # 547 that had no terminal bb

datalist <- with(pheno.t,
                    list( N=nrow(pheno.t),
                          n_sp = length(unique(pheno.t$species.fact)),
                          n_site = length(unique(pheno.t$site.n)),
                          bb = tbb,
                          sp = species.fact,
                          chill = chill.n,
                          photo = photo.n,
                          force = force.n,
                          site = site.n
                    ))

# mdl <- stan("stan/bc.bb.inter.stan",
#             data= datalist
#             ,iter=2000, chains=4)
#gives 200 divergent transitions, 41 transitions that exceed max tree depth, chains were not mixed, with low ESS

mdl.t <- stan("stan/bc.bb.ncpphoto.ncpinter.stan",
          data = datalist,
          iter = 4000, chains=4, control = list(adapt_delta = 0.99))


save(mdl.t, file="output/tbb_ncp_dl.Rda")
#load("output/tbb_ncp_dl.Rda")

#####################################################################
#####################################################################
load("output/tbb_ncp_dl.Rda")
sumt <- summary(mdl.t)$summary
sumt[grep("mu_", rownames(sumt)), ]
sumt
ssm <-  as.shinystan(mdl.t)
launch_shinystan(ssm)

## The model no longer has any divergent transitions for the terminal buds!
 pairs(mdl.t, pars=c("mu_force",
                     "mu_photo",
                     "mu_chill",
                     "mu_site",
                     "mu_inter_fp",
                     "mu_inter_fc",
                     "mu_inter_pc",
                     "mu_inter_fs",
                     "mu_inter_ps",
                     "mu_inter_sc")) # this gives a lot of warning messages and not the figure i was hoping/expected

range(sumt[, "n_eff"]) #1326.783 13829.827
range(sumt[, "Rhat"]) #0.9995324 1.0040581
# # now running the same model for the lateral buds
# pheno.50lat <- pheno[, c("latbb50", "chill.n", "force.n", "photo.n", "site.n", "species")]
# pheno.50l <- pheno.50lat[complete.cases(pheno.50lat), ]
# nrow(pheno.50lat) - nrow(pheno.50l)  # a lot of samples did not reach even 50%! 1084
# 
# pheno.50l$species.fact <- as.numeric(as.factor(pheno.50l$species))
# sort(unique(pheno.50l$species.fact))
# 
# datalist <- with(pheno.50l,
#                list( N = nrow(pheno.50l),
#                      n_sp = length(unique(pheno.50l$species.fact)),
#                      n_site = length(unique(pheno.50l$site.n)),
#                      bb = latbb50,
#                      sp = species.fact,
#                      chill = chill.n,
#                      photo = photo.n,
#                      force = force.n,
#                      site = site.n
#                ))
# 
# # mdl <- stan("stan/bc.bb.inter.stan",
# #             data= datalist
# #             ,iter=2000, chains=4)
# #gives 200 divergent transitions, 41 transitions that exceed max tree depth, chains were not mixed, with low ESS
# 
# mdl.50l <- stan("stan/bc.bb.ncpphoto.ncpinter.stan",
#           data= datalist,
#           iter=4000, chains=4, control = list(adapt_delta = 0.99))
# 
# sum50l <- summary(mdl.50l)$summary
# sum50l[grep("mu_",rownames(sum50l)), ]
# sum50l
# ssm <- as.shinystan(mdl.50l)
# launch_shinystan(ssm)
# 
# ## The model no longer has any divergent transitions for the terminal buds!
# #pairs(sm.sum, pars=c("mu_a","mu_force","mu_chill","mu_photo_ncp")) # this gives a lot of warning messages and not the figure i was hoping/expected
# 
# save(sum50l, file="output/tbb_photo_winter_ncp_lateralbud.Rds")

# now running the same model for the lateral buds
pheno.1lat <- pheno[, c("latbb1", "chill.n", "force.n", "photo.n", "site.n", "species")]
pheno.1l <- pheno.1lat[complete.cases(pheno.1lat), ]
nrow(pheno.1lat) - nrow(pheno.1l)  

pheno.1l$species.fact <- as.numeric(as.factor(pheno.1l$species))
sort(unique(pheno.1l$species.fact))

datalist <- with(pheno.1l,
                 list( N = nrow(pheno.1l),
                       n_sp = length(unique(pheno.1l$species.fact)),
                       n_site = length(unique(pheno.1l$site.n)),
                       bb = latbb1,
                       sp = species.fact,
                       chill = chill.n,
                       photo = photo.n,
                       force = force.n,
                       site = site.n
                 ))

# mdl <- stan("stan/bc.bb.inter.stan",
#             data= datalist
#             ,iter=2000, chains=4)
#gives 200 divergent transitions, 41 transitions that exceed max tree depth, chains were not mixed, with low ESS

mdl.1l <- stan("stan/bc.bb.ncpphoto.ncpinter.stan",
                data= datalist,
                iter=4000, chains=4, control = list(adapt_delta = 0.99))

sum1l <- summary(mdl.1l)$summary
sum1l[grep("mu_",rownames(sum1l)), ]
sum1l
ssm <- as.shinystan(mdl.1l)
launch_shinystan(ssm)

## The model no longer has any divergent transitions for the terminal buds!
#pairs(sm.sum, pars=c("mu_a","mu_force","mu_chill","mu_photo_ncp")) # this gives a lot of warning messages and not the figure i was hoping/expected

save(mdl.1l, file="output/lat_ncp_dl.Rds")
#####################################################################
# PPC 

mdl.slopes <- as.data.frame(sm.sum[grep("b", rownames(sm.sum)), c(1,6)]) 
mdl.int <- as.data.frame(sm.sum[grep("a", rownames(sm.sum)), ]) 
mdl.a <- mdl.int[, 1]
mdl.b <- mdl.slopes[, 1]

ggplot() +
  geom_point(data = mdl.slopes, aes(y = row.names(mdl.slopes), x = mean), color = "darkgreen") +
  labs( x = "doy", y = "Species")

#####################################################################

#plot(mdl, pars="a", ci_level=0.5, outer_level=0.5,col="blue")
# PPC based on the vingette from https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html 

ext <- rstan::extract(mdl.t)

# not the most normal 
hist(ext$a)

y <- pheno.t$tbb

y.ext <- ext$ypred_new # I want this to be a matrix, which it is, with one element for each data point in y

ppc_dens_overlay(y, y.ext[1:50, ])

mean(ext$b_force)
mean(ext$b_chill)
mean(ext$b_photo)

######################################################
# plotting code taken from buds-master Pheno Budburst analysis.R

col4fig <- c("mean","sd","25%","50%","75%","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")

# manually to get right order
mu_params <- c("mu_force",
               "mu_photo",
               "mu_chill",
               "mu_site",
               "mu_inter_fp",
               "mu_inter_fc",
               "mu_inter_pc",
               "mu_inter_fs",
               "mu_inter_ps",
               "mu_inter_sc")

meanzt <- sumt[mu_params, col4fig]
meanzl <- suml[mu_params, col4fig]

rownames(meanzt) = c("Forcing",
                     "Photoperiod",
                     "Chilling",
                     "Site",
                     "Forcing x Photoperiod",
                     "Forcing x Chilling",
                     "Photoperiod x Chilling",
                     "Forcing x Site",
                     "Photoperiod x Site",
                     "Site x Chilling"
  )

rownames(meanzl) = c("Forcing",
                     "Photoperiod",
                     "Chilling",
                     "Site",
                     "Forcing x Photoperiod",
                     "Forcing x Chilling",
                     "Photoperiod x Chilling",
                     "Forcing x Site",
                     "Photoperiod x Site",
                     "Site x Chilling"
  )

meanzt.table <- sumt[mu_params, col4table]
row.names(meanzt.table) <- row.names(meanzt)
head(meanzt.table)
write.table(meanzt.table , "output/termMdlEstiDL.csv", sep = ",", row.names = FALSE)

meanzl.table <- suml[mu_params, col4table]
row.names(meanzl.table) <- row.names(meanzl)
head(meanzl.table)
write.table(meanzl.table , "output/lat.mdl.esti.csv", sep = ",", row.names = FALSE)

# Begin by checking to see what cue is most important and whether there are strong correlations between cues:
df.mean.t <- data.frame(bb.force = sumt[grep("b_force", rownames(sumt)), 1],
                          bb.photo = sumt[grep("b_photo_ncp", rownames(sumt)), 1],
                          bb.chill = sumt[grep("b_chill", rownames(sumt)), 1])

df.mean.l <- data.frame(lat.force = suml[grep("b_force", rownames(sumt)), 1],
                        lat.photo = suml[grep("b_photo_ncp", rownames(sumt)), 1],
                        lat.chill = suml[grep("b_chill", rownames(sumt)), 1])

df.mean.t[which(df.mean.t$bb.force > df.mean.t$bb.photo), ] # species 11- rho alb
df.mean.l[which(df.mean.l$lat.force > df.mean.l$lat.photo), ] #none
df.mean.t[which(df.mean.t$bb.chill > df.mean.t$bb.force), ] # 14
# 3, 5,6,8,9,10,12,13,15,17,18,20
df.mean.l[which(df.mean.l$lat.chill > df.mean.l$lat.force), ] # 16
#1,2,5,6,7,8,10,11,12,13,15,16,17,18,19,20

# all correlated
summary(lm(bb.force~bb.photo, data=df.mean.t))
summary(lm(bb.force~bb.chill, data=df.mean.t))
summary(lm(bb.force~bb.photo, data=df.mean.t))
summary(lm(lat.force~lat.photo, data=df.mean.l))
summary(lm(lat.force~lat.chill, data=df.mean.l))
summary(lm(lat.chill~lat.photo, data=df.mean.l))

pdf(file.path( "figures/changes.pheno.pdf"), width = 7, height = 8)
par(mfrow = c(2,1), mar = c(5, 10, 2, 1))
# Upper panel: bud burst
plot(seq(-22, 
         12,
         length.out = nrow(meanzt)), 
     1:nrow(meanzt),
     type = "n",
     xlab = "",
     ylab = "",
     yaxt = "n")

#legend(x = -20, y = 2, bty="n", legend = "a. Budburst", text.font = 2)
#rasterImage(bbpng, -20, 1, -16, 4)

axis(2, at = nrow(meanzt):1, labels = rownames(meanzt), las = 1, cex.axis = 0.8)
points(meanzt[, 'mean'],
       nrow(meanzt):1,
       pch = 16,
       col = "midnightblue")
arrows(meanzt[, "75%"], nrow(meanzt):1, meanzt[, "25%"], nrow(meanzt):1,
       len = 0, col = "black")
abline(v = 0, lty = 3)
# add advance/delay arrows
par(xpd=NA)
arrows(1, 15.5, 6, 15.5, len = 0.1, col = "black")
legend(5, 16.5, legend = "delay", bty = "n", text.font = 1, cex = 0.75)
arrows(-1, 15.5, -6, 15.5, len = 0.1, col = "black")
legend(-12, 16.5, legend = "advance", bty = "n", text.font = 1, cex = 0.75)
legend(-2, 16.5, legend = "0", bty = "n", text.font = 1, cex = 0.75)
par(xpd = FALSE)

par(mar = c(5, 10, 2, 1))
# Lower panel: leaf-out
plot(seq(-22, 
         12, 
         length.out = nrow(meanzl)), 
     1:nrow(meanzl),
     type = "n",
     xlab = "Model estimate change in day of phenological event",
     ylab = "",
     yaxt = "n")

#legend(x = -24, y = 6, bty="n", legend = "b. Leafout", text.font = 2)
#rasterImage(lopng, -20, 1, -14, 4)

axis(2, at = nrow(meanzl):1, labels = rownames(meanzl), las = 1, cex.axis = 0.8)
points(meanzl[,'mean'],
       nrow(meanzl):1,
       pch = 16,
       col = "midnightblue")
arrows(meanzl[,"75%"], nrow(meanzl):1, meanzl[,"25%"], nrow(meanzl):1,
       len = 0, col = "black")
abline(v = 0, lty = 3)

# add advance/delay arrows
# par(xpd=NA)
# arrows(1, 15.5, 6, 15.5, len=0.1, col = "black")
# legend(5, 16.5, legend="delay", bty="n", text.font = 1, cex=0.75)
# arrows(-1, 15.5, -6, 15.5, len=0.1, col = "black")
# legend(-12, 16.5, legend="advance", bty="n", text.font = 1, cex=0.75)
# legend(-2, 16.5, legend="0", bty="n", text.font = 1, cex=0.75)
# par(xpd=FALSE)
dev.off()

# Comparisons of trees vs shrubs:
shrubs = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG", "alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
trees = c("ACEPEN", "ACERUB", "ACESAC", "ALNINC", "BETALL", "BETLEN", "BETPAP", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL", "acegla","betpap", "poptre", "popbal")

treeshrub = levels(dx$sp)
treeshrub[treeshrub %in% shrubs] = 1
treeshrub[treeshrub %in% trees] = 2
treeshrub = as.numeric(treeshrub)
