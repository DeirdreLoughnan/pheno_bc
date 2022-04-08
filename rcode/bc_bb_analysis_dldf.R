## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

#modified Feb 8, 2022 
# summarizing what models have been run and why they were not good

rm(list=ls()) 
options(stringsAsFactors = FALSE)
#library(scales)
#library(arm)
library(rstan)
library(shinystan)
#library(reshape2)
#library(bayesplot)
#library(ggplot2)
#library(RColorBrewer)
library(dplyr)
library(plyr)
library(stringr)

options(mc.cores = parallel::detectCores())


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  
#else{
#  setwd("~/deirdre/") # for midge
#}

#source('rcode/cleaning/pheno_bb_calc.R')
# head(pheno)
# length(unique(pheno$lab2))
# df <- read.csv("input/lday.csv", header=TRUE, na.strings=c("","NA"))
# 
# head(df)
# df.chill <- read.csv("input/chilling_values_eastern.csv")
# df$population <- df$site
# df$population[df$population == "0"] <- "HF"
# df$population[df$population == "1"] <- "SH"
# 
# df.wchill <- merge(df, df.chill, by =c("population","chill"))

dl <- read.csv("input/dl_allbb.csv")

temp <- str_split_fixed(dl$trt, "_", 3); head(temp)
dl$chill<- temp[,1]
dl$photo <- temp[,2]
dl$force <- temp[,3]

dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")

df <- read.csv("input/df_dxb_prepped_data.csv")
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill <- df.wchill[, c("population", "chill","force","photo","lab2", "bb","species", "treatment","Chill_portions","Utah_Model")]
 
# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)

head(pheno)
#table(pheno$species, pheno$first)
#write.csv(pheno, "input/pheno.wchill.midge.csv")
# 
# combined the data has 3197 unique samples
############################################################
# Preping the data for the model
#1. converting species to a factor
# colnames(pheno)[colnames(pheno) == "day"] <- "tbb"
# pheno <- pheno %>% separate(treatment, c("chill", "photo","force")); pheno <- as.data.frame(pheno)
#2. Adding columns of treatments as numeric values
# pheno$chill.n <- pheno$chill
# pheno$chill.n[pheno$chill.n == "HC"] <- "1"
# pheno$chill.n[pheno$chill.n == "LC"] <- "0"
# pheno$chill.n <- as.numeric(pheno$chill.n)
# Trying to run the model with chill portions 
#pheno$Chill_portions <- as.factor(pheno$Chill_portions)

pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.numeric(pheno$site.n)

head(pheno)
#add dummy/ site level effects:
pheno <- pheno %>%
  mutate ( site2 = if_else(site.n == 2, 1, 0),
           site3 = if_else(site.n == 3, 1, 0),
           site4 = if_else(site.n == 4, 1, 0))

# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
#pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term), ] # 3609

pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 49

nrow(pheno.term) - nrow(pheno.t) # 609


datalist.z2 <- with(pheno.t,
                 list( N=nrow(pheno.t),
                       n_sp = length(unique(pheno.t$species.fact)),
                       n_site = length(unique(pheno.t$population)),
                       lday = bb,
                       sp = species.fact,
                       chill1 = chillport.z2,
                       photo = photo.z2,
                       warm = force.z2,
                       site2 = site2.z2,
                       site3 = site3.z2,
                       site4 = site4.z2
                       ))


mdl.t <- stan("stan/df_mdl_4sites_again_allint_ncp.stan",
              data = datalist.z2,
              iter = 4000, warmup =3000, chains=4
              #, control = list(adapt_delta = 0.99)
              )
save(mdl.t, file="output/tbb_4sites_fullystandardized_ncp.Rda")


## The model no longer has any divergent transitions for the terminal buds!
#pairs(sm.sum, pars=c("mu_a","mu_force","mu_chill","mu_photo_ncp")) # this gives a lot of warning messages and not the figure i was hoping/expected

# range(sumt[, "n_eff"])
# range(sumt[, "Rhat"])

#######################################################################

 load("output/final/dl_df_allbb_4sites.Rda")
# 
 sumt <- summary(mdl.t)$summary
# mu <- sumt[grep("mu_", rownames(sumt)), ]

# ssm <-  as.shinystan(mdl.t)
# launch_shinystan(ssm)

range(sumt[, "n_eff"])
range(sumt[, "Rhat"])

summary(mdl.t)$summary[c("mu_a",
                       "mu_b_warm",
                       "mu_b_photo",
                       "mu_b_chill1",
                       "b_site2",
                       "b_site3",
                       "b_site4",
                       "mu_b_inter_wp",
                       "mu_b_inter_wc1",
                       "mu_b_inter_pc1",
                       "mu_b_inter_ws2",
                       "mu_b_inter_ps2",
                       "mu_b_inter_s2c1",
                       "mu_b_inter_ws3",
                       "mu_b_inter_ps3",
                       "mu_b_inter_s3c1",
                       "mu_b_inter_ws4",
                       "mu_b_inter_ps4",
                       "mu_b_inter_s4c1",
                       "sigma_b_warm",
                       "sigma_b_photo",
                       "sigma_b_chill1",
                       "sigma_a",
                       "sigma_b_inter_wp",
                       "sigma_b_inter_wc1",
                       "sigma_b_inter_pc1",
                       "sigma_b_inter_ws2",
                       "sigma_b_inter_ps2",
                       "sigma_b_inter_s2c1",
                       "sigma_b_inter_ws3",
                       "sigma_b_inter_ps3",
                       "sigma_b_inter_s3c1",
                       "sigma_b_inter_ws4",
                       "sigma_b_inter_ps4",
                       "sigma_b_inter_s4c1",
                       "sigma_y"),"mean"]


post <- rstan::extract(mdl.t)

y<-as.numeric(pheno.t$bb)
yrep<-post$y_hat # I want this to be a matrix, which it is, with one element for each data point in y

pdf("yrepvsypred_4site_stand.pdf")
ppc_dens_overlay(y, yrep[1:50, ])
dev.off()

#####################################################################
a_sp = sumt[grep("mu_a", rownames(sumt)), 1]
mu_b_warm = sumt[grep("mu_b_warm", rownames(sumt)), 1]
mu_b_photo = sumt[grep("mu_b_photo", rownames(sumt)), 1]
mu_b_chill1 = sumt[grep("mu_b_chill1", rownames(sumt)), 1]
mu_b_inter_pc1 = sumt[grep("mu_b_inter_pc1", rownames(sumt)), 1]
mu_b_inter_wp = sumt[grep("mu_b_inter_wp", rownames(sumt)), 1]
mu_b_inter_wc1 = sumt[grep("mu_b_inter_wc1", rownames(sumt)), 1]
mu_b_inter_ws2 = sumt[grep("mu_b_inter_ws2", rownames(sumt)), 1]
mu_b_inter_s2c1 = sumt[grep("mu_b_inter_s2c1", rownames(sumt)), 1]
mu_b_inter_ps2 = sumt[grep("mu_b_inter_ps2", rownames(sumt)), 1]
mu_b_inter_ws3 = sumt[grep("mu_b_inter_ws3", rownames(sumt)), 1]
mu_b_inter_s3c1 = sumt[grep("mu_b_inter_s3c1", rownames(sumt)), 1]
mu_b_inter_ps3 = sumt[grep("mu_b_inter_ps3", rownames(sumt)), 1]
mu_b_inter_ws4 = sumt[grep("mu_b_inter_ws4", rownames(sumt)), 1]
mu_b_inter_s4c1 = sumt[grep("mu_b_inter_s4c1", rownames(sumt)), 1]
mu_b_inter_ps4 = sumt[grep("mu_b_inter_ps4", rownames(sumt)), 1]
b_site2 = sumt[grep("b_site2", rownames(sumt)), 1]
b_site3 = sumt[grep("b_site3", rownames(sumt)), 1]
b_site4 = sumt[grep("b_site4", rownames(sumt)), 1]

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# plot the interactions
# Warm X chill

hfData <- subset(pheno, force == "HF" )
lfData <- subset(pheno, force == "LF")

# Make the other parameters constant
hf <- unique(hfData$force.z2)
lf <- unique(lfData$force.z2)
photo <- -0.5041133
siteSM <- 0
chill1 <- c( -0.7642814, -0.4072595, -0.4023109, -0.3493703,  0.2750890,  0.2977055,  0.4308763,  0.5308110,  0.8457874,  0.9457221)


# plot first for the high forcing
bb_hfc = a_sp + b_site2 * siteSM + b_site3 * siteSM + b_site4 * siteSM + mu_b_warm * hf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (hf*photo) +
  mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*siteSM) + mu_b_inter_ws2 * (hf*siteSM) +mu_b_inter_ps2 * (photo*siteSM) +
  mu_b_inter_s3c1 * (chill1*siteSM) + mu_b_inter_ws3 * (hf*siteSM) +mu_b_inter_ps3 * (photo*siteSM) +
  mu_b_inter_s4c1 * (chill1*siteSM) + mu_b_inter_ws4 * (hf*siteSM) +mu_b_inter_ps4 * (photo*siteSM) 

# plot first for the low forcing
bb_lfc = a_sp + + b_site2 * siteSM + b_site3 * siteSM + b_site4 * siteSM + mu_b_warm * lf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (lf*photo) +
  mu_b_inter_wc1 * (lf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*siteSM) + mu_b_inter_ws2 * (lf*siteSM) +mu_b_inter_ps2 * (photo*siteSM) +
  mu_b_inter_s3c1 * (chill1*siteSM) + mu_b_inter_ws3 * (lf*siteSM) +mu_b_inter_ps3 * (photo*siteSM) +
  mu_b_inter_s4c1 * (chill1*siteSM) + mu_b_inter_ws4 * (lf*siteSM) +mu_b_inter_ps4 * (photo*siteSM) 
# 
pdf("figures/dldf_4sites_interactions.pdf", width =10, height = 3.5)
par(mfrow =c (1,3))
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Chill portions", ylab = "Day of budburst")
points(hfData$chillport.z2, hfData$bb, col = "maroon")
points(lfData$chillport.z2, lfData$bb, col = "darkslategray4")
abline(lm(bb_hfc ~ chill1), col = "darkred", lwd = 3)
abline(lm(bb_lfc ~ chill1), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("low forcing"),
                            expression("high forcing")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# warm and site2 
hf <- unique(hfData$force.z2)
lf <- unique(lfData$force.z2)
photo <- -0.5044652 
site2 <- unique(pheno$site2.z2)

chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)

# plot first for the high forcing
bb_hfsite2 = a_sp + b_site2 * site2 + b_site3 * site2 + b_site4 * site2 + mu_b_warm * hf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (hf*photo) +
  mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site2) + mu_b_inter_ws2 * (hf*site2) +mu_b_inter_ps2 * (photo*site2) +
  mu_b_inter_s3c1 * (chill1*site2) + mu_b_inter_ws3 * (hf*site2) +mu_b_inter_ps3 * (photo*site2) +
  mu_b_inter_s4c1 * (chill1*site2) + mu_b_inter_ws4 * (hf*site2) +mu_b_inter_ps4 * (photo*site2) 

# plot first for the low forcing
bb_lfsite2 = a_sp + b_site2 * site2 + b_site3 * site2 + b_site4 * site2  + mu_b_warm * lf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (lf*photo) +
  mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site2) + mu_b_inter_ws2 * (lf*site2) +mu_b_inter_ps2 * (photo*site2) +
  mu_b_inter_s3c1 * (chill1*site2) + mu_b_inter_ws3 * (lf*site2) +mu_b_inter_ps3 * (photo*site2) +
  mu_b_inter_s4c1 * (chill1*site2) + mu_b_inter_ws4 * (lf*site2) +mu_b_inter_ps4 * (photo*site2) 
# 
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Site2", ylab = "Day of budburst")
points(hfData$site2.z2, hfData$bb, col = "maroon")
points(lfData$site2.z2, lfData$bb, col = "darkslategray4")
abline(lm(bb_hfsite2 ~ site2), col = "darkred", lwd = 3)
abline(lm(bb_lfsite2 ~ site2), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("low forcing"),
                            expression("high forcing")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# warm and site3 
hf <- unique(hfData$force.z2)
lf <- unique(lfData$force.z2)
photo <- -0.5044652 
site3 <- unique(pheno$site3.z2)

chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)

# plot first for the high forcing
bb_hfsite3 = a_sp + b_site2 * site3 + b_site3 * site3 + b_site4 * site3 + mu_b_warm * hf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (hf*photo) +
  mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site3) + mu_b_inter_ws2 * (hf*site3) +mu_b_inter_ps2 * (photo*site3) +
  mu_b_inter_s3c1 * (chill1*site3) + mu_b_inter_ws3 * (hf*site3) +mu_b_inter_ps3 * (photo*site3) +
  mu_b_inter_s4c1 * (chill1*site3) + mu_b_inter_ws4 * (hf*site3) +mu_b_inter_ps4 * (photo*site3) 

# plot first for the low forcing
bb_lfsite3 = a_sp + + b_site2 * site3 + b_site3 * site3 + b_site4 * site3 + mu_b_warm * lf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (lf*photo) +
  mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site3) + mu_b_inter_ws2 * (lf*site3) +mu_b_inter_ps2 * (photo*site3) +
  mu_b_inter_s3c1 * (chill1*site3) + mu_b_inter_ws3 * (lf*site3) +mu_b_inter_ps3 * (photo*site3) +
  mu_b_inter_s4c1 * (chill1*site3) + mu_b_inter_ws4 * (lf*site3) +mu_b_inter_ps4 * (photo*site3) 
# 
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "site3", ylab = "Day of budburst")
points(hfData$site3.z2, hfData$bb, col = "maroon")
points(lfData$site3.z2, lfData$bb, col = "darkslategray4")
abline(lm(bb_hfsite3 ~ site3), col = "darkred", lwd = 3)
abline(lm(bb_lfsite3 ~ site3), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("low forcing"),
                            expression("high forcing")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# warm and site4 
hf <- unique(hfData$force.z2)
lf <- unique(lfData$force.z2)
photo <- -0.5044652 
site4 <- unique(pheno$site4.z2)

chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)

# plot first for the high forcing
bb_hfsite4 = a_sp + b_site2 * site4 + b_site3 * site4 + b_site4 * site4  + mu_b_warm * hf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (hf*photo) +
  mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site4) + mu_b_inter_ws2 * (hf*site4) +mu_b_inter_ps2 * (photo*site4) +
  mu_b_inter_s3c1 * (chill1*site4) + mu_b_inter_ws3 * (hf*site4) +mu_b_inter_ps3 * (photo*site4) +
  mu_b_inter_s4c1 * (chill1*site4) + mu_b_inter_ws4 * (hf*site4) +mu_b_inter_ps4 * (photo*site4) 

# plot first for the low forcing
bb_lfsite4 = a_sp + b_site2 * site4 + b_site3 * site4 + b_site4 * site4 + mu_b_warm * lf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (lf*photo) +
  mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site4) + mu_b_inter_ws2 * (lf*site4) +mu_b_inter_ps2 * (photo*site4) +
  mu_b_inter_s3c1 * (chill1*site4) + mu_b_inter_ws3 * (lf*site4) +mu_b_inter_ps3 * (photo*site4) +
  mu_b_inter_s4c1 * (chill1*site4) + mu_b_inter_ws4 * (lf*site4) +mu_b_inter_ps4 * (photo*site4) 
# 
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "site4", ylab = "Day of budburst")
points(hfData$site4.z2, hfData$bb, col = "maroon")
points(lfData$site4.z2, lfData$bb, col = "darkslategray4")
abline(lm(bb_hfsite4 ~ site4), col = "darkred", lwd = 3)
abline(lm(bb_lfsite4 ~ site4), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("low forcing"),
                            expression("high forcing")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# photo and site2 
hpData <- subset(pheno, photo == "HP" )
lpData <- subset(pheno, photo == "LP")

hp <- unique(hpData$photo.z2)
lp <- unique(lpData$photo.z2)
force <- -0.5077181
site2 <- unique(pheno$site2.z2)

chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)

# plot first for the high forcing
bb_hpsite2 = a_sp + b_site2 * site2 + b_site3 * site2 + b_site4 * site2+ mu_b_warm * force + mu_b_photo * hp + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*hp) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (hp*chill1) +
  mu_b_inter_s2c1 * (chill1*site2) + mu_b_inter_ws2 * (force*site2) +mu_b_inter_ps2 * (hp*site2) +
  mu_b_inter_s3c1 * (chill1*site2) + mu_b_inter_ws3 * (force*site2) +mu_b_inter_ps3 * (hp*site2) +
  mu_b_inter_s4c1 * (chill1*site2) + mu_b_inter_ws4 * (force*site2) +mu_b_inter_ps4 * (hp*site2) 

# plot first for the low forcing
bb_lpsite2 = a_sp + b_site2 * site2 + b_site3 * site2+ b_site4 * site2+ mu_b_warm * force + mu_b_photo * lp + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*lp) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (lp*chill1) +
  mu_b_inter_s2c1 * (chill1*site2) + mu_b_inter_ws2 * (force*site2) +mu_b_inter_ps2 * (lp*site2) +
  mu_b_inter_s3c1 * (chill1*site2) + mu_b_inter_ws3 * (force*site2) +mu_b_inter_ps3 * (lp*site2) +
  mu_b_inter_s4c1 * (chill1*site2) + mu_b_inter_ws4 * (force*site2) +mu_b_inter_ps4 * (lp*site2) 
# 
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Site2", ylab = "Day of budburst")
points(lpData$site2.z2, lpData$bb, col = "darkslategray4")
points(hpData$site2.z2, hpData$bb, col = "maroon")
abline(lm(bb_hpsite2 ~ site2), col = "darkred", lwd = 3)
abline(lm(bb_lpsite2 ~ site2), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("low photoperiod"),
                            expression("high photoperiod")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#\
hpData <- subset(pheno, photo == "HP" )
lpData <- subset(pheno, photo == "LP")

hp <- unique(hpData$photo.z2)
lp <- unique(lpData$photo.z2)
force <- -0.5077181
site3 <- unique(pheno$site3.z2)

chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)

# plot first for the high forcing
bb_hpsite3 = a_sp + b_site2 * site3 + b_site3 * site3 + b_site4 * site3 + mu_b_warm * force + mu_b_photo * hp + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*hp) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (hp*chill1) +
  mu_b_inter_s2c1 * (chill1*site3) + mu_b_inter_ws2 * (force*site3) +mu_b_inter_ps2 * (hp*site3) +
  mu_b_inter_s3c1 * (chill1*site3) + mu_b_inter_ws3 * (force*site3) +mu_b_inter_ps3 * (hp*site3) +
  mu_b_inter_s4c1 * (chill1*site3) + mu_b_inter_ws4 * (force*site3) +mu_b_inter_ps4 * (hp*site3) 

# plot first for the low forcing
bb_lpsite3 = a_sp + b_site2 * site3 + b_site3 * site3 + b_site4 * site3 + mu_b_warm * force + mu_b_photo * lp + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*lp) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (lp*chill1) +
  mu_b_inter_s2c1 * (chill1*site3) + mu_b_inter_ws2 * (force*site3) +mu_b_inter_ps2 * (lp*site3) +
  mu_b_inter_s3c1 * (chill1*site3) + mu_b_inter_ws3 * (force*site3) +mu_b_inter_ps3 * (lp*site3) +
  mu_b_inter_s4c1 * (chill1*site3) + mu_b_inter_ws4 * (force*site3) +mu_b_inter_ps4 * (lp*site3) 

plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "site3", ylab = "Day of budburst")
points(hpData$site3.z2, hpData$bb, col = "maroon")
points(lpData$site3.z2, lpData$bb, col = "darkslategray4")
abline(lm(bb_hpsite3 ~ site3), col = "darkred", lwd = 3)
abline(lm(bb_lpsite3 ~ site3), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("low photoperiod"),
                            expression("high photoperiod")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# warm and site4 
hp <- unique(hpData$photo.z2)
lp <- unique(lpData$photo.z2)
force <- -0.5077181
site4 <- unique(pheno$site4.z2)

chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)

# plot first for the high forcing
bb_hpsite4 = a_sp + b_site2 * site4 + b_site3 * site4 + b_site4 * site4 + mu_b_warm * force + mu_b_photo * hp + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*hp) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (hp*chill1) +
  mu_b_inter_s2c1 * (chill1*site4) + mu_b_inter_ws2 * (force*site4) +mu_b_inter_ps2 * (hp*site4) +
  mu_b_inter_s3c1 * (chill1*site4) + mu_b_inter_ws3 * (force*site4) +mu_b_inter_ps3 * (hp*site4) +
  mu_b_inter_s4c1 * (chill1*site4) + mu_b_inter_ws4 * (force*site4) +mu_b_inter_ps4 * (hp*site4) 

# plot first for the low forcing
bb_lpsite4 = a_sp + + b_site2 * site4 + b_site3 * site4 + b_site4 * site4 + mu_b_warm * force + mu_b_photo * lp + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*lp) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (lp*chill1) +
  mu_b_inter_s2c1 * (chill1*site4) + mu_b_inter_ws2 * (force*site4) +mu_b_inter_ps2 * (lp*site4) +
  mu_b_inter_s3c1 * (chill1*site4) + mu_b_inter_ws3 * (force*site4) +mu_b_inter_ps3 * (lp*site4) +
  mu_b_inter_s4c1 * (chill1*site4) + mu_b_inter_ws4 * (force*site4) +mu_b_inter_ps4 * (lp*site4) 

plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "site4", ylab = "Day of budburst")
points(hpData$site4.z2, hpData$bb, col = "maroon")
points(lpData$site4.z2, lpData$bb, col = "darkslategray4")
abline(lm(bb_hpsite4 ~ site4), col = "darkred", lwd = 3)
abline(lm(bb_lpsite4 ~ site4), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("low photoperiod"),
                            expression("high photoperiod")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# plot the site chilling interactions 

site2_0Data <- subset(pheno, site2 == "0" )
site2_1Data <- subset(pheno, site2 == "1")

site2_0 <- unique(site2_0Data$site2.z2)
site2_1 <- unique(site2_1Data$site2.z2)

# Make the other parameters constant
photo <- -0.5041133
force <- -0.5080665 # zero forcing
chill1 <- c( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)


# plot first site2_0
bb_csite2_0 = a_sp + b_site2 * site2_0 + b_site3 * site2_0 + b_site4 * site2_0 + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*photo) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
    mu_b_inter_s2c1 * (chill1*site2_0) + mu_b_inter_ws2 * (force*site2_0) +mu_b_inter_ps2 * (photo*site2_0) +
    mu_b_inter_s3c1 * (chill1*site2_0) + mu_b_inter_ws3 * (force*site2_0) +mu_b_inter_ps3 * (photo*site2_0) +
    mu_b_inter_s4c1 * (chill1*site2_0) + mu_b_inter_ws4 * (force*site2_0) +mu_b_inter_ps4 * (photo*site2_0) 
  

# plot first for the low forcing
bb_csite2_1 = a_sp + b_site2 * site2_1 + b_site3 * site2_1 + b_site4 * site2_1 + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*photo) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site2_1) + mu_b_inter_ws2 * (force*site2_1) +mu_b_inter_ps2 * (photo*site2_1) +
  mu_b_inter_s3c1 * (chill1*site2_1) + mu_b_inter_ws3 * (force*site2_1) +mu_b_inter_ps3 * (photo*site2_1) +
  mu_b_inter_s4c1 * (chill1*site2_1) + mu_b_inter_ws4 * (force*site2_1) +mu_b_inter_ps4 * (photo*site2_1) 
# 
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Chill portions", ylab = "Day of budburst")
points(site2_0Data$chillport.z2, site2_0Data$bb, col = "maroon")
points(site2_1Data$chillport.z2, site2_1Data$bb, col = "darkslategray4")
abline(lm(bb_csite2_0  ~ chill1), col = "darkred", lwd = 3)
abline(lm(bb_csite2_1 ~ chill1), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("site 2_1"),
                            expression("site 2_0")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# plot the site chilling interactions 

site3_0Data <- subset(pheno, site3 == "0" )
site3_1Data <- subset(pheno, site3 == "1")

site3_0 <- unique(site3_0Data$site3.z2)
site3_1 <- unique(site3_1Data$site3.z2)

# Make the other parameters constant
photo <- -0.5041133
force <- -0.5080665 # zero forcing
chill1 <- c( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)


# plot first site3_0
bb_csite3_0 = a_sp + b_site2 * site3_0 + b_site3 * site3_0 + b_site4 * site3_0 + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*photo) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site3_0) + mu_b_inter_ws2 * (force*site3_0) +mu_b_inter_ps2 * (photo*site3_0) +
  mu_b_inter_s3c1 * (chill1*site3_0) + mu_b_inter_ws3 * (force*site3_0) +mu_b_inter_ps3 * (photo*site3_0) +
  mu_b_inter_s4c1 * (chill1*site3_0) + mu_b_inter_ws4 * (force*site3_0) +mu_b_inter_ps4 * (photo*site3_0) 


# plot first for the low forcing
bb_csite3_1 = a_sp + b_site2 * site3_1 + b_site3 * site3_1 + b_site4 * site3_1 + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*photo) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site3_1) + mu_b_inter_ws2 * (force*site3_1) +mu_b_inter_ps2 * (photo*site3_1) +
  mu_b_inter_s3c1 * (chill1*site3_1) + mu_b_inter_ws3 * (force*site3_1) +mu_b_inter_ps3 * (photo*site3_1) +
  mu_b_inter_s4c1 * (chill1*site3_1) + mu_b_inter_ws4 * (force*site3_1) +mu_b_inter_ps4 * (photo*site3_1) 
# 
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Chill portions", ylab = "Day of budburst")
points(site3_0Data$chillport.z2, site3_0Data$bb, col = "maroon")
points(site3_1Data$chillport.z2, site3_1Data$bb, col = "darkslategray4")
abline(lm(bb_csite3_0  ~ chill1), col = "darkred", lwd = 3)
abline(lm(bb_csite3_1 ~ chill1), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("site 2_1"),
                            expression("site 2_0")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# plot the site chilling interactions 

site4_0Data <- subset(pheno, site4 == "0" )
site4_1Data <- subset(pheno, site4 == "1")

site4_0 <- unique(site4_0Data$site4.z2)
site4_1 <- unique(site4_1Data$site4.z2)

# Make the other parameters constant
photo <- -0.5041133
force <- -0.5080665 # zero forcing
chill1 <- c( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)


# plot first site4_0
bb_csite4_0 = a_sp + b_site2 * site4_0 + b_site3 * site4_0 + b_site4 * site4_0 + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*photo) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site4_0) + mu_b_inter_ws2 * (force*site4_0) +mu_b_inter_ps2 * (photo*site4_0) +
  mu_b_inter_s3c1 * (chill1*site4_0) + mu_b_inter_ws3 * (force*site4_0) +mu_b_inter_ps3 * (photo*site4_0) +
  mu_b_inter_s4c1 * (chill1*site4_0) + mu_b_inter_ws4 * (force*site4_0) +mu_b_inter_ps4 * (photo*site4_0) 


# plot first for the low forcing
bb_csite4_1 = a_sp + b_site2 * site4_1 + b_site3 * site4_1 + b_site4 * site4_1 + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
  mu_b_inter_wp * (force*photo) +
  mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_s2c1 * (chill1*site4_1) + mu_b_inter_ws2 * (force*site4_1) +mu_b_inter_ps2 * (photo*site4_1) +
  mu_b_inter_s3c1 * (chill1*site4_1) + mu_b_inter_ws3 * (force*site4_1) +mu_b_inter_ps3 * (photo*site4_1) +
  mu_b_inter_s4c1 * (chill1*site4_1) + mu_b_inter_ws4 * (force*site4_1) +mu_b_inter_ps4 * (photo*site4_1) 
# 
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Chill portions", ylab = "Day of budburst")
points(site4_0Data$chillport.z2, site4_0Data$bb, col = "maroon")
points(site4_1Data$chillport.z2, site4_1Data$bb, col = "darkslategray4")
abline(lm(bb_csite4_0  ~ chill1), col = "darkred", lwd = 3)
abline(lm(bb_csite4_1 ~ chill1), col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("site 4_1"),
                            expression("site 4_0")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")

# # Begin by checking to see what cue is most important and whether there are strong correlations between cues:
# df.mean.t <- data.frame(bb.force = sumt[grep("b_force", rownames(sumt)), 1],
#                           bb.photo = sumt[grep("b_photo_ncp", rownames(sumt)), 1],
#                           bb.chill = sumt[grep("b_chill", rownames(sumt)), 1])
# 
# # df.mean.l <- data.frame(lat.force = sum1l[grep("b_force", rownames(sumt)), 1],
# #                         lat.photo = sum1l[grep("b_photo_ncp", rownames(sumt)), 1],
# #                         lat.chill = sum1l[grep("b_chill", rownames(sumt)), 1])
# 
# df.mean.t[which(df.mean.t$bb.force > df.mean.t$bb.photo), ] # species 11- rho alb
# # df.mean.l[which(df.mean.l$lat.force > df.mean.l$lat.photo), ] #none
# df.mean.t[which(df.mean.t$bb.chill > df.mean.t$bb.force), ] # 14
# # 3, 5,6,8,9,10,12,13,15,17,18,20
# # df.mean.l[which(df.mean.l$lat.chill > df.mean.l$lat.force), ] # 16
# #1,2,5,6,7,8,10,11,12,13,15,16,17,18,19,20
# 
# # all correlated
# summary(lm(bb.force~bb.photo, data=df.mean.t))
# summary(lm(bb.force~bb.chill, data=df.mean.t))
# summary(lm(bb.force~bb.photo, data=df.mean.t))
# # summary(lm(lat.force~lat.photo, data=df.mean.l))
# # summary(lm(lat.force~lat.chill, data=df.mean.l))
# # summary(lm(lat.chill~lat.photo, data=df.mean.l))
# 
# pdf(file.path( "figures/changes.pheno.pdf"), width = 7, height = 8)
# par(mfrow = c(2,1), mar = c(5, 10, 2, 1))
# # Upper panel: bud burst
# plot(seq(-22, 
#          12,
#          length.out = nrow(meanzt)), 
#      1:nrow(meanzt),
#      type = "n",
#      xlab = "",
#      ylab = "",
#      yaxt = "n")
# 
# #legend(x = -20, y = 2, bty="n", legend = "a. Budburst", text.font = 2)
# #rasterImage(bbpng, -20, 1, -16, 4)
# 
# axis(2, at = nrow(meanzt):1, labels = rownames(meanzt), las = 1, cex.axis = 0.8)
# points(meanzt[, 'mean'],
#        nrow(meanzt):1,
#        pch = 16,
#        col = "midnightblue")
# arrows(meanzt[, "75%"], nrow(meanzt):1, meanzt[, "25%"], nrow(meanzt):1,
#        len = 0, col = "black")
# abline(v = 0, lty = 3)
# # add advance/delay arrows
# par(xpd=NA)
# arrows(1, 15.5, 6, 15.5, len = 0.1, col = "black")
# legend(5, 16.5, legend = "delay", bty = "n", text.font = 1, cex = 0.75)
# arrows(-1, 15.5, -6, 15.5, len = 0.1, col = "black")
# legend(-12, 16.5, legend = "advance", bty = "n", text.font = 1, cex = 0.75)
# legend(-2, 16.5, legend = "0", bty = "n", text.font = 1, cex = 0.75)
# par(xpd = FALSE)
# 
# # par(mar = c(5, 10, 2, 1))
# # # Lower panel: leaf-out
# # plot(seq(-22, 
# #          12, 
# #          length.out = nrow(meanzl)), 
# #      1:nrow(meanzl),
# #      type = "n",
# #      xlab = "Model estimate change in day of phenological event",
# #      ylab = "",
# #      yaxt = "n")
# 
# #legend(x = -24, y = 6, bty="n", legend = "b. Leafout", text.font = 2)
# #rasterImage(lopng, -20, 1, -14, 4)
# 
# # axis(2, at = nrow(meanzl):1, labels = rownames(meanzl), las = 1, cex.axis = 0.8)
# # points(meanzl[,'mean'],
# #        nrow(meanzl):1,
# #        pch = 16,
# #        col = "midnightblue")
# # arrows(meanzl[,"75%"], nrow(meanzl):1, meanzl[,"25%"], nrow(meanzl):1,
# #        len = 0, col = "black")
# # abline(v = 0, lty = 3)
# 
# # add advance/delay arrows
# # par(xpd=NA)
# # arrows(1, 15.5, 6, 15.5, len=0.1, col = "black")
# # legend(5, 16.5, legend="delay", bty="n", text.font = 1, cex=0.75)
# # arrows(-1, 15.5, -6, 15.5, len=0.1, col = "black")
# # legend(-12, 16.5, legend="advance", bty="n", text.font = 1, cex=0.75)
# # legend(-2, 16.5, legend="0", bty="n", text.font = 1, cex=0.75)
# # par(xpd=FALSE)
# dev.off()
# 
# # Comparisons of trees vs shrubs:
# shrubs = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG", "alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
# trees = c("ACEPEN", "ACERUB", "ACESAC", "ALNINC", "BETALL", "BETLEN", "BETPAP", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL", "acegla","betpap", "poptre", "popbal")
# 
# treeshrub = levels(dx$sp)
# treeshrub[treeshrub %in% shrubs] = 1
# treeshrub[treeshrub %in% trees] = 2
# treeshrub = as.numeric(treeshrub)
