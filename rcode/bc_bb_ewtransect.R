## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

# Adding standardization: following the methods outlined by Andrew Gelman in his paper

# changing from site to transect, with transect added as a 0/1 dummy variable, west = 0, east = 1

#Feb 28: model consistently has 7 divergent transitions, these persisted when I increased the iterations to 3000 warmup, 4000 iterations, and increased the adapt delta to 0.99 - from the log posterior vs parameter plots, the worst are sigma_chill with them all butting up against one side, sigma_b_inter_pc -- I am going to try ncp for chilling as well as photo and see if that helps

# March 1: the ncp for chilling produced 79 divergent transitions. The fact that this made it worse got me thinking that maybe the issue is the model has more complexity then it needs, so I removed the ncp on photo

# the model with no ncp on photo or chilling produced dt, running it for 5000 warmup and 6000 iterations produced only 1 dt but the rhat were bad at 1.13

# next I tried running no ncp on photo and chill, but with a bigger adapt delta, 4000:3000 iterations, there were no dt, but the rhat was really bad at 1.39 and the max tree depth exceeded

#March 3: running for 6000:5000 iterations and a greater adapt_delta (0.99) did not help; still 1 dt and poor rhat 1.17

# Mdl with 4000:3000 and an increased tree depth and adapt delat produced no dt, but chains did not mix and the rhat were still bad 1.29

# mdl with cues and site standardized produced 5 divergent transitions and low ess - tried increasing the warmup back up to 4000:300 - 11 divergent transitions - classic banana shape for chilling

# mdl with chilling, forcing, photoperiod, interactions all ncp - 4000:3000 and no increased adapt delta
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
library(stringr)
library(plyr)
library(phytools)
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

dl <- read.csv("input/dl_allbb.csv")

temp <- str_split_fixed(dl$trt, "_", 3); head(temp)
dl$chill<- temp[,1]
dl$photo <- temp[,2]
dl$force <- temp[,3]

dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")
dl.wchill$transect <- "west"

df <- read.csv("input/df_dxb_prepped_data.csv")
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill <- df.wchill[, c("population", "chill","force","photo","lab2", "bb","species", "treatment","Chill_portions","Utah_Model")]
df.wchill$transect <- "east"

# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)
#pheno <- dl

#pheno$first <- ifelse(pheno$tbb < pheno$latbb1,"t", ifelse (pheno$tbb == pheno$latbb1,"tl", "l"))

head(pheno)
#table(pheno$species, pheno$first)
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

pheno$transect.z2 <- (pheno$transect.n-mean(pheno$transect.n,na.rm=TRUE))/(sd(pheno$transect.n,na.rm=TRUE)*2)

#going to split it into analyses of terminal bb and lateral bb
# Starting with the terminal buds:
pheno.term <- pheno[,c("bb", "force.n", "photo.n", "transect.n", "species", "lab2","Utah_Model","Chill_portions","force.z2", "photo.z2", "chillport.z2", "utah.z2","transect.z2")]

pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 49 

# now get the phylogeny and pair it with species names:
spInfo <- read.csv("input/species_list.csv")

head(pheno.t)
# change the df to lowercase:

pheno.t <- merge(pheno.t, spInfo, by = "species")

tree <- read.tree("input/SBphylo_phenobc.tre")
head(tree$tip.label)
length(tree$tip.label) #47

tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
tree$tip.label[tree$tip.label=="Nyssa_sylvatica"] <- "Alnus_viridis"

phymatch <- data.frame(species.name = tree$tip.label, sppnum = c(1:length(tree$tip.label)))

d <- merge(pheno.t, phymatch, by="species.name")

d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)
#nspecies <- 21
cophen_tree <- cophenetic(tree)
vcv_tree <- vcv(tree, cor = TRUE)

phypriors <- list( 
  a_z_prior_mu = 4, # true value
  a_z_prior_sigma = 1,
  lam_interceptsa_prior_alpha = 4, # 
  lam_interceptsa_prior_beta = 6, # 
  sigma_interceptsa_prior_mu = 0.2, # true value
  sigma_interceptsa_prior_sigma = 0.2
)

simu_inits <- function(chain_id) {
  a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
  return(append(list(
    a = a_z.temp),
    phypriors))
}

# head(pheno.t)
datalist_z <- list( N=nrow(pheno.t),
                    Nsite = nrow(pheno.t),
                    n_sp = length(unique(pheno.t$species.fact)),
                    n_site = length(unique(pheno.t$site.n)),
                    lday = pheno.t$bb,
                    sp = pheno.t$species.fact,
                    chill1 = pheno.t$Chill_portions,
                    photo = pheno.t$photo.z2,
                    warm = pheno.t$force.z2,
                    site = pheno.t$transect.z2,
                    Vphy = vcv_tree)

# mdl <- stan("stan/lday_site_sp_chill_inter_poola_ncpwp_2chill_sitefixedstnd_standardized.stan",
#             data = datalist_z,
#             #include = FALSE, pars = c("ypred_new","y_hat"),
#             iter = 4000, warmup = 3000, chains=4
#            , control = list(adapt_delta= 0.99))
# save(mdl, file="output/tbb_cport_stnd_stndsite_ew_adapt.Rda")

mdl.ewphylo <- stan("stan/bc_bb_2sites_standardized_phylogeny_newpriors.stan",
            data = datalist_z,
            #include = FALSE, pars = c("ypred_new","y_hat"),
            iter = 6000, warmup = 4000, chains=4
            , control = list(adapt_delta= 0.99))
save(mdl.ewphylo, file="output/ew_phylo_output_newpriors.Rda")
# 

#lday_site_sp_chill_inter_poola_ncpwp_2chill_sitefixed_standardized: 11 div trans, low ess  

# 2. try running with a greater adapt delta 
#######################################################################

# load("output/final/ew_allncp.Rda")
load("output/final/ew_phylo_output_newpriors.Rda")
# # # # # #
ssm <-  as.shinystan(mdl.ewphylo)
launch_shinystan(ssm)
# # #
sum <- summary(mdl.ewphylo)$summary 
range(sum[, "n_eff"]) # low
range(sum[, "Rhat"])


# summary(mdl)$summary[c("a_z",
#                        "lam_interceptsa",
#                        "sigma_interceptsa",
#                        "mu_b_warm",
#                        "mu_b_photo",
#                        "mu_b_chill1",
#                        "b_site",
#                        "mu_b_inter_wp",
#                        "mu_b_inter_wc1",
#                        "mu_b_inter_pc1",
#                        "mu_b_inter_ws",
#                        "mu_b_inter_ps",
#                        "mu_b_inter_sc1",
#                        "sigma_b_warm",
#                        "sigma_b_photo",
#                        "sigma_b_chill1",
#                        #"sigma_a",
#                        "sigma_b_inter_wp",
#                        "sigma_b_inter_wc1",
#                        "sigma_b_inter_pc1",
#                        "sigma_b_inter_ws",
#                        "sigma_b_inter_ps",
#                        "sigma_b_inter_sc1",
#                        "sigma_y"),"mean"]

summary(mdl.ewphylo)$summary[c("a_z",
                           "lam_interceptsa",
                           "sigma_interceptsa",
                       "mu_b_warm",
                       "mu_b_photo",
                       "mu_b_chill1",
                       "b_site",
                       "mu_b_inter_wp",
                       "mu_b_inter_wc1",
                       "mu_b_inter_pc1",
                       "mu_b_inter_ws",
                       "mu_b_inter_ps",
                       "mu_b_inter_sc1",
                       "sigma_b_warm",
                       "sigma_b_photo",
                       "sigma_b_chill1",
                       #"sigma_a",
                       "sigma_b_inter_wp",
                       "sigma_b_inter_wc1",
                       "sigma_b_inter_pc1",
                       "sigma_b_inter_ws",
                       "sigma_b_inter_ps",
                       "sigma_b_inter_sc1",
                       "sigma_y"),"mean"]

fit <- mdl.ewphylo
y_rep <- as.matrix(fit, pars = "y_hat")

y <- pheno.t$bb
#ppc_hist(y, y_rep[1:8, ], binwidth = 1)

pdf("yvsypred_ew_phylo.pdf")
ppc_dens_overlay(y, y_rep[1:100, ])
dev.off()


# pairs(mdl, pars = c("mu_a", "mu_force","mu_chill","mu_photo","mu_site","mu_inter_fp", "mu_inter_fs","mu_inter_ps","mu_inter_fc","mu_inter_pc","mu_inter_sc", "lp__"))
# 
# pairs(mdl, pars = c("sigma_a", "sigma_force","sigma_chill","sigma_photo","sigma_site","sigma_b_inter_fp","sigma_b_inter_fs","sigma_b_inter_ps","sigma_b_inter_fc","sigma_b_inter_pc","sigma_b_inter_sc", "sigma_y", "lp__"))

# summary(mdl)$summary[c("mu_a", "mu_b_warm", "mu_b_photo", "mu_b_chill1", "b_site",
#                             "mu_b_inter_wp","mu_b_inter_wc1","mu_b_inter_pc1","mu_b_inter_ws","mu_b_inter_ps","mu_b_inter_sc1",
#                             "sigma_b_warm","sigma_b_photo", "sigma_b_chill1","sigma_a", "sigma_b_inter_wp", "sigma_b_inter_wc1", "sigma_b_inter_pc1","sigma_b_inter_ws","sigma_b_inter_ps", "sigma_b_inter_sc1", "sigma_y"),c("mean","25%","75%")]
# 
# ######################################################################
# # let's take a closer look at the interactions and plot the model output against the raw data:\
# pheno$transect <- as.factor(pheno$transect)
# 
# #pred <- sum[grep("y_hat", rownames(sum)), 1]
# a_sp = sum[grep("mu_a", rownames(sum)), 1]
# mu_b_warm = sum[grep("mu_b_warm", rownames(sum)), 1]
# mu_b_photo = sum[grep("mu_b_photo", rownames(sum)), 1]
# mu_b_chill1 = sum[grep("mu_b_chill1", rownames(sum)), 1]
# mu_b_inter_ws = sum[grep("mu_b_inter_ws", rownames(sum)), 1]
# mu_b_inter_sc1 = sum[grep("mu_b_inter_sc1", rownames(sum)), 1]
# mu_b_inter_ps = sum[grep("mu_b_inter_ps", rownames(sum)), 1]
# mu_b_inter_pc1 = sum[grep("mu_b_inter_pc1", rownames(sum)), 1]
# mu_b_inter_wp = sum[grep("mu_b_inter_wp", rownames(sum)), 1]
# mu_b_inter_wc1 = sum[grep("mu_b_inter_wc1", rownames(sum)), 1]
# b_site = sum[grep("b_site", rownames(sum)), 1]
# 
# ######################################################################
# ## transect and chilling:
# west <- subset(pheno.t, transect.n == "0")
# east <- subset(pheno.t, transect.n == "1")
# 
#   force <- -0.5080665 # zero forcing
#   photo <- -0.5044652 # zero photo
#   chill1 <- c( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)
#   west.sites <- unique(west$transect.z2)
#   east.sites <- unique(east$transect.z2)
#   
# # plot first for the west coast
# bb_westc = a_sp + b_site * west.sites + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
#   mu_b_inter_wp * (force*photo) + mu_b_inter_ws * (force*west.sites) +mu_b_inter_ps * (photo*west.sites) +
#   mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_sc1 * (chill1*west.sites)
# 
# # plot first for the east coast
# bb_eastc = a_sp + b_site *  east.sites + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
#   mu_b_inter_wp * (force*photo) +mu_b_inter_ws * (force* east.sites) +mu_b_inter_ps * (photo* east.sites) +
#   mu_b_inter_wc1 * (force*chill1) +mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_sc1 * (chill1* east.sites)
# 
# pdf("figures/ew_interactions.pdf", width = 10, height = 3.5)
# par(mfrow = c(1,3))
# plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Chill portions", ylab = "Day of budburst")
# points(west$chillport.z2,west$bb, col = "maroon")
# points(east$chillport.z2,east$bb, col = "darkslategray4")
# abline(lm(bb_westc ~ chill1), pch =19, col = "darkred", lwd = 3)
# abline(lm(bb_eastc ~ chill1), pch =19, col = "darkslategray", lwd = 3)
# 
# legend("topleft",legend = c(expression("eastern transect"),
#                             expression("western transect")),
#        col = c("darkslategray","maroon"),
#        inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")
# ######################################################################
# ######################################################################
# ## transect and forcing:
# west <- subset(pheno.t, transect.n == "0")
# east <- subset(pheno.t, transect.n == "1")
# 
# #force <- unique(pheno.t$force.z2)
# force <- seq(1:50)
# photo <- -0.5044652 # zero photo
# chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)
# west.sites <- unique(west$transect.z2)
# east.sites <- unique(east$transect.z2)
# 
# # plot first for the west coast
# bb_westf = a_sp + b_site * west.sites + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
#   mu_b_inter_wp * (force*photo) + mu_b_inter_ws * (force*west.sites) +mu_b_inter_ps * (photo*west.sites) +
#   mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_sc1 * (chill1*west.sites)
# 
# # plot first for the east coast
# bb_eastf = a_sp + b_site *  east.sites + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
#   mu_b_inter_wp * (force*photo) +mu_b_inter_ws * (force* east.sites) +mu_b_inter_ps * (photo* east.sites) +
#   mu_b_inter_wc1 * (force*chill1) +mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_sc1 * (chill1* east.sites)
# 
# plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Forcing", ylab = "Day of budburst")
# points(west$force.z2,west$bb, col = "maroon")
# points(east$force.z2,east$bb, col = "darkslategray4")
# abline(lm(bb_westf ~ force), pch =19, col = "darkred", lwd = 3)
# abline(lm(bb_eastf ~ force), pch =19, col = "darkslategray", lwd = 3)
# 
# legend("topleft",legend = c(expression("eastern transect"),
#                             expression("western transect")),
#        col = c("darkslategray","maroon"),
#        inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")
# ######################################################################
# ######################################################################
# ## transect and photoperiod:
# west <- subset(pheno.t, transect.n == "0")
# east <- subset(pheno.t, transect.n == "1")
# 
# #force <- unique(pheno.t$force.z2)
# force <- -0.5080665 # zero forcing
# photo <- seq(1:50)
# chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)
# west.sites <- unique(west$transect.z2)
# east.sites <- unique(east$transect.z2)
# 
# # plot first for the west coast
# bb_westp = a_sp + b_site * west.sites + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
#   mu_b_inter_wp * (force*photo) + mu_b_inter_ws * (force*west.sites) +mu_b_inter_ps * (photo*west.sites) +
#   mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_sc1 * (chill1*west.sites)
# 
# # plot first for the east coast
# bb_eastp = a_sp + b_site *  east.sites + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1 + 
#   mu_b_inter_wp * (force*photo) +mu_b_inter_ws * (force* east.sites) +mu_b_inter_ps * (photo* east.sites) +
#   mu_b_inter_wc1 * (force*chill1) +mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_sc1 * (chill1* east.sites)
# 
# plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Forcing", ylab = "Day of budburst")
# points(west$photo.z2,west$bb, col = "maroon")
# points(east$photo.z2,east$bb, col = "darkslategray4")
# abline(lm(bb_westp ~ photo), pch =19, col = "darkred", lwd = 3)
# abline(lm(bb_eastp ~ photo), pch =19, col = "darkslategray", lwd = 3)
# 
# legend("topleft",legend = c(expression("eastern transect"),
#                             expression("western transect")),
#        col = c("darkslategray","maroon"),
#        inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")
# dev.off()

post <- rstan::extract(mdl.ewphylo)
# histograms
h1 <- hist(rnorm(1000, 4,10))
h2 <- hist(post$a_z)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(0, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rbeta(1000, 1,1))
h2 <- hist(post$lam_interceptsa)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(0, 1))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 1,15))
h2 <- hist(post$sigma_interceptsa)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-50, 50))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,35))
h2 <- hist(post$mu_b_chill1)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-100, 100))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,5))
h2 <- hist(post$b_site)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-100, 100))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,50))
h2 <- hist(post$b_warm_ncp)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-100, 100))
plot(h1, col=rgb(1,0,1,1/4), add = T)

h1 <- hist(rnorm(1000, 0,10))
h2 <- hist(post$sigma_b_inter_wc1)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(-100, 100))
plot(h1, col=rgb(1,0,1,1/4), add = T)
