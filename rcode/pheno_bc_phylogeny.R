# Started Jan 11, 2022 by Deirdre

# the aim of this code is to plot the effects of different cues on the phylogeny and to 
# run a model with phylogeny on the intercept

rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

set.seed(194842)

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc")
} else if(length(grep("Lizzie", getwd())>0)) {
  setwd("~/Documents/git/teaching/stan")
} else if(length(grep("faith", getwd()) > 0)){
  setwd("/home/faith/Documents/mnt/UBC/otherPeople/Deirdre")
}  else {
  setwd("/home/deirdre/pheno_bc") # for midge
}

require(ape)
require(geiger)
require(phytools)
library(rstan)
require(shinystan)
require(bayesplot)
require(tidybayes)
require(truncnorm)
library(ggplot2)
library(dplyr)
library(plyr)
require(lme4)

df <- read.csv("input/day.of.bb.DFlynn.chill0.csv", header=TRUE, na.strings=c("","NA"))
head(df)
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))

dl <- read.csv("input/day.of.bb.DL.csv", header=TRUE, na.strings=c("","NA"))
head(dl)
dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")

# mergeing the my data with DF
# pheno <- rbind.fill(dl.wchill, df.wchill)
pheno <- dl.wchill

pheno$first <- ifelse(pheno$tbb < pheno$latbb1,"t", ifelse (pheno$tbb == pheno$latbb1,"tl", "l"))

head(pheno)
table(pheno$species, pheno$first)
#write.csv(pheno, "input/pheno.w5chill.csv")

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

# pg 153 of statistical rethinking: we have multiple categories, so we want to use index values as integers
# Lizzie pointed out that if you just convert it to integers, it is in alphabetical order, but to keep track of the sites, I am having smithers be 1
pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "0"
pheno$site.n[pheno$site.n == "mp"] <- "1"
pheno$site.n[pheno$site.n == "HF"] <- "2"
pheno$site.n[pheno$site.n == "SH"] <- "3"
pheno$site.n <- as.integer(pheno$site.n)

pheno$Chill_portions <- as.numeric(pheno$Chill_portions)

# z scoring values, but since some are binary I will use 2SD as per the Gelmen et al paper
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)
pheno$utah.z2 <- (pheno$Utah_Model-mean(pheno$Utah_Model,na.rm=TRUE))/(sd(pheno$Utah_Model,na.rm=TRUE)*2)

#going to split it into analyses of terminal bb and lateral bb
# Starting with the terminal buds:
pheno.term <- pheno[,c("tbb", "force.n", "photo.n", "site.n", "species", "lab2","Utah_Model","Chill_portions","force.z2", "photo.z2", "chillport.z2", "utah.z2")]

pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 

pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 

nrow(pheno.term) - nrow(pheno.t) # That had no terminal bb, 609

#creating dummy var
pheno.t <- pheno.t %>%
  mutate ( d2 = if_else(site.n == 2, 1, 0),
           d3 = if_else(site.n == 3, 1, 0),
           d4 = if_else(site.n == 4, 1, 0))

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
#building the required tree & dataframe 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
nsp = 8
nind = 10

param <- list(a_z = 4, # root value intercept
              lam_interceptsa = 0.4, # lambda intercept
              sigma_interceptsa = 0.2 # rate of evolution intercept
)
# Set priors
phypriors <- list(
  mu_grand_prior_mu = 50,
  mu_grand_prior_sigma = 5,
  a_z_prior_mu = 4, # true value
  a_z_prior_sigma = 1,
  lam_interceptsa_prior_alpha = 4, # 
  lam_interceptsa_prior_beta = 6, # 
  sigma_interceptsa_prior_mu = 0.2, # true value
  sigma_interceptsa_prior_sigma = 0.2,
  #b_zf_prior_mu = 0.6, # true value
  #b_zf_prior_sigma = 1,
  # lam_interceptsbf_prior_alpha = 7, #
  # lam_interceptsbf_prior_beta = 3, # 
  # sigma_interceptsbf_prior_mu = 0.1, # true value
  # sigma_interceptsbf_prior_sigma = 0.1,
  sigma_y_mu_prior = 0.01,
  sigma_y_mu_sigma = 1,
  mu_grand_prior_mu = 50,
  mu_grand_prior_sigma = 5,
  mu_force_prior_mu = 0,
  mu_force_prior_sigma = 35,
  mu_chill_prior_mu = 0,
  mu_chill_prior_sigma = 35,
  mu_photo_prior_mu = 0,
  mu_photo_prior_sigma = 35,
  sigma_force_prior_mu = 1,
  sigma_force_prior_sigma = 5,
  sigma_chill_prior_mu = 1,
  sigma_chill_prior_sigma = 5,
  sigma_photo_prior_mu = 1,
  sigma_photo_prior_sigma = 5,
  sigma_a_prior_mu = 5,
  sigma_a_prior_sigma = 5)

spetree <- pbtree(n=nsp, nsim=1, b=1, complete=FALSE,scale=1)
spetree$tip.label <- paste("s", 1:nsp, sep="")

# Generate intercept
scaledtree_intercept <- rescale(spetree, model = "lambda", param[["lam_interceptsa"]])         
intercepts <- fastBM(scaledtree_intercept, a = param[["a_z"]], mu = 0, sig2 = param[["sigma_interceptsa"]] ^ 2)
# Generate bf slope
# scaledtree_bf <- rescale(spetree, model = "lambda", param[["lam_interceptsbf"]])         
# slopes_bf <- fastBM(scaledtree_bf, a = param[["b_zf"]], mu = 0, sig2 = param[["sigma_interceptsbf"]] ^ 2)

# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
#building the dataframe 
# <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
#Set number of data groups
nsite = 4 #number of sites with trait data
nChill <- 6 #This many chilling treatment levels 
nPhoto <- 2
nForce <- 2
nPhenCombo <- nChill*nPhoto*nForce

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

#interactions:
mu_fp <- 3
mu_cp <- 2
mu_fc <- 1

sigma_fp <- 0.5
sigma_cp <- 0.5
sigma_fc <- 0.5

mu_fsite <- 1
mu_csite <- 5
mu_psite  <- 3

sigma_fsite <- 0.5
sigma_csite <- 0.5
sigma_psite <- 0.5

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

fake$alpha.pheno.sp <- rep(intercepts, each =  nsite*nPhenCombo*nind)

alpha.force.sp <- rnorm(nsp, mu_force, sigma_force)
fake$alpha.force.sp <- rep(alpha.force.sp, each =  nsite*nPhenCombo*nind)

alpha.photo.sp <- rnorm(nsp, mu_photo, sigma_photo)
fake$alpha.photo.sp <- rep(alpha.photo.sp, each = nsite*nPhenCombo*nind)

alpha.chill.sp <- rnorm(nsp, mu_chill, sigma_chill)
fake$alpha.chill.sp <- rep(alpha.chill.sp, each =  nsite*nPhenCombo*nind)

alpha.site <- rnorm(nsite, b_site, sigma_site)
fake$alpha.site <- rep(rep(alpha.site, times = nPhenCombo, each = nPhenCombo*nind))

# adding the interaction effects:
alpha.fc <- rnorm(nsp, mu_fc, sigma_fc)
fake$alpha.fc <- rep(alpha.fc, each =  nsite*nForce * nPhoto * nind)

alpha.cp <- rnorm(nsp, mu_cp, sigma_cp)
fake$alpha.cp <- rep(alpha.cp, each =  nsite*nChill * nPhoto *nind)

alpha.fp <- rnorm(nsp, mu_fp, sigma_fp)
fake$alpha.fp <- rep(alpha.fp, each =  nsite*nForce * nPhoto *nind)

#adding interactions between the cues and the 
alpha.fsite2 <- rnorm(nsite, mu_fsite, sigma_fsite)
fake$alpha.fsite2 <- rep(alpha.fsite2, each =  nPhenCombo*nind)

alpha.csite2 <- rnorm(nsite, mu_csite, sigma_csite)
fake$alpha.csite2 <- rep(alpha.csite2, each =  nPhenCombo*nind)

alpha.psite2 <- rnorm(nsite, mu_psite, sigma_psite)
fake$alpha.psite2 <- rep(alpha.psite2, each =  nPhenCombo*nind)

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
fake$chill<- rep(rep(rnorm(n = nChill, mean = 0, sd = 10), each = nForce*nPhoto), times = nsp*nsite*nind)


# for dummy variable test data
fake$bb <-  mu_grand + alpha.site[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
  fake$gen.var

fake$bb_int <-  mu_grand + alpha.site[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 + fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) +
  fake$gen.var

fake$bb_int_site <-  mu_grand + alpha.site[1] + alpha.csite2[1] + alpha.fsite2[1] + alpha.psite2[1] + fake$alpha.pheno.sp + fake$alpha.force.sp *  fake$warm + 
  fake$alpha.chill.sp * fake$chill + fake$alpha.photo.sp * fake$photo + 
  fake$alpha.site*fake$d2  + fake$alpha.site*fake$d3  + fake$alpha.site*fake$d4 +
  fake$alpha.fc * (fake$warm*fake$chill) + fake$alpha.cp * (fake$photo*fake$chill) + fake$alpha.fp * (fake$photo*fake$warm) +
  fake$alpha.fsite2 * (fake$d2 * fake$warm) + fake$alpha.fsite2 * (fake$d3 * fake$warm) + fake$alpha.fsite2 * (fake$d4 * fake$warm) +
  fake$alpha.csite2 * (fake$d2 * fake$chill) + fake$alpha.csite2 * (fake$d3 * fake$chill) + fake$alpha.csite2 * (fake$d4 * fake$chill) +
  fake$alpha.psite2 * (fake$d2 * fake$photo) + fake$alpha.psite2 * (fake$d3 * fake$photo) + fake$alpha.psite2 * (fake$d4 * fake$photo) +
  fake$gen.var


# check if works with lmer or brms

summary(lm(bb ~  warm + photo + chill + d2 + d3 + d4, data = fake)) # 

summary(lm(bb_int ~  warm + photo + chill + d2 + d3 + d4 + warm * chill + chill * photo + photo * warm, data = fake)) # 

summary(lm(bb_int_site ~  warm + photo + chill + d2 + d3 + d4 +
             warm * chill + chill * photo + photo*warm +
             d2 * warm + d3 * warm + d4 * warm + d2 * chill + d3 * chill + d4 * chill + d2 * photo + d3 * photo + d4 * photo
           , data = fake)) # 



mu_grand + alpha.site[1]
alpha.site
# The values are very close between trail data and lmer output

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
                  site4 = fake$d4,
                  Vphy = vcv(spetree, corr = TRUE),
                  phypriors)


# Function for generating "good" initial values
simu_inits <- function(chain_id) {
  a_z.temp <- rnorm(n = nsp, mean = param[["a_z"]], sd = 1)
  #b_z.temp <- rnorm(n = nspecies, mean = param[["b_zf"]], sd = 1)
  return(append(list(a = a_z.temp),
                param))
}

mdl.phylo <- stan("stan/test_model_phylogeny.stan",
                 data = append(list(N=nrow(fake),
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
                                    site4 = fake$d4,
                                    Vphy = vcv(spetree, corr = TRUE)),phypriors),
                 init = simu_inits,
                 iter = 4000,
                 warmup = 2000,
                 chains = 4,
                 seed = 62921
)
save(mdl.phylo, file = "output/simPhylo.Rda")

###########################################
load("output/simPhylo.Rda")
# # 
ssm <-  as.shinystan(mdl.phylo)
launch_shinystan(ssm)
# #
 get_variables(mdl.phylo)
summary(mdl.phylo)$summary[c("mu_grand","mu_force", "mu_chill","mu_photo","a_z","sigma_force","sigma_chill","sigma_photo","sigma_y", "lam_interceptsa","sigma_interceptsa"),"mean"]
#    

###################################################
# plotting the output from the western only model onto the phylogeny:

load("output/tbb_ncp_chillportions_zsc_dl.Rda")

tree <- read.tree("input/SBphylo_phenobc.tre")

sumer <- summary(mdl.t)$summary

betaFSp <- data.frame(sumer[grep("b\\_force", rownames(sumer)), "mean"])
colnames(betaFSp)[colnames(betaFSp) == "sumer.grep..b.._force...rownames.sumer.....mean.."] <- "betaFSp"

betaCSp <- data.frame(sumer[grep("b\\_chill", rownames(sumer)), "mean"])
colnames(betaCSp)[colnames(betaCSp) == "sumer.grep..b.._chill...rownames.sumer.....mean.."] <- "betaCSp"

betaPSp <- data.frame(sumer[grep("b\\_photo", rownames(sumer)), "mean"]); 
betaPSp <- as.data.frame(betaPSp[c(20:38),])
colnames(betaPSp)[colnames(betaPSp) == "betaPSp[c(20:38), ]"] <- "betaPSp"

intSp <- data.frame(sumer[grep("a\\_sp", rownames(sumer)), "mean"])
colnames(intSp)[colnames(intSp) == "sumer.grep..a.._sp...rownames.sumer.....mean.."] <- "intSp"

phylo.dat <- pheno.t[,c("species","species.fact")]d
phylo.dat <- phylo.dat[order(phylo.dat$species.fact),]

phylo.dat$count <- 1
phylo.dat <- aggregate(phylo.dat["count"], phylo.dat[c("species", "species.fact")], FUN = sum);phylo.dat <- phylo.dat[,c("species","species.fact")]

dat.slope <- cbind(phylo.dat, betaFSp, betaCSp, betaPSp, intSp)
head(dat.slope)

## B) generate a comparative.data object merging data and phylogeny
databbslopesphy = comparative.data(tree, dat.slope,names.col="species",
                                   na.omit=TRUE,vcv=TRUE, warn.dropped = TRUE)

phyloplot = databbslopesphy$phy
x = databbslopesphy$data$betaForceSpMean
y = databbslopesphy$data$betaChillSpMean
z = databbslopesphy$data$betaPhotoSpMean
names(x) = names(y) = names(z) = databbslopesphy$phy$tip.label

pdf(paste("figures/force", files[i], ".pdf", sep = ""), width = 10)
par(mfrow=c(1,3))
force <- contMap(phyloplot, x, lwd = 2.5, outline = F,fsize = c(0.8,1))
chill <- contMap(phyloplot, y, lwd = 2.5, outline = F,fsize = c(0.8,1))
photo <- contMap(phyloplot, z, lwd = 2.5, outline = F,fsize = c(0.8,1))
dev.off()

