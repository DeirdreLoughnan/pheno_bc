## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

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
library(phylotools)

options(mc.cores = parallel::detectCores())

rm(list=ls()) 
options(stringsAsFactors = FALSE)

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc")
} else if(length(grep("Lizzie", getwd())>0)) {
  setwd("/home/faith/Documents/mnt/UBC/otherPeople/Deirdre")
}  else {
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

# mergeing the my data with DF
pheno <- dl.wchill
#pheno <- dl
############################################################

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
pheno$site.n <- as.numeric(pheno$site.n)

head(pheno)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
pheno.term <- pheno[,c("bb", "force.n", "photo.n", "site.n", "species", "lab2","Utah_Model","Chill_portions")]

#pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 
pheno.t <- pheno.term
pheno.t$force.z2 <- (pheno.t$force.n-mean(pheno.t$force.n,na.rm=TRUE))/(sd(pheno.t$force.n,na.rm=TRUE)*2)
pheno.t$photo.z2 <- (pheno.t$photo.n-mean(pheno.t$photo.n,na.rm=TRUE))/(sd(pheno.t$photo.n,na.rm=TRUE)*2)
pheno.t$chillport.z2 <- (pheno.t$Chill_portions-mean(pheno.t$Chill_portions,na.rm=TRUE))/(sd(pheno.t$Chill_portions,na.rm=TRUE)*2)
pheno.t$site.z2 <- (pheno.t$site.n-mean(pheno.t$site.n,na.rm=TRUE))/(sd(pheno.t$site.n,na.rm=TRUE)*2)

pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 19, 30 species, 47 with chill0 47 

nrow(pheno.term) - nrow(pheno.t) # 547 that had no terminal bb, 609
str(pheno.t)

pheno.t <- pheno.t[complete.cases(pheno.t$bb), ]

# now get the phylogeny and pair it with species names:
spInfo <- read.csv("input/species_list.csv")
pheno.t <- merge(pheno.t, spInfo, by = "species")

westSp <- unique(pheno.t$species.name)

tree <- read.tree("input/SBphylo_phenobc.tre")
head(tree$tip.label)
length(tree$tip.label) #47
tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
tree$tip.label[tree$tip.label=="Nyssa_sylvatica"] <- "Alnus_viridis"

treeW <- keep.tip(tree, westSp)

phymatchW <- data.frame(species.name = treeW$tip.label, sppnum = c(1:length(treeW$tip.label)))

d <- merge(pheno.t, phymatchW, by="species.name")

d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)
#nspecies <- 21
cophen_tree <- cophenetic(treeW)
vcv_tree <- vcv(treeW, cor = TRUE)

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


datalist <- with(pheno.t,
                 list( N=nrow(pheno.t),
                       n_sp = length(unique(pheno.t$species.fact)),
                       n_site = length(unique(pheno.t$site.n)),
                       lday = bb,
                       sp = species.fact,
                       chill1 = chillport.z2,
                       photo = photo.z2,
                       warm = force.z2,
                       site = site.z2,
                       Vphy = vcv_tree
                 ))

mdl.phylo <- stan("stan/bc_bb_2sites_standardized_phylogeny.stan",
                  data = datalist,
                  iter = 2000, warmup = 1000, chains=4
                  #, control = list(adapt_delta = 0.99)
)
save(mdl.phylo, file="output/dl_phylo_output.Rda")

load("output/dl_phylo_output.Rda")
ssm <-  as.shinystan(mdl.phylo)
launch_shinystan(ssm)
# # #
sum <- summary(mdl.phylo)$summary 
range(sum[, "n_eff"])
range(sum[, "Rhat"])

pairs(mdl.phylo, pars = c("a_z","b_z", "lam_interceptsa","sigma_interceptsa", "lam_interceptsb","sigma_interceptsb", "lp__"))

col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")
# manually to get right order
mu_params <- c("a_z",
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
               "sigma_y"
)


meanz_phylo <- sum[mu_params, col4table]
meanz_phylo


post <- rstan::extract(mdl.phylo)
h1 <- hist(rbeta(1000, 4,1))
h2 <- hist(post$mu_a)
plot(h2, col=rgb(0,0,1,1/4), xlim = c(0, 400))
plot(h1, col=rgb(1,0,1,1/4), add = T)
