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
#library(shinystan)
library(reshape2)
library(bayesplot)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(stringr)
library(phytools)

options(mc.cores = parallel::detectCores())


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
 setwd("~/Documents/github/pheno_bc") # for midge
}


dl <- read.csv("input/dl_allbb_mini.csv")

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

#head(pheno)
# combined the data has 3197 unique samples
############################################################
# Preping the data for the model
# remove chill 1.5 - chill 2
pheno <- subset(pheno, chill != "chill2")

pheno$force.n <- pheno$force

# west mean temp:
# (12*20+12*10)/24 #15
# (12*15+12*5)/24 #10
pheno$force.n[pheno$force.n == "HF" & pheno$population == "mp"] <- "15"
pheno$force.n[pheno$force.n == "HF" & pheno$population == "sm"] <- "15"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "mp"] <- "10"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "sm"] <- "10"

pheno$force.n[pheno$force.n == "HF" & pheno$population == "HF"] <- "13.33"
pheno$force.n[pheno$force.n == "HF" & pheno$population == "SH"] <- "13.33"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "HF"] <- "8.33"
pheno$force.n[pheno$force.n == "LF" & pheno$population == "SH"] <- "8.33"

pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "12"
pheno$photo.n[pheno$photo.n == "LP"] <- "8"
pheno$photo.n <- as.numeric(pheno$photo.n)

pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.numeric(pheno$site.n)

#head(pheno)
#add dummy/ site level effects:
pheno <- pheno %>%
  mutate ( site2 = if_else(site.n == 2, 1, 0),
           site3 = if_else(site.n == 3, 1, 0),
           site4 = if_else(site.n == 4, 1, 0))

pheno$site2 <- as.numeric(pheno$site2)
pheno$site3 <- as.numeric(pheno$site3)
pheno$site4 <- as.numeric(pheno$site4)

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
pheno.term <- pheno[,c("bb", "force.n", "force.z2", "photo.n","photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site2","site3","site4","site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term), ] # 3609

pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 49 

# now get the phylogeny and pair it with species names:
spInfo <- read.csv("input/species_list.csv")
pheno.t <- merge(pheno.t, spInfo, by = "species")

tree <- read.tree("input/SBphylo_phenobc.tre")
head(tree$tip.label)
length(tree$tip.label) #47

tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
tree$tip.label[tree$tip.label=="Alnus_alnobetula"] <- "Alnus_viridis"
tree$tip.label[tree$tip.label== "Fagus_grandifolia_var._caroliniana"] <- "Fagus_grandifolia"
tree$tip.label[tree$tip.label== "Spiraea_alba_var._latifolia"] <- "Spiraea_alba"

phymatch <- data.frame(species.name = tree$tip.label, sppnum = c(1:length(tree$tip.label)))

d <- merge(pheno.t, phymatch, by="species.name")

d <- d[order(d$sppnum),]
nspecies <- max(d$sppnum)
#nspecies <- 21
cophen_tree <- cophenetic(tree)
vcv_tree <- vcv(tree, cor = TRUE)

simu_inits <- function(chain_id) {
  a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
  return(append(list(
    a = a_z.temp),
    phypriors))
}

unique(pheno.t$force.n)

datalist.z2 <- with(pheno.t,
                 list( N= nrow(pheno.t),
                       n_sp = length(unique(pheno.t$species.fact)),
                       n_site = length(unique(pheno.t$population)),
                       lday = bb,
                       sp = species.fact,
                       chill1 = chillport.z2,
                       photo = photo.z2,
                       warm = force.z2,
                       site2 = site2.z2,
                       site3 = site3.z2,
                       site4 = site4.z2,
                       Vphy = vcv_tree
                       ))

mdl.3 <- stan("stan/df_mdl_4sites_again_allint_ncp_phylogeny_natural_newPrior.stan",
                    data = datalist.z2,
                    iter = 4000, warmup = 3000, chains=4,
                    include = FALSE, pars = c("ypred_new","y_hat")
                    #, control = list(adapt_delta = 0.99)
)


save(mdl.3, file="output/bb_phylo.Rda")

load("output/bb_phylo_contphotothermo_2zscoredMay13.Rda")
load("output/bb_phylo_triple.Rda")
sum <- summary(mdl.3)$summary
sumo <- summary(mdl.2)$summary

col4table <- c("mean","2.5%","50%","97.5%")

mu_params_4 <- c("a_z",
                 "lam_interceptsa",
                 "sigma_interceptsa",
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
                 "sigma_b_inter_wp",
                 "sigma_b_inter_wc1",
                 "sigma_b_inter_pc1",
                 "sigma_b_inter_s2c1",
                 "sigma_b_inter_ws2",
                 "sigma_b_inter_ps2",
                 "sigma_b_inter_s3c1",
                 "sigma_b_inter_ws3",
                 "sigma_b_inter_ps3",
                 "sigma_b_inter_s4c1",
                 "sigma_b_inter_ws4",
                 "sigma_b_inter_ps4",
                 "sigma_y")
meanz1 <- sum[mu_params_4, col4table]
round(meanz1,2)

old <- data.frame(rstan::extract(mdl.2))
new <- data.frame(rstan::extract(mdl.3))

old <- old[, mu_params_4]
new <- new[, mu_params_4]

plot(new[,"b_site2"] ~ old[,"b_site2"])
plot(new[,"sigma_b_inter_s4c1"] ~ old[,"sigma_b_inter_s4c1"])

b_chill <- data.frame(sum[grep("b_chill1\\[", rownames(sum)), c("mean")])
b_chill_o <- data.frame(sumo[grep("b_chill1\\[", rownames(sumo)), c("mean","2.5%", "97.5%", "n_eff", "Rhat")])

plot(sum[grep("a_sp\\[", rownames(sum)), c("mean")] ~ sumo[grep("a_sp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "a_sp");abline(0, 1)
plot(sum[grep("b_chill1\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_chill1\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "b_chill");abline(0, 1)
plot(sum[grep("b_force_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_force_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "b_force");abline(0, 1)
plot(sum[grep("b_photo_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_photo_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "b_photo");abline(0, 1)
plot(sum[grep("b_inter_wp_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_wp_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn force x photo");abline(0, 1)
plot(sum[grep("b_inter_wc1_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_wc1_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn force x chill");abline(0, 1)
plot(sum[grep("b_inter_pc1_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_pc1_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn chill x photo");abline(0, 1)
plot(sum[grep("b_inter_s2c1_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_s2c1_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn mp x chill");abline(0, 1)
plot(sum[grep("b_inter_ws2_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_ws2_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn mp x force");abline(0, 1)
plot(sum[grep("b_inter_ps2_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_ps2_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn mp x photo");abline(0, 1)
plot(sum[grep("b_inter_s3c1_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_s3c1_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn hf x chill");abline(0, 1)
plot(sum[grep("b_inter_ws3_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_ws3_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn hf x force");abline(0, 1)
plot(sum[grep("b_inter_ps3_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_ps3_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn hf x photo");abline(0, 1)
plot(sum[grep("b_inter_s4c1_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_s4c1_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn sh x chill");abline(0, 1)
plot(sum[grep("b_inter_ws4_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_ws4_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn sh x force");abline(0, 1)
plot(sum[grep("b_inter_ps4_ncp\\[", rownames(sum)), c("mean")] ~ sumo[grep("b_inter_ps4_ncp\\[", rownames(sumo)), c("mean")], xlab = "Constrained model", ylab ="Wide prior model", main = "intrxn sh x photo");abline(0, 1)
