# code started Sept 7, 2022 to analyze the budburst date of the first lateral bud across the entire plant

# running the new, 4 site model

rm(list=ls()) 
options(stringsAsFactors = FALSE)
library(rstan)
library(shinystan)
library(dplyr)
library(plyr)
library(stringr)
library(phytools)

options(mc.cores = parallel::detectCores())


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  

if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc")
} else if(length(grep("Lizzie", getwd())>0)) {
  setwd("~/Documents/git/projects/others/deirdre/synchrony")
} else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

# dl <- read.csv("input/day.of.bb.DL.csv", header=TRUE, na.strings=c("","NA"))
# head(dl)
# dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")
# 
# dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
# dl.wchill$lab3 <- dl.wchill$lab2f
# dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")
# 
# # load in the df data
# df <- read.csv("input/day.of.bb.DFlynn.chill0.csv", header=TRUE, na.strings=c("","NA"))
# head(df)
# df.chill <- read.csv("input/chilling_values_eastern.csv")
# df.wchill <- merge(df, df.chill, by =c("population","chill"))
# 
# pheno <- rbind.fill(dl.wchill, df.wchill)
# 
# pheno$first <- ifelse(pheno$tbb < pheno$latbb1,"t", ifelse (pheno$tbb == pheno$latbb1,"tl", "l"))
# 
# head(pheno)
# #table(pheno$species, pheno$first)
# ##################################################################################
# # Prep the data for the model
# pheno$force.n <- pheno$force
# pheno$force.n[pheno$force.n == "HF"] <- "1"
# pheno$force.n[pheno$force.n == "LF"] <- "0"
# pheno$force.n <- as.numeric(pheno$force.n)
# 
# pheno$photo.n <- pheno$photo
# pheno$photo.n[pheno$photo.n == "HP"] <- "1"
# pheno$photo.n[pheno$photo.n == "LP"] <- "0"
# pheno$photo.n <- as.numeric(pheno$photo.n)
# 
# pheno$site.n <- pheno$population
# pheno$site.n[pheno$site.n == "sm"] <- "1"
# pheno$site.n[pheno$site.n == "mp"] <- "2"
# pheno$site.n[pheno$site.n == "HF"] <- "3"
# pheno$site.n[pheno$site.n == "SH"] <- "4"
# pheno$site.n <- as.numeric(pheno$site.n)
# 
# head(pheno)
# #add dummy/ site level effects:
# pheno <- pheno %>%
#   mutate ( site2 = if_else(site.n == 2, 1, 0),
#            site3 = if_else(site.n == 3, 1, 0),
#            site4 = if_else(site.n == 4, 1, 0))
# 
# # standardize the 0/1 and standardize sites? 
# pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
# pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
# pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)
# 
# pheno$site2.z2 <- (pheno$site2-mean(pheno$site2,na.rm=TRUE))/(sd(pheno$site2,na.rm=TRUE)*2)
# pheno$site3.z2 <- (pheno$site3-mean(pheno$site3,na.rm=TRUE))/(sd(pheno$site3,na.rm=TRUE)*2)
# pheno$site4.z2 <- (pheno$site4-mean(pheno$site4,na.rm=TRUE))/(sd(pheno$site4,na.rm=TRUE)*2)
# 
# #going to split it into analysis of terminal bb and lateral bb
# # Starting with the terminal buds:
# #pheno.lat <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
# pheno.lat <- pheno[,c("latbb1", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site2.z2", "site3.z2","site4.z2")]
# #pheno.l <- pheno.lat[complete.cases(pheno.lat), ] # 3609
# 
# pheno.l <- pheno.lat[complete.cases(pheno.lat$latbb1), ] #3640
# pheno.l$species <- tolower(pheno.l$species)
# pheno.l$species.fact <- as.numeric(as.factor(pheno.l$species))
# sort(unique(pheno.l$species.fact)) # 49 
# 
# # now get the phylogeny and pair it with species names:
# spInfo <- read.csv("input/species_list.csv")
# 
# head(pheno.l)
# # change the df to lowercase:
# 
# pheno.l <- merge(pheno.l, spInfo, by = "species")
# 
# tree <- read.tree("input/SBphylo_phenobc.tre")
# head(tree$tip.label)
# length(tree$tip.label) #47
# 
# tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
# tree$tip.label[tree$tip.label=="Nyssa_sylvatica"] <- "Alnus_viridis"
# tree$tip.label[tree$tip.label== "Fagus_grandifolia_var._caroliniana"] <- "Fagus_grandifolia"
# tree$tip.label[tree$tip.label== "Spiraea_alba_var._latifolia"] <- "Spiraea_alba"
# 
# phymatch <- data.frame(species.name = tree$tip.label, sppnum = c(1:length(tree$tip.label)))
# 
# d <- merge(pheno.l, phymatch, by="species.name")
# 
# d <- d[order(d$sppnum),]
# nspecies <- max(d$sppnum)
# #nspecies <- 21
# cophen_tree <- cophenetic(tree)
# vcv_tree <- vcv(tree, cor = TRUE)
# 
# phypriors <- list( 
#   a_z_prior_mu = 4, # true value
#   a_z_prior_sigma = 1,
#   lam_interceptsa_prior_alpha = 4, # 
#   lam_interceptsa_prior_beta = 6, # 
#   sigma_interceptsa_prior_mu = 0.2, # true value
#   sigma_interceptsa_prior_sigma = 0.2
# )
# 
# simu_inits <- function(chain_id) {
#   a_z.temp <- rnorm(n = nspecies, mean = phypriors[["a_z_prior_mu"]], sd = phypriors[["a_z_prior_sigma"]])
#   return(append(list(
#     a = a_z.temp),
#     phypriors))
# }
# 
# 
# 
# datalist.z2 <- with(pheno.l,
#                     list( N=nrow(pheno.l),
#                           n_sp = length(unique(pheno.l$species.fact)),
#                           n_site = length(unique(pheno.l$population)),
#                           lday = latbb1,
#                           sp = species.fact,
#                           chill1 = chillport.z2,
#                           photo = photo.z2,
#                           warm = force.z2,
#                           site2 = site2.z2,
#                           site3 = site3.z2,
#                           site4 = site4.z2,
#                           Vphy = vcv_tree
#                     ))


# mdl.l <- stan("stan/df_mdl_4sites_again_allint_ncp.stan",
#               data = datalist.z2,
#               iter = 4000, warmup =3000, chains=4
#               #, control = list(adapt_delta = 0.99)
#               )
# save(mdl.l, file="output/tbb_4sites_fullystandardized_ncp.Rda")

# mdl.4phylo <- stan("stan/df_mdl_4sites_again_allint_ncp_phylogeny.stan",
#                    data = datalist.z2,
#                    iter = 6000, warmup =3000, chains=4,
#                    include = FALSE, pars = c("ypred_new","y_hat")
#                    #, control = list(adapt_delta = 0.99)
# )
# save(mdl.4phylo, file="output/bb_4sites_phylo.Rda")


############################################################
############################################################

# Plot the model with just my data comparing 1lat to 50 lat:
dl <- read.csv("input/day.of.bb.DL.csv", header=TRUE, na.strings=c("","NA"))

dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")

# mergeing the my data with DF
pheno <- dl.wchill

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
pheno.lat <- pheno[,c("latbb1", "latbb50", "force.n", "photo.n", "site.n", "species", "lab2", "lab3","Utah_Model","Chill_portions")]

#pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 

pheno.lat$force.z2 <- (pheno.lat$force.n-mean(pheno.lat$force.n,na.rm=TRUE))/(sd(pheno.lat$force.n,na.rm=TRUE)*2)
pheno.lat$photo.z2 <- (pheno.lat$photo.n-mean(pheno.lat$photo.n,na.rm=TRUE))/(sd(pheno.lat$photo.n,na.rm=TRUE)*2)
pheno.lat$chillport.z2 <- (pheno.lat$Chill_portions-mean(pheno.lat$Chill_portions,na.rm=TRUE))/(sd(pheno.lat$Chill_portions,na.rm=TRUE)*2)
pheno.lat$site.z2 <- (pheno.lat$site.n-mean(pheno.lat$site.n,na.rm=TRUE))/(sd(pheno.lat$site.n,na.rm=TRUE)*2)

pheno.lat$species.fact <- as.numeric(as.factor(pheno.lat$species))
sort(unique(pheno.lat$species.fact)) # 19, 30 species, 47 with chill0 47 

nrow(pheno.lat) - nrow(pheno.lat) # 547 that had no terminal bb, 609
str(pheno.lat)

pheno.1lat <- pheno.lat[complete.cases(pheno.lat$latbb1), ]
pheno.50lat <- pheno.lat[complete.cases(pheno.lat$latbb50), ]

# now get the phylogeny and pair it with species names:
spInfo <- read.csv("input/species_list.csv")
pheno1Lat <- merge(pheno.1lat, spInfo, by = "species")
pheno50Lat <- merge(pheno.50lat, spInfo, by = "species")

westSp <- unique(pheno1Lat$species.name)

tree <- read.tree("input/SBphylo_phenobc.tre")
head(tree$tip.label)
length(tree$tip.label) #47
tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
tree$tip.label[tree$tip.label=="Nyssa_sylvatica"] <- "Alnus_viridis"

treeW <- keep.tip(tree, westSp)

phymatchW <- data.frame(species.name = treeW$tip.label, sppnum = c(1:length(treeW$tip.label)))

d1 <- merge(pheno1Lat, phymatchW, by="species.name")

d1 <- d1[order(d1$sppnum),]
nspecies1 <- max(d1$sppnum)
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


datalist <- with(pheno1Lat,
                 list( N=nrow(pheno1Lat),
                       n_sp = length(unique(pheno1Lat$species.fact)),
                       n_site = length(unique(pheno1Lat$site.n)),
                       lday = latbb1,
                       sp = species.fact,
                       chill1 = chillport.z2,
                       photo = photo.z2,
                       warm = force.z2,
                       site = site.z2,
                       Vphy = vcv_tree
                 ))

mdlLat1 <- stan("stan/bc_bb_2sites_standardized_phylogeny.stan",
                  data = datalist,
                  iter = 2000, warmup = 1000, chains=4
                  #, control = list(adapt_delta = 0.99)
)
save(mdlLat1, file="output/dl_phylo_lat1.Rda")


##########################################################
d50 <- merge(pheno50Lat, phymatchW, by="species.name")

d50 <- d50[order(d50$sppnum),]
nspecies50 <- max(d50$sppnum)
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


datalist <- with(pheno50Lat,
                 list( N=nrow(pheno50Lat),
                       n_sp = length(unique(pheno50Lat$species.fact)),
                       n_site = length(unique(pheno50Lat$site.n)),
                       lday = latbb50,
                       sp = species.fact,
                       chill1 = chillport.z2,
                       photo = photo.z2,
                       warm = force.z2,
                       site = site.z2,
                       Vphy = vcv_tree
                 ))

mdlLat50 <- stan("stan/bc_bb_2sites_standardized_phylogeny.stan",
                data = datalist,
                iter = 2000, warmup = 1000, chains=4
                #, control = list(adapt_delta = 0.99)
)
save(mdlLat50, file="output/dl_phylo_lat50.Rda")



## Now compare how the 1% vs the 50% does:

load("output/final/dl_phylo_lat1.Rda")
sum1 <- summary(mdl.1l)$summary 

range(sum1[, "n_eff"])
range(sum1[, "Rhat"])

# sum1Out <- summary(mdl.1l)$summary[c("mu_grand", "b_z","lam_interceptsa", "lam_interceptsb","sigma_interceptsa", "sigma_interceptsb", "sigma_y"),"mean"]

col4fig <- c("mean","sd","25%","50%","75%","n_eff","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")

mu_params_1l <- c(
  #"a_z","lam_interceptsa",
  "mu_b_warm", "mu_b_photo", "mu_b_chill1", "b_site", "mu_b_inter_wp",
  "mu_b_inter_wc1","mu_b_inter_pc1", "mu_b_inter_ws",
  "mu_b_inter_ps","mu_b_inter_sc1")

meanz1l <- sum1[mu_params_1l, col4fig]

rownames(meanz1l) = c( #"Root trait intercept","Lambda",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Manning Park",
  "Forcing x photoperiod",
  "Forcing x chilling",
  "Photoperiod x chilling",
  "Forcing x Manning Park",
  "Photoperiod x Manning Park",
  "Chilling x Manning Park"
  
)


range(sum1[, "n_eff"])
range(sum1[, "Rhat"])

# sum1Out <- summary(mdl.1l)$summary[c("mu_grand", "b_z","lam_interceptsa", "lam_interceptsb","sigma_interceptsa", "sigma_interceptsb", "sigma_y"),"mean"]
# 

#################################################
load("..//output/final/dl_phylo_lat50.Rda")

sum50 <- summary(mdlLat50)$summary

col4fig <- c("mean","sd","25%","50%","75%","n_eff","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","n_eff","Rhat")

mu_params_50l <- c(
  #"a_z","lam_interceptsa",
  "mu_b_warm", "mu_b_photo", "mu_b_chill1", "b_site", "mu_b_inter_wp",
  "mu_b_inter_wc1","mu_b_inter_pc1", "mu_b_inter_ws",
  "mu_b_inter_ps","mu_b_inter_sc1")

meanz50l <- sum50[mu_params_50l, col4fig]

rownames(meanz50l) = c( #"Root trait intercept","Lambda",
  "Forcing",
  "Photoperiod",
  "Chilling",
  "Manning Park",
  "Forcing x photoperiod",
  "Forcing x chilling",
  "Photoperiod x chilling",
  "Forcing x Manning Park",
  "Photoperiod x Manning Park",
  "Chilling x Manning Park"
  
)

pdf("figures/lat1_lat50_muplot.pdf", width = 7, height = 5)
par(mfrow = c(1,1), mar = c(5, 10, 2, 2))
# Upper panel: bud burst
plot(seq(-15, 
         15,
         length.out = nrow(meanz50l)), 
     1:nrow(meanz50l),
     type = "n",
     xlab = "Estimated change in budburst day",
     ylab = "",
     yaxt = "n"
     )

#legend(x = -20, y = 2, bty="n", legend = "a. Budburst", text.font = 2)
#rasterImage(bbpng, -20, 1, -16, 4)

axis(2, at = nrow(meanz50l):1, labels = rownames(meanz50l), las = 1, cex.axis = 0.8)
points(meanz50l[, 'mean'],
       nrow(meanz50l):1,
       pch = 16,
       col = "#eb8055ff",
       cex = 1.5)
arrows(meanz50l[, "75%"], nrow(meanz50l):1, meanz50l[, "25%"], nrow(meanz50l):1,
       len = 0, col = "#eb8055ff")
abline(v = 0, lty = 3)
########################################
points(meanz1l[, 'mean'],
       nrow(meanz1l):1,
       pch = 16,
       col = "cyan4",
       cex = 1.5)
arrows(meanz1l[, "75%"], nrow(meanz1l):1, meanz1l[, "25%"], nrow(meanz1l):1,
       len = 0, col = "cyan4")
abline(v = 0, lty = 3)

legend("topright",legend = c( "1st lateral budurst", "50% lateral budburst"),
      col = c("black", "black", "black", "black","black", "black","black"),
       pt.bg = c("#eb8055ff", "cyan4"),
       inset = 0.02, pch = c(21, 21), cex = 1, bty = "n")
dev.off()
# add advance/delay arrows
par(xpd=NA)
# arrows(1, 15.5, 6, 15.5, len = 0.1, col = "black")
# legend(5, 16.5, legend = "delay", bty = "n", text.font = 1, cex = 0.75)
# arrows(-1, 15.5, -6, 15.5, len = 0.1, col = "black")
# legend(-12, 16.5, legend = "advance", bty = "n", text.font = 1, cex = 0.75)
# legend(-2, 16.5, legend = "0", bty = "n", text.font = 1, cex = 0.75)
# par(xpd = FALSE)