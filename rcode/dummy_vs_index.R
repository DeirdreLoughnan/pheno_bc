## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

# Adding standardization: following the methods outlined by Andrew Gelman in his paper

# To properly add site to the model I am going to use the indexing approach used on pg 153 of statistical rethinking, modeling code was drafted by Lizzie

#library(scales)
library(arm)
library(rstan)
library(shinystan)
#library(reshape2)
library(bayesplot)
library(ggplot2)
#library(RColorBrewer)
library(dplyr)
library(plyr)

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
pheno <- rbind.fill(dl.wchill, df.wchill)
#pheno <- dl

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
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.integer(pheno$site.n)

pheno$Chill_portions <- as.numeric(pheno$Chill_portions)


pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)
pheno$utah.z2 <- (pheno$Utah_Model-mean(pheno$Utah_Model,na.rm=TRUE))/(sd(pheno$Utah_Model,na.rm=TRUE)*2)

#going to split it into analyses of terminal bb and lateral bb
# Starting with the terminal buds:
pheno.term <- pheno[,c("tbb", "force.n", "photo.n", "site.n", "species", "lab2","Utah_Model","Chill_portions","force.z2", "photo.z2", "chillport.z2", "utah.z2")]

pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 

nrow(pheno.term) - nrow(pheno.t) # That had no terminal bb, 609

#creating dummy var
pheno.t <- pheno.t %>%
    mutate ( d2 = if_else(site.n == 2, 1, 0),
             d3 = if_else(site.n == 3, 1, 0),
             d4 = if_else(site.n == 4, 1, 0))
head(pheno.t)
str(pheno.t$photo.z2)

datalist_z <- list( N=nrow(pheno.t),
                    Nsite = nrow(pheno.t),
                    n_sp = length(unique(pheno.t$species.fact)),
                    n_site = length(unique(pheno.t$site.n)),
                    bb = pheno.t$tbb,
                    sp = pheno.t$species.fact,
                    chill = pheno.t$chillport.z2,
                    photo = pheno.t$photo.z2,
                    force = pheno.t$force.z2,
                    site = pheno.t$site.n,
                    site2 = pheno.t$d2,
                    site3 = pheno.t$d3,
                    site4 = pheno.t$d4)

save(mdl.t, file="output/tbb_cport_stnd_index.Rds")


mdl.i <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_index.stan",
              data = datalist_z,
              iter = 4000, chains=4)
save(mdl.i, file="output/tbb_cport_stnd_index.Rds")

mdl.d <- stan("stan/bc_bb_ncpphoto_ncpinter_standardize_dummy.stan",
              data = datalist_z,
              iter = 4000, chains=1)
save(mdl.d, file="output/tbb_cport_stnd_dummy.Rds")
#######################################################################

load("output/tbb_cport_stnd_dummy.Rds")

ssm <-  as.shinystan(mdl.d)
launch_shinystan(ssm)

#dummy variables:
sumd <- summary(mdl.d)$summary
bforce <- sumd[grep("b_force", rownames(sumd)), "mean"]; bforce
site <- sumd[grep("b_site", rownames(sumd)), "mean"]; site

post <- rstan::extract(mdl.d)

y<-as.numeric(pheno.t$tbb)
yrep<- post$y_hat 
ppc_dens_overlay(y, yrep[1:50, ])

range(sumd[, "n_eff"])
range(sumd[, "Rhat"])

#pairs(mdl.t, pars = c("mu_a","mu_force","mu_chill","mu_photo", "mu_site","mu_inter_fp", "mu_inter_fc"))

########################################################

post3 <- as.matrix(mdl.d, par = c("mu_force", "mu_chill","mu_photo", "mu_inter_fp", "mu_inter_fc"))

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")

mcmc_areas(post3,
           pars = c("mu_force", "mu_chill", "mu_photo"),
           prob = 0.8) + plot_title

col4fig <- c("mean","sd","25%","50%","75%","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")

# manually to get right order
mu_params <- c("mu_grand",
               "mu_a",
               "a_site2",
               "a_site3",
               "a_site4",
               "mu_force",
               "mu_photo",
               "mu_chill",
               "mu_inter_fp",
               "mu_inter_fc",
               "mu_inter_pc")

mu_params <-   c("mu_grand",
                 "mu_a",
                 "a_site2",
                  "a_site3",
                  "a_site4",
                  "mu_force",
                  "mu_photo",
                  "mu_chill",
                  "mu_inter_fp",
                  "mu_inter_fc",
                  "mu_inter_pc")
meanzt <- sumd[mu_params, col4fig]

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

meanzt.table <- sumd[mu_params, col4table]
row.names(meanzt.table) <- row.names(meanzt)
meanzt.table
#write.table(meanzt.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)

################################################################################
load("output/tbb_cport_stnd_index.Rds")

#dummy variables:
sumi <- summary(mdl.i)$summary
bforce <- sumi[grep("b_force", rownames(sumi)), "mean"]; bforce
site <- sumi[grep("a_site", rownames(sumi)), "mean"]; site

post <- rstan::extract(mdl.i)

y<-as.numeric(pheno.t$tbb)
yrep<- post$y_hat 
ppc_dens_overlay(y, yrep[1:50, ])

range(sumi[, "n_eff"])
range(sumi[, "Rhat"])

#pairs(mdl.t, pars = c("mu_a","mu_force","mu_chill","mu_photo", "mu_site","mu_inter_fp", "mu_inter_fc"))

########################################################

post3 <- as.matrix(mdl.i, par = c("mu_force", "mu_chill","mu_photo", "mu_inter_fp", "mu_inter_fc"))

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")

mcmc_areas(post3,
           pars = c("mu_force", "mu_chill", "mu_photo"),
           prob = 0.8) + plot_title

####################################################################3

col4fig <- c("mean","sd","25%","50%","75%","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")

# manually to get right order
mu_params <- c("mu_grand",
               "mu_a",
              "a_site[1]",
               "a_site[2]",
               "a_site[3]",
               "a_site[4]",
               "mu_force",
               "mu_photo",
               "mu_chill",
               "mu_site",
               "mu_inter_fp",
               "mu_inter_fc",
               "mu_inter_pc")

mu_params <-   c("mu_grand",
                 "mu_a",
                 "a_site[1]",
                 "a_site[2]",
                 "a_site[3]",
                 "a_site[4]",
                 "mu_force",
                 "mu_photo",
                 "mu_chill",
                 "mu_inter_fp",
                 "mu_inter_fc",
                 "mu_inter_pc")
meanzt <- sumi[mu_params, col4fig]

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

meanzt.table <- sumi[mu_params, col4table]
row.names(meanzt.table) <- row.names(meanzt)
meanzt.table
#write.table(meanzt.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)

