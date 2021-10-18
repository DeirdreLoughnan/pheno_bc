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
#write.csv(pheno, "input/pheno.wchill.midge.csv")

############################################################

pheno$force.n <- pheno$force
pheno$force.n[pheno$force.n == "HF"] <- "1"
pheno$force.n[pheno$force.n == "LF"] <- "0"
pheno$force.n <- as.numeric(pheno$force.n)

pheno$photo.n <- pheno$photo
pheno$photo.n[pheno$photo.n == "HP"] <- "1"
pheno$photo.n[pheno$photo.n == "LP"] <- "0"
pheno$photo.n <- as.numeric(pheno$photo.n)

# pg 153 of statistical rethinking: we have multiple categories, so we want to use index values as integers
#pheno$population <- as.factor(pheno$population)
pheno$site.n <- pheno$population
# pheno$site.n <- pheno$population
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "2"
pheno$site.n[pheno$site.n == "HF"] <- "3"
pheno$site.n[pheno$site.n == "SH"] <- "4"
pheno$site.n <- as.integer(pheno$site.n)

pheno$Chill_portions <- as.numeric(pheno$Chill_portions)

# pheno$force.z <- (pheno$force.t-mean(pheno$force.t,na.rm=TRUE))/sd(pheno$force.t,na.rm=TRUE)
# pheno$photo.z <- (pheno$photo.t-mean(pheno$photo.t,na.rm=TRUE))/sd(pheno$photo.t,na.rm=TRUE)
# pheno$chillport.z <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/sd(pheno$Chill_portions,na.rm=TRUE)

pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)
pheno$utah.z2 <- (pheno$Utah_Model-mean(pheno$Utah_Model,na.rm=TRUE))/(sd(pheno$Utah_Model,na.rm=TRUE)*2)

#z-score site as well
#pheno$site.z2 <- (pheno$site.n-mean(pheno$site.n,na.rm=TRUE))/(sd(pheno$site.n,na.rm=TRUE)*2)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
#pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.term <- pheno[,c("tbb", "force.n", "photo.n", "site.n", "species", "lab2","Utah_Model","Chill_portions","force.z2", "photo.z2", "chillport.z2", "utah.z2")]

pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 47 

nrow(pheno.term) - nrow(pheno.t) # That had no terminal bb, 609

# creating dummy var
# pheno.t <- pheno.t %>% 
#     mutate ( d2 = if_else(site.n == 2, 1, 0),
#              d3 = if_else(site.n == 3, 1, 0),
#              d4 = if_else(site.n == 4, 1, 0))
# head(pheno.t)
str(pheno.t$photo.z2)

datalist_z <- list( N=nrow(pheno.t),
                         n_sp = length(unique(pheno.t$species.fact)),
                         n_site = length(unique(pheno.t$site.n)),
                         bb = pheno.t$tbb,
                         sp = pheno.t$species.fact,
                         chill = pheno.t$chillport.z2,
                         photo = pheno.t$photo.z2,
                         force = pheno.t$force.z2,
                         site = pheno.t$site.n)
str(datalist_z)
datalist_z$photo


datalist.zutah <- with(pheno.t,
                   list( N=nrow(pheno.t),
                         n_sp = length(unique(pheno.t$species.fact)),
                         n_site = length(unique(pheno.t$site.n)),
                         bb = tbb,
                         sp = species.fact,
                         chill = utah.z2,
                         photo = photo.z2,
                         force = force.z2
                        # site = site.z2
                   ))

datalist_zutah_index <- with(pheno.t,
                       list( N=nrow(pheno.t),
                             n_sp = length(unique(pheno.t$species.fact)),
                             n_site = length(unique(pheno.t$site.n)),
                             bb = tbb,
                             sp = species.fact,
                             chill = utah.z2,
                             photo = photo.z2,
                             force = force.z2,
                             site = pheno.t$site.n
                       ))
N=nrow(pheno.t)
n_sp = length(unique(pheno.t$species.fact))
n_site = length(unique(pheno.t$site.n))
bb = pheno.t$tbb
sp = pheno.t$species.fact
chill = pheno.t$utah.z2
photo = pheno.t$photo.z2
force = pheno.t$force.z2
site = pheno.t$site.n
data=c("N","n_sp","n_site","bb","sp","chill","photo","force","site")

datalist_zutah_index$force.z2
head(pheno.t$force.z2)
mdl.t <- stan("stan/bc_bb_ncp_standardize_oct18.stan",
               data=c("N","n_sp","n_site","bb","sp","chill","photo","force","site"),
              iter = 4000, chains=1)#, control = list(adapt_delta = 0.99))

mdl.i <- stan("stan/bc.bb.ncpphoto.ncpinter.standardize.old.stan",
              data = datalist_z,
              iter = 4000, chains=1)
              #, control = list(adapt_delta = 0.99))
datalist_z$N
length(datalist_z$site)
save(mdl.t, file="output/tbb_utah_cport_stnd_index.Rds")

# no div trans or any warnings of any kind!
head(pheno.t)
temp <- subset(pheno.t, site.n != 1)
dl.t <- subset(temp, site.n != 1)

datalist.z <- with(dl.t,
                   list( N=nrow(dl.t),
                         n_sp = length(unique(dl.t$species.fact)),
                         n_site = length(unique(dl.t$site.n)),
                         bb = tbb,
                         sp = as.numeric(as.factor(dl.t$species)),
                         chill = chillport.z2,
                         photo = photo.z2,
                         force = force.z2,
                         site = site.z2
                   ))

mdl.t <- stan("stan/bc.bb.ncpphoto.ncpinter.standardize.stan",
              data = datalist.z,
              iter = 4000, chains=4, control = list(adapt_delta = 0.99))
#######################################################################

load("output/tbb_ncp_cport_stnd_index_DL.Rds")
load("output/tbb_utah_cport_stnd_index.Rds")

ssm <-  as.shinystan(mdl.t)
launch_shinystan(ssm)

sumt <- summary(mdl.t)$summary
bforce <- sumt[grep("b_force", rownames(sumt)), "mean"]; bforce

post <- rstan::extract(mdl.t)

y<-as.numeric(pheno.t$tbb)
yrep<- post$y_hat 
ppc_dens_overlay(y, yrep[1:50, ])

range(sumt[, "n_eff"])
range(sumt[, "Rhat"])

#pairs(mdl.t, pars = c("mu_a","mu_force","mu_chill","mu_photo", "mu_site","mu_inter_fp", "mu_inter_fc"))

########################################################

post3 <- as.matrix(mdl.t, par = c("mu_force", "mu_chill","mu_photo", "mu_inter_fp", "mu_inter_fc"))

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")

mcmc_areas(post3,
           pars = c("mu_force", "mu_chill", "mu_photo"),
           prob = 0.8) + plot_title

##############################################################################
# histograms
par(mfrow = c(1,1))
# mu_force
h1 <- hist(rnorm(1000, 0, 50))
h2 <- hist(post$mu_force)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-150,150))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# mu_chill
h1 <- hist(rnorm(1000, 0, 35))
h2 <- hist(post$mu_chill)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# mu_photo
h1 <- hist(rnorm(1000, 0, 35))
h2 <- hist(post$mu_photo)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# mu_site
h1 <- hist(rnorm(1000, 1, 35))
h2 <- hist(post$mu_site)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# mu_inter_fp
h1 <- hist(rnorm(1000, 0, 35))
h2 <- hist(post$mu_inter_fp)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_force
h1 <- hist(rnorm(1000, 0, 10))
h2 <- hist(post$sigma_force)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-50,50))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_chill
h1 <- hist(rnorm(1000, 0, 40))
h2 <- hist(post$sigma_chill)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_photo
h1 <- hist(rnorm(1000, 0, 10))
h2 <- hist(post$sigma_photo)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-50,50))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_site
h1 <- hist(rnorm(1000, 0, 40))
h2 <- hist(post$sigma_site)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_interaction
h1 <- hist(rnorm(1000, 1, 10))
h2 <- hist(post$sigma_b_inter_pc, xlim =c(-10,10))
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

####################################################################3
# How to interpret the z.score values:
mu_force_prcnt <- round((sumt[grep("mu_force", rownames(sumt)), "mean"]/4)/(sd(pheno$force.n)*2)*100*2, digits = 2) 
mu_photo_prcnt <- round((sumt[grep("mu_photo", rownames(sumt)), "mean"]/4)/(sd(pheno$photo.n)*2)*100*2, digits = 2) 
mu_chill_prcnt <- round((sumt[grep("mu_chill", rownames(sumt)), "mean"]/4)/(sd(pheno$Chill_portions)*2)*100*2) 

mu_force_bt <- mean(pheno$force.n) + ((sumt[grep("mu_force", rownames(sumt)), "mean"])*(sd(pheno$force.n)*2))
mu_photo_bt <- mean(pheno$photo.n) + ((sumt[grep("mu_photo", rownames(sumt)), "mean"])*(sd(pheno$photo.n)*2))
mu_chillport_bt <- mean(pheno$Chill_portions) + ((sumt[grep("mu_chill", rownames(sumt)), "mean"])*(sd(pheno$Chill_portions)*2))
mu_site_bt <- mean(pheno$site.n) + ((sumt[grep("mu_site", rownames(sumt)), "mean"])*(sd(pheno$site.n)*2))

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

meanzt.table <- sumt[mu_params, col4table]
row.names(meanzt.table) <- row.names(meanzt)
meanzt.table
#write.table(meanzt.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)

#####################################################################
sumt_photo <- sumt[grep("b_photo", rownames(sumt)), "mean"]

plot(sumt[grep("b_force", rownames(sumt)), "mean"], sumt[grep("b_chill", rownames(sumt)), "mean"])
abline(a = 0, b = 1, lty = "dotted")

plot(sumt[grep("b_force", rownames(sumt)), "mean"], sumt_photo[39:76])
abline(a = 0, b = 1, lty = "dotted")

plot(sumt[grep("b_chill", rownames(sumt)), "mean"], sumt_photo[39:76])
abline(a = 0, b = 1, lty = "dotted")

plot(sumt[grep("b_force", rownames(sumt)), "mean"], sumt[grep("b_chill", rownames(sumt)), "mean"])
abline(a = 0, b = 1, lty = "dotted")

plot(sumt_photo[39:76], sumt)
abline(a = 0, b = 1, lty = "dotted")

plot(sumt[grep("b_chill", rownames(sumt)), "mean"], sumt_photo[39:76])
abline(a = 0, b = 1, lty = "dotted")

boxplot(tbb ~ site.n, pheno.t)
plot(tbb ~ Chill_portions, pheno.t)

plot(tbb ~ force.n, pheno.t)
#####################################################################

# now running the same model for the lateral buds
pheno.1lat <- pheno[, c("latbb1", "Chill_portions", "force.n", "photo.n", "site.n","species", "force.z2", "photo.z2", "chillport.z2")]
pheno.1l <- pheno.1lat[complete.cases(pheno.1lat), ]
nrow(pheno.1lat) - nrow(pheno.1l)  

pheno.1l$species.fact <- as.numeric(as.factor(pheno.1l$species))
sort(unique(pheno.1l$species.fact))

# lat.data <- with(pheno.1l,
#                  list( N = nrow(pheno.1l),
#                        n_sp = length(unique(pheno.1l$species.fact)),
#                        # n_site = length(unique(pheno.1l$site.n)),
#                        bb = latbb1,
#                        sp = species.fact,
#                        chill = as.numeric(Chill_portions),
#                        photo = photo.n,
#                        force = force.n
#                        # site = site.n
#                  ))

lat.data.z <- with(pheno.1l,
                 list( N = nrow(pheno.1l),
                       n_sp = length(unique(pheno.1l$species.fact)),
                       n_site = length(unique(pheno.1l$site.n)),
                       bb = latbb1,
                       sp = species.fact,
                       chill = chillport.z2,
                       photo = photo.z2,
                       force = force.z2,
                       site = site.n
                 ))

# mdl <- stan("stan/bc.bb.inter.stan",
#             data= datalist
#             ,iter=2000, chains=4)
#gives 200 divergent transitions, 41 transitions that exceed max tree depth, chains were not mixed, with low ESS

# mdl.1l <- stan("stan/bc.bb.ncpphoto.ncpinter.newpriors.nosite.stan",
#                 data= lat.data,
#                 iter=4000, chains=4)

mdl.1l <- stan("stan/bc.bb.ncpphoto.ncpinter.standardize.stan",
               data= lat.data.z,
               iter=4000, chains=4)


#save(mdl.1l, file="output/tbb.photo.winter.ncp.lateralbud.nosite.zscore.Rds")
load("output/lateralbud_stbd_Aug31.Rds")

sum1l <- summary(mdl.1l)$summary
sum1l[grep("mu_",rownames(sum1l)), ]

ssm <- as.shinystan(mdl.1l)
launch_shinystan(ssm)

post1l <- rstan::extract(mdl.1l)

y<-as.numeric(pheno.1l$latbb1)
yrep<- post1l$y_hat 
ppc_dens_overlay(y, yrep[1:50, ])

range(sum1l[, "n_eff"])
range(sum1l[, "Rhat"])

pairs(mdl.1l, pars = c("mu_a","mu_force","mu_chill","mu_photo", "mu_site","mu_inter_fp", "mu_inter_fc"))

##############################################################################
# histograms
par(mfrow = c(1,1))
# mu_force
h1 <- hist(rnorm(1000, 0, 50))
h2 <- hist(post1l$mu_force)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-150,150))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# mu_chill
h1 <- hist(rnorm(1000, 0, 35))
h2 <- hist(post1l$mu_chill)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# mu_photo
h1 <- hist(rnorm(1000, 0, 35))
h2 <- hist(post1l$mu_photo)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# mu_site
h1 <- hist(rnorm(1000, 1, 35))
h2 <- hist(post1l$mu_site)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# mu_inter_fp
h1 <- hist(rnorm(1000, 0, 35))
h2 <- hist(post1l$mu_inter_fp)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_force
h1 <- hist(rnorm(1000, 0, 10))
h2 <- hist(post1l$sigma_force)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-50,50))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_chill
h1 <- hist(rnorm(1000, 0, 40))
h2 <- hist(post1l$sigma_chill)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_photo
h1 <- hist(rnorm(1000, 0, 10))
h2 <- hist(post1l$sigma_photo)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-50,50))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_site
h1 <- hist(rnorm(1000, 0, 40))
h2 <- hist(post1l$sigma_site)
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)

# sigma_interaction
h1 <- hist(rnorm(1000, 1, 10))
h2 <- hist(post1l$sigma_b_inter_pc, xlim =c(-10,10))
plot(h2, col=rgb(0,0,1,1/4), xlim =c(-100,100))
plot(h1, col=rgb(1,0,1,1/4), add = TRUE)
####################################################################
# PPC 

# mdl.slopes <- as.data.frame(sm.sum[grep("b", rownames(sm.sum)), c(1,6)]) 
# mdl.int <- as.data.frame(sm.sum[grep("a", rownames(sm.sum)), ]) 
# mdl.a <- mdl.int[, 1]
# mdl.b <- mdl.slopes[, 1]
# 
# ggplot() +
#   geom_point(data = mdl.slopes, aes(y = row.names(mdl.slopes), x = mean), color = "darkgreen") +
#   labs( x = "doy", y = "Species")
# 
# #####################################################################
# 
# #plot(mdl, pars="a", ci_level=0.5, outer_level=0.5,col="blue")
# # PPC based on the vingette from https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html 
# 
# ext <- rstan::extract(mdl)
# 
# # not the most normal 
# hist(ext$a)
# 
# y <- pheno.t$tbb
# 
# y.ext <- ext$ypred_new # I want this to be a matrix, which it is, with one element for each data point in y
# 
# ppc_dens_overlay(y, y.ext[1:50, ])
# 
# mean(ext$b_force)
# mean(ext$b_chill)
# mean(ext$b_photo)

######################################################
# plotting code taken from buds-master Pheno Budburst analysis.R
# load("output/tbb_ncp_termianlbud_dldf.Rds")
# load("output/lat1_ncp_dldf.Rds")

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
meanzl <- sum1l[mu_params, col4fig]

rownames(meanzt) = c("Forcing",
                     "Photoperiod",
                     "Chilling",
                     "Site",
                     "Forcing x Photoperiod",
                     "Forcing x Chilling",
                     "Photoperiod x Chilling",
                     "Forcing x Site",
                     "Photoperiod x Site",
                     "Chilling x Site"
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
                     "Chilling x Site"
  )

meanzt.table <- sumt[mu_params, col4table]
row.names(meanzt.table) <- row.names(meanzt)
meanzt.table
#write.table(meanzt.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)

meanzl.table <- sum1l[mu_params, col4table]
row.names(meanzl.table) <- row.names(meanzl)
meanzl.table
#write.table(meanzl.table , "output/lat.mdl.esti.csv", sep = ",", row.names = FALSE)

# Begin by checking to see what cue is most important and whether there are strong correlations between cues:
df.mean.t <- data.frame(bb.force = sumt[grep("b_force", rownames(sumt)), 1],
                          bb.photo = sumt[grep("b_photo_ncp", rownames(sumt)), 1],
                          bb.chill = sumt[grep("b_chill", rownames(sumt)), 1])

df.mean.l <- data.frame(lat.force = sum1l[grep("b_force", rownames(sumt)), 1],
                        lat.photo = sum1l[grep("b_photo_ncp", rownames(sumt)), 1],
                        lat.chill = sum1l[grep("b_chill", rownames(sumt)), 1])

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

pdf(file.path( "figures/changes.pheno.standardized.pdf"), width = 7, height = 8)
par(mfrow = c(2,1), mar = c(5, 10, 2, 1))
# Upper panel: bud burst
plot(seq(-15, 
         15,
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
plot(seq(-15, 
         15, 
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

## Replicating Flynn Figure 2:

b.force.both <- sumt[grep("b_force", rownames(sumt))]
b.photo.both <- sumt[grep("b_photo", rownames(sumt))]; b.photo.both <- b.photo.both[48:94]
b.chill.both <- sumt[grep("b_chill", rownames(sumt))]

shrubs.both = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG", "alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
trees.both = c("ACEPEN", "ACERUB", "ACESAC", "ALNINC", "BETALL", "BETLEN", "BETPAP", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL", "acegla","betpap", "poptre", "popbal")

pheno.term <- pheno[,c("tbb", "chillport.z2", "force.z2", "photo.z2", "species", "lab2","population")]
pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 


species.both <- sort(unique(pheno.t$species))
species.fact.both <-as.numeric( as.factor(unique(pheno.t$species)))
type.both <- c("tree", "tree", "tree","tree", "shrub", "shrub", "shrub", "tree", "tree", "tree", "tree", "tree",
               "tree", "shrub","shrub","tree","tree", "shrub", "shrub","shrub",  "shrub", "shrub", "shrub", "shrub", "shrub",
               "tree","tree","tree","tree","tree","tree","tree","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub")
both <- data.frame(species.both, species.fact.both, b.force.both, b.photo.both, b.chill.both, type.both)


#pdf(file.path( "figures/chill_vs_force_dldf.pdf"), width = 7, height = 8)
cf.both <- ggplot(both, aes(x= b.chill.both, y = b.force.both, col = type.both)) +
  geom_point()+
  # ylim (-25, 1) +
  # xlim (-30, 0) +
  labs( y = "High forcing", x = "High chilling") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()

#pdf(file.path( "figures/chill_vs_photo_dldf.pdf"), width = 7, height = 8)
cp.both <- ggplot(both, aes(x= b.chill.both, y = b.photo.both, col = type.both)) +
  geom_point() +
  labs (x = "High chilling", y = "Long photoperiod") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#dev.off()
#legend.position = "none"

#pdf(file.path( "figures/force_vs_photo_dldf.pdf"), width = 7, height = 8)
fp.both <- ggplot(both, aes(x= b.force.both, y = b.photo.both, col = type.both)) +
  geom_point() +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1)+
  labs(x = "High forcing", y = "Long photoperiod") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()

## Plotting the day to bb with the cues on the y-axis 
term.bb.both <- ddply(pheno, c("species"), summarize, mean = mean(tbb, na.rm = TRUE), mean.lat = mean(latbb50, na.rm = TRUE))
names(term.bb.both) <- c("species.both", "mean","mean.lat")
term.both <- merge(term.bb.both, both, by = "species.both", all =TRUE)
term.both <- term.both[,c("species.both","mean","b.force.both","b.chill.both","b.photo.both")]
term.both <- term.both[complete.cases(term.both), ] 

tf.both <- ggplot(term.both, aes(y = b.force.both, x= mean,col = type.both)) +
  geom_point() +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  labs(x = "Mean day of budburst", y = "High forcing") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

tc.both <- ggplot(term.both, aes(y = b.chill.both, x= mean,col = type.both)) +
  geom_point() +
  labs(x = "Mean day of budburst", y = "High chilling") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

tp.both <- ggplot(term.both, aes(y = b.photo.both, x= mean,col = type.both)) +
  geom_point() +
  labs(x = "Mean day of budburst", y = "Long photoperiod")+
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##### General boxplots across treatments:
# As per Lizzie's July 7 post: I should look at how linear these relationships are
# Plot raw data (bb~ chill) with the chill effect from the model plotted on top

plot(pheno.t$tbb ~ pheno.t$chillport.z2, pch = 19, col = "darkgreen")
abline(a = sumt[grep("mu_a", rownames(sumt)), "mean"], b = sumt[grep("mu_chill", rownames(sumt)), "mean"])

plot(sumt[grep("sigma_chill", rownames(sumt)), "mean"])
plot(sumt[grep("log", rownames(sumt)), "mean"])

plot(pheno.t$tbb , (sumt[grep("b_chill", rownames(sumt)), "mean"]))
abline(a = 0, b = 1, lty = "dotted")

##############################################################################
# plots of the chill portions
# Do we see large differences across sites, or could it just be east/west
head(pheno.t)

plot
