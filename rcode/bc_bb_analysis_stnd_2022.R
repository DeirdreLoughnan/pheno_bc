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

dl <- read.csv("input/day.of.bb.DL.csv", header=TRUE, na.strings=c("","NA"))
head(dl)
dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dl, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")

# mergeing the my data with DF
pheno <- dl.wchill
#pheno <- dl

pheno$first <- ifelse(pheno$tbb < pheno$latbb1,"t", ifelse (pheno$tbb == pheno$latbb1,"tl", "l"))

head(pheno)
table(pheno$species, pheno$first)
#write.csv(pheno, "input/pheno.wchill.midge.csv")
# 
# combined the data has 3197 unique samples
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
pheno.term <- pheno[,c("tbb", "force.n", "photo.n", "site.n", "species", "lab2","Utah_Model","Chill_portions")]

pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 

pheno.t$force.z2 <- (pheno.t$force.n-mean(pheno.t$force.n,na.rm=TRUE))/(sd(pheno.t$force.n,na.rm=TRUE)*2)
pheno.t$photo.z2 <- (pheno.t$photo.n-mean(pheno.t$photo.n,na.rm=TRUE))/(sd(pheno.t$photo.n,na.rm=TRUE)*2)
pheno.t$chillport.z2 <- (pheno.t$Chill_portions-mean(pheno.t$Chill_portions,na.rm=TRUE))/(sd(pheno.t$Chill_portions,na.rm=TRUE)*2)
pheno.t$site.z2 <- (pheno.t$site.n-mean(pheno.t$site.n,na.rm=TRUE))/(sd(pheno.t$site.n,na.rm=TRUE)*2)

pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 19, 30 species, 47 with chill0 47 

nrow(pheno.term) - nrow(pheno.t) # 547 that had no terminal bb, 609
str(pheno.t)
#pheno.term$Chill_portions <- as.factor(pheno.term$Chill_portions)

datalist <- with(pheno.t,
                 list( N=nrow(pheno.t),
                       n_sp = length(unique(pheno.t$species.fact)),
                       n_site = length(unique(pheno.t$site.n)),
                       lday = tbb,
                       sp = species.fact,
                       chill1 = chillport.z2,
                       photo = photo.z2,
                       warm = force.z2,
                       site = site.z2
                 ))


mdl.t <- stan("stan/lday_site_sp_chill_inter_poola_ncpwp_2chill_sitefixedstnd_standardized.stan",
              data = datalist,
              iter = 2000, warmup = 1000, chains=4
             # , control = list(adapt_delta = 0.99)
              )
save(mdl.t, file="output/tbb_ncpint_ncpwp_2chillport_2xstandardized_2stndsite.Rda")

#1. running the model with no ncp - 1417 divergent transitions after warmup, 4 transitions exceed max tree depth, Rhat up to 1.06- poor chain mixing, low ESS
# the sigmas for forcing and photo look smushed, but chilling looks fine, wp kinds smushed, ws badly, ps looks fine, wc and pc look kinda bad, sc looks fine

#2. Running the model with ncp for all the interactions - 426 divergent transitions after warmup; sigma_warm and sigma_photo both look bad
# tried simply increasing the warmup = 367 divergent transitions and the Rhat 1.07

#3. Running with ncp for forcing and photoperiod too, but with 2000:1000! No warning messages!! 

#4. Just to double check - ran model with only photo ncp (like I had before) - 61 divergent transitions

# 5. Also double checked that running the model with the standardization code in both the R code and the model or just in the R code matter, but it actually didn't
############################################################################
#########################################################################

# standardized in both the datalist and in the transformed parameter block
load("output/final/tbb_ncpint_ncpwp_2chillport_2xstandardized.Rda")

sumt <- summary(mdl.t)$summary
range(sumt[, "n_eff"])
range(sumt[, "Rhat"])
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")
# manually to get right order
mu_params <- c("mu_a",
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
               "sigma_a",
               "sigma_b_inter_wp",
               "sigma_b_inter_wc1",
               "sigma_b_inter_pc1",
               "sigma_b_inter_ws",
               "sigma_b_inter_ps",
               "sigma_b_inter_sc1",
               "sigma_y"
)

meanz_2stnd <- sumt[mu_params, col4table]
meanz_2stnd

#########################################################################

# ssm <-  as.shinystan(mdl.t)
# launch_shinystan(mdl.t)

#####################################################################
#####################################################################

# # now running the same model for the lateral buds
pheno.50lat <- pheno[, c("latbb50", "force.n", "photo.n", "site.n", "species", "lab2","Utah_Model","Chill_portions")]
pheno.50l <- pheno.50lat[complete.cases(pheno.50lat), ]
nrow(pheno.50lat) - nrow(pheno.50l)  # a lot of samples did not reach even 50%! 1084

pheno.50l$species.fact <- as.numeric(as.factor(pheno.50l$species))
sort(unique(pheno.50l$species.fact))

pheno.50l$force.z2 <- (pheno.50l$force.n-mean(pheno.50l$force.n,na.rm=TRUE))/(sd(pheno.50l$force.n,na.rm=TRUE)*2)
pheno.50l$photo.z2 <- (pheno.50l$photo.n-mean(pheno.50l$photo.n,na.rm=TRUE))/(sd(pheno.50l$photo.n,na.rm=TRUE)*2)
pheno.50l$chillport.z2 <- (pheno.50l$Chill_portions-mean(pheno.50l$Chill_portions,na.rm=TRUE))/(sd(pheno.50l$Chill_portions,na.rm=TRUE)*2)
pheno.50l$site.z2 <- (pheno.50l$site.n-mean(pheno.50l$site.n,na.rm=TRUE))/(sd(pheno.50l$site.n,na.rm=TRUE)*2)

datalist.50 <- with(pheno.50l,
               list( N = nrow(pheno.50l),
                     n_sp = length(unique(pheno.50l$species.fact)),
                     n_site = length(unique(pheno.50l$site.n)),
                     lday = latbb50,
                     sp = species.fact,
                     chill1 = chillport.z2,
                     photo = photo.z2,
                     warm = force.z2,
                     site = site.z2
               ))
mdl.lat50 <- stan("stan/lday_site_sp_chill_inter_poola_ncpwp_2chill_sitefixedstnd_standardized.stan", 
              data = datalist.50,
              iter = 4000, warmup = 3000, chains=4
              # , control = list(adapt_delta = 0.99)
)
save(mdl.lat50, file="output/l50_ncpint_ncpwp_2chillport_2xstandardized_stndsite.Rda")

# sum50l <- summary(mdl.50l)$summary
# sum50l[grep("mu_",rownames(sum50l)), ]
# sum50l
# ssm <- as.shinystan(mdl.50l)
# launch_shinystan(ssm)
# 
# #pairs(sm.sum, pars=c("mu_a","mu_force","mu_chill","mu_photo_ncp")) # this gives a lot of warning messages and not the figure i was hoping/expected
# 
# save(sum50l, file="output/tbb_photo_winter_ncp_lateralbud.Rds")

############################################################################################
# now running the same model for the lateral buds
pheno.1lat <- pheno[, c("latbb1", "force.n", "photo.n", "site.n", "species", "lab2","Utah_Model","Chill_portions")]
pheno.1l <- pheno.1lat[complete.cases(pheno.1lat), ]
nrow(pheno.1lat) - nrow(pheno.1l)

pheno.1l$species.fact <- as.numeric(as.factor(pheno.1l$species))
sort(unique(pheno.1l$species.fact))

pheno.1l$force.z2 <- (pheno.1l$force.n-mean(pheno.1l$force.n,na.rm=TRUE))/(sd(pheno.1l$force.n,na.rm=TRUE)*2)
pheno.1l$photo.z2 <- (pheno.1l$photo.n-mean(pheno.1l$photo.n,na.rm=TRUE))/(sd(pheno.1l$photo.n,na.rm=TRUE)*2)
pheno.1l$chillport.z2 <- (pheno.1l$Chill_portions-mean(pheno.1l$Chill_portions,na.rm=TRUE))/(sd(pheno.1l$Chill_portions,na.rm=TRUE)*2)
pheno.1l$site.z2 <- (pheno.1l$site.n-mean(pheno.1l$site.n,na.rm=TRUE))/(sd(pheno.1l$site.n,na.rm=TRUE)*2)

datalist.1 <- with(pheno.1l,
                 list( N = nrow(pheno.1l),
                       n_sp = length(unique(pheno.1l$species.fact)),
                       n_site = length(unique(pheno.1l$site.n)),
                       lday = latbb1,
                       sp = species.fact,
                       chill1 = chillport.z2,
                       photo = photo.z2,
                       warm = force.z2,
                       site = site.z2
                 ))

# mdl <- stan("stan/bc.bb.inter.stan",
#             data= datalist
#             ,iter=2000, chains=4)
#gives 200 divergent transitions, 41 transitions that exceed max tree depth, chains were not mixed, with low ESS

mdl.1l <- stan("stan/lday_site_sp_chill_inter_poola_ncpwp_2chill_sitefixedstnd_standardized.stan", 
                            data = datalist.1,
                            iter = 2000, warmup = 1000, chains=4
                            # , control = list(adapt_delta = 0.99)
)
save(mdl.1l, file="output/l1_ncpint_ncpwp_2chillport_2xstandardized_stndsite.Rda")

sum1l <- summary(mdl.1l)$summary
sum1l[grep("mu_",rownames(sum1l)), ]
sum1l
ssm <- as.shinystan(mdl.1l)
launch_shinystan(ssm)

############################################################################
#looking at the timing of first budburst regardless of terminal or lateral
dlall <- read.csv("input/dl_allbb.csv")

temp <- str_split_fixed(dlall$trt, "_", 3); head(temp)
dlall$chill<- temp[,1]
dlall$photo <- temp[,2]
dlall$force <- temp[,3]

dl.chill <- read.csv("input/chilling_values_Hope_Smithers.csv")

dl.wchill <- merge(dlall, dl.chill, by = c("population","chill"))
dl.wchill$lab3 <- dl.wchill$lab2
dl.wchill$lab2 <- paste(dl.wchill$species, dl.wchill$population, dl.wchill$rep, sep = "_")

# mergeing the my data with DF
pheno <- dl.wchill

head(pheno)

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
pheno$site.n <- as.numeric(pheno$site.n)

head(pheno)
#
# standardize the 0/1 and standardize sites? 
pheno$force.z2 <- (pheno$force.n-mean(pheno$force.n,na.rm=TRUE))/(sd(pheno$force.n,na.rm=TRUE)*2)
pheno$photo.z2 <- (pheno$photo.n-mean(pheno$photo.n,na.rm=TRUE))/(sd(pheno$photo.n,na.rm=TRUE)*2)
pheno$chillport.z2 <- (pheno$Chill_portions-mean(pheno$Chill_portions,na.rm=TRUE))/(sd(pheno$Chill_portions,na.rm=TRUE)*2)

pheno$site.z2 <- (pheno$site-mean(pheno$site,na.rm=TRUE))/(sd(pheno$site,na.rm=TRUE)*2)
pheno.all<- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site.z2")]
pheno.all <- pheno.all[complete.cases(pheno.all), ] # 3609

pheno.all$species.fact <- as.numeric(as.factor(pheno.all$species))
sort(unique(pheno.all$species.fact)) # 19, 30 species, 47 with chill0 47 

datalist.all <- with(pheno.all,
                   list( N = nrow(pheno.all),
                         n_sp = length(unique(pheno.all$species.fact)),
                         n_site = length(unique(pheno.all$site.z2)),
                         lday = bb,
                         sp = species.fact,
                         chill1 = chillport.z2,
                         photo = photo.z2,
                         warm = force.z2,
                         site = site.z2
                   ))


mdl.all <- stan("stan/lday_site_sp_chill_inter_poola_ncpwp_2chill_sitefixedstnd_standardized.stan", 
               data = datalist.all,
               iter = 2000, warmup = 1000, chains=4
               # , control = list(adapt_delta = 0.99)
)
save(mdl.all, file="output/allbuds_ncpint_ncpwp_2chillport_2xstandardized_stndsite.Rda")
## The model no longer has any divergent transitions for the terminal buds!
#pairs(sm.sum, pars=c("mu_a","mu_force","mu_chill","mu_photo_ncp")) # this gives a lot of warning messages and not the figure i was hoping/expected
suma <- summary(mdl.all)$summary
range(sumt[, "n_eff"])
range(sumt[, "Rhat"])
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")
# manually to get right order
mu_params <- c("mu_a",
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
               "sigma_a",
               "sigma_b_inter_wp",
               "sigma_b_inter_wc1",
               "sigma_b_inter_pc1",
               "sigma_b_inter_ws",
               "sigma_b_inter_ps",
               "sigma_b_inter_sc1",
               "sigma_y"
)

meanz_all <- suma[mu_params, col4table]
meanz_all
####################################################################
### PPC

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
ext <- rstan::extract(mdl.t)
# 
# # not the most normal 
# hist(ext$a)
# 
y <- pheno.t$tbb

y.ext <- ext$y_hat # I want this to be a matrix, which it is, with one element for each data point in y

ppc_dens_overlay(y, y.ext[1:150, ])
# 
# mean(ext$b_force)
# mean(ext$b_chill)
# mean(ext$b_photo)


# # Begin by checking to see what cue is most important and whether there are strong correlations between cues:
# df.mean.t <- data.frame(bb.force = sumt[grep("b_force", rownames(sumt)), 1],
#                           bb.photo = sumt[grep("b_photo_ncp", rownames(sumt)), 1],
#                           bb.chill = sumt[grep("b_chill", rownames(sumt)), 1])
# 
# df.mean.l <- data.frame(lat.force = sum1l[grep("b_force", rownames(sumt)), 1],
#                         lat.photo = sum1l[grep("b_photo_ncp", rownames(sumt)), 1],
#                         lat.chill = sum1l[grep("b_chill", rownames(sumt)), 1])
# 
# df.mean.t[which(df.mean.t$bb.force > df.mean.t$bb.photo), ] # species 11- rho alb
# df.mean.l[which(df.mean.l$lat.force > df.mean.l$lat.photo), ] #none
# df.mean.t[which(df.mean.t$bb.chill > df.mean.t$bb.force), ] # 14
# # 3, 5,6,8,9,10,12,13,15,17,18,20
# df.mean.l[which(df.mean.l$lat.chill > df.mean.l$lat.force), ] # 16
# #1,2,5,6,7,8,10,11,12,13,15,16,17,18,19,20
# 
# # all correlated
# summary(lm(bb.force~bb.photo, data=df.mean.t))
# summary(lm(bb.force~bb.chill, data=df.mean.t))
# summary(lm(bb.force~bb.photo, data=df.mean.t))
# summary(lm(lat.force~lat.photo, data=df.mean.l))
# summary(lm(lat.force~lat.chill, data=df.mean.l))
# summary(lm(lat.chill~lat.photo, data=df.mean.l))
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
# par(mar = c(5, 10, 2, 1))
# # Lower panel: leaf-out
# plot(seq(-22, 
#          12, 
#          length.out = nrow(meanzl)), 
#      1:nrow(meanzl),
#      type = "n",
#      xlab = "Model estimate change in day of phenological event",
#      ylab = "",
#      yaxt = "n")
# 
# #legend(x = -24, y = 6, bty="n", legend = "b. Leafout", text.font = 2)
# #rasterImage(lopng, -20, 1, -14, 4)
# 
# axis(2, at = nrow(meanzl):1, labels = rownames(meanzl), las = 1, cex.axis = 0.8)
# points(meanzl[,'mean'],
#        nrow(meanzl):1,
#        pch = 16,
#        col = "midnightblue")
# arrows(meanzl[,"75%"], nrow(meanzl):1, meanzl[,"25%"], nrow(meanzl):1,
#        len = 0, col = "black")
# abline(v = 0, lty = 3)
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
