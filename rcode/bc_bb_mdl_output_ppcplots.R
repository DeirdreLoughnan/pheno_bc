## Started Nov 2020 by D. Loughnan 

# Budburst experiment for bc species in 2019
# Code largely based off of budchill code written by D. Flynn and Lizzie --> budchill_analysis.R

library(scales)
library(arm)
library(rstan)
library(shinystan)
library(reshape2)
library(bayesplot)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plyr)

options(mc.cores = parallel::detectCores())

rm(list=ls()) 
options(stringsAsFactors = FALSE)

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
pheno <- rbind.fill(dl, df)
head(pheno)
#pheno <- dl
# because I only had two chilling treatments, I am removing the DF zero chill
pheno <- subset(pheno, chill != "chill0")
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
pheno$site.n[pheno$site.n == "sm"] <- "1"
pheno$site.n[pheno$site.n == "mp"] <- "0"
pheno$site.n[pheno$site.n == "HF"] <- "2"
pheno$site.n[pheno$site.n == "SH"] <- "3"
pheno$site.n <- as.numeric(pheno$site.n)

head(pheno)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 

pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 30 species 

## Load the model output:
load("output/tbb_ncp_dl.Rds")
load("output/latbb50_ncp_dl.Rds")

sumt <- summary(mdl.t)$summary
suml <- summary(mdl.l)$summary

load("output/tbb_ncp_dldf.Rds")
sumtdf <- summary(mdl.t.df)$summary

## The model no longer has any divergent transitions for the terminal buds!
#pairs(sm.sum, pars=c("mu_a","mu_force","mu_chill","mu_photo_ncp")) # this gives a lot of warning messages and not the figure i was hoping/expected

range(sumt[, "n_eff"])
range(sumt[, "Rhat"])

range(suml[, "n_eff"])
range(suml[, "Rhat"])

range(sumtdf[, "n_eff"])
range(sumtdf[, "Rhat"])
#####################################################################
# PPC 

mdl.slopes <- as.data.frame(sumt[grep("b", rownames(sumt)), c(1,6)]) 
mdl.int <- as.data.frame(sumt[grep("a", rownames(sumt)), ]) 
mdl.a <- mdl.int[, 1]
mdl.b <- mdl.slopes[, 1]

# ggplot() +
#   geom_point(data = mdl.slopes, aes(y = row.names(mdl.slopes), x = mean), color = "darkgreen") +
#   labs( x = "doy", y = "Species")

#####################################################################

#plot(mdl, pars="a", ci_level=0.5, outer_level=0.5,col="blue")
# PPC based on the vingette from https://cran.r-project.org/web/packages/bayesplot/vignettes/graphical-ppcs.html 

ext.t <- rstan::extract(mdl.t)
ext.l <- rstan::extract(mdl.l)
ext.dldf <- rstan::extract(mdl.t.df)

# not the most normal 
hist(ext.t$a)

hist(ext.t$b_force)
hist(ext.t$b_chill)
hist(ext.t$b_photo)

y <- pheno.t$tbb

y.ext <- ext.dldf$ypred_new # I want this to be a matrix, which it is, with one element for each data point in y

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
meanztdf <- sumtdf[mu_params, col4fig]

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

rownames(meanztdf) = c("Forcing",
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
#write.table(meanzt.table , "output/term_mdl_esti_dldf.csv", sep = ",", row.names = FALSE)

meanzl.table <- suml[mu_params, col4table]
row.names(meanzl.table) <- row.names(meanzl)
head(meanzl.table)
#write.table(meanzl.table , "output/lat.mdl.esti.csv", sep = ",", row.names = FALSE)

meanztdf.table <- sumtdf[mu_params, col4table]
row.names(meanztdf.table) <- row.names(meanztdf)
head(meanztdf.table)
#write.table(meanzt.table , "output/term_mdl_esti_dldf.csv", sep = ",", row.names = FALSE)

#
# meanzl1.table <- suml1[mu_params, col4table]
# row.names(meanzl1.table) <- row.names(meanzl)
# head(meanzl1.table)
# write.table(meanzl.table , "output/lat.mdl.esti.csv", sep = ",", row.names = FALSE)

# Begin by checking to see what cue is most important and whether there are strong correlations between cues:
df.mean.t <- data.frame(bb.force = sumt[grep("b_force", rownames(sumt)), 1],
                          bb.photo = sumt[grep("b_photo_ncp", rownames(sumt)), 1],
                          bb.chill = sumt[grep("b_chill", rownames(sumt)), 1])

df.mean.l <- data.frame(lat.force = suml[grep("b_force", rownames(sumt)), 1],
                        lat.photo = suml[grep("b_photo_ncp", rownames(sumt)), 1],
                        lat.chill = suml[grep("b_chill", rownames(sumt)), 1])

df.mean.tdf <- data.frame(bb.force = sumtdf[grep("b_force", rownames(sumtdf)), 1],
                        bb.photo = sumtdf[grep("b_photo_ncp", rownames(sumtdf)), 1],
                        bb.chill = sumtdf[grep("b_chill", rownames(sumtdf)), 1])

df.mean.t[which(df.mean.t$bb.force > df.mean.t$bb.photo), ] # none
df.mean.l[which(df.mean.l$lat.force > df.mean.l$lat.photo), ] #none
df.mean.t[which(df.mean.t$bb.chill > df.mean.t$bb.force), ] # 20
df.mean.l[which(df.mean.l$lat.chill > df.mean.l$lat.force), ] # 29

# all correlated
summary(lm(bb.force~bb.photo, data=df.mean.t))# ns
summary(lm(bb.force~bb.chill, data=df.mean.t)) #s
summary(lm(bb.force~bb.photo, data=df.mean.t)) #ns
summary(lm(lat.force~lat.photo, data=df.mean.l))#ns
summary(lm(lat.force~lat.chill, data=df.mean.l)) # s
summary(lm(lat.chill~lat.photo, data=df.mean.l)) #ns

#pdf(file.path( "figures/changes_pheno_dldf.pdf"), width = 7, height = 8)
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

## Plotting both my and the df terminal bud data:
par(mfrow = c(1,1))
# Upper panel: bud burst
plot(seq(-22, 
         12,
         length.out = nrow(meanztdf)), 
     1:nrow(meanztdf),
     type = "n",
     xlab = "",
     ylab = "",
     yaxt = "n")

#legend(x = -20, y = 2, bty="n", legend = "a. Budburst", text.font = 2)
#rasterImage(bbpng, -20, 1, -16, 4)

axis(2, at = nrow(meanztdf):1, labels = rownames(meanztdf), las = 1, cex.axis = 0.8)
points(meanztdf[, 'mean'],
       nrow(meanztdf):1,
       pch = 16,
       col = "midnightblue")
arrows(meanztdf[, "75%"], nrow(meanztdf):1, meanztdf[, "25%"], nrow(meanztdf):1,
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

# add advance/delay arrows
# par(xpd=NA)
# arrows(1, 15.5, 6, 15.5, len=0.1, col = "black")
# legend(5, 16.5, legend="delay", bty="n", text.font = 1, cex=0.75)
# arrows(-1, 15.5, -6, 15.5, len=0.1, col = "black")
# legend(-12, 16.5, legend="advance", bty="n", text.font = 1, cex=0.75)
# legend(-2, 16.5, legend="0", bty="n", text.font = 1, cex=0.75)
# par(xpd=FALSE)
#dev.off()

##########################################################################################################
## Replicating Flynn Figure 2:

b_force <- sumt[grep("b_force", rownames(sumt))]
b_photo <- sumt[grep("b_photo", rownames(sumt))]; b_photo <- b_photo[20:38]
b_chill <- sumt[grep("b_chill", rownames(sumt))]

# Comparisons of trees vs shrubs:
shrubs = c("alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
trees = c("acegla","betpap", "poptre", "popbal")

sp.temp <- c("alninc","alnvir","amelan", "corsto","loninv", "rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","vacmem","vibedu","acegla","betpap", "poptre", "popbal")
species.fact <- as.numeric(as.factor(sort(sp.temp)))
species <- sort(unique(sp.temp))
type <- c("tree", "shrub", "shrub", "shrub", "tree", "shrub", "shrub",  "tree", "tree","shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub")

# if the df data is included:
# shrubs = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG", "alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
# trees = c("ACEPEN", "ACERUB", "ACESAC", "ALNINC", "BETALL", "BETLEN", "BETPAP", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL", "acegla","betpap", "poptre", "popbal")
# 
# species.fact <- sort(unique(pheno.t$species.fact))
# species <- sort(unique(pheno.t$species))
# type <- c("tree","tree", "tree","tree", "shrub", "shrub", "shrub", "tree", "tree", "tree", "shrub","tree" , "shrub", "shrub", "tree", "tree", "tree", "tree", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub", "shrub")
test <- data.frame(species, species.fact, b_force, b_photo, b_chill, type)
head(test)

pdf(file.path( "figures/chill_vs_force_dldf.pdf"), width = 7, height = 8)
ggplot(test, aes(x= b_chill, y = b_force, col = type)) +
   geom_point() +
  ylim (-25, 1) +
  xlim (-30, 0) +
  geom_text(aes(label=species),hjust=0.5, vjust= 1) +
  theme_bw()
dev.off()

pdf(file.path( "figures/chill_vs_photo_dldf.pdf"), width = 7, height = 8)
ggplot(test, aes(x= b_chill, y = b_photo, col = type)) +
  geom_point() +
  ylim (-3.5, 1) +
  xlim (-30, 0) +
  geom_text(aes(label=species),hjust=0.5, vjust= 1) +
  theme_bw()
  geom_text(aes(label=species),hjust=0.5, vjust= 1)
dev.off()

pdf(file.path( "figures/force_vs_photo_dldf.pdf"), width = 7, height = 8)
ggplot(test, aes(x= b_force, y = b_photo, col = type)) +
  geom_point() +
  geom_text(aes(label=species),hjust=0.5, vjust= 1)+
  ylim (-3.5, 0.5) +
  xlim (-22, 0) +
  theme_bw()
dev.off()

# treeshrub = levels(pheno$species)
# treeshrub[treeshrub %in% shrubs] = 1
# treeshrub[treeshrub %in% trees] = 2
# treeshrub = as.numeric(treeshrub)
# par(mar=rep(1,4))
# layout(matrix(c(1, 2, 3, # use layout instead of par(mfrow for more control of where labels end up
#                 4, 5, 6,
#                 7, 8, 9),ncol = 3, byrow = T),
#        widths = c(1, 4, 4),
#        heights = c(4, 4, 1))
# plotblank = function(){plot(1:10, type="n",bty="n",xaxt="n",yaxt="n",ylab="",xlab="")}
# 
# plotblank() 
# text(5,5, "Budburst \n Change (days) due to 5° warming", font = 2, srt = 90) # \n\n add two line breaks
# 
# plot( "b_photo", "b_warm",
#          #  ylab = "Advance due to 5° warming", 
#          # xlab = "Advance due to 4 hr longer photoperiod", 
#          ylim = c(-27, 0.5),
#          xlim = c(-16, 0.5),
#          #  xaxt="n", 
#          group = treeshrub,
#          data = sumt)
# 
# legend("topleft", bty = "n", inset = 0.035, legend = "A.", text.font=2)
# 
# legend("bottomright",
#        pch = "+",
#        col = colz,
#        legend = c("Shrubs","Trees"),
#        inset = 0.02, 
#        bg = 'white')
# 
# plotlet("b_chill1", "b_warm", 
#         # ylab = "Advance due to 5° warming", 
#         #  xlab = "Advance due to 30d 4° chilling", 
#         ylim = c(-27, 0.5),
#         xlim = c(-28, -8),
#         yaxt="n",
#         # xaxt="n", 
#         group = treeshrub,
#         data = sumerb)
# axis(2, seq(0, -25, by = -5), labels = FALSE)
# legend("topleft", bty = "n", inset = 0.035, legend = "B.", text.font=2)
# 
# plotblank()
# text(5,5, "Leafout \n Change (days) due to 5° warming", font = 2, srt = 90)
# 
# plotlet("b_photo", "b_warm", 
#         #    ylab = "Advance due to 5° warming", 
#         #     xlab = "Advance due to 4 hr longer photoperiod", 
#         ylim = c(-27, 0.5),
#         xlim = c(-16, 0.5),
#         group = treeshrub,
#         data = sumerl)
# legend("topleft", bty = "n", inset = 0.035, legend = "C.", text.font=2)
# plotlet("b_chill1", "b_warm", 
#         #   ylab = "Advance due to 5° warming", 
#         #   xlab = "Advance due to 30d 4° chilling", 
#         ylim = c(-27, 0.5),
#         xlim = c(-28, -8),
#         yaxt="n",
#         group = treeshrub,
#         data = sumerl)
# axis(2, seq(0, -25, by = -5), labels = FALSE)
# legend("topleft", bty = "n", inset = 0.035, legend = "D.", text.font=2)
# plotblank()
# 
# plotblank()
# text(5.5, 5, "Change (days) due to 4 hr longer photoperiod", font = 2, pos = 3)
# 
# plotblank()
# text(5.5, 5, "Change (days) due to 30d 4° chilling", font = 2, pos = 3)
# 
# #dev.off();#system(paste("open", file.path(figpath, "Fig2_4panel.pdf"), "-a /Applications/Preview.app"))
# 
# 



