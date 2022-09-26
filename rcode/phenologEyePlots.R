# started Sept 9, 2022 by Deirdre

# aim of this code is to create eye plot to compare the east to west and shrubs to trees

# Started May 3, 2022 by deirde

# the aim of this code is to generate the model output for my phenology ms
rm(list=ls()) 
options(stringsAsFactors = FALSE)

require(rstan)
require(ggbiplot)
require(gridExtra)
require(ggplot2)
require(bayesplot)
require(dplyr)
library(reshape2)
library(viridis)
library(bayesplot)
library(tidybayes)# for arranging plots 
library(patchwork) # another way of arranging plots 
# library(rethinking)
require(cowplot)
require(plotrix)
require(stringr)


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  

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

head(pheno)
# combined the data has 3197 unique samples
############################################################
# Preping the data for the model

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

pheno$transect.n <- pheno$transect
pheno$transect.n[pheno$transect.n == "east"] <- "1"
pheno$transect.n[pheno$transect.n == "west"] <- "0"
pheno$transect.n <- as.numeric(pheno$transect.n)

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

pheno$transect.n <- pheno$transect
pheno$transect.n[pheno$transect.n == "east"] <- "1"
pheno$transect.n[pheno$transect.n == "west"] <- "0"
pheno$transect.n <- as.numeric(pheno$transect.n)

#going to split it into analysis of terminal bb and lateral bb
# Starting with the terminal buds:
#pheno.term <- pheno[,c("tbb", "chill.n", "force.n", "photo.n", "site.n", "species", "lab2")]
pheno.term <- pheno[,c("bb", "force.z2", "photo.z2", "population", "species", "lab2","Utah_Model","Chill_portions","chillport.z2", "site2.z2", "site3.z2","site4.z2")]
pheno.t <- pheno.term[complete.cases(pheno.term), ] # 3609

pheno.t <- pheno.term[complete.cases(pheno.term$bb), ] # 1780 rows data 
pheno.t$species <- tolower(pheno.t$species)
pheno.t$species.fact <- as.numeric(as.factor(pheno.t$species))
sort(unique(pheno.t$species.fact)) # 49 

# now get the phylogeny and pair it with species names:
spInfo <- read.csv("input/species_list.csv")
head(spInfo)
head(pheno.t)
spInfo$transect[spInfo$transect == "e"] <- "Eastern"
spInfo$transect[spInfo$transect == "w"] <- "Western"
spInfo$transect[spInfo$transect == "ew"] <- "Both"

# change the df to lowercase:

pheno.t <- merge(pheno.t, spInfo, by = "species")

# now load the stan output:

load("output/final/bb_4sites_phylo.Rda")
sum <- summary(mdl.4phylo)$summary 

fit <- rstan::extract(mdl.4phylo)

chillB <- data.frame(fit$b_chill1)
chillBInterSite2 <- data.frame(fit$b_chill1+fit$b_inter_s2c1)
chillBInterSite3 <- data.frame(fit$b_chill1+fit$b_inter_s3c1)
chillBInterSite4 <- data.frame(fit$b_chill1+fit$b_inter_s4c1)

colnames(chillB) <- unique(pheno.t$species.name)

# site2 <- data.frame(fit[,c("b_site2","b_inter_s2c1")])
# 
# site3 <- data.frame(fit$b_site3)
# 
# site4 <- data.frame(fit$b_site4)
# 
# site <- cbind(site2, site3, site4)

colnames(site) <- c("site2","site3", "site4")
longsite <- melt(site)

chillBsite2 <- chillB + site2$fit.b_site2
chillBsite3 <- chillB + site3$fit.b_site3
chillBsite4 <- chillB + site4$fit.b_site4

chillBsite2Inter <- chillBInterSite2 + site2$fit.b_site2
chillBsite3Inter <- chillBInterSite3 + site3$fit.b_site3
chillBsite4Inter <- chillBInterSite4 + site4$fit.b_site4

longChill <- melt(chillB)
colnames(longChill) <- c("species.name","chillNoSite")

longChillSite2 <- melt(chillBsite2)
colnames(longChillSite2) <- c("species.name","chillSite2")

longChillSite3 <- melt(chillBsite3)
colnames(longChillSite3) <- c("species.name","chillSite3")

longChillSite4 <- melt(chillBsite4)
colnames(longChillSite4) <- c("species.name","chillSite4")

longChillSite2Inter <- melt(chillBsite2Inter)
colnames(longChillSite2Inter) <- c("species.name","chillSite2Inter")

longChillSite3Inter <- melt(chillBsite3Inter)
colnames(longChillSite3Inter) <- c("species.name","chillSite3Inter")

longChillSite4Inter <- melt(chillBsite4Inter)
colnames(longChillSite4Inter) <- c("species.name","chillSite4Inter")

# add in the needed factors

longChillInfo <- merge(longChill, spInfo, by = "species.name") # slow to run
#longChillInfo <- merge(longChillInfo, longChillSite2, by = "species.name") # slow to run

longChillInfo <- cbind(longChillInfo, longChillSite2[,2], longChillSite3[,2], longChillSite4[,2])

longChillInterInfo <- cbind(longChillInfo, longChillSite2Inter[,2], longChillSite3Inter[,2], longChillSite4Inter[,2])

head(longChillInfo)
colnames(longChillInfo) <- c("species.name","value", "species","type", "transect","Site2", "Site3", "Site4")

colnames(longChillInterInfo) <- c("species.name","value", "species","type", "transect","Site2", "Site3", "Site4","Site2Inter", "Site3Inter", "Site4Inter")

head(longChillInfo)
head(longChillInterInfo)

# Now make the eye plots:
unique(longChillInfo$transect)

longChillInfo$cue <- "Chilling"
longChillInterInfo$cue <- "Chilling"

###################################################################
# Photoperiod:
photoB <- data.frame(fit$b_photo)
colnames(photoB) <- unique(pheno.t$species.name)

longPhoto <- melt(photoB)

colnames(longPhoto) <- c("species.name","value")

photoBsite2 <- photoB + site2$fit.b_site2
photoBsite3 <- photoB + site3$fit.b_site3
photoBsite4 <- photoB + site4$fit.b_site4

photoBInterSite2 <- data.frame(fit$b_photo+fit$b_inter_s2c1)
photoBInterSite3 <- data.frame(fit$b_photo+fit$b_inter_s3c1)
photoBInterSite4 <- data.frame(fit$b_photo+fit$b_inter_s4c1)

photoBsite2Inter <- photoBInterSite2 + site2$fit.b_site2
photoBsite3Inter <- photoBInterSite3 + site3$fit.b_site3
photoBsite4Inter <- photoBInterSite4 + site4$fit.b_site4


longPhoto <- melt(photoB)
colnames(longPhoto) <- c("species.name","photoNoSite")

longPhotoSite2 <- melt(photoBsite2)
colnames(longPhotoSite2) <- c("species.name","photoSite2")

longPhotoSite3 <- melt(photoBsite3)
colnames(longPhotoSite3) <- c("species.name","photoSite3")

longPhotoSite4 <- melt(photoBsite4)
colnames(longPhotoSite4) <- c("species.name","photoSite4")

longPhotoSite2Inter <- melt(photoBsite2Inter)
colnames(longPhotoSite2Inter) <- c("species.name","photoSite2Inter")

longPhotoSite3Inter <- melt(photoBsite3Inter)
colnames(longPhotoSite3Inter) <- c("species.name","photoSite3Inter")

longPhotoSite4Inter <- melt(photoBsite4Inter)
colnames(longPhotoSite4Inter) <- c("species.name","photoSite4Inter")

# add in the needed factors

longPhotoInfo <- merge(longPhoto, spInfo, by = "species.name") # slow to run
#longPhotoInfo <- merge(longPhotoInfo, longPhotoSite2, by = "species.name") # slow to run

longPhotoInfo <- cbind(longPhotoInfo, longPhotoSite2[,2], longPhotoSite3[,2], longPhotoSite4[,2])

longPhotoInterInfo <- cbind(longPhotoInfo, longPhotoSite2Inter[,2], longPhotoSite3Inter[,2], longPhotoSite4Inter[,2])

head(longPhotoInfo)
colnames(longPhotoInfo) <- c("species.name","value", "species","type", "transect","Site2", "Site3", "Site4")

colnames(longPhotoInterInfo) <- c("species.name","value", "species","type", "transect","Site2", "Site3", "Site4","Site2Inter", "Site3Inter", "Site4Inter")

head(longPhotoInfo)
head(longPhotoInterInfo)

colnames(longPhotoInfo) <- c("species.name","value", "species","type", "transect","Site2", "Site3", "Site4")


#longPhotoInfo <- merge(longPhoto, spInfo, by = "species.name") # slow to run

longPhotoInfo$cue <- "Photoperiod"
longPhotoInterInfo$cue <- "Photoperiod"

###################################################################
# Forcing:
forceB <- data.frame(fit$b_warm)
colnames(forceB) <- unique(pheno.t$species.name)

longForce <- melt(forceB)

colnames(longForce) <- c("species.name","value")

forceBsite2 <- forceB + site2$fit.b_site2
forceBsite3 <- forceB + site3$fit.b_site3
forceBsite4 <- forceB + site4$fit.b_site4

forceBInterSite2 <- data.frame(fit$b_warm+fit$b_inter_s2c1)
forceBInterSite3 <- data.frame(fit$b_warm+fit$b_inter_s3c1)
forceBInterSite4 <- data.frame(fit$b_warm+fit$b_inter_s4c1)

forceBsite2Inter <- forceBInterSite2 + site2$fit.b_site2
forceBsite3Inter <- forceBInterSite3 + site3$fit.b_site3
forceBsite4Inter <- forceBInterSite4 + site4$fit.b_site4

longForce <- melt(forceB)
colnames(longForce) <- c("species.name","forceNoSite")

longForceSite2 <- melt(forceBsite2)
colnames(longForceSite2) <- c("species.name","forceSite2")

longForceSite3 <- melt(forceBsite3)
colnames(longForceSite3) <- c("species.name","forceSite3")

longForceSite4 <- melt(forceBsite4)
colnames(longForceSite4) <- c("species.name","forceSite4")

longForceSite2Inter <- melt(forceBsite2Inter)
colnames(longForceSite2Inter) <- c("species.name","forceSite2Inter")

longForceSite3Inter <- melt(forceBsite3Inter)
colnames(longForceSite3Inter) <- c("species.name","forceSite3Inter")

longForceSite4Inter <- melt(forceBsite4Inter)
colnames(longForceSite4Inter) <- c("species.name","forceSite4Inter")

# add in the needed factors

longForceInfo <- merge(longForce, spInfo, by = "species.name") # slow to run

longForceInfo <- cbind(longForceInfo, longForceSite2[,2], longForceSite3[,2], longForceSite4[,2])

longForceInterInfo <- cbind(longForceInfo, longForceSite2Inter[,2], longForceSite3Inter[,2], longForceSite4Inter[,2])

head(longForceInfo)
colnames(longForceInfo) <- c("species.name","value", "species","type", "transect","Site2", "Site3", "Site4")

colnames(longForceInterInfo) <- c("species.name","value", "species","type", "transect","Site2", "Site3", "Site4","Site2Inter", "Site3Inter", "Site4Inter")

head(longForceInfo)
head(longForceInterInfo)

colnames(longForceInfo) <- c("species.name","value", "species","type", "transect","Site2", "Site3", "Site4")

head(longForceInfo)
# Now make the eye plots:
unique(longForceInfo$transect)

longForceInfo$cue <- "Forcing"
longForceInterInfo$cue <- "Forcing"

###################################################################
longCues <- rbind(longForceInfo, longChillInfo, longPhotoInfo)

cueEW <- ggplot(data = longCues, aes(x = transect, y = value)) + 
  stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = cue), position = position_dodge(0.9)) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Transect", y = "Cue response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
  scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))


pdf("figures/cueEW.pdf", width = 8, height = 5)
cueEW
dev.off()
# + scale_x_discrete(labels = c("primary consumer \n(n = 556)", "primary producer \n(n = 472)", "secondary consumer \n(n = 247)", "tertiary consumer \n(n = 6)"))

#cues  by architecture
cueST <- ggplot(data = longCues, aes(x = type, y = value)) + 
  stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = cue), position = position_dodge(0.9))+
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Plant type", y = "Cue response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 10, label = "a)", cex =5) + scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))

pdf("figures/cueST.pdf", width = 8, height = 5)
cueST
dev.off()

# Cues by species
chillSp <- ggplot() + 
  stat_eye(data = longChillInfo, aes(x = species.name, y = value), .width = c(.90, .5), cex = 0.75, fill = "#cc6a70ff")+
  theme_classic() +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Transect", y = "Chilling response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "a)", cex =5) 
forceSp <- ggplot() + 
  stat_eye(data = longForceInfo, aes(x = species.name, y = value), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff")+
  theme_classic() +  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Transect", y = "Forcing response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "a)", cex =5) 

photoSp <- ggplot() + 
  stat_eye(data = longPhotoInfo, aes(x = species.name, y = value), .width = c(.90, .5), cex = 0.75, fill = "cyan4")+
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Transect", y = "Photoperiod response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "a)", cex =5) 
pdf("figures/cueSp.pdf", height =12, width = 12)
grid.arrange(chillSp,forceSp, photoSp, nrow = 3)
dev.off()

##### Site specific plots:
#site 2
cueEWSite2 <- ggplot(data = longCues, aes(x = cue, y = Site2)) + 
  stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = transect), position = position_dodge(0.9)) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Transect", y = "Cue response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
  scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))

cueEWSite3 <- ggplot(data = longCues, aes(x = cue, y = Site3)) + 
  stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = transect), position = position_dodge(0.9)) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Transect", y = "Cue response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
  scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))

cueEWSite4 <- ggplot(data = longCues, aes(x = cue, y = Site4)) + 
  stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = transect), position = position_dodge(0.9)) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Transect", y = "Cue response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
  scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))


ggplot() + 
  stat_eye(data = longCues, aes(x = cue, y = Site2, fill = "#cc6a70ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  stat_eye(data = longCues, aes(x = cue, y = Site3, fill = "#f9b641ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  stat_eye(data = longCues, aes(x = cue, y = Site4, fill = "cyan4"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Cue", y = "Cue response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
  scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))


###########################################
longCuesInter <- rbind(longForceInterInfo, longChillInterInfo, longPhotoInterInfo)

ggplot() + 
  stat_eye(data = longCuesInter, aes(x = cue, y = Site2Inter, fill = "#cc6a70ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  stat_eye(data = longCuesInter, aes(x = cue, y = Site3Inter, fill = "#f9b641ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  stat_eye(data = longCuesInter, aes(x = cue, y = Site4Inter, fill = "cyan4"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Cue", y = "Cue response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
  scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))
