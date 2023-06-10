# started Sept 9, 2022 by Deirdre

# aim of this code is to create eye plot to compare the east to west and shrubs to trees

# Started May 3, 2022 by deirde

# the aim of this code is to generate the model output for my phenology ms
rm(list=ls()) 
options(stringsAsFactors = FALSE)
require(stringr)
require(rstan)
require(plyr)
require(ggbiplot)
require(gridExtra)
require(ggplot2)
library(reshape2)
library(viridis)
#library(tidybayes)# for arranging plots 
library(patchwork) # another way of arranging plots 
# library(rethinking)
require(cowplot)
require(plotrix)

require(ggpubr)
require(dplyr)
require(bayesplot)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
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
dl.wchill$transect <- "west"

df <- read.csv("input/df_dxb_prepped_data.csv")
df.chill <- read.csv("input/chilling_values_eastern.csv")
df.wchill <- merge(df, df.chill, by =c("population","chill"))
df.wchill <- df.wchill[, c("population", "chill","force","photo","lab2", "bb","species", "treatment","Chill_portions","Utah_Model")]
df.wchill$transect <- "east"

# mergeing the my data with DF
pheno <- rbind.fill(dl.wchill, df.wchill)

#head(pheno)
# combined the data has 3197 unique samples
############################################################
# Preping the data for the model
pheno <- subset(pheno, chill != "chill2")

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

# pheno$site2 <- ifelse(pheno$site.n == "2", "1", pheno$site.n)
# pheno$site2 <- ifelse(pheno$site.n == c("1"), "0", pheno$site2)
# pheno$site2 <- ifelse(pheno$site.n == c("3"), "0", pheno$site2)
# pheno$site2 <- ifelse(pheno$site.n == c("4"), "0", pheno$site2)
# 
# pheno$site3 <- ifelse(pheno$site.n == "3", "1", pheno$site.n)
# pheno$site3 <- ifelse(pheno$site.n == c("1"), "0", pheno$site3)
# pheno$site3 <- ifelse(pheno$site.n == c("2"), "0", pheno$site3)
# pheno$site3 <- ifelse(pheno$site.n == c("4"), "0", pheno$site3)
# 
# pheno$site4 <- ifelse(pheno$site.n == "4", "1", pheno$site.n)
# pheno$site4 <- ifelse(pheno$site.n == c("1"), "0", pheno$site4)
# pheno$site4 <- ifelse(pheno$site.n == c("2"), "0", pheno$site4)
# pheno$site4 <- ifelse(pheno$site.n == c("3"), "0", pheno$site4)

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
sort(unique(pheno.t$species.fact)) # 47

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

#load("output/final/bb_4sites_phylo.Rda")

#load("output/bb_4sites_phylo_mini.Rda")
#load("output/bb_4sites_phylo_contphotothermo_zscored_.Rda")
load("output/bb_phylo_contphotothermo_2zscoredMay13.Rda")
sum <- summary(mdl.2z)$summary
fit <- rstan::extract(mdl.2z)

# sum <- summary(mdl.4phyloMini)$summary 
# sum <- summary(mdl.4phyloContWP)$summary 
# fit <- rstan::extract(mdl.4phyloMini)
# fit <- rstan::extract(mdl.4phyloContWp)

chillB <- data.frame(fit$b_chill1)
chillBInterSite2 <- data.frame(fit$b_chill1+fit$b_inter_s2c1)
chillBInterSite3 <- data.frame(fit$b_chill1+fit$b_inter_s3c1)
chillBInterSite4 <- data.frame(fit$b_chill1+fit$b_inter_s4c1)

colnames(chillB) <- unique(pheno.t$species.name)

site2 <- data.frame(fit$b_site2)
site3 <- data.frame(fit$b_site3)
site4 <- data.frame(fit$b_site4)

site <- cbind(site2, site3, site4)

colnames(site) <- c("site2","site3", "site4")
longsite <- melt(site)

chillBsite1 <- chillB
chillBsite2 <- chillB + site2$fit.b_site2
chillBsite3 <- chillB + site3$fit.b_site3
chillBsite4 <- chillB + site4$fit.b_site4

chillBsite2Inter <- chillBInterSite2 + site2$fit.b_site2
chillBsite3Inter <- chillBInterSite3 + site3$fit.b_site3
chillBsite4Inter <- chillBInterSite4 + site4$fit.b_site4

longChill <- melt(chillB)
colnames(longChill) <- c("species.name","chillNoSite")

longChillSite1 <- melt(chillB)
colnames(longChillSite1) <- c("species.name","chillSite")
longChillSite1$site <- "Smithers"

longChillSite2 <- melt(chillBsite2)
colnames(longChillSite2) <- c("species.name","chillSite")
longChillSite2$site <- "Manning Park"

longChillSite3 <- melt(chillBsite3)
colnames(longChillSite3) <- c("species.name","chillSite")
longChillSite3$site <- "Harvard Forest"

longChillSite4 <- melt(chillBsite4)
colnames(longChillSite4) <- c("species.name","chillSite")
longChillSite4$site <- "St. Hippolyte"

longChillSite2Inter <- melt(chillBsite2Inter)
colnames(longChillSite2Inter) <- c("species.name","chillSiteInter")
longChillSite2Inter$site <- "Manning Park"

longChillSite3Inter <- melt(chillBsite3Inter)
colnames(longChillSite3Inter) <- c("species.name","chillSiteInter")
longChillSite3Inter$site <- "Harvard Forest"

longChillSite4Inter <- melt(chillBsite4Inter)
colnames(longChillSite4Inter) <- c("species.name","chillSiteInter")
longChillSite4Inter$site <- "St. Hippolyte"

# add in the needed factors

longChillSite <- rbind(longChillSite1, longChillSite2, longChillSite3, longChillSite4)

longChillSiteInter <- rbind(longChillSite2Inter, longChillSite3Inter, longChillSite4Inter)

longChillInfo <- merge(longChill, spInfo, by = "species.name") # slow to run
#longChillInfo <- merge(longChillInfo, longChillSite2, by = "species.name") # slow to run

longChillInfo <- cbind(longChillInfo, longChillSite1[,2], longChillSite2[,2], longChillSite3[,2], longChillSite4[,2])

longChillInterInfo <- cbind(longChillInfo, longChillSite2Inter[,2], longChillSite3Inter[,2], longChillSite4Inter[,2])

longChillInfo$cue <- "Chilling"
longChillInterInfo$cue <- "Chilling"
head(longChillInfo)
colnames(longChillInfo) <- c("species.name","value", "species","type", "transect","Site1","Site2", "Site3", "Site4", "cue")

colnames(longChillInterInfo) <- c("species.name","value", "species","type", "transect","Site1","Site2", "Site3", "Site4","cue","Site2Inter", "Site3Inter", "Site4Inter")

head(longChillInfo)
head(longChillInterInfo)

# Now make the eye plots:
unique(longChillInfo$transect)



###################################################################
# Photoperiod:
photoB <- data.frame(fit$b_photo)
colnames(photoB) <- unique(pheno.t$species.name)

longPhoto <- melt(photoB)

colnames(longPhoto) <- c("species.name","value")

photoBsite1 <- photoB 
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

longPhotoSite1 <- melt(photoBsite1)
colnames(longPhotoSite1) <- c("species.name","photoSite")
longPhotoSite1$site <- "Smithers"

longPhotoSite2 <- melt(photoBsite2)
colnames(longPhotoSite2) <- c("species.name","photoSite")
longPhotoSite2$site <- "Manning Park"

longPhotoSite3 <- melt(photoBsite3)
colnames(longPhotoSite3) <- c("species.name","photoSite")
longPhotoSite3$site <- "Harvard Forest"
  
longPhotoSite4 <- melt(photoBsite4)
colnames(longPhotoSite4) <- c("species.name","photoSite")
longPhotoSite4$site <- "St. Hippolyte"

longPhotoSite2Inter <- melt(photoBsite2Inter)
colnames(longPhotoSite2Inter) <- c("species.name","photoSiteInter")
longPhotoSite2Inter$site <- "Manning Park"
  
longPhotoSite3Inter <- melt(photoBsite3Inter)
colnames(longPhotoSite3Inter) <- c("species.name","photoSiteInter")
longPhotoSite3Inter$site <- "Harvard Forest"
  
longPhotoSite4Inter <- melt(photoBsite4Inter)
colnames(longPhotoSite4Inter) <- c("species.name","photoSiteInter")
longPhotoSite4Inter$site <- "St. Hippolyte"

# add in the needed factors
longPhotoSite <- rbind(longPhotoSite1, longPhotoSite2, longPhotoSite3, longPhotoSite4)
longPhotoSiteInter <- rbind(longPhotoSite2Inter, longPhotoSite3Inter, longPhotoSite4Inter)

longPhotoInfo <- merge(longPhoto, spInfo, by = "species.name") # slow to run
#longPhotoInfo <- merge(longPhotoInfo, longPhotoSite2, by = "species.name") # slow to run

longPhotoInfo <- cbind(longPhotoInfo, longPhotoSite1[,2], longPhotoSite2[,2], longPhotoSite3[,2], longPhotoSite4[,2])

longPhotoInterInfo <- cbind(longPhotoInfo, longPhotoSite2Inter[,2], longPhotoSite3Inter[,2], longPhotoSite4Inter[,2])


longPhotoInfo$cue <- "Photoperiod"
longPhotoInterInfo$cue <- "Photoperiod"

head(longPhotoInfo)
colnames(longPhotoInfo) <- c("species.name","value", "species","type", "transect","Site1","Site2", "Site3", "Site4","cue")

colnames(longPhotoInterInfo) <- c("species.name","value", "species","type", "transect","Site1","Site2", "Site3", "Site4","Site2Inter", "Site3Inter", "Site4Inter","cue")

head(longPhotoInfo)
head(longPhotoInterInfo)



#longPhotoInfo <- merge(longPhoto, spInfo, by = "species.name") # slow to run


###################################################################
# Forcing:
forceB <- data.frame(fit$b_warm)
colnames(forceB) <- unique(pheno.t$species.name)

longForce <- melt(forceB)

colnames(longForce) <- c("species.name","value")

forceBsite1 <- forceB 
forceBsite2 <- forceB + site2$fit.b_site2
forceBsite3 <- forceB + site3$fit.b_site3
forceBsite4 <- forceB + site4$fit.b_site4

forceBInterSite2 <- data.frame(forceB+fit$b_inter_s2c1)
forceBInterSite3 <- data.frame(forceB+fit$b_inter_s3c1)
forceBInterSite4 <- data.frame(forceB+fit$b_inter_s4c1)

forceBsite2Inter <- forceBInterSite2 + site2$fit.b_site2
forceBsite3Inter <- forceBInterSite3 + site3$fit.b_site3
forceBsite4Inter <- forceBInterSite4 + site4$fit.b_site4

longForce <- melt(forceB)
colnames(longForce) <- c("species.name","forceNoSite")

longForceSite1 <- melt(forceBsite1)
colnames(longForceSite1) <- c("species.name","forceSite")
longForceSite1$site <- "Smithers"

longForceSite2 <- melt(forceBsite2)
colnames(longForceSite2) <- c("species.name","forceSite")
longForceSite2$site <- "Manning Park"

longForceSite3 <- melt(forceBsite3)
colnames(longForceSite3) <- c("species.name","forceSite")
longForceSite3$site <- "Harvard Forest"

longForceSite4 <- melt(forceBsite4)
colnames(longForceSite4) <- c("species.name","forceSite")
longForceSite4$site <- "St. Hippolyte"

longForceSite2Inter <- melt(forceBsite2Inter)
colnames(longForceSite2Inter) <- c("species.name","forceSiteInter")
longForceSite2Inter$site <- "Manning Park"

longForceSite3Inter <- melt(forceBsite3Inter)
colnames(longForceSite3Inter) <- c("species.name","forceSiteInter")
longForceSite3Inter$site <- "Harvard Forest"

longForceSite4Inter <- melt(forceBsite4Inter)
colnames(longForceSite4Inter) <- c("species.name","forceSiteInter")
longForceSite4Inter$site <- "St. Hippolyte"

# add in the needed factors
longForceSite <- rbind(longForceSite1, longForceSite2, longForceSite3, longForceSite4)
longForceSiteInter <- rbind(longForceSite2Inter, longForceSite3Inter, longForceSite4Inter)

longForceInfo <- merge(longForce, spInfo, by = "species.name") # slow to run

longForceInfo <- cbind(longForceInfo, longForceSite1[,2], longForceSite2[,2], longForceSite3[,2], longForceSite4[,2])

longForceInterInfo <- cbind(longForceInfo, longForceSite2Inter[,2], longForceSite3Inter[,2], longForceSite4Inter[,2])

longForceInfo$cue <- "Forcing"
longForceInterInfo$cue <- "Forcing"

head(longForceInfo)
colnames(longForceInfo) <- c("species.name","value", "species","type", "transect","Site1","Site2", "Site3", "Site4","cue")

colnames(longForceInterInfo) <- c("species.name","value", "species","type", "transect","Site1","Site2", "Site3", "Site4","Site2Inter", "Site3Inter", "Site4Inter","cue")

head(longForceInfo)
head(longForceInterInfo)


head(longForceInfo)
# Now make the eye plots:
unique(longForceInfo$transect)


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

cueST2 <- ggplot(data = longCues, aes(x = cue, y = value)) + 
  stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = type), position = position_dodge(0.9))+
  theme_classic()+
 # ylim(-40,5) +
  theme(axis.text.x = element_text( size=12,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Plant type", y = "Cue response", main = NA)+ 
  theme(legend.title = element_blank()) +  scale_fill_manual(values = c("#cc6a70ff","cyan4"))

# pdf("figures/cueST.pdf", width = 8, height = 5)
# cueST
# dev.off()

pdf("figures/cueST28Chill.pdf", width = 8, height = 5)
cueST2
dev.off()

#head(meanz4)
# Cues by species

# want them to be ordered by tree 
tree <- read.tree("input/SBphylo_phenobc.tre")
length(tree$tip.label) #47

tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
tree$tip.label[tree$tip.label=="Rhamnus_arguta"] <- "Rhamnus_frangula"
tree$tip.label[tree$tip.label=="Alnus_alnobetula"] <- "Alnus_viridis"
tree$tip.label[tree$tip.label== "Fagus_grandifolia_var._caroliniana"] <- "Fagus_grandifolia"
tree$tip.label[tree$tip.label== "Spiraea_alba_var._latifolia"] <- "Spiraea_alba"


spOrder <- tree$tip.label
spOrder <- gsub("_"," ", spOrder)
longChillInfo$mean <- rowMeans(longChillInfo[,c("Site1","Site2","Site3","Site4")], na.rm=TRUE)
longChillInfo$species.name <- gsub("_"," ", longChillInfo$species.name)

chillSp <- ggplot() + 
  stat_eye(data = longChillInfo, aes(x = factor(species.name, level = spOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "#cc6a70ff") +
  geom_hline(yintercept = sum["mu_b_chill1",1], linetype="dashed") +
  theme_classic() +  
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
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

longForceInfo$mean <- rowMeans(longForceInfo[,c("Site1","Site2","Site3","Site4")], na.rm=TRUE)
longForceInfo$species.name <- gsub("_"," ", longForceInfo$species.name)

forceSp <- ggplot() + 
  stat_eye(data = longForceInfo, aes(x = factor(species.name, level = spOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  geom_hline(yintercept = sum["mu_b_warm",1], linetype="dashed") +
  theme_classic() +  
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
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
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 

longPhotoInfo$mean <- rowMeans(longPhotoInfo[,c("Site1","Site2","Site3","Site4")], na.rm=TRUE)
longPhotoInfo$species.name <- gsub("_"," ", longPhotoInfo$species.name)

photoSp <- ggplot() + 
  stat_eye(data = longPhotoInfo, aes(x = factor(species.name, level = spOrder), y = value), .width = c(.90, .5), cex = 0.75, fill = "cyan4")+
  geom_hline(yintercept = sum["mu_b_photo",1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Photoperiod response", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "c)", cex =5) 

# pdf("figures/cueSp.pdf", height =12, width = 12)
# grid.arrange(chillSp,forceSp, photoSp, nrow = 3)
# dev.off()

pdf("figures/cueSp8Chill.pdf", height =12, width = 12)
plot_grid(chillSp,forceSp, photoSp, nrow = 3, align = "v", rel_heights = c(1/4, 1/4, 1.2/3))
dev.off()
##### Site specific plots:
#site 2
# cueEWSite2 <- ggplot(data = longCues, aes(x = cue, y = Site2)) + 
#   stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = transect), position = position_dodge(0.9)) +
#   theme_classic() +   
#   theme(axis.text.x = element_text( size=10,
#                                     #angle = 78, 
#                                     hjust=1),
#         axis.title=element_text(size=9) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Transect", y = "Cue response", main = NA)+ 
#   scale_color_identity(name = "Model fit",
#                        breaks = c("black"),
#                        labels = c("Model Posterior"),
#                        guide = guide_legend(override.aes = list(
#                          linetype = c(NA),
#                          shape = c(8)))) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
#   scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))
# 
# cueEWSite3 <- ggplot(data = longCues, aes(x = cue, y = Site3)) + 
#   stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = transect), position = position_dodge(0.9)) +
#   theme_classic() +   
#   theme(axis.text.x = element_text( size=10,
#                                     #angle = 78, 
#                                     hjust=1),
#         axis.title=element_text(size=9) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Transect", y = "Cue response", main = NA)+ 
#   scale_color_identity(name = "Model fit",
#                        breaks = c("black"),
#                        labels = c("Model Posterior"),
#                        guide = guide_legend(override.aes = list(
#                          linetype = c(NA),
#                          shape = c(8)))) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
#   scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))
# 
# cueEWSite4 <- ggplot(data = longCues, aes(x = cue, y = Site4)) + 
#   stat_eye( .width = c(.90, .5), cex = 0.75, aes(fill = transect), position = position_dodge(0.9)) +
#   theme_classic() +   
#   theme(axis.text.x = element_text( size=10,
#                                     #angle = 78, 
#                                     hjust=1),
#         axis.title=element_text(size=9) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Transect", y = "Cue response", main = NA)+ 
#   scale_color_identity(name = "Model fit",
#                        breaks = c("black"),
#                        labels = c("Model Posterior"),
#                        guide = guide_legend(override.aes = list(
#                          linetype = c(NA),
#                          shape = c(8)))) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
#   scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))
# 
# 
# ggplot() + 
#   stat_eye(data = longCues, aes(x = cue, y = Site2, fill = "#cc6a70ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
#   stat_eye(data = longCues, aes(x = cue, y = Site3, fill = "#f9b641ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
#   stat_eye(data = longCues, aes(x = cue, y = Site4, fill = "cyan4"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
#   theme_classic() +   
#   theme(axis.text.x = element_text( size=10,
#                                     #angle = 78, 
#                                     hjust=1),
#         axis.title=element_text(size=9) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Cue", y = "Cue response", main = NA)+ 
#   scale_color_identity(name = "Model fit",
#                        breaks = c("black"),
#                        labels = c("Model Posterior"),
#                        guide = guide_legend(override.aes = list(
#                          linetype = c(NA),
#                          shape = c(8)))) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 10, label = "a)", cex =5)  +
#   scale_fill_manual(values = c("#cc6a70ff","#f9b641ff","cyan4"))

# How do you jitter the eye plots?
siteOrder <- c("Smithers","Manning Park", "St. Hippolyte","Harvard Forest")

siteChill <- ggplot() + 
  stat_eye(data = longChillSite, aes(x = factor(site, level = siteOrder), y = chillSite, fill = "#cc6a70ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
    ylim (-60,30) +
    theme_classic() +   
    theme(legend.position = "none",
          axis.title = element_text( size=15),
          axis.text.y = element_text( size=15),
          axis.text.x = element_text( size=15,
          angle = 78, 
          hjust=1)) +  
  annotate("text", x = 0.85, y = 25, label = "a)", cex =8) +
  labs( x = "Site", y = "Chilling response", main = "Site level chilling", cex = 5)+
  scale_fill_manual(values = c("#cc6a70ff"))

siteForce <- ggplot() + 
  stat_eye(data = longForceSite, aes(x = factor(site, level = siteOrder), y = forceSite, fill = "#f9b641ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9), color = "black") +
  ylim (-60,30) +
  theme_classic() +   
  theme(legend.position = "none", 
        axis.title = element_text( size=15),
        axis.text.y = element_text( size=15),
        axis.text.x = element_text( size=15, angle = 78,  hjust=1))+  
  annotate("text", x = 0.85, y = 25, label = "b)", cex =8) +
  labs( x = "Site", y = "Forceing response", main = "Site level forcing", cex = 5)+
  scale_fill_manual(values = c("#f9b641ff"))

sitePhoto <- ggplot() + 
  stat_eye(data = longPhotoSite, aes(x = factor(site, level = siteOrder), y = photoSite, fill = "cyan4"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  ylim (-60,30) +
  theme_classic() +   
  theme(legend.position = "none", 
        axis.title = element_text( size=15),
        axis.text.y = element_text( size=15),
        axis.text.x = element_text( size=15, angle = 78, hjust=1)) +  
  annotate("text", x = 0.85, y = 25, label = "c)", cex =8) +
  labs( x = "Site", y = "Photoperiod response", main = "Site level photoperiod", cex = 5)+
  scale_fill_manual(values = c("cyan4"))

pdf("figures/site4Cue.pdf", width = 15, height = 6)
plot_grid(siteChill, siteForce, sitePhoto, ncol = 3, nrow = 1)
dev.off()


####################################
# Make the dataset long
longForceInfo <- longForceInfo[, c("species.name","species","type","transect","Site1","Site2","Site3","Site4","cue")]
longF <- melt(longForceInfo, id = c("species.name","species","type","transect", "cue"))
names(longF)<- c("species.name","species","type","transect","cue","site","value")

longChillInfo <- longChillInfo[, c("species.name","species","type","transect","Site1","Site2","Site3","Site4","cue")]
longC <- melt(longChillInfo, id = c("species.name","species","type","transect", "cue"))
names(longC)<- c("species.name","species","type","transect","cue","site","value")

longPhotoInfo <- longPhotoInfo[, c("species.name","species","type","transect","Site1","Site2","Site3","Site4","cue")]
longP <- melt(longPhotoInfo, id = c("species.name","species","type","transect", "cue"))
names(longP)<- c("species.name","species","type","transect","cue","site","value")

longest <- rbind(longC, longF, longP)

longest$site <- as.character(longest$site)
longest$site[which(longest$site == "Site1")] <- "Smithers"
longest$site[which(longest$site == "Site2")] <- "Manning Park"
longest$site[which(longest$site == "Site3")] <- "Harvard Forest"
longest$site[which(longest$site == "Site4")] <- "St. Hippolyte"
longest$site <- as.factor(longest$site)

siteOrder <- c("Smithers","Manning Park", "St. Hippolyte","Harvard Forest")


pdf("figures/site4CueGrouped8Chill.pdf", width = 12, height =4)
ggplot() + 
  stat_eye(data = longest, aes(x = as.factor(cue), y = value, fill = factor(site, level = siteOrder)), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
#  ylim (-40,10) +
  theme_classic() +   
  theme(legend.position = "right", 
        legend.title = element_blank(),
        axis.text.x = element_text( size= 16),
        axis.text.y = element_text( size= 12),
        axis.title=element_text(size = 18)) +
  labs( x = "Treatment cue", y = "Cue response", main = NA)+
  scale_fill_manual(values = c("Smithers" = "deepskyblue3",
                              "Manning Park" = "palegreen4", 
                              "St. Hippolyte"="darkorchid3", 
                              "Harvard Forest" = "tomato3"))
dev.off()

# "bb_hpsite3" ="darkred",
# "bb_hpsite4"="darkorchid4",
# "bb_hpsite1"="deepskyblue3",
# "bb_hpsite2"= "forestgreen",
# "bb_lpsite1"="deepskyblue1",
# "bb_lpsite2"="palegreen3",
# "bb_lpsite3"="tomato1",
# "bb_lpsite4"="darkorchid1"
pdf("figures/site4CueGrouped.pdf", width = 12, height =4)
ggplot() + 
  stat_eye(data = longest, aes(x = factor(site, level = siteOrder), y = value, fill = cue, group = cue), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  ylim (-45,10) +
  theme_classic() +   
  theme(legend.position = "right", 
        axis.text.x = element_text( size= 12),
        axis.text.y = element_text( size= 12),
        axis.title=element_text(size = 15)) +
  labs( x = "Site", y = "Cue response", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 0.5, y = 9, label = "a)", cex =5) + scale_fill_manual(values = c("#cc6a70ff","cyan4", "#f9b641ff"))
dev.off()

# + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Plant type", y = "Cue response", main = NA)+ 
#   scale_color_identity(name = "Model fit",
#                        breaks = c("black"),
#                        labels = c("Model Posterior"),
#                        guide = guide_legend(override.aes = list(
#                          linetype = c(NA),
#                          shape = c(8)))) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.75, y = 10, label = "a)", cex =5) + scale_fill_manual(values = c("#cc6a70ff","cyan4", "#f9b641ff"))

# WIth the interaction:
siteChillInter <- ggplot() + 
  stat_eye(data = longChillSiteInter, aes(x = site, y = chillSiteInter, fill = "#cc6a70ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  ylim (-60,50) +
  theme_classic() +   
  theme(legend.position = "none") +
  labs( x = "Site", y = "Chilling response", main = "Site level chilling with site interactions")+
  scale_fill_manual(values = c("#cc6a70ff"))

siteForceInter <- ggplot() + 
  stat_eye(data = longForceSiteInter, aes(x = site, y = forceSiteInter, fill = "#f9b641ff"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9), color = "black") +
  ylim (-60,50) +
  theme_classic() +   
  theme(legend.position = "none") +
  labs( x = "Site", y = "Forceing response", main = "Site level forcing with site interactions")+
  scale_fill_manual(values = c("#f9b641ff"))

sitePhoto <- ggplot() + 
  stat_eye(data = longPhotoSiteInter, aes(x = site, y = photoSiteInter, fill = "cyan4"), .width = c(.90, .5), cex = 0.75, position = position_dodge(0.9)) +
  ylim (-60,50) +
  theme_classic() +   
  theme(legend.position = "none") +
  labs( x = "Site", y = "Photoperiod response", main = "Site level photoperiod with site interactions") +
  scale_fill_manual(values = c("cyan4"))

pdf("figures/siteCueInter.pdf", width = 15, height = 5)
ggarrange(siteChill, siteForce, sitePhoto,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
dev.off()
###########################################
longCuesInter <- rbind(longForceInterInfo, longChillInterInfo, longPhotoInterInfo)

longChillSub <- longChillInterInfo[,c("species.name","species","type","transect","value", "cue", "Site2Inter", "Site3Inter","Site4Inter" )]
longerChill <- melt(longChillSub, id = c("species.name","species","type","transect","value", "cue"))

colnames(longerChill) <-  c("species.name","species","type", "transect","value","cue", "site","siteValue")

ggplot() + 
  stat_eye(data = longerChill, aes(x = site, y = siteValue, fill = "#cc6a70ff"), .width = c(.90, .5), cex = 0.75 
                   #, position = position_jitter(width = 1)
  ) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Cue", y = "Cue response", main = NA) + 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.position = "none") +  annotate("text", x = 0.5, y = 50, label = "a)", cex =5)  +
  scale_fill_manual(values = c("#cc6a70ff"))


###################################################################\
longForceSub <- longForceInterInfo[,c("species.name","species","type","transect","value", "cue", "Site2Inter", "Site3Inter","Site4Inter" )]
longerForce <- melt(longForceSub, id = c("species.name","species","type","transect","value", "cue"))

colnames(longerForce) <-  c("species.name","species","type", "transect","value","cue", "site","siteValue")

ggplot() + 
  stat_eye(data = longerForce, aes(x = site, y = siteValue, fill = "#f9b641ff"), .width = c(.90, .5), cex = 0.75 
           #, position = position_jitter(width = 1)
  ) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Cue", y = "Cue response", main = NA) + 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.position = "none") +  annotate("text", x = 0.5, y = 50, label = "a)", cex =5)  +
  scale_fill_manual(values = c("#f9b641ff"))
#############################################################################
longPhotoSub <- longPhotoInterInfo[,c("species.name","species","type","transect","value", "cue", "Site2Inter", "Site3Inter","Site4Inter" )]
longerPhoto <- melt(longPhotoSub, id = c("species.name","species","type","transect","value", "cue"))

colnames(longerPhoto) <-  c("species.name","species","type", "transect","value","cue", "site","siteValue")

ggplot() + 
  stat_eye(data = longerPhoto, aes(x = site, y = siteValue, fill = "cyan4"), .width = c(.90, .5), cex = 0.75 
           #, position = position_jitter(width = 1)
  ) +
  theme_classic() +   
  theme(axis.text.x = element_text( size=10,
                                    #angle = 78, 
                                    hjust=1),
        axis.title=element_text(size=9) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Cue", y = "Cue response", main = NA) + 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.position = "none") +  annotate("text", x = 0.5, y = 50, label = "a)", cex =5)  +
  scale_fill_manual(values = c("cyan4"))
##################################################################################
ggplot() + 
  ggdist::stat_eye(data = longCuesInter, aes(x = cue, y = Site2Inter, fill = "#cc6a70ff"), .width = c(.90, .5), cex = 0.75 
                   , position = position_dodge(2)) +
  stat_eye(data = longCuesInter, aes(x = cue, y = Site3Inter, fill = "#f9b641ff"), .width = c(.90, .5), cex = 0.75
           , position = position_dodge(2)) +
           
  stat_eye(data = longCuesInter, aes(x = cue, y = Site4Inter, fill = "cyan4"), .width = c(.90, .5), cex = 0.75, position = position_dodge(2)) +
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


### Histograms of shrub vs tree
col1 <- rgb(204 / 255, 102 / 255, 119 / 255, alpha = 0.8)
col2 <- rgb(68 / 255, 170 / 255, 153 / 255, alpha = 0.6)

treeLongC <- subset(longC, type == "tree")
shrubLongC <- subset(longC, type == "shrub")

xmin = -(mean(treeLongC$value)) + quantile(treeLongC$value, c(0.9))
xmax = (mean(treeLongC$value) + quantile(treeLongC$value, c(0.9)))

tempTC <- treeLongC[treeLongC$value < xmin & treeLongC$value > xmax, ]

xmin = -(mean(shrubLongC$value)) + quantile(treeLongC$value, c(0.9))
xmax = (mean(shrubLongC$value) + quantile(treeLongC$value, c(0.9)))
tempSC <- shrubLongC[shrubLongC$value < xmin & shrubLongC$value > xmax, ]


### Forcing 
treeLongF <- subset(longF, type == "tree")
shrubLongF <- subset(longF, type == "shrub")

xmin = -(mean(treeLongF$value)) + quantile(treeLongF$value, c(0.9))
xmax = (mean(treeLongF$value) + quantile(treeLongF$value, c(0.9)))

tempTF <- treeLongF[treeLongF$value < xmin & treeLongF$value > xmax, ]

xmin = -(mean(shrubLongF$value)) + quantile(shrubLongF$value, c(0.9))
xmax = (mean(shrubLongF$value) + quantile(shrubLongF$value, c(0.9)))
tempSF <- shrubLongF[shrubLongF$value < xmin & shrubLongF$value > xmax, ]

### Photoperiod
treeLongP <- subset(longP, type == "tree")
shrubLongP <- subset(longP, type == "shrub")

xmin = -(mean(treeLongP$value)) + quantile(treeLongP$value, c(0.9))
xmax = (mean(treeLongP$value) + quantile(treeLongP$value, c(0.9)))

tempTP <- treeLongP[treeLongP$value < xmin & treeLongP$value > xmax, ]

xmin = -(mean(shrubLongP$value)) + quantile(shrubLongP$value, c(0.9))
xmax = (mean(shrubLongP$value) + quantile(shrubLongP$value, c(0.9)))
tempSP <- shrubLongP[shrubLongP$value < xmin & shrubLongP$value > xmax, ]


pdf("figures/fullHistogramST8Chill.pdf", width =12, height = 6)
par(mfrow = c(1,3))
hist(shrubLongC$value, col = col1, xlab = "Cue Response", main = NA, ylab = "Frequency", cex.main =2, cex.lab = 2, cex.axis = 1.75, ylim = c(0, 125000))
hist(treeLongC$value, col = col2, add = T)
title(main = "Chilling", adj = 0, cex.main = 2)
text(-50, 120000, label = "a)", cex = 2)

hist(shrubLongF$value, col = col1, xlab = "Cue Response", main = NA, ylab = "Frequency", cex.main =2, cex.lab = 2, cex.axis = 1.75, ylim = c(0, 75000))
hist(treeLongF$value, col = col2, add = T)
title(main = "Forcing", adj = 0, cex.main = 2)
text(-25, 72000, label = "b)", cex = 2)

hist(shrubLongP$value, col = col1, xlab = "Cue Response", main = NA, ylab = "Frequency", cex.main =2, cex.lab = 2, cex.axis = 1.75, ylim = c(0, 125000))
hist(treeLongP$value, col = col2, add = T)
title(main = "Photoperiod", adj = 0, cex.main = 2)
text(-12, 120000, label = "c)", cex = 2)

legend("topright", c("Shrubs", "Trees"), col = c( col1, col2), bty = "n", pt.cex =3, pch = 19, cex = 2)

dev.off()

pdf("figures/fullHistogramSTNorm8Chill.pdf", width =12, height = 6)
par(mfrow = c(1,3), mar = c(5.1, 5.1, 4.1, 2.1))
hist(shrubLongC$value, col = col1, xlab = NA, main = NA, ylab = "Normalized frequency", cex.main =2, cex.lab = 2.5, cex.axis = 1.75, prob = T,ylim = c(0, 0.125))
hist(treeLongC$value, col = col2, add = T, prob = T)
title(main = "Chilling", adj = 0, cex.main = 2)
text(-48, 0.123, label = "a)", cex = 2)
abline(v =mean(shrubLongC$value), col = "#CC6677", lwd = 3)
abline(v =mean(treeLongC$value), col = "cyan4", lwd = 3)

par(mar = c(5.1, 4.1, 4.1, 2.1))
hist(shrubLongF$value, col = col1, xlab = "Cue Response", main = NA, ylab = NA, cex.main =2, cex.lab = 2.5, cex.axis = 1.75, prob = T,ylim = c(0, 0.125))
hist(treeLongF$value, col = col2, add = T, prob = T)
title(main = "Forcing", adj = 0, cex.main = 2)
text(-25, 0.123, label = "b)", cex = 2)
abline(v =mean(shrubLongF$value), col = "#CC6677", lwd = 3)
abline(v =mean(treeLongF$value), col = "cyan4", lwd = 3)

hist(shrubLongP$value, col = col1, xlab = NA, main = NA, ylab = NA, cex.main =2, cex.lab = 2, cex.axis = 1.75, prob = T,ylim = c(0, 0.4))
hist(treeLongP$value, col = col2, add = T, prob = T)
title(main = "Photoperiod", adj = 0, cex.main = 2)
text(-12, 0.4, label = "c)", cex = 2)
abline(v =mean(shrubLongP$value), col = "#CC6677", lwd = 3)
abline(v =mean(treeLongP$value), col = "cyan4", lwd = 3)

legend("topright", c("Shrubs", "Trees"), col = c( col1, col2), bty = "n", pt.cex =3, pch = 19, cex = 2)
dev.off()

# now fig just showing 90%
hist(tempSC$value, col = col1, xlab = "Cue Response", main = NA, ylab = "Frequency", cex.main =2, cex.lab = 2, cex.axis = 1.75)
hist(tempTC$value, add = T,col = col2)
title(main = "Chilling", adj = 0, cex.main = 2)

hist(tempT$value, col = col1, xlab = "Cue Response", main = NA, ylab = "Frequency", cex.main =2, cex.lab = 2, cex.axis = 1.75, ylim = c(0, 45000))
hist(tempS$value, add = T,col = col2)
title(main = "Forcing", adj = 0, cex.main = 2)
text(-10, 0.06, label = "e)", cex = 2)
legend("topright", c("Non-interacting", "Interacting"), col = c( col1.sp, col4.sp), bty = "n", pt.cex =2, pch = 19, cex = 1.45, inset = c(-0.03, 0))
dev.off()