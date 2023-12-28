# Started May 3, 2022 by deirde

# the aim of this code is to generate the model output for my phenology ms
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(rstan)
library(shinystan)
#library(reshape2)
#library(bayesplot)
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(phytools)
#library(ggpubr)
library(lattice)
require(cowplot)


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  

#load("output/final/ew_phylo_output_newpriors.Rda")
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

head(pheno)
# combined the data has 3197 unique samples
############################################################
# Preping the data for the model
pheno <- subset(pheno, chill != "chill2")

# pheno$force.n <- pheno$force
# pheno$force.n[pheno$force.n == "HF"] <- "1"
# pheno$force.n[pheno$force.n == "LF"] <- "0"
# pheno$force.n <- as.numeric(pheno$force.n)
pheno$force.n <- pheno$force
# pheno$force.n[pheno$force.n == "HF"] <- "1"
# pheno$force.n[pheno$force.n == "LF"] <- "0"
# pheno$force.n <- as.numeric(pheno$force.n)
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
# pheno <- pheno %>%
#   mutate ( site2 = if_else(site.n == 2, 1, 0),
#            site3 = if_else(site.n == 3, 1, 0),
#            site4 = if_else(site.n == 4, 1, 0))

pheno$site2 <- ifelse(pheno$site.n == "2", "1", pheno$site.n)
pheno$site2 <- ifelse(pheno$site.n == c("1"), "0", pheno$site2)
pheno$site2 <- ifelse(pheno$site.n == c("3"), "0", pheno$site2)
pheno$site2 <- ifelse(pheno$site.n == c("4"), "0", pheno$site2)

pheno$site3 <- ifelse(pheno$site.n == "3", "1", pheno$site.n)
pheno$site3 <- ifelse(pheno$site.n == c("1"), "0", pheno$site3)
pheno$site3 <- ifelse(pheno$site.n == c("2"), "0", pheno$site3)
pheno$site3 <- ifelse(pheno$site.n == c("4"), "0", pheno$site3)

pheno$site4 <- ifelse(pheno$site.n == "4", "1", pheno$site.n)
pheno$site4 <- ifelse(pheno$site.n == c("1"), "0", pheno$site4)
pheno$site4 <- ifelse(pheno$site.n == c("2"), "0", pheno$site4)
pheno$site4 <- ifelse(pheno$site.n == c("3"), "0", pheno$site4)

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
sort(unique(pheno.t$species.fact)) # 49 

# now get the phylogeny and pair it with species names:
spInfo <- read.csv("input/species_list.csv")
head(spInfo)
head(pheno.t)
# change the df to lowercase:

pheno.t <- merge(pheno.t, spInfo, by = "species")

#################################################################
# load("output/final/ew_phylo_output_newpriors_allncp.Rda")
# sumew <- summary(mdl.ewphylo)$summary

#load("output/final/bb_4sites_phylo_mini.Rda")
#load("output/bb_4sites_phylo_contphotothermo_zscored_Apr19.Rda")
load("output/bb_phylo_contphotothermo_2zscoredMay13.Rda")
# sum <- summary(mdl.4phyloContWP)$summary
# post <- rstan::extract(mdl.4phyloContWP)
sum <- summary(mdl.2z)$summary
post <- rstan::extract(mdl.2z)
#############################################
col4fig <- c("mean","sd","25%","50%","75%","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")

mu_params_4 <- c( 
               #"a_z",
               # "lam_interceptsa",
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
                "mu_b_inter_s4c1")

meanz4m <- sum[mu_params_4, col4table]

rownames(meanz4m) = c( 
                      #"Root trait intercept", "Lambda",
                      "Forcing",
                      "Photoperiod",
                      "Chilling",
                      "Manning Park",
                       "Harvard Forest",
                       "St. Hippolyte",
                      "Forcing x photoperiod",
                      "Forcing x chilling",
                      "Photoperiod x chilling",
                      "Forcing x Manning Park",
                      "Photoperiod x Manning Park",
                      "Chilling x Manning Park",
                      "Forcing x Harvard Forest",
                      "Photoperiod x Harvard Forest",
                      "Chilling x Harvard Forest",
                      "Forcing x St. Hippolyte",
                      "Photoperiod x St. Hippolyte",
                      "Chilling x St. Hippolyte"            
)

# meanz4.table <- sum[mu_params_4, col4table]
# row.names(meanz4.table) <- row.names(meanz4)
# head(meanz4.table)
#write.table(meanzew.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)

###########################################################################################
##### 4 site plots ########################################################################
###########################################################################################

pdf(file.path( "figures/changes_pheno_4sitesMay13.pdf"), width = 7, height = 5)
par(mfrow = c(1,1), mar = c(5, 10, 2, 1))
# Upper panel: bud burst
plot(seq(-20, 
         20,
         length.out = nrow(meanz4m)), 
     1:nrow(meanz4m),
     type = "n",
     xlab = "Estimated change in budburst day",
     ylab = "",
     yaxt = "n")

#legend(x = -20, y = 2, bty="n", legend = "a. Budburst", text.font = 2)
#rasterImage(bbpng, -20, 1, -16, 4)

axis(2, at = nrow(meanz4m):1, labels = rownames(meanz4m), las = 1, cex.axis = 0.8)
points(meanz4m[, 'mean'],
       nrow(meanz4m):1,
       pch = 16,
       col = "cyan4",
       cex = 2)
arrows(meanz4m[, "97.5%"], nrow(meanz4m):1, meanz4m[, "2.5%"], nrow(meanz4m):1,
       len = 0, col = "black")
abline(v = 0, lty = 3)

dev.off()

## Replicating Flynn Figure 2:

bForceSp <-  data.frame(sum[grep("^b_warm", rownames(sum)), c("mean","25%", "75%")])
bChillSp <-  data.frame(sum[grep("^b_chill1", rownames(sum)), c("mean","25%", "75%")])
bPhotoSp <-  data.frame(sum[grep("^b_photo\\[", rownames(sum)), c("mean","25%", "75%")])

bForceSp$species <- sort(unique(pheno.t$species))
bChillSp$species <- sort(unique(pheno.t$species))
bPhotoSp$species <- sort(unique(pheno.t$species))


# shrubs.both = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG", "alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
# trees.both = c("ACEPEN", "ACERUB", "ACESAC", "BETALL", "BETLEN", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL", "acegla","betpap", "poptre", "popbal")
head(spInfo)
bForceSp <- merge(bForceSp, spInfo, by = "species", all =T)
bChillSp <- merge(bChillSp, spInfo, by = "species", all =T)
bPhotoSp <- merge(bPhotoSp, spInfo, by = "species", all =T)

names(bForceSp) <- c("species", "bForceMean", "bForce25", "bForce75", "species.name","type","transect")
names(bChillSp) <- c("species", "bChillMean", "bChill25", "bChill75", "species.name","type","transect")
names(bPhotoSp) <- c("species", "bPhotoMean", "bPhoto25", "bPhoto75", "species.name","type","transect")

cueOut <- merge(bForceSp, bChillSp, by = c("species", "species.name","type","transect"))
cueOut <- merge(cueOut, bPhotoSp, by = c("species", "species.name","type","transect"))

# pdf(file.path( "figures/cueCompare.pdf"), width = 5, height = 6)
photoForce <-  ggplot(cueOut, aes(x = bPhotoMean, y = bForceMean , col =type, shape = transect)) +
  geom_point(size = 2) +
  #guides(color = "none", size = "none") +
ylim (-17, 1) +
xlim (-7, 0) +
labs( y = "High Forcing", x = "Long Photoperiod") +
geom_text(cueOut, mapping = aes(label= species),hjust=0.5, vjust= 1, show.legend = F) +
  # geom_bar(stat="identity", color="black", 
  #          position=position_dodge()) +
  geom_errorbar(aes(xmin= bPhoto25, xmax = bPhoto75), width= 0) +
  geom_errorbar(aes(ymin= bForce25, ymax = bForce75), width= 0) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none",
      legend.key=element_rect(fill="white"))  +
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))


# color_scheme_set("viridis")
#pdf(file.path( "figures/chill_vs_force_dldf4site.pdf"), width = 5, height = 6)
chillForce <- ggplot(cueOut, aes(x= bChillMean, y = bForceMean, col = type, shape = transect)) +
  geom_point() +
  ylim (-17, 1) +
  xlim (-35, 0) +
  labs( y = "High forcing", x = "High chilling") +
  geom_text(aes(label=species),hjust= 0.5, vjust= 1.5, show.legend = F) +
  geom_errorbar(aes(xmin= bChill25, xmax = bChill75), width= 0) +
  geom_errorbar(aes(ymin= bForce25, ymax = bForce75), width= 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="none",
        legend.key=element_rect(fill="white"))  +
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
#dev.off()

#pdf(file.path( "figures/chill_vs_photo_dldf4site.pdf"), width = 5, height = 6)
chillPhoto <- ggplot(cueOut, aes(y= bChillMean, x = bPhotoMean, col = type, shape = transect)) +
  geom_point() +
  xlim (-7, 1) +
  ylim (-35, 0) +
  labs( x = "Long Photoperiod", y = "High chilling") +
  geom_text(aes(label=species),hjust= 0.5, vjust= 1.5, show.legend = F) +
  geom_errorbar(aes(ymin= bChill25, ymax = bChill75), width= 0) +
  geom_errorbar(aes(xmin= bPhoto25, xmax = bPhoto75), width= 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))

pdf("figures/cueCompare.pdf", width = 12, height = 4)
ggarrange(chillForce, photoForce, chillPhoto,
             labels = c("A", "B", "C"),
             ncol = 3)
dev.off()

## Plotting the day to bb with the cues on the y-axis 
head(pheno)
pheno.t$species <- tolower(pheno.t$species)

phenoInfoMean <- ddply(pheno.t, c("species" ), summarize, mean = mean(bb, na.rm = TRUE),
                      # n = n(),
                       sd = sd(bb, na.rm = T))
                       #se = sd/sqrt(n))
# names(term.bb.both) <- c("species.both", "mean")

term.both <- merge(phenoInfoMean, cueOut, by = c("species"), all =TRUE)
term.both <- term.both[,c("species","mean","bForceMean","bChillMean","bPhotoMean", "type", "transect")]
names(term.both) <- c("species","mean","bForceMean","bChillMean","bPhotoMean", "type", "transect")
#term.both <- term.both[complete.cases(term.both), ] 


tf.both <-  ggplot(term.both, aes(y = bForceMean, x= mean,col = type, shape = transect)) +
  geom_point() +
  ylim (-17, 0) +
  xlim (5, 55) +
  labs(x = "Mean day of budburst", y = "High forcing") +
  geom_text(aes(label=species),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
tf.both
#dev.off()

#pdf(file.path( "figures/chill_dobb_dldf4sites.pdf"), width = 5, height = 6)
tc.both <- ggplot(term.both, aes(y = bChillMean, x= mean,col = type, shape = transect)) +
  geom_point() +
  ylim (-32, -4) +
  xlim (10, 55) +
  labs(x = "Mean day of budburst", y = "Chilling") +
  geom_text(aes(label=species),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
tc.both
#dev.off()

#pdf(file.path( "figures/photo_dobb_dldf.pdf"), width = 5, height = 6)
tp.both <- ggplot(term.both, aes(y = bPhotoMean, x= mean,col = type, shape = transect)) +
  geom_point() +
  ylim (-6.5, -1) +
  xlim (10, 55) +
  labs(x = "Mean day of budburst", y = "Long photoperiod")+
  geom_text(aes(label=species),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"))+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
tp.both

pdf(file.path( "figures/cues_dobb_dldf4sites.pdf"), width = 15, height = 6)

ggarrange (tc.both, tf.both,tp.both,
           labels = c("A", "B", "C"),
           ncol = 3, nrow = 1)

dev.off()

# Let's plot some interactions:
a_sp = mean(sum[grep("a_sp", rownames(sum)), c("mean")])
a_sp5 <- round(quantile(post$a_sp, c(0.05)),1)
a_sp95 <- round(quantile(post$a_sp, c(0.95)),1)
a_sp25 <- round(quantile(post$a_sp, c(0.25)),1)
a_sp75 <- round(quantile(post$a_sp, c(0.75)),1)
a_sp <- cbind(a_sp, a_sp5,a_sp95, a_sp25,a_sp75)

#a_z = (sum[grep("a_z", rownames(sum)), c("mean")])
mu_b_warm = sum[grep("mu_b_warm", rownames(sum)), c("mean")]
mu_b_warm5 <- round(quantile(post$mu_b_warm, c(0.05)),1)
mu_b_warm95 <- round(quantile(post$mu_b_warm, c(0.95)),1)
mu_b_warm25 <- round(quantile(post$mu_b_warm, c(0.25)),1)
mu_b_warm75 <- round(quantile(post$mu_b_warm, c(0.75)),1)
mu_b_warm <- (cbind(mu_b_warm, mu_b_warm5,mu_b_warm95, mu_b_warm25,mu_b_warm75))

mu_b_photo = sum[grep("mu_b_photo", rownames(sum)), c("mean")]
mu_b_photo5 <- round(quantile(post$mu_b_photo, c(0.05)),1)
mu_b_photo95 <- round(quantile(post$mu_b_photo, c(0.95)),1)
mu_b_photo25 <- round(quantile(post$mu_b_photo, c(0.25)),1)
mu_b_photo75 <- round(quantile(post$mu_b_photo, c(0.75)),1)
mu_b_photo <- (cbind(mu_b_photo, mu_b_photo5,mu_b_photo95, mu_b_photo25,mu_b_photo75))

mu_b_chill1 = sum[grep("mu_b_chill1", rownames(sum)), c("mean")]
mu_b_chill5 <- round(quantile(post$mu_b_chill, c(0.05)),1)
mu_b_chill95 <- round(quantile(post$mu_b_chill, c(0.95)),1)
mu_b_chill25 <- round(quantile(post$mu_b_chill, c(0.25)),1)
mu_b_chill75 <- round(quantile(post$mu_b_chill, c(0.75)),1)
mu_b_chill1 <- (cbind(mu_b_chill1, mu_b_chill5,mu_b_chill95, mu_b_chill25,mu_b_chill75))

mu_b_inter_pc1 = sum[grep("mu_b_inter_pc1", rownames(sum)), c("mean")]
mu_b_inter_pc15 <- round(quantile(post$mu_b_inter_pc1, c(0.05)),1)
mu_b_inter_pc195 <- round(quantile(post$mu_b_inter_pc1, c(0.95)),1)
mu_b_inter_pc125 <- round(quantile(post$mu_b_inter_pc1, c(0.25)),1)
mu_b_inter_pc175 <- round(quantile(post$mu_b_inter_pc1, c(0.75)),1)
mu_b_inter_pc1 <- (cbind(mu_b_inter_pc1, mu_b_inter_pc15,mu_b_inter_pc195, mu_b_inter_pc125,mu_b_inter_pc175))

mu_b_inter_wp = sum[grep("mu_b_inter_wp", rownames(sum)), c("mean")]
mu_b_inter_wp5 <- round(quantile(post$mu_b_inter_wp, c(0.05)),1)
mu_b_inter_wp95 <- round(quantile(post$mu_b_inter_wp, c(0.95)),1)
mu_b_inter_wp25 <- round(quantile(post$mu_b_inter_wp, c(0.25)),1)
mu_b_inter_wp75 <- round(quantile(post$mu_b_inter_wp, c(0.75)),1)
mu_b_inter_wp <- (cbind(mu_b_inter_wp, mu_b_inter_wp5,mu_b_inter_wp95, mu_b_inter_wp25,mu_b_inter_wp75))

mu_b_inter_wc1 = sum[grep("mu_b_inter_wc1", rownames(sum)), c("mean")]
mu_b_inter_wc15 <- round(quantile(post$mu_b_inter_wc1, c(0.05)),1)
mu_b_inter_wc195 <- round(quantile(post$mu_b_inter_wc1, c(0.95)),1)
mu_b_inter_wc125 <- round(quantile(post$mu_b_inter_wc1, c(0.25)),1)
mu_b_inter_wc175 <- round(quantile(post$mu_b_inter_wc1, c(0.75)),1)
mu_b_inter_wc1 <- (cbind(mu_b_inter_wc1, mu_b_inter_wc15,mu_b_inter_wc195, mu_b_inter_wc125,mu_b_inter_wc175))

mu_b_inter_ws2 = sum[grep("mu_b_inter_ws2", rownames(sum)), c("mean")]
mu_b_inter_ws25 <- round(quantile(post$mu_b_inter_ws2, c(0.05)),1)
mu_b_inter_ws295 <- round(quantile(post$mu_b_inter_ws2, c(0.95)),1)
mu_b_inter_ws225 <- round(quantile(post$mu_b_inter_ws2, c(0.25)),1)
mu_b_inter_ws275 <- round(quantile(post$mu_b_inter_ws2, c(0.75)),1)
mu_b_inter_ws2 <- (cbind(mu_b_inter_ws2, mu_b_inter_ws25,mu_b_inter_ws295, mu_b_inter_ws225,mu_b_inter_ws275))

mu_b_inter_s2c1 = sum[grep("mu_b_inter_s2c1", rownames(sum)), c("mean")]
mu_b_inter_s2c15 <- round(quantile(post$mu_b_inter_s2c1, c(0.05)),1)
mu_b_inter_s2c195 <- round(quantile(post$mu_b_inter_s2c1, c(0.95)),1)
mu_b_inter_s2c125 <- round(quantile(post$mu_b_inter_s2c1, c(0.25)),1)
mu_b_inter_s2c175 <- round(quantile(post$mu_b_inter_s2c1, c(0.75)),1)
mu_b_inter_s2c1 <- (cbind(mu_b_inter_s2c1, mu_b_inter_s2c15,mu_b_inter_s2c195, mu_b_inter_s2c125,mu_b_inter_s2c175))

mu_b_inter_ps2 = sum[grep("mu_b_inter_ps2", rownames(sum)), c("mean")]
mu_b_inter_ps25 <- round(quantile(post$mu_b_inter_ps2, c(0.05)),1)
mu_b_inter_ps295 <- round(quantile(post$mu_b_inter_ps2, c(0.95)),1)
mu_b_inter_ps225 <- round(quantile(post$mu_b_inter_ps2, c(0.25)),1)
mu_b_inter_ps275 <- round(quantile(post$mu_b_inter_ps2, c(0.75)),1)
mu_b_inter_ps2 <- (cbind(mu_b_inter_ps2, mu_b_inter_ps25,mu_b_inter_ps295, mu_b_inter_ps225,mu_b_inter_ps275))

mu_b_inter_ws3 = sum[grep("mu_b_inter_ws3", rownames(sum)), c("mean")]
mu_b_inter_ws35 <- round(quantile(post$mu_b_inter_ws3, c(0.05)),1)
mu_b_inter_ws395 <- round(quantile(post$mu_b_inter_ws3, c(0.95)),1)
mu_b_inter_ws325 <- round(quantile(post$mu_b_inter_ws3, c(0.25)),1)
mu_b_inter_ws375 <- round(quantile(post$mu_b_inter_ws3, c(0.75)),1)
mu_b_inter_ws3 <- (cbind(mu_b_inter_ws3, mu_b_inter_ws35,mu_b_inter_ws395, mu_b_inter_ws325,mu_b_inter_ws375))

mu_b_inter_s3c1 = sum[grep("mu_b_inter_s3c1", rownames(sum)), c("mean")]
mu_b_inter_s3c15 <- round(quantile(post$mu_b_inter_s3c1, c(0.05)),1)
mu_b_inter_s3c195 <- round(quantile(post$mu_b_inter_s3c1, c(0.95)),1)
mu_b_inter_s3c125 <- round(quantile(post$mu_b_inter_s3c1, c(0.25)),1)
mu_b_inter_s3c175 <- round(quantile(post$mu_b_inter_s3c1, c(0.75)),1)
mu_b_inter_s3c1 <- (cbind(mu_b_inter_s3c1, mu_b_inter_s3c15,mu_b_inter_s3c195, mu_b_inter_s3c125,mu_b_inter_s3c175))

mu_b_inter_ps3 = sum[grep("mu_b_inter_ps3", rownames(sum)), c("mean")]
mu_b_inter_ps35 <- round(quantile(post$mu_b_inter_ps3, c(0.05)),1)
mu_b_inter_ps395 <- round(quantile(post$mu_b_inter_ps3, c(0.95)),1)
mu_b_inter_ps325 <- round(quantile(post$mu_b_inter_ps3, c(0.25)),1)
mu_b_inter_ps375 <- round(quantile(post$mu_b_inter_ps3, c(0.75)),1)
mu_b_inter_ps3 <- (cbind(mu_b_inter_ps3, mu_b_inter_ps35,mu_b_inter_ps395, mu_b_inter_ps325,mu_b_inter_ps375))

mu_b_inter_ws4 = sum[grep("mu_b_inter_ws4", rownames(sum)), c("mean")]
mu_b_inter_ws45 <- round(quantile(post$mu_b_inter_ws4, c(0.05)),1)
mu_b_inter_ws495 <- round(quantile(post$mu_b_inter_ws4, c(0.95)),1)
mu_b_inter_ws425 <- round(quantile(post$mu_b_inter_ws4, c(0.25)),1)
mu_b_inter_ws475 <- round(quantile(post$mu_b_inter_ws4, c(0.75)),1)
mu_b_inter_ws4 <- (cbind(mu_b_inter_ws4, mu_b_inter_ws45,mu_b_inter_ws495, mu_b_inter_ws425,mu_b_inter_ws475))

mu_b_inter_s4c1 = sum[grep("mu_b_inter_s4c1", rownames(sum)), c("mean")]
mu_b_inter_s4c15 <- round(quantile(post$mu_b_inter_s4c1, c(0.05)),1)
mu_b_inter_s4c195 <- round(quantile(post$mu_b_inter_s4c1, c(0.95)),1)
mu_b_inter_s4c125 <- round(quantile(post$mu_b_inter_s4c1, c(0.25)),1)
mu_b_inter_s4c175 <- round(quantile(post$mu_b_inter_s4c1, c(0.75)),1)
mu_b_inter_s4c1 <- (cbind(mu_b_inter_s4c1, mu_b_inter_s4c15,mu_b_inter_s4c195, mu_b_inter_s4c125,mu_b_inter_s4c175))

mu_b_inter_ps4 = sum[grep("mu_b_inter_ps4", rownames(sum)), c("mean")]
mu_b_inter_ps45 <- round(quantile(post$mu_b_inter_ps4, c(0.05)),1)
mu_b_inter_ps495 <- round(quantile(post$mu_b_inter_ps4, c(0.95)),1)
mu_b_inter_ps425 <- round(quantile(post$mu_b_inter_ps4, c(0.25)),1)
mu_b_inter_ps475 <- round(quantile(post$mu_b_inter_ps4, c(0.75)),1)
mu_b_inter_ps4 <- (cbind(mu_b_inter_ps4, mu_b_inter_ps45,mu_b_inter_ps495, mu_b_inter_ps425,mu_b_inter_ps475))

b_site2 = sum[grep("b_site2", rownames(sum)), c("mean")]
b_site25 <- round(quantile(post$b_site2, c(0.05)),1)
b_site295 <- round(quantile(post$b_site2, c(0.95)),1)
b_site225 <- round(quantile(post$b_site2, c(0.25)),1)
b_site275 <- round(quantile(post$b_site2, c(0.75)),1)
b_site2 <- (cbind(b_site2, b_site25,b_site295, b_site225,b_site275))

b_site3 = sum[grep("b_site3", rownames(sum)), c("mean")]
b_site35 <- round(quantile(post$b_site3, c(0.05)),1)
b_site395 <- round(quantile(post$b_site3, c(0.95)),1)
b_site325 <- round(quantile(post$b_site3, c(0.25)),1)
b_site375 <- round(quantile(post$b_site3, c(0.75)),1)
b_site3 <- (cbind(b_site3, b_site35,b_site395, b_site325,b_site375))

b_site4 = sum[grep("b_site4", rownames(sum)), c("mean")]
b_site45 <- round(quantile(post$b_site4, c(0.05)),1)
b_site495 <- round(quantile(post$b_site4, c(0.95)),1)
b_site425 <- round(quantile(post$b_site4, c(0.25)),1)
b_site475 <- round(quantile(post$b_site4, c(0.75)),1)
b_site4 <- (cbind(b_site4, b_site45,b_site495, b_site425,b_site475))

b_warm = sum[grep("b_warm\\[", rownames(sum)), c("mean")]
mu_b_warm5 <- round(quantile(post$mu_b_warm, c(0.05)),1)
mu_b_warm95 <- round(quantile(post$mu_b_warm, c(0.95)),1)
mu_b_warm25 <- round(quantile(post$mu_b_warm, c(0.25)),1)
mu_b_warm75 <- round(quantile(post$mu_b_warm, c(0.75)),1)
mu_b_warm <- (cbind(mu_b_warm, mu_b_warm5,mu_b_warm95, mu_b_warm25,mu_b_warm75))

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# # plot the interactions
# # Warm X chill
# 
hfData <- subset(pheno, force == "HF" )
lfData <- subset(pheno, force == "LF")

# Make the other parameters constant
hf <- unique(hfData$force.z2)
lf <- unique(lfData$force.z2)

photo <- -0.5033863
siteSM <- 0
chill1 <- seq( -2, 2, length.out =30)
lf <- seq( -0.6723508, -0.3568628,  length.out = 30)
hf <- seq( 0.2722241, 0.5877121, length.out =30)
  # c( -2, -1, -0.7642814, -0.4072595, -0.4023109, -0.3493703,  0.2750890,  0.2977055,  0.4308763,  0.5308110,  0.8457874,  0.9457221, 1,2)


# plot first for the high forcing
bb_hfc = a_sp[1] + b_site2[1] * siteSM + b_site3[1] * siteSM + b_site4[1] * siteSM + mu_b_warm[1] * hf + 
  mu_b_photo[1] * photo + mu_b_chill1[1] * chill1 +
  mu_b_inter_wp[1] * (hf*photo) +
  mu_b_inter_wc1[1] * (hf*chill1) + mu_b_inter_pc1[1] * (photo*chill1) +
  mu_b_inter_s2c1[1] * (chill1*siteSM) + mu_b_inter_ws2[1] * (hf*siteSM) +mu_b_inter_ps2[1] * (photo*siteSM) +
  mu_b_inter_s3c1[1] * (chill1*siteSM) + mu_b_inter_ws3[1] * (hf*siteSM) +mu_b_inter_ps3[1] * (photo*siteSM) +
  mu_b_inter_s4c1[1] * (chill1*siteSM) + mu_b_inter_ws4[1] * (hf*siteSM) +mu_b_inter_ps4[1] * (photo*siteSM)

bb_hfc5 = a_sp[2] + b_site2[2] * siteSM + b_site3[2] * siteSM + b_site4[2] * siteSM + mu_b_warm[2] * hf + 
  mu_b_photo[2] * photo + mu_b_chill1[2] * chill1 +
  mu_b_inter_wp[2] * (hf*photo) +
  mu_b_inter_wc1[2] * (hf*chill1) + mu_b_inter_pc1[2] * (photo*chill1) +
  mu_b_inter_s2c1[2] * (chill1*siteSM) + mu_b_inter_ws2[2] * (hf*siteSM) +mu_b_inter_ps2[2] * (photo*siteSM) +
  mu_b_inter_s3c1[2] * (chill1*siteSM) + mu_b_inter_ws3[2] * (hf*siteSM) +mu_b_inter_ps3[2] * (photo*siteSM) +
  mu_b_inter_s4c1[2] * (chill1*siteSM) + mu_b_inter_ws4[2] * (hf*siteSM) +mu_b_inter_ps4[2] * (photo*siteSM)

bb_hfc95 = a_sp[3] + b_site2[3] * siteSM + b_site3[3] * siteSM + b_site4[3] * siteSM + mu_b_warm[3] * hf + 
  mu_b_photo[3] * photo + mu_b_chill1[3] * chill1 +
  mu_b_inter_wp[3] * (hf*photo) +
  mu_b_inter_wc1[3] * (hf*chill1) + mu_b_inter_pc1[3] * (photo*chill1) +
  mu_b_inter_s2c1[3] * (chill1*siteSM) + mu_b_inter_ws2[3] * (hf*siteSM) +mu_b_inter_ps2[3] * (photo*siteSM) +
  mu_b_inter_s3c1[3] * (chill1*siteSM) + mu_b_inter_ws3[3] * (hf*siteSM) +mu_b_inter_ps3[3] * (photo*siteSM) +
  mu_b_inter_s4c1[3] * (chill1*siteSM) + mu_b_inter_ws4[3] * (hf*siteSM) +mu_b_inter_ps4[3] * (photo*siteSM)

# plot first for the low forcing
bb_lfc = a_sp[1]  + b_site2[1] * siteSM + b_site3[1] * siteSM + b_site4[1] * siteSM + mu_b_warm[1] * lf + mu_b_photo[1] * photo + mu_b_chill1[1] * chill1 +
  mu_b_inter_wp[1] * (lf*photo) +
  mu_b_inter_wc1[1] * (lf*chill1) + mu_b_inter_pc1[1] * (photo*chill1) +
  mu_b_inter_s2c1[1] * (chill1*siteSM) + mu_b_inter_ws2[1] * (lf*siteSM) +mu_b_inter_ps2[1] * (photo*siteSM) +
  mu_b_inter_s3c1[1] * (chill1*siteSM) + mu_b_inter_ws3[1] * (lf*siteSM) +mu_b_inter_ps3[1] * (photo*siteSM) +
  mu_b_inter_s4c1[1] * (chill1*siteSM) + mu_b_inter_ws4[1] * (lf*siteSM) +mu_b_inter_ps4[1] * (photo*siteSM)

bb_lfc5 = a_sp[2]  + b_site2[2] * siteSM + b_site3[2] * siteSM + b_site4[2] * siteSM + mu_b_warm[2] * lf + mu_b_photo[2] * photo + mu_b_chill1[2] * chill1 +
  mu_b_inter_wp[2] * (lf*photo) +
  mu_b_inter_wc1[2] * (lf*chill1) + mu_b_inter_pc1[2] * (photo*chill1) +
  mu_b_inter_s2c1[2] * (chill1*siteSM) + mu_b_inter_ws2[2] * (lf*siteSM) +mu_b_inter_ps2[2] * (photo*siteSM) +
  mu_b_inter_s3c1[2] * (chill1*siteSM) + mu_b_inter_ws3[2] * (lf*siteSM) +mu_b_inter_ps3[2] * (photo*siteSM) +
  mu_b_inter_s4c1[2] * (chill1*siteSM) + mu_b_inter_ws4[2] * (lf*siteSM) +mu_b_inter_ps4[2] * (photo*siteSM)

bb_lfc95 = a_sp[3]  + b_site2[3] * siteSM + b_site3[3] * siteSM + b_site4[3] * siteSM + mu_b_warm[3] * lf + mu_b_photo[3] * photo + mu_b_chill1[3] * chill1 +
  mu_b_inter_wp[3] * (lf*photo) +
  mu_b_inter_wc1[3] * (lf*chill1) + mu_b_inter_pc1[3] * (photo*chill1) +
  mu_b_inter_s2c1[3] * (chill1*siteSM) + mu_b_inter_ws2[3] * (lf*siteSM) +mu_b_inter_ps2[3] * (photo*siteSM) +
  mu_b_inter_s3c1[3] * (chill1*siteSM) + mu_b_inter_ws3[3] * (lf*siteSM) +mu_b_inter_ps3[3] * (photo*siteSM) +
  mu_b_inter_s4c1[3] * (chill1*siteSM) + mu_b_inter_ws4[3] * (lf*siteSM) +mu_b_inter_ps4[3] * (photo*siteSM)
#"#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"
# pdf("figures/chill_forcing_4sites_interactions_mini.pdf", width =12, height = 12)
# par(mfrow =c (2,2), mar = c(5.1, 5.1, 4.1, 2.1))
# plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Z-scored chill portions", ylab = "Day of budburst", cex.lab = 2)
# # points(hfData$chillport.z2, hfData$bb, bg = "#f9b641ff", pch =21, cex = 2.5)
# # points(lfData$chillport.z2, lfData$bb,  bg = "cyan4", pch = 21, cex = 2.5)
# abline(lm(bb_hfc ~ chill1), col = "#f9b641ff", lwd = 3)
# abline(lm(bb_lfc ~ chill1), col = "cyan4", lwd = 3)
# 
# legend("topright",legend = c(expression("low forcing"),
#                             expression("high forcing")),
#        col = c("cyan4","#f9b641ff"),
#        inset = 0.02, pch = c(19, 19 ),  cex = 1.25, bty = "n")
#dev.off()

intCF <- data.frame(bb_hfc = c(bb_hfc),bb_hfc5 = c(bb_hfc5),bb_hfc95 = c(bb_hfc95), bb_lfc = c(bb_lfc), bb_lfc5 = c(bb_lfc5), bb_lfc95 = c(bb_lfc95), chill = c(chill1))

intrxnCF <- ggplot(intCF, aes(x= chill, group =1)) +
  geom_line(aes(y = bb_hfc, col = "#2F5061"), size =1.5) +
  geom_line(aes(y = bb_lfc, col = "#4297A0"), size = 1.5) +
  xlim (-1.5,2) + 
  ylim (-2,50) +
  xlab("Chilling (standardized chill portions)") + ylab("Estimated day of budburst") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20))+
  theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
  #scale_fill_manual( labels = c("Low force", "High force")) +
  scale_color_manual(values = c("sienna4","tan2"), labels = c("High forcing", "Low forcing"), name = "") +
  #scale_colour_discrete(labels=c("High forcing","Low forcing"), name = "") +
  theme(legend.title = element_blank())# +  annotate("text", x = -5.4, y = 125, label = "a)", cex = 10) 
intrxnCF


intrxnCF <- ggplot(intCF, aes(x= chill, group =1)) +
  geom_line(aes(y = bb_hfc, col = "#2F5061"), size =1.5) +
  geom_line(aes(y = bb_lfc, col = "#4297A0"), size = 1.5) + 
   geom_ribbon(data = intCF, aes(ymin = bb_hfc5, ymax = bb_hfc95, x= chill), alpha = 0.2, fill = "sienna4") +
  geom_ribbon(data = intCF, aes(ymin = bb_lfc5, ymax = bb_lfc95, x= chill), alpha = 0.2, fill = "tan2") + 
  #scale_color_manual(values = c("cyan4", "red"), labels = c("High forcing", "Low forcing","High forcing", "Low forcing"), name = "") +
  xlab("Chilling (standardized chill portions) ") + ylab("Estimated day of budburst") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20))+
  theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
  #scale_fill_manual( labels = c("Low force", "High force")) +
  scale_color_manual(values = c("sienna4","tan2"), labels = c("High forcing", "Low forcing"), name = "") +
  #scale_colour_discrete(labels=c("High forcing","Low forcing"), name = "") +
  theme(legend.title = element_blank())# +  annotate("text", x = -5.4, y = 125, label = "a)", cex = 10) 


#"#E57F84","hotpink4"
# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# warm and site4
hfWest <- subset(hfData, population == "sm")
hfW <- unique(hfWest$force.z2)
lfWest <- subset(lfData, population == "sm")
lfW <- unique(lfWest$force.z2)

hfEast <- subset(hfData, population == "SH")
hfE <- unique(hfEast$force.z2)
lfEast <- subset(lfData, population == "SH")
lfE <- unique(lfEast$force.z2)

photo <- -0.5033863

site2 <- unique(pheno$site2.z2)
site3 <- unique(pheno$site3.z2)
site4 <- unique(pheno$site4.z2)

chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)

lcData <- subset(pheno, chill == "LC")
lcSmithers <- subset(lcData, population == "sm")
lcsm <- unique(lcSmithers$chillport.z2)

lcManning <- subset(lcData, population == "mp")
lcmp <- unique(lcManning$chillport.z2)

lcDataE <- subset(pheno, chill == "chill0")
lcStHipp <- subset(lcDataE, population == "SH")
lcsh <- unique(lcStHipp$chillport.z2)

lcHarv <- subset(lcDataE, population == "HF")
lchf <- unique(lcHarv$chillport.z2)

# St. Hippolyte plot first for the high forcing - site = 1 is the second value (site4[2])
bb_hfsite4 = a_sp[1:5]   + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[2]  + mu_b_warm[1:5] * hfE + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lcsh +
  mu_b_inter_wp[1:5] * (hfE*photo) +
  mu_b_inter_wc1[1:5] * (hfE*lcsh) + mu_b_inter_pc1[1:5] * (photo*lcsh) +
  mu_b_inter_s2c1[1:5] * (lcsh*site2[2]) + mu_b_inter_ws2[1:5] * (hfE*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lcsh*site3[1]) + mu_b_inter_ws3[1:5] * (hfE*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcsh*site4[2]) + mu_b_inter_ws4[1:5] * (hfE*site4[2]) +mu_b_inter_ps4[1:5] * (photo*site4[2])

# plot first for the low forcing
bb_lfsite4 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[2] + mu_b_warm[1:5] * lfE + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lcsh +
  mu_b_inter_wp[1:5] * (lfE*photo) +
  mu_b_inter_wc1[1:5] * (lfE*lcsh) + mu_b_inter_pc1[1:5] * (photo*lcsh) +
  mu_b_inter_s2c1[1:5] * (lcsh*site2[2]) + mu_b_inter_ws2[1:5] * (lfE*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lcsh*site3[1]) + mu_b_inter_ws3[1:5] * (lfE*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcsh*site4[2]) + mu_b_inter_ws4[1:5] * (lfE*site4[2]) +mu_b_inter_ps4[1:5] * (photo*site4[2])
#

# Manning park trends - site2[1]
bb_hfsite2 = a_sp[1:5] + b_site2[1:5] * site2[1] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * hfW + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lcmp +
  mu_b_inter_wp[1:5] * (hfW*photo) +
  mu_b_inter_wc1[1:5] * (hfW*lcmp) + mu_b_inter_pc1[1:5] * (photo*lcmp) +
  mu_b_inter_s2c1[1:5] * (lcmp*site2[1]) + mu_b_inter_ws2[1:5] * (hfW*site2[1]) +mu_b_inter_ps2[1:5] * (photo*site2[1]) +
  mu_b_inter_s3c1[1:5] * (lcmp*site3[1]) + mu_b_inter_ws3[1:5] * (hfW*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcmp*site4[1]) + mu_b_inter_ws4[1:5] * (hfW*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])

# plot first for the low forcing
bb_lfsite2 = a_sp[1:5] + b_site2[1:5] * site2[1] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1] + mu_b_warm[1:5] * lfW + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lcmp +
  mu_b_inter_wp[1:5] * (lfW*photo) +
  mu_b_inter_wc1[1:5] * (lfW*lcmp) + mu_b_inter_pc1[1:5] * (photo*lcmp) +
  mu_b_inter_s2c1[1:5] * (lcmp*site2[1]) + mu_b_inter_ws2[1:5] * (lfW*site2[1]) +mu_b_inter_ps2[1:5] * (photo*site2[1]) +
  mu_b_inter_s3c1[1:5] * (lcmp*site3[1]) + mu_b_inter_ws3[1:5] * (lfW*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcmp*site4[1]) + mu_b_inter_ws4[1:5] * (lfW*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])

# Smithers trends - site2[2]
bb_hfsite1 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * hfW + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lcsm +
  mu_b_inter_wp[1:5] * (hfW*photo) +
  mu_b_inter_wc1[1:5] * (hfW*lcsm) + mu_b_inter_pc1[1:5] * (photo*lcsm) +
  mu_b_inter_s2c1[1:5] * (lcsm*site2[2]) + mu_b_inter_ws2[1:5] * (hfW*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lcsm*site3[1]) + mu_b_inter_ws3[1:5] * (hfW*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcsm*site4[1]) + mu_b_inter_ws4[1:5] * (hfW*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])

# plot first for the low forcing
bb_lfsite1 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1] + mu_b_warm[1:5] * lfW + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lcsm +
  mu_b_inter_wp[1:5] * (lfW*photo) +
  mu_b_inter_wc1[1:5] * (lfW*lcsm) + mu_b_inter_pc1[1:5] * (photo*lcsm) +
  mu_b_inter_s2c1[1:5] * (lcsm*site2[2]) + mu_b_inter_ws2[1:5] * (lfW*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lcsm*site3[1]) + mu_b_inter_ws3[1:5] * (lfW*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcsm*site4[1]) + mu_b_inter_ws4[1:5] * (lfW*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])

# HarvardForest trends
bb_hfsite3 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[2] + b_site4[1:5] * site4[1] + mu_b_warm[1:5]* hfE + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lchf +
  mu_b_inter_wp[1:5] * (hfE*photo) +
  mu_b_inter_wc1[1:5] * (hfE*lchf) + mu_b_inter_pc1[1:5] * (photo*lchf) +
  mu_b_inter_s2c1[1:5] * (lchf*site2[2]) + mu_b_inter_ws2[1:5] * (hfE*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lchf*site3[2]) + mu_b_inter_ws3[1:5] * (hfE*site3[2]) +mu_b_inter_ps3[1:5] * (photo*site3[2]) +
  mu_b_inter_s4c1[1:5] * (lchf*site4[1]) + mu_b_inter_ws4[1:5] * (hfE*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])

# plot first for the low forcing
bb_lfsite3 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[2] + b_site4[1:5] * site4[1]+ mu_b_warm[1:5] * lfE + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lchf +
  mu_b_inter_wp[1:5] * (lfE*photo) +
  mu_b_inter_wc1[1:5]* (lfE*lchf) + mu_b_inter_pc1[1:5] * (photo*lchf) +
  mu_b_inter_s2c1[1:5] * (lchf*site2[2]) + mu_b_inter_ws2[1:5] * (lfE*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lchf*site3[2]) + mu_b_inter_ws3[1:5] * (lfE*site3[2]) +mu_b_inter_ps3[1:5] * (photo*site3[2]) +
  mu_b_inter_s4c1[1:5] * (lchf*site4[1]) + mu_b_inter_ws4[1:5] * (lfE*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])


temp <- rbind(bb_lfsite1,bb_lfsite2,bb_lfsite3,bb_lfsite4,bb_hfsite1,bb_hfsite2,bb_hfsite3,bb_hfsite4)
siteForce <- data.frame(temp, site = c("Smithers","Manning park", "Harvard forest","St.Hippolyte","Smithers","Manning park", "Harvard forest","St.Hippolyte"), force = c("low forcing","low forcing","low forcing","low forcing", "high forcing","high forcing", "high forcing", "high forcing"))
names(siteForce) <- c("mean", "lower1", "upper1","lower2","upper2","site","force")


siteForce <- siteForce[order(siteForce$site),]
siteForce$temp <- rownames(siteForce)
siteForce$tempSite <- paste(siteForce$force, siteForce$site, sep = " ")
siteOrder <- c("Smithers", "Manning park", "St.Hippolyte", "Harvard forest")

siteFPoint <- ggplot()  +
  geom_pointrange(siteForce, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=lower1, ymax=upper1, col = temp, shape = tempSite), size =0.55, position=position_dodge(width=0.5)) +
  geom_linerange( siteForce, mapping = aes(x = factor(site, level = siteOrder), ymin = lower2, ymax=upper2, col = temp), size =1.25, position=position_dodge(width=0.5)) +
  xlab("Population") + ylab("Estimated day of budburst") +
  ylim(0,80) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20),axis.text.x = element_text( size= 15,angle = 55, hjust=1)) +
  # scale_color_manual(values = c("cyan4", "#CC6677"), labels = c("High forcing", "Low forcing"), name = "") +
  scale_color_manual(values = c(
    "bb_hfsite1"="forestgreen",
    "bb_lfsite1"="palegreen3",
    "bb_hfsite2"= "forestgreen",
    "bb_lfsite2"="palegreen3",
    "bb_hfsite4"="darkorchid4",
    "bb_lfsite4"="darkorchid1",
    "bb_hfsite3" ="darkorchid4",
    "bb_lfsite3"="darkorchid1"
  ), labels = c("Smithers - High forcing", "Smithers - Low forcing","Manning Park - High forcing", "Manning Park - Low forcing", "St. Hippolyte - High forcing", "St. Hippolyte - Low forcing","Harvard Forest - High forcing", "Harvard Forest - Low forcing"), name = "") +theme(legend.position = "none") +  annotate("text", x = 1.15, y = 75, label = "a) Forcing", cex = 9) + scale_shape_manual(values = c("low forcing Smithers" = 0,  "high forcing Smithers" = 15,  "low forcing Manning park" = 0, "high forcing Manning park" = 15,  "low forcing Harvard forest" = 2, "high forcing Harvard forest" = 17,  "high forcing St.Hippolyte" = 17, "low forcing St.Hippolyte" = 2 ), breaks = c("low forcing Smithers", "high forcing Manning park"), label = c("Western", "Eastern")) +
scale_x_discrete(breaks = c("Smithers", "Manning park", "St.Hippolyte", "Harvard forest"),
                       #limits = c(0,100),
                       labels = c(expression("Smithers (54.8"*~degree*N*")"), 
                                  expression("Manning park (49.1"*~degree*N*")"),
                                  expression("St.Hippolyte (45.9"*~degree*N*")"),
                                  expression("Harvard forest (42.5"*~degree*N*")")
                       )
  )
siteFPoint
  # theme(legend.key=element_blank(), legend.position=c(.83,.85),legend.text = element_text(size = 15)) +
  # theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 50, label = "b)", cex = 10) 

# siteForce$force_pop <- paste(siteForce$force, siteForce$site, sep = "_")
# siteForce$col <- "black"
# siteForce$col[which(siteForce$force_pop == "18")] <- 1978
# 
# 
# ggplot() +
#   geom_pointrange(siteForce, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=X2.5., ymax=X97.5., col = force_pop),
#                   position=position_dodge(width=0.5), size =1.25) +
#   xlab("Population") + ylab("Estimated day of budburst") +
#   ylim(0,50) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
#   scale_color_manual(values = c( "#D48C70", 
#                                 "#29A0B1",
#                                 "#98D7C2",
#                                 "#167D7F",
#                                 "#116530",
#                                 
#                                 "#C85250",
#                                 "#EBBBB0", 
#                                 "#E98980"), labels = c("High forcing", "Low forcing"), name = "") +
#   theme(legend.position = "none") + 
#   #theme(legend.key=element_blank(), legend.position=c(.83,.85),legend.text = element_text(size = 15)) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 50, label = "b)", cex = 10) 

#98D7C2, #29A0B1
#167D7F
#116530
#18A558

#D48C70 #C85250
#EBBBB0 #E98980

### Site 4 and chill:
hcData <- subset(pheno, chill == "HC" )
lcData <- subset(pheno, chill == "LC")
c1Data <- subset(pheno, chill == "chill1" )
c2Data <- subset(pheno, chill == "chill2")
c0Data <- subset(pheno, chill == "chill0" )


hc <- unique(hcData$chillport.z2)
lc <- unique(lcData$chillport.z2)
# c2 <- unique(c2Data$chillport.z2)
c1 <- unique(c1Data$chillport.z2)
c0 <- unique(c0Data$chillport.z2)

lcData <- subset(pheno, chill == "LC")
lcSmithers <- subset(lcData, population == "sm")
lcsm <- unique(lcSmithers$chillport.z2)

lcManning <- subset(lcData, population == "mp")
lcmp <- unique(lcManning$chillport.z2)

lcDataE <- subset(pheno, chill == "chill0")
lcStHipp <- subset(lcDataE, population == "SH")
lcsh <- unique(lcStHipp$chillport.z2)

lcHarv <- subset(lcDataE, population == "HF")
lchf <- unique(lcHarv$chillport.z2)

# High chill #############
hcData <- subset(pheno, chill == "HC")
hcSmithers <- subset(hcData, population == "sm")
hcsm <- unique(hcSmithers$chillport.z2)

hcManning <- subset(hcData, population == "mp")
hcmp <- unique(hcManning$chillport.z2)

hcDataE <- subset(pheno, chill == "chill1")
hcStHipp <- subset(hcDataE, population == "SH")
hcsh <- unique(hcStHipp$chillport.z2)

hcHarv <- subset(hcDataE, population == "HF")
hchf <- unique(hcHarv$chillport.z2)

photo <- -0.5033863

site4 <- unique(pheno$site4.z2)
site3 <- unique(pheno$site3.z2)
site2 <- unique(pheno$site2.z2)

#Site 4
# plot first for the high chill
# bb_hc2site4 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[2]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * c2[2] +
#   mu_b_inter_wp[1:5] * (force*photo) +
#   mu_b_inter_wc1[1:5] * (force*c2[2]) + mu_b_inter_pc1[1:5] * (photo*c2[2]) +
#   mu_b_inter_s2c1[1:5] * (c2[2]*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
#   mu_b_inter_s3c1[1:5] * (c2[2]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
#   mu_b_inter_s4c1[1:5] * (c2[2]*site4[2]) + mu_b_inter_ws4[1:5] * (force*site4[2]) +mu_b_inter_ps4[1:5] * (photo*site4[2])
force <- lfE

bb_hc1site4 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[2]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * c1[2] +
  mu_b_inter_wp[1:5] * (force*photo) +
  mu_b_inter_wc1[1:5] * (force*c1[2]) + mu_b_inter_pc1[1:5] * (photo*c1[2]) +
  mu_b_inter_s2c1[1:5] * (c1[2]*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (c1[2]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (c1[2]*site4[2]) + mu_b_inter_ws4[1:5] * (force*site4[2]) +mu_b_inter_ps4[1:5] * (photo*site4[2])


bb_lc0site4 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[2]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * c0[2] +
  mu_b_inter_wp[1:5] * (force*photo) +
  mu_b_inter_wc1[1:5] * (force*c0[2]) + mu_b_inter_pc1[1:5] * (photo*c0[2]) +
  mu_b_inter_s2c1[1:5] * (c0[2]*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (c0[2]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (c0[2]*site4[2]) + mu_b_inter_ws4[1:5] * (force*site4[2]) +mu_b_inter_ps4[1:5] * (photo*site4[2])

# site 3: Harvard forest

# bb_hc2site3 = a_sp[1:5]  + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[2] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * c2[1] +
#   mu_b_inter_wp[1:5] * (force*photo) +
#   mu_b_inter_wc1[1:5] * (force*c2[1]) + mu_b_inter_pc1[1:5] * (photo*c2[1]) +
#   mu_b_inter_s2c1[1:5] * (c2[1]*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
#   mu_b_inter_s3c1[1:5] * (c2[1]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
#   mu_b_inter_s4c1[1:5] * (c2[1]*site4[2]) + mu_b_inter_ws4[1:5] * (force*site4[2]) +mu_b_inter_ps4[1:5] * (photo*site4[2])

bb_hc1site3 = a_sp[1:5]  + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[2] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * c1[1] +
  mu_b_inter_wp[1:5] * (force*photo) +
  mu_b_inter_wc1[1:5] * (force*c1[1]) + mu_b_inter_pc1[1:5] * (photo*c1[1]) +
  mu_b_inter_s2c1[1:5] * (c1[1]*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (c1[1]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (c1[1]*site4[2]) + mu_b_inter_ws4[1:5] * (force*site4[2]) +mu_b_inter_ps4[1:5] * (photo*site4[2])


bb_lc0site3 = a_sp[1:5]  + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[2] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * c0[1] +
  mu_b_inter_wp[1:5] * (force*photo) +
  mu_b_inter_wc1[1:5] * (force*c0[1]) + mu_b_inter_pc1[1:5] * (photo*c0[1]) +
  mu_b_inter_s2c1[1:5] * (c0[1]*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) + mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (c0[1]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) + mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (c0[1]*site4[2]) + mu_b_inter_ws4[1:5] * (force*site4[2]) + mu_b_inter_ps4[1:5] * (photo*site4[2])

# now for the western sites:

#Site 2 Manning park
force <- lfW

# plot first for the high chill
bb_hcsite2 = a_sp[1:5]  + b_site2[1:5] * site2[1] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * hc[1] +
  mu_b_inter_wp[1:5] * (force*photo) +
  mu_b_inter_wc1[1:5]* (force*hc[1]) + mu_b_inter_pc1[1:5] * (photo*hc[1]) +
  mu_b_inter_s2c1[1:5] * (hc[1]*site2[1]) + mu_b_inter_ws2[1:5] * (force*site2[1]) +mu_b_inter_ps2[1:5] * (photo*site2[1]) +
  mu_b_inter_s3c1[1:5] * (hc[1]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (hc[1]*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])

bb_lcsite2 = a_sp[1:5] + b_site2[1:5] * site2[1] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lc[1] +
  mu_b_inter_wp[1:5] * (force*photo) +
  mu_b_inter_wc1[1:5] * (force*lc[1]) + mu_b_inter_pc1[1:5] * (photo*lc[1]) +
  mu_b_inter_s2c1[1:5] * (lc[1]*site2[1]) + mu_b_inter_ws2[1:5] * (force*site2[1]) +mu_b_inter_ps2[1:5] * (photo*site2[1]) +
  mu_b_inter_s3c1[1:5] * (lc[1]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lc[1]*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])

#smithers
bb_hcsite1 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * hc[2] +
  mu_b_inter_wp[1:5] * (force*photo) +
  mu_b_inter_wc1[1:5] * (force*hc[2]) + mu_b_inter_pc1[1:5] * (photo*hc[2]) +
  mu_b_inter_s2c1[1:5] * (hc[2]*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) + mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (hc[2]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) + mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (hc[2]*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) + mu_b_inter_ps4[1:5] * (photo*site4[1])

bb_lcsite1 = a_sp[1:5]+ b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * photo + mu_b_chill1[1:5] * lc[2] +
  mu_b_inter_wp[1:5] * (force*photo) +
  mu_b_inter_wc1[1:5] * (force*lc[2]) + mu_b_inter_pc1[1:5] * (photo*lc[2]) +
  mu_b_inter_s2c1[1:5] * (lc[2]*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (photo*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lc[2]*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (photo*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lc[2]*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) +mu_b_inter_ps4[1:5] * (photo*site4[1])

tempChill <- rbind(bb_hc1site4,bb_lc0site4,bb_hc1site3,bb_lc0site3,bb_hcsite2,bb_lcsite2,bb_hcsite1,bb_lcsite1) #bb_hc2site4,bb_hc2site3,)
siteChill <- data.frame(tempChill , chill = c("Eastern - High Chill", "Eastern - Low Chill","Eastern - High Chill", "Eastern - Low Chill","Western - High Chill","Western - Low Chill", "Western - High Chill","Western - Low Chill"), site = c("St.Hippolyte","St.Hippolyte", "Harvard forest", "Harvard forest","Manning park","Manning park","Smithers","Smithers"), col = c("darkorchid4","darkorchid3","darkforestgreen","springgreen4","darkred","indianred3","deepskyblue3","deepskyblue1"))
                        #col = c("#593d9cff", "cyan4","#593d9cff", "cyan4","#CC6677","#f9b641ff","#CC6677","#f9b641ff"))
names(siteChill) <- c("mean", "lower1", "upper1", "lower2", "upper2","chill","site","col")
siteChill <- siteChill[order(siteChill$site),]
siteChill$temp <- rownames(siteChill)

siteOrder <- c("Smithers", "Manning park", "Harvard forest", "St.Hippolyte")

siteChill$transect <- siteChill$site
siteChill$transect[siteChill$transect == "Smithers"] <- "Western"
siteChill$transect[siteChill$transect == "Manning park"] <- "Western"
siteChill$transect[siteChill$transect == "Harvard forest"] <- "Eastern"
siteChill$transect[siteChill$transect == "St.Hippolyte"] <- "Eastern"

siteChill$tranSite <- paste(siteChill$transect, siteChill$site, sep = " ")
siteChill$tempSite <- paste(siteChill$chill, siteChill$site, sep = " ")

siteCPoint <- ggplot() +
  geom_pointrange(siteChill, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=lower1, ymax=upper1, col = temp, shape = tempSite),
                  position=position_dodge(width=0.5), size = 0.55) +
  geom_linerange( siteChill, mapping = aes(x = factor(site, level = siteOrder), ymin = lower2, ymax=upper2, col = temp), size =1.25, position=position_dodge(width=0.5)) + 
  ylim(0,80) +
  xlab("Population") + ylab("Estimated day of budburst") +
  #ylim(0,50) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20),axis.text.x = element_text( size= 15,angle = 55, hjust=1), legend.position = "none") + scale_color_manual(values = c(
    "bb_hcsite1"="forestgreen",
    "bb_lcsite1"="palegreen3",
    "bb_hcsite2"= "forestgreen",
    "bb_lcsite2"="palegreen3",
    "bb_hc1site4"="darkorchid4",
    "bb_lc0site4"="darkorchid1",
    "bb_hc1site3" ="darkorchid4",
    "bb_lc0site3"="darkorchid1"
  ), labels = c("Smithers - High  Chill", "Smithers - Low  Chill","Manning Park - High  Chill", "Manning Park - Low  Chill", "St. Hippolyte - High  Chill", "St. Hippolyte - Low  Chill","Harvard Forest - High  Chill", "Harvard Forest - Low Chill")) +  annotate("text", x = 1.15, y = 75, label = "b) Chilling", cex = 9) + scale_shape_manual(values = c("Western - Low Chill Smithers" = 0,  "Western - High Chill Smithers" = 15,  "Western - Low Chill Manning park" = 0, "Western - High Chill Manning park" = 15,  "Eastern - Low Chill Harvard forest" = 2, "Eastern - High Chill Harvard forest" = 17,  "Eastern - High Chill St.Hippolyte" = 17, "Eastern - Low Chill St.Hippolyte" = 2 ), breaks = c("Western - Low Chill Smithers", "Western - High Chill Manning park"), label = c("Western", "Eastern"))+
  scale_x_discrete(breaks = c("Smithers", "Manning park", "St.Hippolyte", "Harvard forest"),
                   #limits = c(0,100),
                   labels = c(expression("Smithers (54.8"*~degree*N*")"), 
                              expression("Manning park (49.1"*~degree*N*")"),
                              expression("St.Hippolyte (45.9"*~degree*N*")"),
                              expression("Harvard forest (42.5"*~degree*N*")")
                   )
  )
siteCPoint
  # ggplot() +
  # geom_pointrange(siteChill, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=lower, ymax=upper, col = col),
  #                 position=position_dodge(width=0.5), size = 1.25) +
  # #ylim (0,50) +
  # xlab("Population") + ylab("Estimated day of budburst") +
  # #ylim(0,50) +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"),
  #       axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
  # scale_color_manual(values = c("darkorchid4",
  #                               "darkorchid3",
  #                               "forestgreen",
  #                               "springgreen4",
  #                               "darkred",
  #                               "indianred3",
  #                               "deepskyblue3",
  #                               "deepskyblue1"), labels = c("High chill - East",
  #                                                        "Low chill - East",
  #                                                       # "No chill - East",
  #                                                        "High chill - West",
  #                                                        "Low chill - West"), name = "") +
  # theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
  # theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 50, label = "c)", cex = 10) 


# cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# remove the low chill:

# tempChillNoHigh <- rbind(bb_hc2site4,bb_lc0site4,bb_hc2site3,bb_lc0site3,bb_hcsite2,bb_lcsite2,bb_hcsite1,bb_lcsite1)
# siteChillNoHigh <- data.frame(tempChillNoHigh , chill = c("Eastern - High Chill","Eastern - No Chill","Eastern - High Chill", "Eastern - No Chill", "Western - High Chill","Western - Low Chill", "Western - High Chill","Western - Low Chill"), site = c("St.Hippolyte","St.Hippolyte", "Harvard forest",  "Harvard forest","Manning park","Manning park","Smithers","Smithers"), col = c("#593d9cff", "cyan4","#593d9cff", "cyan4","#CC6677","#f9b641ff","#CC6677","#f9b641ff"))
# 
# siteChillNoHigh <- siteChillNoHigh[order(siteChillNoHigh$site),]
# 
# siteOrder <- c("Smithers", "Manning park", "Harvard forest", "St.Hippolyte")
# 
# siteCPointNH <- ggplot() +
#   geom_pointrange(siteChillNoHigh, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=X2.5., ymax=X97.5., col = chill),
#                   position=position_dodge(width=0.5), size = 1.25) +
#   xlab("Population") + ylab("Estimated day of budburst") +
#   ylim(0,50) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
#   scale_color_manual(values = c("cyan4", 
#                                 "orchid4", 
#                                 "#CC6677",
#                                 "#f9b641ff"), labels = c("High chill - East",
#                                                          "Low chill - East",
#                                                          "High chill - West",
#                                                          "Low chill - West"), name = "") +
#   theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 50, label = "c)", cex = 10) 


# plot(hcData$site4.z2, hcData$bb, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Site", ylab = "Day of budburst", cex.lab = 2)
# # points(jitter(hcData$site4.z2, 10), hcData$bb, bg = "#f9b641ff", pch =21, cex = 2.5)
# # points(jitter(lcData$site4.z2, 8), lcData$bb, bg = "cyan4", pch = 21, cex = 2.5)
# # points((c0Data$site4.z2), c0Data$bb, bg = "#a65c85ff", pch =21, cex = 2.5)
# # points(c1Data$site4.z2 ,c1Data$bb, bg = "#7e4e90ff", pch =21, cex = 2.5)
# # points(c2Data$site4.z2, c2Data$bb, bg = "#cc6a70ff", pch =21, cex = 2.5)
# 
# # abline(lm(bb_hcsite4 ~ site4), col = "#f9b641ff", lwd = 3)
# # abline(lm(bb_lcsite4 ~ site4), col = "cyan4", lwd = 3)
# abline(lm(bb_c0site4 ~ site4), col = "cyan4", lwd = 3, lty = 1)
# abline(lm(bb_c1site4 ~ site4), col = "cyan4", lwd = 3, lty = 2)
# abline(lm(bb_c2site4 ~ site4), col = "cyan4", lwd = 3, lty = 3)
# 
# abline(lm(bb_c0site3 ~ site3), col = "#f9b641ff", lwd = 3, lty = 1)
# abline(lm(bb_c1site3 ~ site3), col = "#f9b641ff", lwd = 3, lty = 2)
# abline(lm(bb_c2site3 ~ site3), col = "#f9b641ff", lwd = 3, lty = 3)
# 
# abline(lm(bb_hcsite2 ~ site2), col = "#a65c85ff", lwd = 3, lty = 6)
# abline(lm(bb_lcsite2 ~ site2), col = "#a65c85ff", lwd = 3, lty = 1)
# 
# legend("topleft",legend = c(expression("St. Hippolyte"),
#                             expression("Harvard Forest"),
#                             expression("Manning Park"),
#                             expression("c0"),
#                             expression("c1"),
#                             expression("c2"),
#                              expression("hc"),
#                              expression("lc")),
#        col = c("cyan4","#f9b641ff","#a65c85ff", "black", "black", "black", "black", "black" ),
#        inset = 0.02, lty = c(1,1,1,1,2, 1,6 ), lwd =2,  cex = 1.1, bty = "n")


### Photo and site 4
hpData <- subset(pheno, photo == "HP" )
lpData <- subset(pheno, photo == "LP")

hp <- unique(hpData$photo.z2)
lp <- unique(lpData$photo.z2)

site4 <- unique(pheno$site4.z2)
site3 <- unique(pheno$site3.z2)
site2 <- unique(pheno$site2.z2)

#Site 4
# plot first for the high chill
force <- hfE
bb_hpsite4 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[2]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * hp + mu_b_chill1[1:5] *lcsh +
  mu_b_inter_wp[1:5] * (force*hp) +
  mu_b_inter_wc1[1:5] * (force*lcsh) + mu_b_inter_pc1[1:5] * (hp*lcsh) +
  mu_b_inter_s2c1[1:5] * (lcsh*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (hp*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lcsh*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (hp*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcsh*site4[2]) + mu_b_inter_ws4[1:5] * (force*site4[2]) +mu_b_inter_ps4[1:5] * (hp*site4[2])

bb_lpsite4 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[2]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * lp + mu_b_chill1[1:5] *lcsh +
  mu_b_inter_wp[1:5] * (force*lp) +
  mu_b_inter_wc1[1:5] * (force*lcsh) + mu_b_inter_pc1[1:5] * (lp*lcsh) +
  mu_b_inter_s2c1[1:5] * (lcsh*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) + mu_b_inter_ps2[1:5] * (lp*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lcsh*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) + mu_b_inter_ps3[1:5] * (lp*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcsh*site4[2]) + mu_b_inter_ws4[1:5] * (force*site4[2]) + mu_b_inter_ps4[1:5] * (lp*site4[2])

# harvard forest
bb_hpsite3 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[2] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * hp + mu_b_chill1[1:5] *lchf +
  mu_b_inter_wp[1:5] * (force*hp) +
  mu_b_inter_wc1[1:5] * (force*lchf) + mu_b_inter_pc1[1:5] * (hp*lchf) +
  mu_b_inter_s2c1[1:5] * (lchf*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (hp*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lchf*site3[2]) + mu_b_inter_ws3[1:5] * (force*site3[2]) +mu_b_inter_ps3[1:5] * (hp*site3[2]) +
  mu_b_inter_s4c1[1:5] * (lchf*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) +mu_b_inter_ps4[1:5] * (hp*site4[1])

bb_lpsite3 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[2] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * lp + mu_b_chill1[1:5] *lchf +
  mu_b_inter_wp[1:5] * (force*lp) +
  mu_b_inter_wc1[1:5] * (force*lchf) + mu_b_inter_pc1[1:5] * (lp*lchf) +
  mu_b_inter_s2c1[1:5] * (lchf*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (lp*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lchf*site3[2]) + mu_b_inter_ws3[1:5] * (force*site3[2]) +mu_b_inter_ps3[1:5] * (lp*site3[2]) +
  mu_b_inter_s4c1[1:5] * (lchf*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) +mu_b_inter_ps4[1:5] * (lp*site4[1])

# Manning park
force <- lfW
bb_hpsite2 = a_sp[1:5] + b_site2[1:5] * site2[1] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * hp + mu_b_chill1[1:5] *lcmp +
  mu_b_inter_wp[1:5] * (force*hp) +
  mu_b_inter_wc1[1:5] * (force*lcmp) + mu_b_inter_pc1[1:5] * (hp*lcmp) +
  mu_b_inter_s2c1[1:5] * (lcmp*site2[1]) + mu_b_inter_ws2[1:5] * (force*site2[1]) + mu_b_inter_ps2[1:5] * (hp*site2[1]) +
  mu_b_inter_s3c1[1:5] * (lcmp*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) + mu_b_inter_ps3[1:5] * (hp*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcmp*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) + mu_b_inter_ps4[1:5] * (hp*site4[1])

bb_lpsite2 = a_sp[1:5] + b_site2[1:5] * site2[1] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * lp + mu_b_chill1[1:5] *lcmp +
  mu_b_inter_wp[1:5] * (force*lp) +
  mu_b_inter_wc1[1:5] * (force*lcmp) + mu_b_inter_pc1[1:5] * (lp*lcmp) +
  mu_b_inter_s2c1[1:5] * (lcmp*site2[1]) + mu_b_inter_ws2[1:5] * (force*site2[1]) + mu_b_inter_ps2[1:5] * (lp*site2[1]) +
  mu_b_inter_s3c1[1:5] * (lcmp*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) + mu_b_inter_ps3[1:5] * (lp*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcmp*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) + mu_b_inter_ps4[1:5] * (lp*site4[1])

# Smithers
bb_hpsite1 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * hp + mu_b_chill1[1:5] *lcsm +
  mu_b_inter_wp[1:5] * (force*hp) +
  mu_b_inter_wc1[1:5] * (force*lcsm) + mu_b_inter_pc1[1:5] * (hp*lcsm) +
  mu_b_inter_s2c1[1:5] * (lcsm*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (hp*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lcsm*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (hp*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcsm*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) +mu_b_inter_ps4[1:5] * (hp*site4[1])

bb_lpsite1 = a_sp[1:5] + b_site2[1:5] * site2[2] + b_site3[1:5] * site3[1] + b_site4[1:5] * site4[1]  + mu_b_warm[1:5] * force + mu_b_photo[1:5] * lp + mu_b_chill1[1:5] *lcsm +
  mu_b_inter_wp[1:5] * (force*lp) +
  mu_b_inter_wc1[1:5] * (force*lcsm) + mu_b_inter_pc1[1:5] * (lp*lcsm) +
  mu_b_inter_s2c1[1:5] * (lcsm*site2[2]) + mu_b_inter_ws2[1:5] * (force*site2[2]) +mu_b_inter_ps2[1:5] * (lp*site2[2]) +
  mu_b_inter_s3c1[1:5] * (lcsm*site3[1]) + mu_b_inter_ws3[1:5] * (force*site3[1]) +mu_b_inter_ps3[1:5] * (lp*site3[1]) +
  mu_b_inter_s4c1[1:5] * (lcsm*site4[1]) + mu_b_inter_ws4[1:5] * (force*site4[1]) +mu_b_inter_ps4[1:5] * (lp*site4[1])

tempPhoto <- rbind(bb_hpsite4,bb_hpsite3,bb_hpsite2,bb_hpsite1,bb_lpsite4,bb_lpsite3,bb_lpsite2,bb_lpsite1)

sitePhoto <- data.frame(tempPhoto, photo = c("High photoperiod","High photoperiod","High photoperiod","High photoperiod","Low photoperiod","Low photoperiod","Low photoperiod","Low photoperiod"), site = c("St.Hippolyte", "Harvard forest","Manning park","Smithers","St.Hippolyte", "Harvard forest","Manning park","Smithers"))
names(sitePhoto) <- c("mean", "lower1", "upper1", "lower2", "upper2","photo","site")

siteOrder <- c("Smithers","Manning park", "St.Hippolyte", "Harvard forest")
sitePhoto <- sitePhoto[order(sitePhoto$site),]
sitePhoto$temp <- rownames(sitePhoto)

sitePhoto$transect <- sitePhoto$site
sitePhoto$transect[sitePhoto$transect == "Smithers"] <- "Western"
sitePhoto$transect[sitePhoto$transect == "Manning park"] <- "Western"
sitePhoto$transect[sitePhoto$transect == "Harvard forest"] <- "Eastern"
sitePhoto$transect[sitePhoto$transect == "St.Hippolyte"] <- "Eastern"

sitePhoto$tranSite <- paste(sitePhoto$transect, sitePhoto$site, sep = " ")
sitePhoto$tempSite <- paste(sitePhoto$photo, sitePhoto$site, sep = " ")

sitePPoint <- ggplot() +
  geom_pointrange(sitePhoto, mapping = 
                    aes(x = factor(site, level = siteOrder), y = mean, ymin=lower1, ymax=upper1, 
                        col = temp, shape = tempSite),
                  position=position_dodge(width=0.5), size =0.55) + 
  ylim(0,80) +
  geom_linerange( sitePhoto, mapping = aes(x = factor(site, level = siteOrder), ymin = lower2, ymax=upper2, col = temp), size =1.25, position=position_dodge(width=0.5))+
  xlab("Population") + ylab("Estimated day of budburst")+
  scale_color_manual(values = c("bb_hpsite1"="forestgreen",
                                "bb_lpsite1"="palegreen3",
                                "bb_hpsite2"= "forestgreen",
                                "bb_lpsite2"="palegreen3",
                                "bb_hpsite3" ="darkorchid4",
                                "bb_lpsite3"="darkorchid1",
                                "bb_hpsite4"="darkorchid4",
                                "bb_lpsite4"="darkorchid1"),
                     breaks = c("bb_hpsite1","bb_lpsite1"), label = c("High", "Low"), name = "") +
  scale_shape_manual(values = c("Low photoperiod Smithers" = 0,
                                "High photoperiod Smithers" = 15,
                                "Low photoperiod Manning park" = 0,
                                "High photoperiod Manning park" = 15,
                                "Low photoperiod Harvard forest" = 2,
                                "High photoperiod Harvard forest" = 17,
                                "High photoperiod St.Hippolyte" = 17,
                                "Low photoperiod St.Hippolyte" = 2),
  breaks = c("High photoperiod Manning park", "High photoperiod St.Hippolyte"), label = c("Western", "Eastern"), name = "") +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 15), axis.title = element_text(size = 20), axis.text.x = element_text( size= 15,angle = 55, hjust=1), legend.key=element_rect(fill="white"),legend.text=element_text(size=25))+  annotate("text", x = 1.55, y = 75, label = "c) Photoperiod", cex = 9) +
  scale_x_discrete(breaks = c("Smithers", "Manning park", "St.Hippolyte", "Harvard forest"),
                   #limits = c(0,100),
                   labels = c(expression("Smithers (54.8"*~degree*N*")"), 
                              expression("Manning park (49.1"*~degree*N*")"),
                              expression("St.Hippolyte (45.9"*~degree*N*")"),
                              expression("Harvard forest (42.5"*~degree*N*")")
                   )
  )
sitePPoint

# sitePPointNL <- ggplot() +
#   geom_pointrange(sitePhoto, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=lower1, ymax=upper1, col = temp),
#                   position=position_dodge(width=0.5), size =.55) +
#   geom_linerange( sitePhoto, mapping = aes(x = factor(site, level = siteOrder), ymin = lower2, ymax=upper2, col = temp), size =1.25, position=position_dodge(width=0.5))+
#   xlab("Population") + ylab("Estimated day of budburst") +
#    ylim(0,80) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         axis.text = element_text(size = 15), axis.title = element_text(size = 20), axis.text.x = element_text( size= 15,angle = 55, hjust=1))+
#   scale_color_manual(values = c(
#     "bb_hpsite1"="deepskyblue3",
#     "bb_lpsite1"="deepskyblue1",
#     "bb_hpsite2"= "forestgreen",
#     "bb_lpsite2"="palegreen3",
#     "bb_hpsite3" ="darkred",
#     "bb_lpsite3"="tomato1",
#     "bb_hpsite4"="darkorchid4",
#     "bb_lpsite4"="darkorchid1"), labels = c("Smithers - High cue", "Smithers - Low cue","Manning Park - High cue", "Manning Park - Low cue", "St. Hippolyte - High cue", "St. Hippolyte - Low cue","Harvard Forest - High cue", "Harvard Forest - Low cue")) + theme(legend.position = "none")  +  annotate("text", x = 1.75, y = 80, label = "c) Photoperiod", cex = 10) 

  # theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
  # theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 50, label = "d)", cex = 10) 


# pdf("figures/intrxnPlotsApril19.pdf", height =12, width = 15)
# plot_grid(intrxnCF, siteFPoint, siteCPoint, sitePPoint, nrow = 2, ncol = 2, align = "v")
# dev.off()

l <- get_legend(sitePPoint) 

pdf("figures/intrxnForceChill.pdf", height = 2, width = 2)
intrxnCF
dev.off()

pdf("figures/intrxnPlots8Chill3Lat.pdf", height =8, width = 20)
plot_grid( siteFPoint, siteCPoint, sitePPoint , ncol = 3, nrow =1,align = "v", rel_widths = c(1,1,1.25))
dev.off()

pdf("figures/intrxnPlots8Chill3.pdf", height =5, width = 20)
plot_grid( siteFPoint, siteCPoint, sitePPointNL, l , ncol = 4, nrow =1,align = "v")
dev.off()


##################################

# not a bar plot - make it points and lines 


# plot first for the high forcing
# bb_hpsite4 = a_sp + b_site2 * 0 + b_site3 * 0 + b_site4 * site4  + mu_b_warm * hp + mu_b_photo * photo + mu_b_chill1 * chill1 +
#   mu_b_inter_wp * (hp*photo) +
#   mu_b_inter_wc1 * (hp*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site4) + mu_b_inter_ws2 * (hp*site4) +mu_b_inter_ps2 * (photo*site4) +
#   mu_b_inter_s3c1 * (chill1*site4) + mu_b_inter_ws3 * (hp*site4) +mu_b_inter_ps3 * (photo*site4) +
#   mu_b_inter_s4c1 * (chill1*site4) + mu_b_inter_ws4 * (hp*site4) +mu_b_inter_ps4 * (photo*site4)
# 
# # plot first for the low forcing
# bb_lpsite4 = a_sp + b_site2 * 0 + b_site3 * 0 + b_site4 * site4 + mu_b_warm * lp + mu_b_photo * photo + mu_b_chill1 * chill1 +
#   mu_b_inter_wp * (lp*photo) +
#   mu_b_inter_wc1 * (lp*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site4) + mu_b_inter_ws2 * (lp*site4) +mu_b_inter_ps2 * (photo*site4) +
#   mu_b_inter_s3c1 * (chill1*site4) + mu_b_inter_ws3 * (lp*site4) +mu_b_inter_ps3 * (photo*site4) +
#   mu_b_inter_s4c1 * (chill1*site4) + mu_b_inter_ws4 * (lp*site4) +mu_b_inter_ps4 * (photo*site4)
# #
# # site 3
# bb_hpsite3 = a_sp + b_site2 * 0 + b_site3 * site3 + b_site4 * 0  + mu_b_warm * hp + mu_b_photo * photo + mu_b_chill1 * chill1 +
#   mu_b_inter_wp * (hp*photo) +
#   mu_b_inter_wc1 * (hp*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site3) + mu_b_inter_ws2 * (hp*site3) +mu_b_inter_ps2 * (photo*site3) +
#   mu_b_inter_s3c1 * (chill1*site3) + mu_b_inter_ws3 * (hp*site3) +mu_b_inter_ps3 * (photo*site3) +
#   mu_b_inter_s4c1 * (chill1*site3) + mu_b_inter_ws4 * (hp*site3) +mu_b_inter_ps4 * (photo*site3)
# 
# # plot first for the low forcing
# bb_lpsite3 = a_sp + b_site2 * 0 + b_site3 * site3 + b_site4 * 0 + mu_b_warm * lp + mu_b_photo * photo + mu_b_chill1 * chill1 +
#   mu_b_inter_wp * (lp*photo) +
#   mu_b_inter_wc1 * (lp*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site3) + mu_b_inter_ws2 * (lp*site3) +mu_b_inter_ps2 * (photo*site3) +
#   mu_b_inter_s3c1 * (chill1*site3) + mu_b_inter_ws3 * (lp*site3) +mu_b_inter_ps3 * (photo*site3) +
#   mu_b_inter_s4c1 * (chill1*site3) + mu_b_inter_ws4 * (lp*site3) +mu_b_inter_ps4 * (photo*site3)
# #
# #Site 2
# bb_hpsite2 = a_sp + b_site2 * site2 + b_site3 * 0 + b_site4 * 0 + mu_b_warm * hp + mu_b_photo * photo + mu_b_chill1 * chill1 +
#   mu_b_inter_wp * (hp*photo) +
#   mu_b_inter_wc1 * (hp*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site2) + mu_b_inter_ws2 * (hp*site2) +mu_b_inter_ps2 * (photo*site2) +
#   mu_b_inter_s3c1 * (chill1*site2) + mu_b_inter_ws3 * (hp*site2) +mu_b_inter_ps3 * (photo*site2) +
#   mu_b_inter_s4c1 * (chill1*site2) + mu_b_inter_ws4 * (hp*site2) +mu_b_inter_ps4 * (photo*site2)
# 
# # plot first for the low forcing
# bb_lpsite2 = a_sp + b_site2 * site2 + b_site3 * 0 + b_site4 * 0 + mu_b_warm * lp + mu_b_photo * photo + mu_b_chill1 * chill1 +
#   mu_b_inter_wp * (lp*photo) +
#   mu_b_inter_wc1 * (lp*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site2) + mu_b_inter_ws2 * (lp*site2) +mu_b_inter_ps2 * (photo*site2) +
#   mu_b_inter_s3c1 * (chill1*site2) + mu_b_inter_ws3 * (lp*site2) +mu_b_inter_ps3 * (photo*site2) +
#   mu_b_inter_s4c1 * (chill1*site2) + mu_b_inter_ws4 * (lp*site2) +mu_b_inter_ps4 * (photo*site2)
# #
# plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Site", ylab = "Day of budburst", cex.lab = 2)
# # points(hpData$site4.z2, hpData$bb, bg = "#f9b641ff", pch =21, cex = 2.5)
# # points(lpData$site4.z2, lpData$bb, bg = "cyan4", pch = 21, cex = 2.5)
# 
# abline(lm(bb_hpsite2 ~ site2), col = "#a65c85ff", lwd = 3)
# abline(lm(bb_lpsite2 ~ site2), col = "#a65c85ff", lwd = 3, lty =2)
# 
# abline(lm(bb_hpsite4 ~ site4), col = "cyan4", lwd = 3)
# abline(lm(bb_lpsite4 ~ site4), col = "cyan4", lwd = 3, lty =2)
# 
# abline(lm(bb_hpsite3 ~ site3), col = "#f9b641ff", lwd = 3)
# abline(lm(bb_lpsite3 ~ site3), col = "#f9b641ff", lwd = 3, lty =2)
# 
# legend("topleft",legend = c(expression("St. Hippolyte"),
#                             expression("Harvard Forest"),
#                             expression("Manning Park"),
#                             expression("Long photoperiod"),
#                             expression("Short photoperiod")),
#        col = c("cyan4","#f9b641ff", "#a65c85ff", "black", "black" ),
#        inset = 0.02, lty = c(1,1,1,1,2 ), lwd =2,  cex = 1.1, bty = "n")
# 

dev.off()
##### EW Figures ##########################################################################

mu_params_ew_plot <- c(
                  "mu_b_warm",
                  "mu_b_photo",
                  "mu_b_chill",
                  "b_site",
                  "mu_b_inter_wp",
                  "mu_b_inter_wc1",
                  "mu_b_inter_pc1",
                  "mu_b_inter_ws",
                  "mu_b_inter_ps",
                  "mu_b_inter_sc1"
)
meanzew <- sumew[mu_params_ew_plot, col4fig]

rownames(meanzew) = c(
                      "Forcing",
                      "Photoperiod",
                      "Chilling",
                      "Site",
                      "Forcing x photoperiod",
                      "Forcing x chilling",
                      "Photoperiod x chilling",
                      "Forcing x site",
                      "Photoperiod x site",
                      "Chilling x site")


pdf(file.path( "figures/changes_pheno_ewtemp.pdf"), width = 7, height = 5)
par(mfrow = c(1,1), mar = c(2, 10, 2, 2))
# Upper panel: bud burst
plot(seq(-15, 
         15,
         length.out = nrow(meanzew)), 
     1:nrow(meanzew),
     type = "n",
     xlab = "",
     ylab = "",
     yaxt = "n")

#legend(x = -20, y = 2, bty="n", legend = "a. Budburst", text.font = 2)
#rasterImage(bbpng, -20, 1, -16, 4)

axis(2, at = nrow(meanzew):1, labels = rownames(meanzew), las = 1, cex.axis = 0.8)
points(meanzew[, 'mean'],
       nrow(meanzew):1,
       pch = 16,
       col = "Darkgreen",
       cex = 1.5)
arrows(meanzew[, "75%"], nrow(meanzew):1, meanzew[, "25%"], nrow(meanzew):1,
       len = 0, col = "black")
abline(v = 0, lty = 3)
# add advance/delay arrows
par(xpd=NA)
# arrows(1, 15.5, 6, 15.5, len = 0.1, col = "black")
# legend(5, 16.5, legend = "delay", bty = "n", text.font = 1, cex = 0.75)
# arrows(-1, 15.5, -6, 15.5, len = 0.1, col = "black")
# legend(-12, 16.5, legend = "advance", bty = "n", text.font = 1, cex = 0.75)
# legend(-2, 16.5, legend = "0", bty = "n", text.font = 1, cex = 0.75)
# par(xpd = FALSE)
dev.off()

pdf(file.path( "figures/changes.pheno.4sites.pdf"), width = 7, height = 8)
par(mfrow = c(1,1), mar = c(5, 10, 2, 1))
# Upper panel: bud burst
plot(seq(-15, 
         15,
         length.out = nrow(meanz4)), 
     1:nrow(meanz4),
     type = "n",
     xlab = "",
     ylab = "",
     yaxt = "n")

#legend(x = -20, y = 2, bty="n", legend = "a. Budburst", text.font = 2)
#rasterImage(bbpng, -20, 1, -16, 4)

axis(2, at = nrow(meanz4):1, labels = rownames(meanz4), las = 1, cex.axis = 0.8)
points(meanz4[, 'mean'],
       nrow(meanz4):1,
       pch = 16,
       col = "midnightblue")
arrows(meanz4[, "75%"], nrow(meanz4):1, meanz4[, "25%"], nrow(meanz4):1,
       len = 0, col = "black")
abline(v = 0, lty = 3)
# add advance/delay arrows
# par(xpd=NA)
# arrows(1, 15.5, 6, 15.5, len = 0.1, col = "black")
# legend(5, 16.5, legend = "delay", bty = "n", text.font = 1, cex = 0.75)
# arrows(-1, 15.5, -6, 15.5, len = 0.1, col = "black")
# legend(-12, 16.5, legend = "advance", bty = "n", text.font = 1, cex = 0.75)
# legend(-2, 16.5, legend = "0", bty = "n", text.font = 1, cex = 0.75)
# par(xpd = FALSE)
dev.off()

## Replicating Flynn Figure 2:

# b.force.both <- sumt[grep("^b_warm", rownames(sumt))]
# b.photo.both <- sumt[grep("^b_photo\\[", rownames(sumt))]
# b.chill.both <- sumt[grep("^b_chill1", rownames(sumt))]

b.force.both.ew <- sumew[grep("^b_warm\\[", rownames(sumew))]
b.photo.both.ew <- sumew[grep("^b_photo\\[", rownames(sumew))]
b.chill.both.ew <- sumew[grep("^b_chill1", rownames(sumew))]

shrubs.both = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG", "alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
trees.both = c("ACEPEN", "ACERUB", "ACESAC", "BETALL", "BETLEN", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL", "acegla","betpap", "poptre", "popbal")

pheno.term <- pheno[,c("bb", "chillport.z2", "force.z2", "photo.z2", "species", "lab2","transect")]


#pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 


species.both <- sort(unique(tolower(pheno.t$species)))
species.fact.both <-as.numeric( as.factor(unique(tolower(pheno.t$species))))
Type <- c("tree", "tree", "tree","tree", "shrub", "shrub", "shrub", "tree", "tree", "tree", "tree", "shrub","shrub","shrub",
        "tree", "shrub","shrub","shrub", "shrub", "shrub","shrub",  "shrub", "shrub", "tree","tree", "tree", "tree", "tree",
        "tree","tree","tree","shrub","shrub","shrub","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub", "shrub"
        ,"shrub", "shrub", "shrub", "shrub","shrub")
Transect <- c("western","eastern","eastern","eastern","both","western","western","eastern","eastern","eastern","both","eastern","western","eastern","eastern","eastern","eastern","eastern","eastern","western","eastern","western","eastern","both","eastern","western","eastern","eastern","eastern","eastern","eastern","western","eastern","western","western","western","western","western","western","western","western","western","western","eastern","eastern","western","eastern")

#both <- data.frame(species.both, species.fact.both, b.force.both, b.photo.both, b.chill.both, type.both)

both.ew <- data.frame(species.both, Transect, species.fact.both, b.force.both.ew, b.photo.both.ew, b.chill.both.ew, Type)
both.ew$Type[both.ew$species.both == "faggra"] <- "tree"

pdf(file.path( "figures/chill_vs_force_dldfew.pdf"), width = 5, height = 6)
#cf.both <- 
  # ggplot(both, aes(x= b.chill.both, y = b.force.both, col = type.both)) +
  # geom_point() +
  # #ylim (-25, 1) +
  # #xlim (-40, 0) +
  # labs( y = "High forcing", x = "High chilling") +
  # geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()
  color_scheme_set("viridis")
  ggplot(both.ew, aes(x= b.chill.both.ew, y = b.force.both.ew, col = Type, shape = Transect)) +
    geom_point() +
    labs( y = "High forcing", x = "High chilling") +
    geom_text(aes(label=species.both),hjust= 0.5, vjust= 1.5, show.legend = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
dev.off()
  
  
pdf(file.path( "figures/chill_vs_photo_dldfew.pdf"), width = 5, height = 6)
#cp.both <- 
  # ggplot(both, aes(x= b.chill.both, y = b.photo.both, col = Type) +
  # geom_point() +
  # #ylim (-5.5, 1) +
  # #xlim (-30, 0) +
  # labs (x = "High chilling", y = "Long photoperiod", color = type.both, shape = Transect) +
  # geom_text(aes(label=species.both), hjust=0.5, vjust= 1) +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"))

#dev.off()
  
  ggplot(both.ew, aes(x= b.chill.both.ew, y = b.photo.both.ew, col = Type, shape = Transect)) +
    geom_point() +
    ylim (-5.5, -2) +
    xlim (-30,-2) +
    labs (x = "High chilling", y = "Long photoperiod") +
    geom_text(aes(label=species.both),hjust=0.5, vjust= 2, show.legend = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
#legend.position = "none"
  dev.off()

#pdf(file.path( "figures/force_vs_photo_dldf.pdf"), width = 7, height = 8)
fp.both <- ggplot(both.ew, aes(x= b.force.both.ew, y = b.photo.both.ew, col = Type)) +
  geom_point() +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1, show.legend = F)+
  ylim (-5.5, -2) +
  xlim (-22, 0.5) +
  labs(x = "High forcing", y = "Long photoperiod") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
#dev.off()

## Plotting the day to bb with the cues on the y-axis 
pheno$species <- tolower(pheno$species)
term.bb.both <- ddply(pheno, c("species"), summarize, mean = mean(bb, na.rm = TRUE))
names(term.bb.both) <- c("species.both", "mean")

term.both <- merge(term.bb.both, both.ew, by = "species.both", all =TRUE)
term.both <- term.both[,c("species.both","mean","b.force.both.ew","b.chill.both.ew","b.photo.both.ew")]
term.both <- term.both[complete.cases(term.both), ] 

pdf(file.path( "figures/force_dobb_dldfew.pdf"), width = 5, height = 6)
tf.both <-  ggplot(term.both, aes(y = b.force.both.ew, x= mean,col = Type, shape = Transect)) +
  geom_point() +
  ylim (-13, -4) +
  xlim (10, 55) +
  labs(x = "Mean day of budburst", y = "High forcing") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
tf.both
dev.off()

pdf(file.path( "figures/chill_dobb_dldfew.pdf"), width = 5, height = 6)
tc.both <- ggplot(term.both, aes(y = b.chill.both.ew, x= mean,col = Type, shape = Transect)) +
  geom_point() +
  ylim (-27, -4) +
  xlim (10, 55) +
  labs(x = "Mean day of budburst", y = "High chilling") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
tc.both
dev.off()

pdf(file.path( "figures/photo_dobb_dldfew.pdf"), width = 5, height = 6)
tp.both <- ggplot(term.both, aes(y = b.photo.both.ew, x= mean,col = Type, shape = Transect)) +
  geom_point() +
  ylim (-5.5, -2) +
  xlim (10, 55) +
  labs(x = "Mean day of budburst", y = "Long photoperiod")+
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff","#eb8055ff","#f9b641ff","#a65c85ff","#7e4e90ff", "cyan4"))
tp.both
dev.off()
##### General boxplots across treatments:
# As per Lizzie's July 7 post: I should look at how linear these relationships are
# Plot raw data (bb~ chill) with the chill effect from the model plotted on top

plot(pheno.t$tbb ~ pheno.t$chillport.z2, pch = 19, col = "darkgreen")
abline(a = sumt[grep("mu_a", rownames(sumt)), "mean"], b = sumt[grep("mu_chill", rownames(sumt)), "mean"])

plot(sumt[grep("sigma_chill", rownames(sumt)), "mean"])
plot(sumt[grep("log", rownames(sumt)), "mean"])


## Do we see differences across shrubs and trees:
tree <- subset(both.ew, Type == "tree")
treeMean <- colMeans(tree[,c(4,5,6)])
#b.force.both.ew b.photo.both.ew b.chill.both.ew 
#-9.157388       -3.566382      -14.261067 

shrub <- subset(both.ew, Type == "shrub")
shrubMean <- colMeans(shrub[,c(4,5,6)])

#b.force.both.ew b.photo.both.ew b.chill.both.ew 
#-8.628631       -3.361155      -12.952969 
# par(mar =c (6,5,1,1))
# pdf(file="figures/dltrt_boxplot.pdf")
# west <- boxplot(dl$tbb ~ dl$treatment, las =2, xlab ="", ylab = "Day of terminal bb of western spp." )
# dev.off()

# pdf(file="figures/dftrt_boxplot.pdf")
# east <- boxplot(df$tbb ~ df$treatment, las =2, xlab ="", ylab = "Day of terminal bb of eastern spp.")
# dev.off()
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
# plotblank = function(){plot(1:10, type="n",bty="n",xaxt="n",yaxt="n",y="",x="")}
# 
# plotblank() 
# text(5,5, "Budburst \n Change (days) due to 5° warming", font = 2, srt = 90) # \n\n add two line breaks
# 
# plot( "b.photo", "b_warm",
#          #  y = "Advance due to 5° warming", 
#          # x = "Advance due to 4 hr longer photoperiod", 
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
# plotlet("b.chill1", "b_warm", 
#         # y = "Advance due to 5° warming", 
#         #  x = "Advance due to 30d 4° chilling", 
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
# plotlet("b.photo", "b_warm", 
#         #    y = "Advance due to 5° warming", 
#         #     x = "Advance due to 4 hr longer photoperiod", 
#         ylim = c(-27, 0.5),
#         xlim = c(-16, 0.5),
#         group = treeshrub,
#         data = sumerl)
# legend("topleft", bty = "n", inset = 0.035, legend = "C.", text.font=2)
# plotlet("b.chill1", "b_warm", 
#         #   y = "Advance due to 5° warming", 
#         #   x = "Advance due to 30d 4° chilling", 
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

a_sp = sumew[grep("a_sp", rownames(sumew)), 1]
mu_b_warm = sumew[grep("mu_b_warm", rownames(sumew)), 1]
mu_b_photo = sumew[grep("mu_b_photo", rownames(sumew)), 1]
mu_b_chill1 = sumew[grep("mu_b_chill", rownames(sumew)), 1]
mu_b_inter_ws = sumew[grep("mu_b_inter_ws", rownames(sumew)), 1]
mu_b_inter_sc1 = sumew[grep("mu_b_inter_sc1", rownames(sumew)), 1]
mu_b_inter_ps = sumew[grep("mu_b_inter_ps", rownames(sumew)), 1]
mu_b_inter_pc1 = sumew[grep("mu_b_inter_pc1", rownames(sumew)), 1]
mu_b_inter_wp = sumew[grep("mu_b_inter_wp", rownames(sumew)), 1]
mu_b_inter_wc1 = sumew[grep("mu_b_inter_wc1", rownames(sumew)), 1]
b_site = sumew[grep("b_site", rownames(sumew)), 1]

## transect and chilling:
west <- subset(pheno, transect.n == "0")
east <- subset(pheno, transect.n == "1")

force <- -0.5080665 # zero forcing
photo <- -0.5044652 # zero photo
chill1 <- seq( -1, 1, by = 0.01)
west.sites <- unique(west$transect.n)
east.sites <- unique(east$transect.n)

# plot first for the west coast
bb_westc <- matrix(NA, ncol = 1, nrow = length(chill1))
for(k in 1:length(chill1)){
  bb_westc[k] <- a_sp + b_site * west.sites + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1[k] +
    mu_b_inter_wp * (force*photo) + mu_b_inter_ws * (force*west.sites) +mu_b_inter_ps * (photo*west.sites) +
    mu_b_inter_wc1 * (force*chill1[k]) + mu_b_inter_pc1 * (photo*chill1[k]) +
    mu_b_inter_sc1 * (chill1[k]*west.sites)
}
# plot first for the east coast
bb_eastc <- matrix(NA, ncol = 1, nrow = length(chill1))
for(k in 1:length(chill1)){
  bb_eastc [k] <- a_sp + b_site *  east.sites + mu_b_warm * force + mu_b_photo * photo + mu_b_chill1 * chill1[k] +
    mu_b_inter_wp * (force*photo) +mu_b_inter_ws * (force* east.sites) +mu_b_inter_ps * (photo* east.sites) +
    mu_b_inter_wc1 * (force*chill1[k]) +mu_b_inter_pc1 * (photo*chill1[k]) +
    mu_b_inter_sc1 * (chill1[k]* east.sites)
}
pdf("figures/ew_interactions_chillSite.pdf", width = 5, height = 5)
par(mfrow = c(1,1))
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Chill portions", ylab = "Day of budburst")
# points(west$chillport.z2,west$bb, col = "maroon")
# points(east$chillport.z2,east$bb, col = "darkslategray4")

points(bb_westc ~ chill1, type = "l", col = "darkred", lwd = 3)
points(bb_eastc ~ chill1, type = "l", col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("eastern transect"),
                            expression("western transect")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(19, 19  ),  cex = 1, bty = "n")
dev.off()
######################################################################
######################################################################
## transect and forcing:
west <- subset(pheno.t, transect.n == "0")
east <- subset(pheno.t, transect.n == "1")

#force <- unique(pheno.t$force.z2)
force <- seq(-1, 1, by = 0.01)
photo <- -0.5044652 # zero photo
chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)
west.sites <- unique(west$transect.z2)
east.sites <- unique(east$transect.z2)

# plot first for the west coast

bb_westf <- matrix(NA, ncol = 1, nrow = length(force))
for(k in 1:length(force)){
  bb_westf[k] <- a_sp + b_site * west.sites + mu_b_warm * force[k] + mu_b_photo * photo + mu_b_chill1 * chill1 +
    mu_b_inter_wp * (force[k]*photo) + mu_b_inter_ws * (force[k]*west.sites) +mu_b_inter_ps * (photo*west.sites) +
    mu_b_inter_wc1 * (force[k]*chill1) + mu_b_inter_pc1 * (photo*chill1) +
    mu_b_inter_sc1 * (chill1*west.sites)
}

# plot first for the east coast
bb_eastf <- matrix(NA, ncol = 1, nrow = length(chill1))
for(k in 1:length(force)){
  bb_eastf[k] <- a_sp + b_site *  east.sites + mu_b_warm * force[k] + mu_b_photo * photo + mu_b_chill1 * chill1 +
    mu_b_inter_wp * (force[k]*photo) +mu_b_inter_ws * (force[k]* east.sites) +mu_b_inter_ps * (photo* east.sites) +
    mu_b_inter_wc1 * (force[k]*chill1) +mu_b_inter_pc1 * (photo*chill1) +
    mu_b_inter_sc1 * (chill1* east.sites)
}

pdf(file="figures/ew_interaction_siteForce.pdf", width = 5, height = 5)
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Standardized forcing", ylab = "Day of budburst")
points(west$force.z2,west$bb, col = "maroon")
points(east$force.z2,east$bb, col = "darkslategray4")

points(bb_westf ~force, type = "l", col = "darkred", lwd = 3)
points(bb_eastf ~force, type = "l", col = "darkslategray", lwd = 3)

# abline(lm(bb_westf ~ force), pch =19, col = "darkred", lwd = 3)
# abline(lm(bb_eastf ~ force), pch =19, col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("eastern transect"),
                            expression("western transect")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(19,19 ),  cex = 1, bty = "n")
dev.off()
# ######################################################################
# ######################################################################
## transect and photoperiod:
west <- subset(pheno.t, transect.n == "0")
east <- subset(pheno.t, transect.n == "1")

#force <- unique(pheno.t$force.z2)
force <- -0.5080665 # zero forcing
photo <- seq(-1,1, by = 0.01)
chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)
west.sites <- unique(west$transect.z2)
east.sites <- unique(east$transect.z2)

# plot first for the west coast
bb_westp <- matrix(NA, ncol = 1, nrow = length(photo))
for(k in 1:length(photo)){
  bb_westp [k] <- a_sp + b_site * west.sites + mu_b_warm * force + mu_b_photo * photo[k] + mu_b_chill1 * chill1 +
    mu_b_inter_wp * (force*photo[k]) + mu_b_inter_ws * (force*west.sites) +mu_b_inter_ps * (photo[k]*west.sites) +
    mu_b_inter_wc1 * (force*chill1) + mu_b_inter_pc1 * (photo*chill1) +
    mu_b_inter_sc1 * (chill1*west.sites)
}
# plot first for the east coast
bb_eastp <- matrix(NA, ncol = 1, nrow = length(photo))
for(k in 1:length(photo)){
  bb_eastp[k] <- a_sp + b_site *  east.sites + mu_b_warm * force + mu_b_photo * photo[k] + mu_b_chill1 * chill1 +
    mu_b_inter_wp * (force*photo[k]) +mu_b_inter_ws * (force* east.sites) +mu_b_inter_ps * (photo[k] * east.sites) +
    mu_b_inter_wc1 * (force*chill1) +mu_b_inter_pc1 * (photo[k] * chill1) +
    mu_b_inter_sc1 * (chill1* east.sites)
}
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Forcing", ylab = "Day of budburst")
points(west$photo.z2,west$bb, col = "maroon")
points(east$photo.z2,east$bb, col = "darkslategray4")

points(bb_westp ~ photo, type = "l", col = "darkred", lwd = 3)
points(bb_eastp ~ photo, type = "l", col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("eastern transect"),
                            expression("western transect")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(19,19 ),  cex = 0.75, bty = "n")
#####################################################################
## transect and chilling:
lfData <- subset(pheno.t, force.n == "0")
hfData <- subset(pheno.t, force.n == "1")

lf <- unique(lfData$force.z2) # zero forcing
hf <- unique(hfData$force.z2) # zero forcing

photo <- -0.5044652 # zero photo
chill1 <- seq( -1, 1, by = 0.01)
sites <- mean(unique(pheno.t$transect.z2))


# plot first for the west coast
bb_hf_c <- matrix(NA, ncol = 1, nrow = length(chill1))
for(k in 1:length(chill1)){
  bb_hf_c[k] <- a_sp + b_site * sites + mu_b_warm * hf + mu_b_photo * photo + mu_b_chill1 * chill1[k] +
    mu_b_inter_wp * (hf*photo) + mu_b_inter_ws * (hf*sites) +mu_b_inter_ps * (photo*sites) +
    mu_b_inter_wc1 * (hf*chill1[k]) + mu_b_inter_pc1 * (photo*chill1[k]) +
    mu_b_inter_sc1 * (chill1[k]*sites)
}

# plot first for the east coast
bb_lf_c <- matrix(NA, ncol = 1, nrow = length(chill1))
for(k in 1:length(chill1)){
  bb_lf_c [k] <- a_sp + b_site * sites + mu_b_warm * lf + mu_b_photo * photo + mu_b_chill1 * chill1[k] +
    mu_b_inter_wp * (lf*photo) + mu_b_inter_ws * (lf*sites) +mu_b_inter_ps * (photo*sites) +
    mu_b_inter_wc1 * (lf*chill1[k]) + mu_b_inter_pc1 * (photo*chill1[k]) +
    mu_b_inter_sc1 * (chill1[k]*sites)
}
#pdf("figures/ew_interactions.pdf", width = 10, height = 3.5)
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Chill portions", ylab = "Day of budburst")
points(lfData$chillport.z2,  lfData$bb, col = "maroon")
points(hfData$chillport.z2, hfData$bb, col = "darkslategray4")

points(bb_lf_c ~ chill1, type = "l", col = "darkred", lwd = 3)
points(bb_hf_c ~ chill1, type = "l", col = "darkslategray", lwd = 3)

legend("topleft",legend = c(expression("eastern transect"),
                            expression("western transect")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")
dev.off()

# Forcing X chilling

hfData <- subset(pheno, force == "HF" )
lfData <- subset(pheno, force == "LF")

# Make the other parameters constant
hf <- unique(hfData$force.z2)
lf <- unique(lfData$force.z2)
photo <- -0.5041133
siteSM <- 0
chill1 <- c( -0.7642814, -0.4072595, -0.4023109, -0.3493703,  0.2750890,  0.2977055,  0.4308763,  0.5308110,  0.8457874,  0.9457221)


# plot first for the high forcing
bb_hfc = a_sp + b_site * siteSM + mu_b_warm * hf + mu_b_photo * photo + mu_b_chill1 * chill1 +
  mu_b_inter_wp * (hf*photo) +
  mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_sc1 * (chill1*siteSM) + mu_b_inter_ws * (hf*siteSM) +mu_b_inter_ps * (photo*siteSM) 

# plot first for the low forcing
bb_lfc = a_sp + b_site * siteSM + mu_b_warm * lf + mu_b_photo * photo + mu_b_chill1 * chill1 +
  mu_b_inter_wp * (lf*photo) +
  mu_b_inter_wc1 * (lf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
  mu_b_inter_sc1 * (chill1*siteSM) + mu_b_inter_ws * (lf*siteSM) +mu_b_inter_ps * (photo*siteSM) 

pdf("figures/chill_forcing_4sites_interactions.pdf", width =5, height = 5)
par(mfrow =c (1,1))
plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "Standardized chill portions", ylab = "Day of budburst")
points(hfData$chillport.z2, hfData$bb, col = "maroon")
points(lfData$chillport.z2, lfData$bb, col = "darkslategray4")
ablineclip(lm(bb_hfc ~ chill1), col = "darkred", lwd = 3, x1 =-1, x2 = 1)
ablineclip(lm(bb_lfc ~ chill1), col = "darkslategray", lwd = 3, x1 =-1, x2 = 1)

legend("topright",legend = c(expression("low forcing"),
                             expression("high forcing")),
       col = c("darkslategray","maroon"),
       inset = 0.02, pch = c(19, 19 ),  cex = 1, bty = "n")
dev.off()

############ Old figures ##########
siteF <- ggplot(siteForce, aes(x = factor(site, level = siteOrder), y = value, fill = force)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_discrete(name="Forcing level",
                      labels=c("high forcing"="High forcing", "low forcing" = "Low forcing")) +
  xlab("Population") + ylab("Estimated day of budburst") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15)) +
  scale_fill_manual(values = c("cyan4", "#CC6677"), labels = c("High forcing", "Low forcing"), name = "") + theme(legend.key=element_blank(), legend.position=c(.9,.85)) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 60, label = "b)", cex =5) 

siteC <- ggplot(siteChill, aes(x = factor(site, level = siteOrder), y = value, fill = col)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_discrete(name="Forcing level",
                      labels=c("high forcing"="High forcing", "low forcing" = "Low forcing")) +
  xlab("Population") + ylab("Estimated day of budburst") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 15), axis.title = element_text(size = 15)) +
  scale_fill_manual(values = c(
    "cyan4", 
    "salmon2",
    "orchid4", 
    "#CC6677",
    "#f9b641ff"), labels = c("High chill - East",
                             "High chill - West",
                             "Low chill - West",
                             "Low chill - East",
                             "No chill - East"), name = "") + theme(legend.key=element_blank(), legend.position=c(.9,.85)) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 40, label = "c)", cex =5) 

siteP <- ggplot(sitePhoto, aes(x = factor(site, level = siteOrder), y = value, fill = photo)) +
  geom_bar(stat="identity", position="dodge") +
  xlab("Population") + ylab("Estimated day of budburst") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 15)) +
  scale_fill_manual(values = c("cyan4", "#CC6677"), labels = c("Long photoperiod", "Short photoperiod"), name = "") + theme(legend.key=element_blank(), legend.position=c(.8,.85)) +
  theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 60, label = "d)", cex =5) 



# warm and site3
# hf <- unique(hfData$force.z2)
# lf <- unique(lfData$force.z2)
# photo <- -0.5044652
# site2 <- unique(pheno$site2)
# site3 <- unique(pheno$site3)
# 
# chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)
# 
# # plot first for the high forcing
# bb_hfsite3 = a_sp + b_site2 * site3 + b_site3 * site3 + b_site4 * site3 + mu_b_warm * hf + mu_b_photo * photo + mu_b_chill1 * chill1 +
#   mu_b_inter_wp * (hf*photo) +
#   mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site2) + mu_b_inter_ws2 * (hf*site2) +mu_b_inter_ps2 * (photo*site2) +
#   mu_b_inter_s3c1 * (chill1*site2) + mu_b_inter_ws3 * (hf*site2) +mu_b_inter_ps3 * (photo*site2) +
#   mu_b_inter_s4c1 * (chill1*site2) + mu_b_inter_ws4 * (hf*site2) +mu_b_inter_ps4 * (photo*site2)
# 
# # plot first for the low forcing
# bb_lfsite3 = a_sp + b_site2 * site3 + b_site3 * site3 + b_site4 * site3  + mu_b_warm * lf + mu_b_photo * photo + mu_b_chill1 * chill1 +
#   mu_b_inter_wp * (lf*photo) +
#   mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site2) + mu_b_inter_ws2 * (lf*site2) +mu_b_inter_ps2 * (photo*site2) +
#   mu_b_inter_s3c1 * (chill1*site2) + mu_b_inter_ws3 * (lf*site2) +mu_b_inter_ps3 * (photo*site2) +
#   mu_b_inter_s4c1 * (chill1*site2) + mu_b_inter_ws4 * (lf*site2) +mu_b_inter_ps4 * (photo*site2)
# #
# pdf("figures/Site3_forcing_4sites_interactions.pdf", width =5, height = 5)
# plot(0, type = "n",  xlim = c(-0.5,1.5), ylim = c(-5,90), xlab = "Harvard Forest", ylab = "Day of budburst")
# points(hfData$site3, hfData$bb, bg = "#f9b641ff", pch =21)
# points(lfData$site3, lfData$bb, bg = "cyan4", pch =21)
# abline(lm(bb_hfsite3 ~ site2), bg = "#f9b641ff", pch =21)
# abline(lm(bb_lfsite3 ~ site2), col = "cyan4", lwd = 3)
# 
# #legend("topleft",legend = c(expression("low forcing"),
# #                            expression("high forcing")),
# #       col = c("darkslategray","maroon"),
# #       inset = 0.02, pch = c(19, 19 ),  cex = 1, bty = "n")
# dev.off()

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# # warm and site3 
# hf <- unique(hfData$force.z2)
# lf <- unique(lfData$force.z2)
# photo <- -0.5044652 
# site3 <- unique(pheno$site3.z2)
# 
# chill1 <- mean( -0.3482404,  0.9462697,  0.8463799, -0.7629649,  0.5315452,  0.4316554,0.2985445, -0.4011572,  0.2759381, -0.4061035)
# 
# # plot first for the high forcing
# bb_hfsite3 = a_sp + b_site2 * site3 + b_site3 * site3 + b_site4 * site3 + mu_b_warm * hf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
#   mu_b_inter_wp * (hf*photo) +
#   mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site3) + mu_b_inter_ws2 * (hf*site3) +mu_b_inter_ps2 * (photo*site3) +
#   mu_b_inter_s3c1 * (chill1*site3) + mu_b_inter_ws3 * (hf*site3) +mu_b_inter_ps3 * (photo*site3) +
#   mu_b_inter_s4c1 * (chill1*site3) + mu_b_inter_ws4 * (hf*site3) +mu_b_inter_ps4 * (photo*site3) 
# 
# # plot first for the low forcing
# bb_lfsite3 = a_sp + + b_site2 * site3 + b_site3 * site3 + b_site4 * site3 + mu_b_warm * lf + mu_b_photo * photo + mu_b_chill1 * chill1 + 
#   mu_b_inter_wp * (lf*photo) +
#   mu_b_inter_wc1 * (hf*chill1) + mu_b_inter_pc1 * (photo*chill1) +
#   mu_b_inter_s2c1 * (chill1*site3) + mu_b_inter_ws2 * (lf*site3) +mu_b_inter_ps2 * (photo*site3) +
#   mu_b_inter_s3c1 * (chill1*site3) + mu_b_inter_ws3 * (lf*site3) +mu_b_inter_ps3 * (photo*site3) +
#   mu_b_inter_s4c1 * (chill1*site3) + mu_b_inter_ws4 * (lf*site3) +mu_b_inter_ps4 * (photo*site3) 
# # 
# plot(0, type = "n",  xlim = c(-1,1), ylim = c(-5,90), xlab = "site3", ylab = "Day of budburst")
# points(hfData$site3.z2, hfData$bb, bg = "#f9b641ff", pch =21)
# points(lfData$site3.z2, lfData$bb, bg = "cyan4", pch =21)
# abline(lm(bb_hfsite3 ~ site3), col = "#f9b641ff", lwd = 3)
# abline(lm(bb_lfsite3 ~ site3), col = "cyan4", lwd = 3)
# 
# legend("topleft",legend = c(expression("low forcing"),
#                             expression("high forcing")),
#        col = c("darkslategray","maroon"),
#        inset = 0.02, pch = c(21,21 ),  cex = 0.75, bty = "n")
# 

##### New figure with fewer colours:

siteForce$transect <- siteForce$site
siteForce$transect[siteForce$transect == "Smithers"] <- "Western"
siteForce$transect[siteForce$transect == "Manning park"] <- "Western"
siteForce$transect[siteForce$transect == "Harvard forest"] <- "Eastern"
siteForce$transect[siteForce$transect == "St.Hippolyte"] <- "Eastern"

siteForce$tranSite <- paste(siteForce$transect, siteForce$site, sep = " ")

siteFPoint <- ggplot()  +
  geom_pointrange(siteForce, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=lower1, ymax=upper1, col = temp, shape = transect), size =0.55, position=position_dodge(width=0.5)) +
  geom_linerange( siteForce, mapping = aes(x = factor(site, level = siteOrder), ymin = lower2, ymax=upper2, col = temp), size =1.25, position=position_dodge(width=0.5)) +
  xlab("Population") + ylab("Estimated day of budburst") +
  #ylim(0,50) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20),axis.text.x = element_text( size= 15,angle = 55, hjust=1)) +
  # scale_color_manual(values = c("cyan4", "#CC6677"), labels = c("High forcing", "Low forcing"), name = "") +
  scale_color_manual(values = c(
    "bb_hfsite1"="forestgreen",
    "bb_lfsite1"="palegreen3",
    "bb_hfsite2"= "forestgreen",
    "bb_lfsite2"="palegreen3",
    "bb_hfsite4"="darkorchid4",
    "bb_lfsite4"="darkorchid1",
    "bb_hfsite3" ="darkorchid4",
    "bb_lfsite3"="darkorchid1"
  ), labels = c("Smithers - High forcing", "Smithers - Low forcing","Manning Park - High forcing", "Manning Park - Low forcing", "St. Hippolyte - High forcing", "St. Hippolyte - Low forcing","Harvard Forest - High forcing", "Harvard Forest - Low forcing"), name = "") +theme(legend.position = "none") +  annotate("text", x = 1.5, y = 75, label = "a) Forcing", cex = 10) 
siteFPoint
# theme(legend.key=element_blank(), legend.position=c(.83,.85),legend.text = element_text(size = 15)) +
# theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 50, label = "b)", cex = 10) 
siteFPoint

siteChill$transect <- siteChill$site
siteChill$transect[siteChill$transect == "Smithers"] <- "Western"
siteChill$transect[siteChill$transect == "Manning park"] <- "Western"
siteChill$transect[siteChill$transect == "Harvard forest"] <- "Eastern"
siteChill$transect[siteChill$transect == "St.Hippolyte"] <- "Eastern"

siteChill$tranSite <- paste(siteChill$transect, siteChill$site, sep = " ")


siteCPoint <- ggplot() +
  geom_pointrange(siteChill, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=lower1, ymax=upper1, col = temp, shape = transect),
                  position=position_dodge(width=0.5), size = 0.55) +
  geom_linerange( siteChill, mapping = aes(x = factor(site, level = siteOrder), ymin = lower2, ymax=upper2, col = temp), size =1.25, position=position_dodge(width=0.5))+
  #ylim (0,50) +
  xlab("Population") + ylab("Estimated day of budburst") +
  #ylim(0,50) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20),axis.text.x = element_text( size= 15,angle = 55, hjust=1)) +
  scale_color_manual(values = c(
    "bb_hcsite1"="forestgreen",
    "bb_lcsite1"="palegreen3",
    "bb_hcsite2"= "forestgreen",
    "bb_lcsite2"="palegreen3",
    "bb_hc1site4"="darkorchid4",
    "bb_lc0site4"="darkorchid1",
    "bb_hc1site3" ="darkorchid4",
    "bb_lc0site3"="darkorchid1"
  ), labels = c("Smithers - High  chill", "Smithers - Low  chill","Manning Park - High  chill", "Manning Park - Low  chill", "St. Hippolyte - High  chill", "St. Hippolyte - Low  chill","Harvard Forest - High  chill", "Harvard Forest - Low chill")) + theme(legend.position = "none") +  annotate("text", x = 1.5, y = 80, label = "b) Chilling", cex = 10) 
siteCPoint

sitePhoto$transect <- sitePhoto$site
sitePhoto$transect[sitePhoto$transect == "Smithers"] <- "Western"
sitePhoto$transect[sitePhoto$transect == "Manning park"] <- "Western"
sitePhoto$transect[sitePhoto$transect == "Harvard forest"] <- "Eastern"
sitePhoto$transect[sitePhoto$transect == "St.Hippolyte"] <- "Eastern"

sitePhoto$tranSite <- paste(sitePhoto$transect, sitePhoto$site, sep = " ")

sitePPoint <- ggplot() +
  geom_pointrange(sitePhoto, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=lower1, ymax=upper1, col = tranSite, shape = transect),
                  position=position_dodge(width=0.5), size =0.55) +
  geom_linerange( sitePhoto, mapping = aes(x = factor(site, level = siteOrder), ymin = lower2, ymax=upper2, col =tranSite), size =1.25, position=position_dodge(width=0.5))+
  xlab("Population") + ylab("Estimated day of budburst") +
  # ylim(0,50) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20), axis.text.x = element_text( size= 15,angle = 55, hjust=1))+
  scale_color_manual(values = c(
    "Western Smithers"="forestgreen",
    "Western Manning Park"="palegreen3",
    "Eastern Harvard forest"= "darkorchid1",
    "Eastern St.Hippolyte"="darkorchid4"), labels = c("Western - High cue", "Western - Low cue","Easterb - High cue", "Eastern - Low cue")) + theme(legend.title = element_blank(), legend.key=element_blank(), legend.text = element_text(size=20))  +  annotate("text", x = 1.5, y = 80, label = "c)", cex = 10) 
sitePPoint

sitePPointNL <- ggplot() +
  geom_pointrange(sitePhoto, mapping = aes(x = factor(site, level = siteOrder), y = mean, ymin=lower1, ymax=upper1, col = temp, shape = transect),
                  position=position_dodge(width=0.5), size =.55) +
  geom_linerange( sitePhoto, mapping = aes(x = factor(site, level = siteOrder), ymin = lower2, ymax=upper2, col = temp), size =1.25, position=position_dodge(width=0.5))+
  xlab("Population") + ylab("Estimated day of budburst") +
  # ylim(0,50) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20), axis.text.x = element_text( size= 15,angle = 55, hjust=1))+
  scale_color_manual(values = c(
    "bb_hpsite1"="forestgreen",
    "bb_lpsite1"="palegreen3",
    "bb_hpsite2"= "forestgreen",
    "bb_lpsite2"="palegreen3",
    "bb_hpsite3" ="darkorchid4",
    "bb_lpsite3"="darkorchid1",
    "bb_hpsite4"="darkorchid4",
    "bb_lpsite4"="darkorchid1"), labels = c("Smithers - High cue", "Smithers - Low cue","Manning Park - High cue", "Manning Park - Low cue", "St. Hippolyte - High cue", "St. Hippolyte - Low cue","Harvard Forest - High cue", "Harvard Forest - Low cue")) + theme(legend.position = "none")  +  annotate("text", x = 1.75, y = 80, label = "c) Photoperiod", cex = 10) 

# theme(legend.key=element_blank(), legend.position=c(.8,.85),legend.text = element_text(size = 15)) +
# theme(legend.title = element_blank()) +  annotate("text", x = 0.6, y = 50, label = "d)", cex = 10) 


# pdf("figures/intrxnPlotsApril19.pdf", height =12, width = 15)
# plot_grid(intrxnCF, siteFPoint, siteCPoint, sitePPoint, nrow = 2, ncol = 2, align = "v")
# dev.off()

l <- get_legend(sitePPoint) 

pdf("figures/intrxnForceChill.pdf", height = 5, width = 5)
intrxnCF
dev.off()

pdf("figures/intrxnPlots8Chill3_color2.pdf", height =5, width = 20)
plot_grid( siteFPoint, siteCPoint, sitePPointNL, l , ncol = 4, nrow =1,align = "v")
dev.off()


## Need a nicer phylogy - colour tips based on intercept values:
tree <- read.tree("input/SBphylo_phenobc.tre")
head(tree$tip.label)
tree <- read.tree("input/SBphylo_phenobc.tre")
head(tree$tip.label)
length(tree$tip.label) #47
tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
tree$tip.label[tree$tip.label=="Alnus_alnobetula"] <- "Alnus_viridis"
tree$tip.label[tree$tip.label== "Fagus_grandifolia_var._caroliniana"] <- "Fagus_grandifolia"
tree$tip.label[tree$tip.label== "Spiraea_alba_var._latifolia"] <- "Spiraea_alba"
tree$tip.label[tree$tip.label== "Rhamnus_arguta"] <- "Rhamnus_frangula"


# first getting the intercept values paired with the species names:
spInfo <- read.csv("input/species_list.csv")


spFact <- spInfo$species.name
rownames(spInfo) <- spFact

intercepts <- data.frame(sum[grep("a_sp\\[", rownames(sum)), c("mean","25%")])
phylo.dat <- spInfo

dat.int <- cbind(phylo.dat, intercepts)

# tree needs to be rooted:
namesphy <- tree$tip.label
tree$root.edge <- 0
root(tree, outgroup = "Acer_glabrum")

is.rooted(tree)
tree$node.label<-NULL

dataPhy = comparative.data(tree, dat.int, names.col = "species.name", na.omit = T,
  vcv = T, warn.dropped = T)

phyloplot = dataPhy$phy
x = dataPhy$data$mean
names(x)=dataPhy$phy$tip.label

slope <- contMap(tree, x, plot = T)

slopeCol <- setMap(slope, colors=c("blue","purple","red"))
h<-max(nodeHeights(slopeCol$tree))

pdf("figures/phyloIntColor.pdf", height = 9, width = 7)
plot(slopeCol,legend = F, lwd=3, ylim=c(1-0.09*(Ntip(slopeCol$tree)),Ntip(slopeCol$tree)))

add.color.bar(60, slopeCol$cols, title = "Intercept (days)", subtitle="", digits = 1,lims = c(10,50),  prompt = F,x=0.2*h, y = -2)
dev.off()
# LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
# lines(x=rep(-0.2*h+LWD*15/2,2),y=c(2,Ntip(slopeCol$tree)-1))
# nticks<-10
# Y<-cbind(round(seq(2,Ntip(obj$tree)-1,length.out=nticks),1),
#   round(seq(2,Ntip(obj$tree)-1,length.out=nticks),1))
# X<-cbind(rep(-0.2*h+LWD*15/2,nticks),
#   rep(-0.2*h+LWD*15/2+0.02*h,nticks))
# for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
# ticks<-seq(round(slopeCol$lims[1],1),round(slopeCol$lims[2],1),length.out=nticks)
# text(x=X[,2],y=Y[,2],round(ticks,3),pos=4,cex=0.8)
# 
# 
# 
# # plot(tree)
# # x <- as.numeric(sum[grep("a_sp\\[", rownames(sum)), c("mean")])
# # 
# # rownames(intercepts) <- spFact
# # 
# # obj<-contMap(tree,x,plot=FALSE)
# # obj
# # plot(obj)
# 
# obj<-setMap(obj,colors=c("blue","purple","red"))
# plot(obj,leg.txt="Intercept)", lims = c(10,50),lwd=3)
# cols=c("blue","purple","red") 
# leg <- 50
# add.color.bar(leg, cols, lims = c(10,50))
# add.simmap.legend(leg = Null, cols, prompt = T)
# data(mammal.tree)
# data(mammal.data)
# ## extract character of interest
# ln.bodyMass<-(setNames(intercepts$mean,
#   rownames(intercepts)))
# ## create "contMap" object
# mammal.contMap<-contMap(tree,
#   ln.bodyMass,plot=T)
# ## change color scheme
# mammal.contMap<-setMap(mammal.contMap,
#   c("white","#FFFFB2","#FECC5C","#FD8D3C",
#     "#E31A1C"))
# plot(mammal.contMap,fsize=c(0.7,0.8),
#   leg.txt="log(body mass)")
# par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default
# 
# 
# 
# 
# # now get the tree:
# tree <- read.tree("input/SBphylo_phenobc.tre")
# head(tree$tip.label)
# length(tree$tip.label) #47
# tree$tip.label[tree$tip.label=="Cornus_asperifolia"] <- "Cornus_stolonifera"
# tree$tip.label[tree$tip.label=="Alnus_alnobetula"] <- "Alnus_viridis"
# tree$tip.label[tree$tip.label== "Fagus_grandifolia_var._caroliniana"] <- "Fagus_grandifolia"
# tree$tip.label[tree$tip.label== "Spiraea_alba_var._latifolia"] <- "Spiraea_alba"
# tree$tip.label[tree$tip.label== "Rhamnus_arguta"] <- "Rhamnus_frangula"

# tree.contMAP <- contMap(tree, intercepts, method = "user")
# plot(tree.contMAP,ftype="off",lwd=c(3,6))
# 
# plot(tree, tip.color = spInfo$mean)
# ?dottree
# dotTree(tree, spInfo$mean)
# # all_tree <- root(all_tree, 'SAL_AB9236AA_AS') # This is an outgroup I picked for the tree.
# # 
# # # Just shrinking some long branches so it's clearer
# # # Don't distort your actual data this way without good reason.
# # all_tree$edge.length[all_tree$edge.length  > 100  ]  <- 100
# # 
# # p1 <- ggtree(all_tree) %<+% info +
# #   geom_tippoint(aes(color=Country)) + # Colour code the tips with country
# #   # Adding in a scale
# #   geom_treescale(x=0, y=45, fontsize=4, linesize=2, offset=2, width=10)
# # 
# # plot(p1)
# # 
# # 
# # 
# # a_sp
# # 
# # ggtree(beast_tree, aes(color=rate)) +
# #   scale_color_continuous(low='darkgreen', high='red') +
# #   theme(legend.position="right")
# # 
# # p1 <- ggtree(tree, aes(color=trait), layout = 'circular', 
# #   ladderize = FALSE, continuous = 'colour', size=2) +
# #   scale_color_gradientn(colours=c("red", 'orange', 'green', 'cyan', 'blue')) +
# #   geom_tiplab(hjust = -.1) + 
# #   xlim(0, 1.2) + 
# #   theme(legend.position = c(.05, .85)) 
# # 
# # ggtree(tree, aes(color=trait), continuous = 'colour', yscale = "trait") + 
# #   scale_color_viridis_c() + theme_minimal()
# # 
# # set.seed(2020)
# # x <- rtree(30)
# # d <- data.frame(label=x$tip.label, var1=abs(rnorm(30)), var2=abs(rnorm(30)))
# # tree <- full_join(x, d, by='label')
# # trs <- list(TREE1 = tree, TREE2 = tree)
# # class(trs) <- 'treedataList'
# # ggtree(trs) + facet_wrap(~.id) + 
# #   geom_tippoint(aes(subset=.id == 'TREE1', colour=var1)) + 
# #   scale_colour_gradient(low='blue', high='red') +  
# #   ggnewscale::new_scale_colour()  + 
# #   geom_tippoint(aes(colour=var2), data=td_filter(.id == "TREE2")) + 
# #   scale_colour_viridis_c()