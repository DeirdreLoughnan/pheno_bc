# started Feb 15, 2023 by Deirdre 

# Aim of this code is to make a cool figure that estimates BB from set conditions and them ranks the early to late bb individuals; do the cues responses similarly vary?
rm(list=ls())
options(stringsAsFactors = FALSE)

# the set cues will be: 12h photoperiod, 20C, high chill - 75/10
require(rstan)
library(forcats)
library(ggdist)
library(reshape2)
require(cowplot)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

# load("output/bb_4sites_phylo_mini.Rda")
# sum <- summary(mdl.4phyloMini)$summary 

#load("output/bb_4sites_phylo_contphotothermo_zscored_Apr19.Rda")
#load("output/bb_phylo_contphotothermo_2zscoredMay13.Rda")
load("output/bb_phylo_contphotothermo_2zscored_oct172024_triple.Rda")
sum <- summary(mdl.3)$summary
post <- rstan::extract(mdl.3)

a_sp = (sum[grep("a_sp", rownames(sum)), 1])
b_photo = sum[grep("b_photo\\[", rownames(sum)), 1]
b_chill = sum[grep("b_chill1\\[", rownames(sum)), 1]
b_force = sum[grep("b_warm\\[", rownames(sum)), 1]

a_sp5 <- vector()
for(i in 1:ncol(post$a_sp)){
  quantU <- round(quantile(post$a_sp[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  a_sp5 <- rbind(a_sp5, quantU)
}
colnames(a_sp5) <- c("Int5","Int95","Int25","Int75")

b_chill5 <- vector()
for(i in 1:ncol(post$b_chill1)){
  quantU <- round(quantile(post$b_chill1[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_chill5 <- rbind(b_chill5, quantU)
}
colnames(b_chill5) <- c("chill5","chill95","chill25","chill75")

b_force5 <- vector()
for(i in 1:ncol(post$b_warm)){
  quantU <- round(quantile(post$b_warm[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_force5 <- rbind(b_force5, quantU)
}
colnames(b_force5) <- c("force5","force95","force25","force75")

b_photo5 <- vector()
for(i in 1:ncol(post$b_photo)){
  quantU <- round(quantile(post$b_photo[,i], c(0.05, 0.95, 0.25, 0.75)),1)
  b_photo5 <- rbind(b_photo5, quantU)
}
colnames(b_photo5) <- c("photo5","photo95","photo25","photo75")

#If we are using the old model, we will use the z-scored values for the parameters
photo <- -0.5033863 #8 h photo
siteSM <- 0
force <- -0.3568628 #5/15 C trt
chill <- -0.3546922 # low chill

m <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(m)){
      m[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM + 
        post$b_warm[it,sp] * force + post$b_photo[it, sp] * photo + post$b_chill[it,sp] * chill +
        post$b_inter_wp[it,sp] * (force*photo) + post$b_inter_wc1[it,sp] * (force*chill) + post$b_inter_pc1[it,sp] * (photo*chill) +
        post$b_inter_s2c1[it,sp] * (chill*siteSM) + post$b_inter_ws2[it,sp] * (force*siteSM) + post$b_inter_ps2[it,sp] * (photo*siteSM) +
        post$b_inter_s3c1[it,sp] * (chill*siteSM) + post$b_inter_ws3[it,sp] * (force*siteSM) + post$b_inter_ps3[it,sp] * (photo*siteSM) +
        post$b_inter_s4c1[it,sp] * (chill*siteSM) + post$b_inter_ws4[it,sp] * (force*siteSM) + post$b_inter_ps4[it,sp] * (photo*siteSM)
}
}

photoHigh <- 0.4965051 #8 h photo
siteSM <- 0
forceHigh <- 0.5877121 
chillHigh <- 0.3660412 # high chill for Smithers

mHigh <- matrix(nrow = 1000, ncol = 47)

for(sp in 1:47){
  for (it in 1:nrow(mHigh)){
    mHigh[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM + 
      post$b_warm[it,sp] * forceHigh + post$b_photo[it, sp] * photoHigh + post$b_chill[it,sp] * chillHigh +
      post$b_inter_wp[it,sp] * (forceHigh*photoHigh) + post$b_inter_wc1[it,sp] * (forceHigh*chillHigh) + post$b_inter_pc1[it,sp] * (photoHigh*chillHigh) +
      post$b_inter_s2c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws2[it,sp] * (forceHigh*siteSM) + post$b_inter_ps2[it,sp] * (photoHigh*siteSM) +
      post$b_inter_s3c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws3[it,sp] * (forceHigh*siteSM) + post$b_inter_ps3[it,sp] * (photoHigh*siteSM) +
      post$b_inter_s4c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws4[it,sp] * (forceHigh*siteSM) + post$b_inter_ps4[it,sp] * (photoHigh*siteSM)
  }
}


spInfo <- read.csv("input/species_list.csv")

spInfo <- spInfo[order(spInfo$species),]
spInfo$meanBB <- colMeans(m)
colnames(m) <- spInfo$species.name

spInfo$meanBBHigh <- colMeans(mHigh)
colnames(mHigh) <- spInfo$species.name

spInfo$Int <- a_sp
spInfo <- cbind(spInfo, a_sp5,b_force5, b_chill5,b_photo5)

spInfo$force <- b_force
spInfo$chill <- b_chill
spInfo$photo <- b_photo

quantile595 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.05, 0.95))
  return(returnQuanilte)
}

quantile75.25 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.75, 0.25))
  return(returnQuanilte)
}

bb_quan <- apply(m, 2, quantile595)
bb_t <- t(bb_quan)
bb_df <- data.frame(bb_t)
colnames(bb_df)[colnames(bb_df) == "X5."] <- "bb5"
colnames(bb_df)[colnames(bb_df) == "X95."] <- "bb95"

bb_quan75.25 <- apply(m, 2, quantile75.25)
bb_t75.25 <- t(bb_quan75.25)
bb_df75.25 <- data.frame(bb_t75.25)
colnames(bb_df75.25)[colnames(bb_df75.25) == "X75."] <- "bb75"
colnames(bb_df75.25)[colnames(bb_df75.25) == "X25."] <- "bb25"
colnames(bb_df)[colnames(bb_df) == "X25."] <- "bb25"
colnames(bb_df)[colnames(bb_df) == "X75."] <- "bb75"

bb_quanHigh <- apply(mHigh, 2, quantile595)
bb_tHigh <- t(bb_quanHigh)
bb_dfHigh <- data.frame(bb_tHigh)
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X5."] <- "bb5High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X95."] <- "bb95High"

bb_quan75.25High <- apply(mHigh, 2, quantile75.25)
bb_t75.25High <- t(bb_quan75.25High)
bb_df75.25High <- data.frame(bb_t75.25High)
colnames(bb_df75.25High)[colnames(bb_df75.25High) == "X75."] <- "bb75High"
colnames(bb_df75.25High)[colnames(bb_df75.25High) == "X25."] <- "bb25High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X25."] <- "bb25High"
colnames(bb_dfHigh)[colnames(bb_dfHigh) == "X75."] <- "bb75High"

spInfo <- cbind(spInfo, bb_df)
spInfo <- cbind(spInfo, bb_df75.25)
spInfo$value <- spInfo$meanBB

spInfo <- cbind(spInfo, bb_dfHigh)
spInfo <- cbind(spInfo, bb_df75.25High)
spInfo$valueHigh <- spInfo$meanBBHigh

m <- data.frame(m)

long <- melt(m)
names(long) <- c("species.name", "valueLow")

mHigh <- data.frame(mHigh)

longHigh <- melt(mHigh)
names(longHigh) <- c("species.name", "valueHigh")

long <- cbind(long, longHigh[,2])

long <- merge(long,spInfo, by = "species.name")

spOrderData <- spInfo[order(spInfo$meanBB),]
spOrder <- as.factor(spOrderData$species.name)

long <- long[order(long$species),]

# longPhotoInfo$mean <- rowMeans(longPhotoInfo[,c("Site1","Site2","Site3","Site4")], na.rm=TRUE)

bChill <- data.frame(post$b_chill1[1:1000,])
colnames(bChill) <- (spInfo$species.name)
longChill <- melt(bChill)
names(longChill) <- c("species.name", "chill")

long <- cbind(long, longChill$chill)

# Add forcing
bForce <- data.frame(post$b_warm[1:1000,])
colnames(bForce) <- (spInfo$species.name)
longForce <- melt(bForce)
names(longForce) <- c("species.name", "force")

long <- cbind(long, longForce$force)

# photoperiod
bPhoto <- data.frame(post$b_photo[1:1000,])
colnames(bPhoto) <- (spInfo$species.name)
longPhoto <- melt(bPhoto)
names(longPhoto) <- c("species.name", "photo")

long <- cbind(long, longPhoto$photo)

# intercept
aSp <- data.frame(post$a_sp[1:1000,])
colnames(aSp) <- (spInfo$species.name)
longInt <- melt(aSp)
names(longInt) <- c("species.name", "int")

long <- cbind(long, longInt$int)

data <- long[order(long$meanBB),]
# 
# data$species.name <- factor(data$species.name, levels=unique(data$species.name) )
#data <- transform(data, variable=reorder(species.name, -meanBB) ) 

names(data) <- c("species.name","valueLow", "valueHigh","species","type","transect","meanBB","meanBBHigh", "Int","Int5","Int95","Int25","Int75","force5","force95","force25","force75","chill5", "chill95", "chill25", "chill75","photo5", "photo95", "photo25", "photo75","spMeanForce", "spMeanChill", "spMeanPhoto","bb5","bb95","bb75","bb25", "spacing","bb5High","bb95High","bb75High","bb25High","valueHigh2","chill", "force","photo","intercept")

#####################################################################
## Eastern spp only #################################################

east <- subset(spInfo, transect != "west")
eastSp <- unique(east$species.name)

dataEast <- data[data$species.name %in% eastSp, ]

overlappingE <- c("aromel","betlen", "betpap", "lyolig","faggra", "betall", "prupen","poptre","rhafra")
spMiniE <- east[!east$species %in% overlappingE,]
spTopE <- east[east$species %in% overlappingE,]

#####################################################################
## Western spp only #################################################

west <- subset(spInfo, transect != "east")
westSp <- unique(west$species.name)

dataWest <- data[data$species.name %in% westSp, ]

overlappingW <- c("spialb","betpap","popbal","rhoalb","alninc")
spMiniW <- west[!west$species %in% overlappingW,]
spTopW <- west[west$species %in% overlappingW,]

meanPtW <- aggregate(dataWest[c("meanBB", "meanBBHigh", "Int")], dataWest[c("species.name","type","transect")], FUN = mean)
names(meanPtW) <- c("species.name","type","transect","Budburst","Buddy","Intercept")

dotW <- ggplot(meanPtW) +
  geom_point(aes(y= Budburst, x = Budburst, colour = "Budburst"), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, colour = "Intercept"), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtW, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) + ylim(-1.5,65) +
  theme(axis.text.x = element_text( size=15,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=20), legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
  labs( x = "Species ordered by mean predicted budburst date", y = "Estimated parameter (days/standardized units)", main = NA) +
  theme(legend.title = element_blank(), legend.text = element_text(size =25), legend.position = "top") +  annotate("text", x = 18, y = 60, label = "b)      Western transect", cex =8) +
  annotate("text", x = spTopW[1,5], y = -1, label = spTopW[1,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -1, label = spTopW[2,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -1, label = spTopW[3,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -1, label = spTopW[4,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -1, label = spTopW[5,2], cex = 5, angle = 78) +
  scale_color_manual(values = c("cyan4", "firebrick4"))

dotWBw <- ggplot(meanPtW) +
  geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst"), size = 5) +
  geom_point(aes(y= Buddy, x = Budburst, shape = "Buddy"), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept"), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtW, col = "black") +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Buddy), data = meanPtW, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) + ylim(-1.5,65) +
  theme(axis.text.x = element_text( size=15,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=20), legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,70)) +
  labs( x = "Species ordered by mean predicted budburst date", y = "Day of budburst (days/standardized units)", main = NA) +
  theme(legend.position = "none") +  #annotate("text", x = 21.5, y = 60, label = "b)      Western transect", cex =8) +
  annotate("text", x = spTopW[1,5], y = -1, label = spTopW[1,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -1, label = spTopW[2,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -1, label = spTopW[3,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -1, label = spTopW[4,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -1, label = spTopW[5,2], cex = 5, angle = 78) +
  scale_shape_manual(values = c(1,2,16)) +
  annotate("text", x = 68, y = 62.5, label = "Low cues", cex = 7) +  
  annotate("text", x = 68, y = 46.5, label = "Intercept", cex = 7) +
  annotate("text", x = 68, y = 29.5, label = "High cues", cex = 7) + 
  geom_segment(aes(x = 63, y = 62.5, xend = 64.5 , yend = 62.5)) +
  geom_segment(aes(x = 63, y = 46.5, xend = 64.5 , yend = 46.5)) +
  geom_segment(aes(x = 63, y = 29.5, xend = 64.5 , yend = 29.5))

# +  annotate("text", x = 62, y = 39, label = "Intercept", cex = 3) +
#   geom_segment(aes(x = 62, y = 60, xend = 59 , yend = 60)) +
#                  geom_segment(aes(x = 62, y = 39, xend = 59 , yend = 39))
dotWBw

ggplot(meanPtW) +
  geom_point(aes(y= Budburst, x = Intercept, shape = "Budburst"), size = 5) +
  geom_point(aes(y= Buddy, x = Intercept, shape = "Buddy"), size = 5) +
  geom_point(aes(y= Intercept, x = Intercept, shape = "Intercept"), size = 5) +
  geom_segment(aes(x = Intercept, y = Intercept, xend = Intercept, yend = Budburst), data = meanPtW, col = "black") +
  geom_segment(aes(x = Intercept, y = Intercept, xend = Intercept, yend = Buddy), data = meanPtW, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) + ylim(-1.5,65) +
  theme(axis.text.x = element_text( size=15,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=20), legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$Int, labels = spMiniW$species) +
  labs( x = "Species ordered by mean predicted budburst date", y = "Day of budburst (days/standardized units)", main = NA) +
  theme(legend.position = "none") +  annotate("text", x = 18, y = 60, label = "b)      Western transect", cex =8) +
  annotate("text", x = spTopW[1,7], y = -1, label = spTopW[1,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[2,7], y = -1, label = spTopW[2,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[3,7], y = -1, label = spTopW[3,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[4,7], y = -1, label = spTopW[4,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[5,7], y = -1, label = spTopW[5,2], cex = 5, angle = 78) +

  scale_shape_manual(values = c(1,2,16)) +
  annotate("text", x = 66, y = 60.5, label = "High cues", cex = 7) +  

  annotate("text", x = 66, y = 44.5, label = "Intercept", cex = 7) +
  annotate("text", x = 66, y = 27.5, label = "High cues", cex = 7) + 
  geom_segment(aes(x = 61, y = 60.5, xend = 62.5 , yend = 60.5)) +
  geom_segment(aes(x = 61, y = 44.5, xend = 62.5 , yend = 44.5)) +
  geom_segment(aes(x = 61, y = 27.5, xend = 62.5 , yend = 27.5))

#### Eastern plot
meanPtE <- aggregate(dataEast[c("meanBB", "meanBBHigh","Int")], dataEast[c("species.name","type","transect")], FUN = mean)
names(meanPtE) <- c("species.name","type","transect","Budburst", "Buddy","Intercept")

dotE <- ggplot(meanPtE) +
  geom_point(aes(y= Budburst, x = Budburst, colour = "Budburst"), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, colour = "Intercept"), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtE, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  ylim(-1.5,65) + theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=20)) +
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(24,64)) +
  labs( x = "Species ordered by mean predicted budburst date", y = "Estimated parameter (days/standardized units)", main = NA) +
  theme(legend.title = element_blank(), legend.text = element_text(size =25), legend.position = "top") +  annotate("text", x = 29, y = 60, label = "a)      Eastern transect", cex =8) +
  annotate("text", x = spTopE[1,5], y = -1, label = spTopE[1,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -1, label = spTopE[2,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -1, label = spTopE[3,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[4,5], y = -1, label = spTopE[4,2], cex = 5, angle = 78) +
 # annotate("text", x = 38, y = 0, label = spTopE[4,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -1, label = spTopE[5,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[6,5], y = -1, label = spTopE[6,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[7,5], y = -1, label = spTopE[7,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[8,5], y = -1, label = spTopE[8,2], cex = 5, angle = 78) +
  scale_color_manual(values = c("forestgreen","orange3"))


dotEBw <- ggplot(meanPtE) +
  geom_point(aes(y= Budburst, x = Budburst, shape = "Budburst"), size = 5) +
  geom_point(aes(y= Buddy, x = Budburst, shape = "Buddy"), size = 5) +
  geom_point(aes(y= Intercept, x = Budburst, shape = "Intercept"), size = 5) +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Budburst), data = meanPtE, col = "black") +
  geom_segment(aes(x = Budburst, y = Intercept, xend = Budburst, yend = Buddy), data = meanPtE, col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white")) +
  ylim(-1.5,65) + theme(axis.text.x = element_text( size=15, angle = 78,  hjust=1),
                        axis.text.y=element_text(size = 15),
                        axis.title=element_text(size=20)) +
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(24,70)) +
  labs( x = "", y = "Day of budburst (days/standardized units)", main = NA) +
  theme(legend.position = "none") +  #annotate("text", x = 29, y = 60, label = "Eastern transect", cex =8) +#annotate("text", x = 29, y = 60, label = "a)      Eastern transect", cex =8) +
  annotate("text", x = spTopE[1,5], y = -1, label = spTopE[1,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -1, label = spTopE[2,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -1, label = spTopE[3,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[4,5], y = -1, label = spTopE[4,2], cex = 5, angle = 78) +
  # annotate("text", x = 38, y = 0, label = spTopE[4,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -1, label = spTopE[5,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[6,5], y = -1, label = spTopE[6,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[7,5], y = -1, label = spTopE[7,2], cex = 5, angle = 78) +
  annotate("text", x = spTopE[8,5], y = -1, label = spTopE[8,2], cex = 5, angle = 78) +
  scale_shape_manual(values = c(1,2,16)) +       
  annotate("text", x = 69, y = 65, label = "Low cues", cex = 7) +  
  annotate("text", x = 69, y = 50, label = "Intercept", cex = 7) +
  annotate("text", x = 69, y = 35, label = "High cues", cex = 7) + 
  geom_segment(aes(x = 65, y = 65, xend = 66 , yend = 65)) +
  geom_segment(aes(x = 65, y = 50, xend = 66.5 , yend = 50)) +
  geom_segment(aes(x = 65, y = 35, xend = 66 , yend = 35))
dotEBw


# pdf("figures/dotBetaAlphaLongColor.pdf", width = 12, height = 16)
# plot_grid(dotE, dotW, nrow = 2, ncol = 1, align = "v")
# dev.off()
pdf("figures/dotBetaAlphaLongBW_LHa.pdf", width = 12, height = 8)
dotEBw
dev.off()

pdf("figures/dotBetaAlphaLongBW_LHb.pdf", width = 12, height = 8)
dotWBw
dev.off()

pdf("figures/dotBetaAlphaLongBW_LH_April20.pdf", width = 12, height = 16)
plot_grid(dotEBw, dotWBw, nrow = 2, ncol = 1, align = "v")
dev.off()
#############################################
# old plots
####### Old plot not spaced out ####################3

overlapping <- c("aromel","prupen","vacmyr","spialb","ilemuc", "vibcas","betpap","betlen","symalb","acerub","lyolig","rhopri")
spMini <- spInfo[!spInfo$species %in% overlapping,]
spTop <- spInfo[spInfo$species %in% overlapping,]

#########################################################
###### Remake fig but with dots and bars ################

chillPtW <- aggregate(dataWest[c("meanBB", "Int","chill","force","photo","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")], dataWest[c("species.name","type","transect")], FUN = mean)
# names(chillPtW) <- c("species.name","type","transect","Budburst","Intercept","Int5", "Int95", "Int25", "Int75","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")


chilldotW <- ggplot(chillPtW, aes(y= chill, x = meanBB, color = transect), size = 7) + geom_point(size =7, shape = 15) +
  geom_errorbar(aes(ymin= chill5, ymax = chill95,xmin= meanBB, xmax = meanBB,color = transect), width= 0, size =0.5) +
  geom_errorbar(aes(ymin= chill25, ymax = chill75,xmin= meanBB, xmax = meanBB,color = transect), width= 0, size =0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-40,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size= 17),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
  labs( x = "", y = 
         ""# "Days per standardized chill portion"
        , main = NA) +
  theme(legend.title = element_blank()) +  #annotate("text", x = 15, y = 2, label = "b)", cex =10) +
  annotate("text", x = 38, y = 3, label = "Western transect", cex =10) +
  annotate("text", x = spTopW[1,5], y = -37, label = spTopW[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -37, label = spTopW[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -37, label = spTopW[3,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -37, label = spTopW[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -37, label = spTopW[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[6,5], y = -37, label = spTopW[6,2], cex = 6, angle = 78) +
  scale_fill_manual(values = c("darkred","#cc6a70ff","#cc6a70ff")) + 
  scale_color_manual(values = c("darkred","#cc6a70ff","#cc6a70ff")) + 
  geom_segment(aes(x = 15, y = 5, xend = 15 , yend = -2),
  arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 15, y = -3, label = "Earlier", cex =5) 
chilldotW

pdf("figures/dotCFPEWSizeColorUnitsShapeb.pdf", width = 10, height = 7)
chilldotW
dev.off()

forcedotW <- ggplot(chillPtW,aes(y= force, x = meanBB, color = transect), size = 7) +
  geom_point(size = 7, shape = 15) +
  geom_errorbar(aes(ymin= force5, ymax = force95,xmin= meanBB, xmax = meanBB, color = transect), width= 0, size =0.5) +
  geom_errorbar(aes(ymin= force25, ymax = force75,xmin= meanBB, xmax = meanBB, color = transect), width= 0, size =1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-25,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
  labs( x = "", 
        y = ""# "Days per standardized forcing" )+ # expression("Forcing response (days/"*~degree*C*")")
  , main = NA) +
  theme(legend.title = element_blank()) +  #annotate("text", x = 15, y = 2, label = "d)", cex =10) +
  # annotate("text", x = 38, y = 10, label = "Western transect", cex =5) +
  annotate("text", x = spTopW[1,5], y = -23, label = spTopW[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -23, label = spTopW[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -23, label = spTopW[3,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -23, label = spTopW[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -23, label = spTopW[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[6,5], y = -23, label = spTopW[6,2], cex = 6, angle = 78) + 
  scale_fill_manual(values = c("tan4","#f9b641ff","#f9b641ff")) + 
  scale_color_manual(values = c("tan4","#f9b641ff","#f9b641ff")) + 
  geom_segment(aes(x = 15, y = 5, xend = 15 , yend = 1),
    arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 15, y = 0, label = "Earlier", cex =5) 
forcedotW

pdf("figures/dotCFPEWSizeColorUnitsShaped.pdf", width = 10, height = 7)
forcedotW
dev.off()

photodotW <- ggplot(chillPtW,aes(y= photo, x = meanBB, color = transect), size = 7) +
  geom_point(size = 7, shape =15) +
  geom_errorbar(aes(ymin= photo5, ymax = photo95,xmin= meanBB, xmax = meanBB, color = transect), width= 0, size = 0.5) +
  geom_errorbar(aes(ymin= photo25, ymax = photo75,xmin= meanBB, xmax = meanBB, color = transect), width= 0, size =1.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-10,1) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=21),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,65)) +
  labs( x = "Species ordered by predicted budburst date", y = "",#"Days per standardized photoperiod",
        main = NA) +
  theme(legend.title = element_blank()) +  #annotate("text", x = 15, y = 0, label = "f)", cex =10) +
  # annotate("text", x = 38, y = 10, label = "Western transect", cex =5) +
  annotate("text", x = spTopW[1,5], y = -9.3, label = spTopW[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -9.3, label = spTopW[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -9.3, label = spTopW[3,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -9.3, label = spTopW[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -9.3, label = spTopW[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[6,5], y = -9.3, label = spTopW[6,2], cex = 6, angle = 78)+ 
  scale_fill_manual(values = c("navy","cyan4", "cyan4")) + 
  scale_color_manual(values = c("navy","cyan4", "cyan4")) + 
  geom_segment(aes(x = 15, y = 0.5, xend =15 , yend = -0.7),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 15, y = -1, label = "Earlier", cex =5) 
photodotW

pdf("figures/dotCFPEWSizeColorUnitsShapef.pdf", width = 10, height = 7)
photodotW
dev.off()

chillPtE <- aggregate(dataEast[c("meanBB", "Int","chill","force","photo","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")], dataEast[c("species.name","type","transect")], FUN = mean)
# names(chillPtE) <- c("species.name","type","transect","Budburst","Intercept","Chilling","Forcing","Photoperiod","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")


chilldotE <- ggplot(chillPtE,aes(y= chill, x = meanBB, color = transect), size = 7) +
  geom_point(size = 7,  shape  = 17) +
  geom_errorbar(aes(ymin= chill5, ymax = chill95,xmin= meanBB, xmax = meanBB, color = transect), width= 0, size = 0.5) + geom_errorbar(aes(ymin= chill25, ymax = chill75,xmin= meanBB, xmax = meanBB,color = transect), width= 0, size = 1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-45,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=21),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
  labs( x = "", y = "Chill response (days/standardized unit)"#"Days per standardized chill portion"
        , main = NA) +
  theme(legend.title = element_blank()) +  #annotate("text", x =23, y = 2, label = "a)", cex =10) +
  annotate("text", x = 45, y = 0, label = "Eastern transect", cex = 10) +
  annotate("text", x = spTopE[1,5], y = -41.5, label = spTopE[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -41.5, label = spTopE[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -41.5, label = spTopE[3,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[4,5], y = -41.5, label = spTopE[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -41.5, label = spTopE[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[6,5], y = -41.5, label = spTopE[6,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[7,5], y = -41.5, label = spTopE[7,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[8,5], y = -41.5, label = spTopE[8,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[9,5], y = -41.5, label = spTopE[9,2], cex = 6, angle = 78) +
  scale_fill_manual(values = c("#cc6a70ff","darkred","#cc6a70ff")) + 
  scale_color_manual(values = c("#cc6a70ff","darkred","#cc6a70ff")) + 
  geom_segment(aes(x = 23, y = 5, xend = 23 , yend = -2),arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 23, y = -3, label = "Earlier", cex =5) 
chilldotE

pdf("figures/dotCFPEWSizeColorUnitsShapea.pdf", width = 10, height = 7)
chilldotE
dev.off()

forcedotE <- ggplot(chillPtE,aes(y= force, x = meanBB, colour = transect), size = 7) +
  geom_point(size = 7, shape  = 17) +
  geom_errorbar(aes(ymin= force5, ymax = force95,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 0.5) +
  geom_errorbar(aes(ymin= force25, ymax = force75,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +  ylim(-25,5) +
  theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=21),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
  labs( x = "", y = "Forcing response (days/standardized unit)"#"Days per standardized forcing")+#expression("Forcing response (days/"*~degree*C*")")
  , main = NA) +
  theme(legend.title = element_blank()) +  #annotate("text", x = 23, y = 2, label = "c)", cex =10) +
  annotate("text", x = spTopE[1,5], y = -23.5, label = spTopE[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -23.5, label = spTopE[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -23.5, label = spTopE[3,2], cex = 6, angle = 78) +
  annotate("text", x =  spTopE[4,5], y = -23.5, label = spTopE[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -23.5, label = spTopE[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[6,5], y = -23.5, label = spTopE[6,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[7,5], y = -23.5, label = spTopE[7,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[8,5], y = -23.5, label = spTopE[8,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[9,5], y = -23.5, label = spTopE[9,2], cex = 6, angle = 78) +
  scale_fill_manual(values = c("#f9b641ff","tan4","#f9b641ff")) + 
  scale_color_manual(values = c("#f9b641ff","tan4","#f9b641ff")) + 
  geom_segment(aes(x = 23, y = 5, xend =23 , yend = 1),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 23, y = 0, label = "Earlier", cex =5) 
forcedotE

pdf("figures/dotCFPEWSizeColorUnitsShapec.pdf", width = 10, height = 7)
forcedotE
dev.off()

photodotE <- ggplot(chillPtE,aes(y= photo, x = meanBB, colour = transect), size = 7) +
  geom_point(size = 7, shape  = 17) +
  geom_errorbar(aes(ymin= photo5, ymax = photo95,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 0.5) +
  geom_errorbar(aes(ymin= photo25, ymax = photo75,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 1.5) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +  ylim(-10,1) +
  theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=  21),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
  labs( x = "Species ordered by predicted budburst date", y = "Photoperiod response (days/standardized unit)"#"Days per standardized photoperiod"
        , main = NA) +
  theme(legend.title = element_blank()) +  #annotate("text", x = 23, y = 0, label = "e)", cex =10) + 
  annotate("text", x = spTopE[1,5], y = -9.3, label = spTopE[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -9.3, label = spTopE[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -9.3, label = spTopE[3,2], cex = 6, angle = 78) +
  annotate("text", x =  spTopE[4,5], y = -9.3, label = spTopE[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -9.3, label = spTopE[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[6,5], y = -9.3, label = spTopE[6,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[7,5], y = -9.3, label = spTopE[7,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[8,5], y = -9.3, label = spTopE[8,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[9,5], y = -9.3, label = spTopE[9,2], cex = 6, angle = 78) +
  scale_fill_manual(values = c("cyan4", "navy","cyan4")) + 
  scale_color_manual(values = c("cyan4","navy", "cyan4")) + 
  geom_segment(aes(x = 23, y = 0.5, xend =23 , yend = -0.7),
    arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 23, y = -1, label = "Earlier", cex =5) 
photodotE

pdf("figures/dotCFPEWSizeColorUnitsShapee.pdf", width = 10, height = 7)
photodotE
dev.off()

pdf("figures/dotCFPEWSizeColorUnitsShape.pdf", width = 20, height = 20)
plot_grid(chilldotE, chilldotW, forcedotE, forcedotW, photodotE, photodotW, nrow = 3, ncol = 2, align = "v")
dev.off()

## How does species rank change with high low intercept?

rank <- spInfo[,c("species.name","species","type","transect","meanBB","meanBBHigh","Int")]
rank <- rank[order(rank$Int),]
rank$rankInt <- seq(1:nrow(rank))

rank <- rank[order(rank$meanBB),]
rank$rankLowC <- seq(1:nrow(rank))

rank <- rank[order(rank$meanBBHigh),]
rank$rankHighC <- seq(1:nrow(rank))

# pdf("figures/rankEstiBB.pdf", width = 9, height =3)
# colTran <- c("maroon","navy","forestgreen")
# par(mfrow = c(1,3), mar = c(5.1, 4.8, 4.1, 2.1))
# plot(rank$rankHighC~rank$rankInt, 
#      col = colTran[factor(rank$transect)], 
#      pch = 19,
#      xlab = "Intercept rank",
#      ylab = "High cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# plot(rank$rankLowC~rank$rankInt, 
#      col = colTran[factor(rank$transect)], 
#      pch = 19,
#      xlab = "Intercept rank",
#      ylab = "Low cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# plot(rank$rankHighC~rank$rankLowC, 
#      col = colTran[factor(rank$transect)], 
#      pch = 19,
#      xlab = "Low cue rank",
#      ylab = "High cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# legend("topleft",legend = c( expression("East"),
#                             expression("West"),
#                             expression("Both")),
#        col = c("maroon","forestgreen","navy"),
#        #pt.bg = c("#042333ff","#cc6a70ff","#593d9cff","#f9b641ff","#13306dff","#efe350ff","#eb8055ff"),
#        pt.bg = c( "maroon","forestgreen","navy"),
#        inset = 0.02, pch = c(21, 21, 21 ), cex = 1.25, bty = "n")
# dev.off()
# 
# colType <- c("maroon","navy")
# par(mfrow = c(1,3), mar = c(5.1, 4.8, 4.1, 2.1))
# plot(rank$rankHighC~rank$rankInt, 
#      col = colType[factor(rank$type)], 
#      pch = 19,
#      xlab = "Intercept rank",
#      ylab = "High cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# plot(rank$rankLowC~rank$rankInt, 
#      col = colType[factor(rank$type)], 
#      pch = 19,
#      xlab = "Intercept rank",
#      ylab = "Low cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)
# 
# plot(rank$rankHighC~rank$rankLowC, 
#      col = colType[factor(rank$type)], 
#      pch = 19,
#      xlab = "Low cue rank",
#      ylab = "High cue rank",
#      cex.lab =1.5,
#      cex =1.5)
# abline(0,1)

### Rank within transect:
rankE <- subset(rank, transect != "west")

rankE <- rankE[order(rankE$Int),]
rankE$rankInt <- seq(1:nrow(rankE))

rankE <- rankE[order(rankE$meanBB),]
rankE$rankLowC <- seq(1:nrow(rankE))

rankE <- rankE[order(rankE$meanBBHigh),]
rankE$rankHighC <- seq(1:nrow(rankE))

pdf("figures/rankstiBBTransect.pdf", width = 12, height = 8)
par(mfrow = c(2,3), mar = c(5.1, 4.8, 4.1, 2.1))
colTran <- c("maroon","goldenrod")
#par(mfrow = c(1,1), mar = c(5.1, 4.8, 4.1, 2.1))
plot(rankE$rankHighC~rankE$rankInt, 
  col = colTran[factor(rankE$transect)], 
  pch = 19,
  xlab = "Intercept rank",
  ylab = "High cue rank",
  main = "",
  cex.lab =1.5,
  xlim = c(0,30),
  ylim = c(0,30),
  cex = 2, cex.lab = 2, cex.axis = 2,
  axes = F)
abline(0,1)
axis(side=1, at=seq(-5,30,5), cex.axis = 1.75)
axis(side=2,at=seq(-5,30,5), cex.axis = 1.75)
title("Eastern Sites", adj = 0, cex.main = 2)
text(0.8,29, label = "(a)", cex = 1.5,cex.main = 2,)
text(26,18, label = expression(italic("Populus")), cex = 1.25)
text(7,14, label = expression(italic("Betula")), cex = 1.25)
text(10,5, label = expression(italic("Alnus")), cex = 1.25)
#text(10,30, label = "Eastern Communities", cex = 1.75)
legend("topleft",legend = c( expression("Eastern"),
  expression("Both")),
  col = c("maroon","goldenrod"),
  #pt.bg = c("#042333ff","#cc6a70ff","#593d9cff","#f9b641ff","#13306dff","#efe350ff","#eb8055ff"),
  pt.bg = c( "maroon","goldenrod"),
  inset = 0.05, pch = c(21, 21), cex = 1.5, bty = "n")
#dev.off()

#pdf("figures/rankstiBBTransectb.pdf", width = 6, height = 6)
#par( mar = c(5.1, 4.8, 4.1, 2.1))
plot(rankE$rankLowC~rankE$rankInt, 
  col = colTran[factor(rankE$transect)], 
  pch = 19,
  xlab = "Intercept rank",
  ylab = "Low cue rank",
  cex.lab =1.5,
  xlim = c(0,30),
  ylim = c(0,30),
  cex = 2, cex.lab = 2, cex.axis = 2,axes = F)
abline(0,1)
axis(side=1, at=seq(-5,30,5), cex.axis = 1.75)
axis(side=2,at=seq(-5,30,5), cex.axis = 1.75)
text(0.8,29, label = "(b)", cex = 1.5)
text(26,19, label = expression(italic("Populus")), cex = 1.25)
text(7,13.5, label = expression(italic("Betula")), cex = 1.25)
text(10,5.5, label = expression(italic("Alnus")), cex = 1.25)
#dev.off()

#pdf("figures/rankstiBBTransectc.pdf", width = 6, height = 6)
#par( mar = c(5.1, 4.8, 4.1, 2.1))
plot(rankE$rankHighC~rankE$rankLowC, 
  col = colTran[factor(rankE$transect)], 
  pch = 19,
  xlim = c(0,30),
  ylim = c(0,30),
  xlab = "Low cue rank",
  ylab = "High cue rank",
  cex.lab =1.5,
  cex = 2, cex.lab = 2, cex.axis = 2, axes = F)
abline(0,1)
axis(side=1, at=seq(-5,30,5), cex.axis = 1.75)
axis(side=2,at=seq(-5,30,5), cex.axis = 1.75)
text(0.8,29, label = "(c)", cex = 1.5)
text(26,18, label = expression(italic("Populus")), cex = 1.25)
text(7,14, label = expression(italic("Betula")), cex = 1.25)
text(10,5, label = expression(italic("Alnus")), cex = 1.25)

# legend("topleft",legend = c( expression("Eastern"),
#   expression("Both")),
#   col = c("maroon","goldenrod"),
#   #pt.bg = c("#042333ff","#cc6a70ff","#593d9cff","#f9b641ff","#13306dff","#efe350ff","#eb8055ff"),
#   pt.bg = c( "maroon","goldenrod"),
#   inset = 0.05, pch = c(21, 21), cex = 1.5, bty = "n")
#dev.off()

rankW <- subset(rank, transect != "east")

rankW <- rankW[order(rankW$Int),]
rankW$rankInt <- seq(1:nrow(rankW))

rankW <- rankW[order(rankW$meanBB),]
rankW$rankLowC <- seq(1:nrow(rankW))

rankW <- rankW[order(rankW$meanBBHigh),]
rankW$rankHighC <- seq(1:nrow(rankW))

colTran <- c("goldenrod", "cyan4")


#pdf("figures/rankstiBBTransectd.pdf", width = 6, height = 6)
#par( mar = c(5.1, 4.8, 4.1, 2.1))
plot(rankW$rankHighC~rankW$rankInt, 
  col = colTran[factor(rankW$transect)], 
  pch = 19,
  xlab = "Intercept rank",
  ylab = "High cue rank",
  main = "",
  cex.lab =2,
  cex = 2, cex.axis = 2, cex.main = 2,
  xlim = c(0,30),
  ylim = c(0,30),axes = F)
abline(0,1)
axis(side=1, at=seq(-5,30,5), cex.axis = 1.75)
axis(side=2,at=seq(-5,30,5), cex.axis = 1.75)
title("Western Sites", adj = 0, cex.main = 2)
text(0.8,29, label = "(d)", cex = 1.5)
text(24,17.5, label = expression(italic("Populus")), cex = 1.25)
text(12.2,20, label = expression(italic("Betula")), cex = 1.25)
text(16.3,11.5, label = expression(italic("Alnus")), cex = 1.25)
#text(10,30, label = "Western Communities", cex = 1.75)

legend("topleft",legend = c( expression("Western"),
  expression("Both")),
  col = c("cyan4","goldenrod"),
  #pt.bg = c("#042333ff","#cc6a70ff","#593d9cff","#f9b641ff","#13306dff","#efe350ff","#eb8055ff"),
  pt.bg = c( "cyan4","goldenrod"),
  inset = 0.05, pch = c(21, 21), cex = 1.5, bty = "n")


#dev.off()

#pdf("figures/rankstiBBTransecte.pdf", width = 6, height = 6)
#par( mar = c(5.1, 4.8, 4.1, 2.1))
plot(rankW$rankLowC~rankW$rankInt, 
  col = colTran[factor(rankW$transect)], 
  pch = 19,
  xlab = "Intercept rank",
  ylab = "Low cue rank",
  cex.lab =1.5,
  cex = 2, cex.lab = 2, cex.axis = 2,
  xlim = c(0,30),
  ylim = c(0,30),axes = F)
abline(0,1)
axis(side=1, at=seq(-5,30,5), cex.axis = 1.75)
axis(side=2,at=seq(-5,30,5), cex.axis = 1.75)
text(0.8,29, label = "(e)", cex = 1.5)
text(24,17.5, label = expression(italic("Populus")), cex = 1.25)
text(19,11.5, label = expression(italic("Betula")), cex = 1.25)
text(15,9, label = expression(italic("Alnus")), cex = 1.25)
#dev.off()

#pdf("figures/rankstiBBTransectf.pdf", width = 6, height = 6)
#par( mar = c(5.1, 4.8, 4.1, 2.1))
plot(rankW$rankHighC~rankW$rankLowC, 
  col = colTran[factor(rankW$transect)], 
  pch = 19,
  xlab = "Low cue rank",
  ylab = "High cue rank",
  cex.lab =1.5,
  cex = 2, cex.lab = 2, cex.axis = 2,
  xlim = c(0,30),
  ylim = c(0,30),axes = F)
abline(0,1)
axis(side=1, at=seq(-5,30,5), cex.axis = 1.75)
axis(side=2,at=seq(-5,30,5), cex.axis = 1.75)
text(0.8,29, label = "(f)", cex = 1.5)
text(24,17.5, label = expression(italic("Populus")), cex = 1.25)
text(7,17.5, label = expression(italic("Betula")), cex = 1.25)
text(10,13, label = expression(italic("Alnus")), cex = 1.25)

dev.off()
