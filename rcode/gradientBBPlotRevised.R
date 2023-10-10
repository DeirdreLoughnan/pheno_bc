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
load("output/bb_phylo_contphotothermo_2zscoredMay13.Rda")
sum <- summary(mdl.2z)$summary
post <- rstan::extract(mdl.2z)


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
head(spInfo)
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

names(data) <- c("species.name","valueLow", "valueHigh","species","type","transect","meanBB","meanBBHigh", "Int","Int5","Int95","Int25","Int75","force5","force95","force25","force75","chill5", "chill95", "chill25", "chill75","photo5", "photo95", "photo25", "photo75","spMeanForce", "spMeanChill", "spMeanPhoto","bb5","bb95","bb75","bb25", "spacing","bb5High","bb95High","bb75High","bb25High","spacingHigh","chill", "force","photo","intercept")

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

# Making lizzie's other suggested plot - 2 connected dots:

meanPtW <- aggregate(dataWest[c("meanBB", "meanBBHigh", "Int")], dataWest[c("species.name","type","transect")], FUN = mean)
names(meanPtW) <- c("species.name","type","transect","Budburst","Buddy","Intercept")

#### Eastern plot
meanPtE <- aggregate(dataEast[c("meanBB", "meanBBHigh","Int")], dataEast[c("species.name","type","transect")], FUN = mean)
names(meanPtE) <- c("species.name","type","transect","Budburst", "Buddy","Intercept")

overlapping <- c("aromel","prupen","vacmyr","spialb","ilemuc", "vibcas","betpap","betlen","symalb","acerub","lyolig","rhopri")
spMini <- spInfo[!spInfo$species %in% overlapping,]
spTop <- spInfo[spInfo$species %in% overlapping,]


#########################################################
###### Remake fig but with dots and bars ################

chillPtW <- aggregate(dataWest[c("meanBB", "Int","chill","force","photo","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")], dataWest[c("species.name","type","transect")], FUN = mean)
# names(chillPtW) <- c("species.name","type","transect","Budburst","Intercept","Int5", "Int95", "Int25", "Int75","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")
chillPtW$colT <- "#cc6a70ff"
both <- c("Betula_papyrifera", "Alnus_incana","Populus_tremuloides")
chillPtWBoth <- chillPtW[chillPtW$species.name %in% both,]

chilldotW <- ggplot(chillPtW,aes(y= chill, x = meanBB), size = 7) + geom_point(aes(colour = transect),size =7, shape = 15)  + 
  geom_errorbar(aes(ymin= chill5, ymax = chill95,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 0.5)+
  geom_errorbar(aes(ymin= chill25, ymax = chill75,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 1.5) + 
  geom_point(aes(colour = transect),size =7, shape = 15) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-40,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size= 17),
        legend.position = "none") + scale_color_manual(values = c("east/west"="darkred",
                                                                  "west"="#cc6a70ff"),
        breaks = c("east/west","west"), label = c("High", "Low")) +
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
  labs( x = "", y = 
          "Chill response (days/standardized unit)"# "Days per standardized chill portion"
        , main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 2, label = "b)", cex =10) +
  annotate("text", x = 38, y = 3, label = "Western transect", cex =10) +
  annotate("text", x = spTopW[1,5], y = -37, label = spTopW[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -37, label = spTopW[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -37, label = spTopW[3,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -37, label = spTopW[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -37, label = spTopW[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[6,5], y = -37, label = spTopW[6,2], cex = 6, angle = 78) +
  scale_fill_manual(values = c("darkred","darkred")) + 
  geom_segment(aes(x = 62, y = 0, xend = 62 , yend = -7),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 62, y = -9, label = "Earlier", cex =5) 
chilldotW

forcedotW <- ggplot(chillPtW,aes(y= force, x = meanBB), size = 7) +
  geom_point(aes(colour = transect),size = 7, shape = 15) +
  geom_errorbar(aes(ymin= force5, ymax = force95,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size =0.5) +
  geom_errorbar(aes(ymin= force25, ymax = force75,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size =1.5) +
  geom_point(aes(colour = transect),size =7, shape = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-25,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none")  + scale_color_manual(values = c("east/west"="darkorange3",
                                                                   "west"="#f9b641ff"),
                                                        breaks = c("east/west","west"), label = c("High", "Low")) + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
  labs( x = "", 
        y = "Forcing response (days/standardized unit)"# "Days per standardized forcing" )+ # expression("Forcing response (days/"*~degree*C*")")
        , main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 2, label = "d)", cex =10) +
  # annotate("text", x = 38, y = 10, label = "Western transect", cex =5) +
  annotate("text", x = spTopW[1,5], y = -23, label = spTopW[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -23, label = spTopW[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -23, label = spTopW[3,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -23, label = spTopW[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -23, label = spTopW[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[6,5], y = -23, label = spTopW[6,2], cex = 6, angle = 78) + 
  geom_segment(aes(x = 62, y = 2, xend = 62 , yend = -3),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 62, y = -5, label = "Earlier", cex =5) 
forcedotW

photodotW <- ggplot(chillPtW,aes(y= photo, x = meanBB, colour = "photo"), size = 7) +
  geom_point(aes(colour = transect),size = 7, shape =17) +
  geom_errorbar(aes(ymin= photo5, ymax = photo95,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 0.5) +
  geom_errorbar(aes(ymin= photo25, ymax = photo75,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size =1.5) +
  geom_point(aes(colour = transect),size =7, shape = 17) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylim(-10,0) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,65)) +
  labs( x = "Species ordered by predicted budburst date", y = "Photoperiod response (days/standardized unit)",#"Days per standardized photoperiod",
        main = NA) + scale_color_manual(values = c("east/west"="navy",
                                                   "west"="cyan4"),
                                        breaks = c("east/west","west"), label = c("High", "Low")) +
  theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 0, label = "f)", cex =10) +
  # annotate("text", x = 38, y = 10, label = "Western transect", cex =5) +
  annotate("text", x = spTopW[1,5], y = -9.3, label = spTopW[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -9.3, label = spTopW[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -9.3, label = spTopW[3,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -9.3, label = spTopW[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -9.3, label = spTopW[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopW[6,5], y = -9.3, label = spTopW[6,2], cex = 6, angle = 78)+ 
  geom_segment(aes(x = 63, y = -0.2, xend =63 , yend = -2.3),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 63, y = -2.7, label = "Earlier", cex =5) 
photodotW

chillPtE <- aggregate(dataEast[c("meanBB", "Int","chill","force","photo","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")], dataEast[c("species.name","type","transect")], FUN = mean)
# names(chillPtE) <- c("species.name","type","transect","Budburst","Intercept","Chilling","Forcing","Photoperiod","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")


chilldotE <- ggplot(chillPtE,aes(y= chill, x = meanBB, colour ="#cc6a70ff"), size = 7) +
  geom_point(aes(colour = transect), size = 7, shape  = 17) +
  geom_errorbar(aes(ymin= chill5, ymax = chill95,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 0.5) + geom_errorbar(aes(ymin= chill25, ymax = chill75,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 1.5) +
  geom_point(aes(colour = transect),size =7, shape = 17) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-45,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none") + scale_color_manual(values = c("east/west"="darkred",
                                                                  "east"="#cc6a70ff"),
                                                       breaks = c("east/west","west"), label = c("High", "Low")) + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
  labs( x = "", y = "Chill response (days/standardized unit)"#"Days per standardized chill portion"
        , main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x =23, y = 2, label = "a)", cex =10) +
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
  scale_fill_manual(values = c("#cc6a70ff","#cc6a70ff"))+
  geom_segment(aes(x = 65, y = -2, xend = 65 , yend = -9),arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 65, y = -11, label = "Earlier", cex =5) 

chilldotE

forcedotE <- ggplot(chillPtE,aes(y= force, x = meanBB, colour = "#f9b641ff"), size = 7) +
  geom_point(aes(colour = transect),size = 7, shape  = 17) +
  geom_errorbar(aes(ymin= force5, ymax = force95,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 0.5) +
  geom_errorbar(aes(ymin= force25, ymax = force75,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 1.5) +
  geom_point(aes(colour = transect),size =7, shape = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +  ylim(-25,5) +
  theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none") + scale_color_manual(values = c("east/west"="darkorange3",
                                                                  "east"="#f9b641ff"),
                                                       breaks = c("east/west","west"), label = c("High", "Low")) + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
  labs( x = "", y = "Forcing response (days/standardized unit)"#"Days per standardized forcing")+#expression("Forcing response (days/"*~degree*C*")")
        , main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 2, label = "c)", cex =10) +
  annotate("text", x = spTopE[1,5], y = -23, label = spTopE[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -23, label = spTopE[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -23, label = spTopE[3,2], cex = 6, angle = 78) +
  annotate("text", x =  spTopE[4,5], y = -23, label = spTopE[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -23, label = spTopE[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[6,5], y = -23, label = spTopE[6,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[7,5], y = -23, label = spTopE[7,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[8,5], y = -23, label = spTopE[8,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[9,5], y = -23, label = spTopE[9,2], cex = 6, angle = 78) +
  scale_fill_manual(values = c("#f9b641ff","#f9b641ff")) +
  geom_segment(aes(x = 65, y = 2, xend = 65 , yend = -3),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 65, y = -5, label = "Earlier", cex =5) 
forcedotE


photodotE <- ggplot(chillPtE,aes(y= photo, x = meanBB, colour = "Photoperiod"), size = 7) +
  geom_point(aes(colour = transect),size = 7, shape  = 17) +
  geom_errorbar(aes(ymin= photo5, ymax = photo95,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 0.5) +
  geom_errorbar(aes(ymin= photo25, ymax = photo75,xmin= meanBB, xmax = meanBB, colour = transect), width= 0, size = 1.5) +
  geom_point(aes(colour = transect),size =7, shape = 17) +
  geom_segment(aes(x = 65, y = -0.2, xend = 65 , yend = -2),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +  ylim(-10,0) +
  theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=  17),
        legend.position = "none") + scale_color_manual(values = c("east/west"="navy",
                                                                  "east"="cyan4"),
                                                       breaks = c("east/west","west"), label = c("High", "Low")) + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
  labs( x = "Species ordered by predicted budburst date", y = "Photoperiod response (days/standardized unit)"#"Days per standardized photoperiod"
        , main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 0, label = "e)", cex =10) +  
  #annotate("text", x = 23, y = -2.3, label = "Later", cex =10) +
  annotate("text", x = 65, y = -2.5, label = "Earlier", cex =5) +
  annotate("text", x = spTopE[1,5], y = -9.3, label = spTopE[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -9.3, label = spTopE[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -9.3, label = spTopE[3,2], cex = 6, angle = 78) +
  annotate("text", x =  spTopE[4,5], y = -9.3, label = spTopE[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -9.3, label = spTopE[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[6,5], y = -9.3, label = spTopE[6,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[7,5], y = -9.3, label = spTopE[7,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[8,5], y = -9.3, label = spTopE[8,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[9,5], y = -9.3, label = spTopE[9,2], cex = 6, angle = 78) +
  scale_fill_manual(values = c("#f9b641ff","#f9b641ff"))
photodotE

pdf("figures/dotCFPEWSizeColorUnitsShapeBoth.pdf", width = 20, height = 16)
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

pdf("figures/rankEstiBB.pdf", width = 9, height =3)
colTran <- c("maroon","navy","forestgreen")
par(mfrow = c(1,3), mar = c(5.1, 4.8, 4.1, 2.1))
plot(rank$rankHighC~rank$rankInt, 
     col = colTran[factor(rank$transect)], 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "High cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)

plot(rank$rankLowC~rank$rankInt, 
     col = colTran[factor(rank$transect)], 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "Low cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)

plot(rank$rankHighC~rank$rankLowC, 
     col = colTran[factor(rank$transect)], 
     pch = 19,
     xlab = "Low cue rank",
     ylab = "High cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)

legend("topleft",legend = c( expression("East"),
                             expression("West"),
                             expression("Both")),
       col = c("maroon","forestgreen","navy"),
       #pt.bg = c("#042333ff","#cc6a70ff","#593d9cff","#f9b641ff","#13306dff","#efe350ff","#eb8055ff"),
       pt.bg = c( "maroon","forestgreen","navy"),
       inset = 0.02, pch = c(21, 21, 21 ), cex = 1.25, bty = "n")
dev.off()

colType <- c("maroon","navy")
par(mfrow = c(1,3), mar = c(5.1, 4.8, 4.1, 2.1))
plot(rank$rankHighC~rank$rankInt, 
     col = colType[factor(rank$type)], 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "High cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)

plot(rank$rankLowC~rank$rankInt, 
     col = colType[factor(rank$type)], 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "Low cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)

plot(rank$rankHighC~rank$rankLowC, 
     col = colType[factor(rank$type)], 
     pch = 19,
     xlab = "Low cue rank",
     ylab = "High cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)

### Rank within transect:
rankE <- subset(rank, transect != "west")

rankE <- rankE[order(rankE$Int),]
rankE$rankInt <- seq(1:nrow(rankE))

rankE <- rankE[order(rankE$meanBB),]
rankE$rankLowC <- seq(1:nrow(rankE))

rankE <- rankE[order(rankE$meanBBHigh),]
rankE$rankHighC <- seq(1:nrow(rankE))

pdf("figures/rankstiBBTransect.pdf", width = 9, height =6)

par(mfrow = c(2,3), mar = c(5.1, 4.8, 4.1, 2.1))
plot(rankE$rankHighC~rankE$rankInt, 
     col = "purple3",  
     pch = 19,
     xlab = "Intercept rank",
     ylab = "High cue rank",
     main = "Eastern transect",
     cex.lab =1.5,
     cex =2, cex.main = 2)
abline(0,1)
text(1.5,27, label = "a)", cex = 1.5)


plot(rankE$rankLowC~rankE$rankInt, 
     col = "purple3",  
     pch = 19,
     xlab = "Intercept rank",
     ylab = "Low cue rank",
     cex.lab =1.5,
     cex =2)
abline(0,1)
text(1.5,27, label = "b)", cex = 1.5)

plot(rankE$rankHighC~rankE$rankLowC, 
     col = "purple3", 
     pch = 19,
     xlab = "Low cue rank",
     ylab = "High cue rank",
     cex.lab =1.5,
     cex =2)
abline(0,1)
text(1.5,27, label = "c)", cex = 1.5)


rankW <- subset(rank, transect != "east")

rankW <- rankW[order(rankW$Int),]
rankW$rankInt <- seq(1:nrow(rankW))

rankW <- rankW[order(rankW$meanBB),]
rankW$rankLowC <- seq(1:nrow(rankW))

rankW <- rankW[order(rankW$meanBBHigh),]
rankW$rankHighC <- seq(1:nrow(rankW))



#par(mfrow = c(2,3), mar = c(5.1, 4.8, 4.1, 2.1))
plot(rankW$rankHighC~rankW$rankInt, 
     col = "forestgreen", 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "High cue rank",main = "Western transect",
     cex.lab =1.5,
     cex =2, cex.main = 2,
     xlim = c(0,25),
     ylim = c(0,25)
)
abline(0,1)
text(0.8,24, label = "d)", cex = 1.5)


plot(rankW$rankLowC~rankW$rankInt, 
     col = "forestgreen", 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "Low cue rank",
     cex.lab =1.5,
     cex =2,
     xlim = c(0,25),
     ylim = c(0,25)
)
abline(0,1)
text(0.8,24, label = "e)", cex = 1.5)

plot(rankW$rankHighC~rankW$rankLowC, 
     col = "forestgreen", 
     pch = 19,
     xlab = "Low cue rank",
     ylab = "High cue rank",
     cex.lab =1.5,
     cex =2,
     xlim = c(0,25),
     ylim = c(0,25)
)
abline(0,1)
text(0.8,24, label = "f)", cex = 1.5)
dev.off()

### How does the ranking change for the species in both transects?
both <- c("Betula_papyrifera", "Alnus_incana","Populus_tremuloides")

rankBe <- rankE[rankE$species.name %in% both,]
rankBe$diffe <- abs(rankBe$rankHighC - rankBe$rankLowC)

rankBw <- rankW[rankW$species.name %in% both,]
rankBw$diffw <- abs(rankBw$rankHighC - rankBw$rankLowC)

rankDiff <- cbind(rankBe, rankBw$diffw)
rankDiff$species.name <- c("Betula papyrifera", "Alnus incana","Populus tremuloides")

pdf("figures/bothRanked.pdf", width = 5, height = 5)
ggplot(rankDiff) +
  geom_point(aes(x= species.name, y = diffe, col = "forestgreen"), size = 7) +
  geom_point (aes(x= species.name, y = rankBw$diffw, col = "purple4"), size =7) +
  scale_color_manual(values = c("forestgreen","purple4"),label = c("Eastern", "Western")) + 
  theme(axis.text.x = element_text( size=12,angle = 15,  hjust=1)) +
  theme(legend.title= element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key=element_rect(fill="white")) +  ylim(0,8) +
  labs( x = "Species name", y = "Absolute differene in rank between high and low cues"#"Days per standardized photoperiod"
        , main = NA) 
dev.off()
