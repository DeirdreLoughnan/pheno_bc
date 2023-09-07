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

# sum <- summary(mdl.4phyloContWP)$summary 
# post <- rstan::extract(mdl.4phyloContWP)

 # a_sp = mean(sum[grep("a_sp", rownames(sum)), 1])
# mu_b_warm = sum[grep("b_warm", rownames(sum)), 1]
# mu_b_photo = sum[grep("mu_b_photo", rownames(sum)), 1]
# mu_b_chill1 = sum[grep("mu_b_chill1", rownames(sum)), 1]
# mu_b_inter_pc1 = sum[grep("mu_b_inter_pc1", rownames(sum)), 1]
# mu_b_inter_wp = sum[grep("mu_b_inter_wp", rownames(sum)), 1]
# mu_b_inter_wc1 = sum[grep("mu_b_inter_wc1", rownames(sum)), 1]
# mu_b_inter_ws2 = sum[grep("mu_b_inter_ws2", rownames(sum)), 1]
# mu_b_inter_s2c1 = sum[grep("mu_b_inter_s2c1", rownames(sum)), 1]
# mu_b_inter_ps2 = sum[grep("mu_b_inter_ps2", rownames(sum)), 1]
# mu_b_inter_ws3 = sum[grep("mu_b_inter_ws3", rownames(sum)), 1]
# mu_b_inter_s3c1 = sum[grep("mu_b_inter_s3c1", rownames(sum)), 1]
# mu_b_inter_ps3 = sum[grep("mu_b_inter_ps3", rownames(sum)), 1]
# mu_b_inter_ws4 = sum[grep("mu_b_inter_ws4", rownames(sum)), 1]
# mu_b_inter_s4c1 = sum[grep("mu_b_inter_s4c1", rownames(sum)), 1]
# mu_b_inter_ps4 = sum[grep("mu_b_inter_ps4", rownames(sum)), 1]

# b_site2 = sum[grep("b_site2", rownames(sum)), 1]
# b_site3 = sum[grep("b_site3", rownames(sum)), 1]
# b_site4 = sum[grep("b_site4", rownames(sum)), 1]

a_sp = (sum[grep("a_sp", rownames(sum)), 1])
# a_sp97.5 = (sum[grep("a_sp", rownames(sum)), "97.5%"])
# a_sp2.5 = (sum[grep("a_sp", rownames(sum)), "2.5%"])

b_photo = sum[grep("b_photo\\[", rownames(sum)), 1]
b_chill = sum[grep("b_chill1\\[", rownames(sum)), 1]
b_force = sum[grep("b_warm\\[", rownames(sum)), 1]

# b_photo25 = sum[grep("b_photo\\[", rownames(sum)), "25%"]
# b_chill25 = sum[grep("b_chill1\\[", rownames(sum)), "25%"]
# b_force25 = sum[grep("b_warm\\[", rownames(sum)), "25%"]
# 
# b_photo75 = sum[grep("b_photo\\[", rownames(sum)), "75%"]
# b_chill75 = sum[grep("b_chill1\\[", rownames(sum)), "75%"]
# b_force75 = sum[grep("b_warm\\[", rownames(sum)), "75%"]
# 
# b_photo2.5 = sum[grep("b_photo\\[", rownames(sum)), "2.5%"]
# b_chill2.5 = sum[grep("b_chill1\\[", rownames(sum)), "2.5%"]
# b_force2.5 = sum[grep("b_warm\\[", rownames(sum)), "2.5%"]
# 
# b_photo97.5 = sum[grep("b_photo\\[", rownames(sum)), "97.5%"]
# b_chill97.5 = sum[grep("b_chill1\\[", rownames(sum)), "97.5%"]
# b_force97.5 = sum[grep("b_warm\\[", rownames(sum)), "97.5%"]

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


# par(mfrow = c(1,3))
# plot(spInfo$meanBB ~ spInfo$chill, xlab = "Chilling response", ylab = "Estimated budburst"); abline(lm(meanBB~chill, spInfo))
# plot(spInfo$meanBB ~ spInfo$force, xlab = "Forcing response", ylab = "Estimated budburst"); abline(lm(meanBB~force, spInfo))
# plot(spInfo$meanBB ~ spInfo$photo, xlab = "Photoperiod response", ylab = "Estimated budburst"); abline(lm(meanBB~photo, spInfo))

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# I think what we want is a loop that goes through each iteration of the posteriors and calculates the bb, but using 20 for forcing, 12 for photoperiod, 75 (75/10 when rescaled), and smithers to start
# 

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

# it <- 2
# sp <- 2
# m[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM + 
#   post$b_warm[it,sp] * force + post$b_photo[it, sp] * photo + post$b_chill[it,sp] * chill +
#   post$b_inter_wp[it,sp] * (force*photo) + post$b_inter_wc1[it,sp] * (force*chill) + post$b_inter_pc1[it,sp] * (photo*chill) +
#   post$b_inter_s2c1[it,sp] * (chill*siteSM) + post$b_inter_ws2[it,sp] * (force*siteSM) + post$b_inter_ps2[it,sp] * (photo*siteSM) +
#   post$b_inter_s3c1[it,sp] * (chill*siteSM) + post$b_inter_ws3[it,sp] * (force*siteSM) + post$b_inter_ps3[it,sp] * (photo*siteSM) +
#   post$b_inter_s4c1[it,sp] * (chill*siteSM) + post$b_inter_ws4[it,sp] * (force*siteSM) + post$b_inter_ps4[it,sp] * (photo*siteSM)
# 
# mHigh[it,sp] <- post$a_sp[it,sp]+ post$b_site2[it] * siteSM + post$b_site3[it] * siteSM + post$b_site4[it] * siteSM + 
#   post$b_warm[it,sp] * forceHigh + post$b_photo[it, sp] * photoHigh + post$b_chill[it,sp] * chillHigh +
#   post$b_inter_wp[it,sp] * (forceHigh*photoHigh) + post$b_inter_wc1[it,sp] * (forceHigh*chillHigh) + post$b_inter_pc1[it,sp] * (photoHigh*chillHigh) +
#   post$b_inter_s2c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws2[it,sp] * (forceHigh*siteSM) + post$b_inter_ps2[it,sp] * (photoHigh*siteSM) +
#   post$b_inter_s3c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws3[it,sp] * (forceHigh*siteSM) + post$b_inter_ps3[it,sp] * (photoHigh*siteSM) +
#   post$b_inter_s4c1[it,sp] * (chillHigh*siteSM) + post$b_inter_ws4[it,sp] * (forceHigh*siteSM) + post$b_inter_ps4[it,sp] * (photoHigh*siteSM)
# now get the order of diff spp bb that I can use to order the figure
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

# spInfo$force2.5 <- b_force2.5
# spInfo$chill2.5 <- b_chill2.5
# spInfo$photo2.5 <- b_photo2.5

# spInfo$force97.5 <- b_force97.5
# spInfo$chill97.5 <- b_chill97.5
# spInfo$photo97.5 <- b_photo97.5

# spInfo$force25 <- b_force25
# spInfo$chill25 <- b_chill25
# spInfo$photo25 <- b_photo25
# 
# spInfo$force75 <- b_force75
# spInfo$chill75 <- b_chill75
# spInfo$photo75 <- b_photo75


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

names(data) <- c("species.name","valueLow", "valueHigh","species","type","transect","meanBB","meanBBHigh", "Int","Int5","Int95","Int25","Int75","force5","force95","force25","force75","chill5", "chill95", "chill25", "chill75","photo5", "photo95", "photo25", "photo75","spMeanForce", "spMeanChill", "spMeanPhoto","bb5","bb95","bb75","bb25", "spacing","bb5High","bb95High","bb75High","bb25High","chill", "force","photo","intercept")

#####################################################################
## Eastern spp only #################################################

east <- subset(spInfo, transect != "west")
eastSp <- unique(east$species.name)

dataEast <- data[data$species.name %in% eastSp, ]

overlappingE <- c("aromel","betlen", "betpap", "lyolig","faggra", "betall", "prupen","poptre","rhafra")
spMiniE <- east[!east$species %in% overlappingE,]
spTopE <- east[east$species %in% overlappingE,]

# bbSpaceE <- ggplot() + 
#   geom_violin(dataEast, mapping=aes(x = meanBB, y = value, group = species.name, width =1, fill = type)) +
#   geom_point(east, mapping= aes(x = meanBB, y = meanBB)) +
#   geom_linerange(data = east,aes(x = meanBB, y = meanBB, ymin=bb2.5, ymax = bb97.5) ) +
#   theme_classic() +  
#  ylim(-3,80) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
#   labs( x = "Species", y = "Estimated budburst", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 80, label = "a)", cex =5) +
#   annotate("text", x =40, y = 80, label = "Eastern transect", cex =5) +
#   annotate("text", x = spTopE[1,5], y = -0.2, label = spTopE[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[2,5], y = -0.2, label = spTopE[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[3,5], y = -0.2, label = spTopE[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[4,5], y = -0.2, label = spTopE[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[5,5], y = -0.2, label = spTopE[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[6,5], y = -0.2, label = spTopE[6,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[7,5], y = -0.2, label = spTopE[7,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[8,5], y = -0.2, label = spTopE[8,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("mediumpurple2","mediumpurple2"))
# 
# # intercept
# aSpaceE<-  ggplot() + 
#   geom_violin(dataEast, mapping=aes(x = meanBB, y = intercept, group = species.name, width =1, fill = type)) +
#   geom_point(east, mapping= aes(x = meanBB, y = Int)) +
#   geom_linerange(data = east,aes(x = meanBB, y = Int, ymin=Int5, ymax = Int95)) +
#   theme_classic() +  
#   ylim(-6,65) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
#   labs( x = "Species", y = "Species intercept", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 60, label = "c)", cex =5) +
#   annotate("text", x = spTopE[1,5], y = -4, label = spTopE[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[2,5], y = -4, label = spTopE[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[3,5], y = -4, label = spTopE[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[4,5], y = -4, label = spTopE[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[5,5], y = -4, label = spTopE[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[6,5], y = -4, label = spTopE[6,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[7,5], y = -4, label = spTopE[7,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[8,5], y = -4, label = spTopE[8,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("darkolivegreen4","darkolivegreen4"))
# 
# chillSpaceE <- ggplot() + 
#   geom_violin(dataEast, mapping=aes(x = meanBB, y = chill, group = species.name, width =1, fill = type)) +
#   geom_point(east, mapping= aes(x = meanBB, y = chill)) +
#   geom_linerange(data = east, aes(x = meanBB, y = chill, ymin=chill2.5, ymax = chill97.5)) +
#   theme_classic() +  
#   ylim(-45,15) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
#   labs( x = "Species", y = "Chilling response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x =23, y = 15, label = "a)", cex =5) +
#   annotate("text", x = 45, y = 15, label = "Eastern transect", cex =5) +
#   annotate("text", x = spTopE[1,5], y = -44.5, label = spTopE[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[2,5], y = -44.5, label = spTopE[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[3,5], y = -44.5, label = spTopE[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[4,5], y = -44.5, label = spTopE[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[5,5], y = -44.5, label = spTopE[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[6,5], y = -44.5, label = spTopE[6,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("#cc6a70ff","#cc6a70ff"))
# 
# forceSpaceE <- ggplot() + 
#   geom_violin(dataEast, mapping=aes(x = meanBB, y = force, group = species.name, width =1, fill = type)) +
#   geom_point(east, mapping= aes(x = meanBB, y = force)) +
#   geom_linerange(data = east,aes(x = meanBB, y = force, ymin=force2.5, ymax = force97.5)) +
#   theme_classic() +  
#   ylim(-25,10) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
#   labs( x = "Species", y = "Forcing response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 10, label = "c)", cex =5) +
#   annotate("text", x = spTopE[1,5], y = -24.5, label = spTopE[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[2,5], y = -24.5, label = spTopE[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[3,5], y = -24.5, label = spTopE[3,2], cex = 3, angle = 78) +
#   annotate("text", x =  spTopE[4,5], y = -24.5, label = spTopE[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[5,5], y = -24.5, label = spTopE[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[6,5], y = -24.5, label = spTopE[6,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[7,5], y = -24.5, label = spTopE[7,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[8,5], y = -24.5, label = spTopE[8,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("#f9b641ff","#f9b641ff"))
# 
# 
# photoSpaceE <- ggplot() + 
#   geom_violin(dataEast, mapping=aes(x = meanBB, y = photo, group = species.name, width =1, fill = type)) +
#   geom_point(east, mapping= aes(x = meanBB, y = photo)) +
#   geom_linerange(data = east,aes(x = meanBB, y = photo, ymin=photo2.5, ymax = photo97.5)) +
#   theme_classic() +  
#   ylim(-25,10) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 10, label = "e)", cex =5) +
#   annotate("text", x = spTopE[1,5], y = -24.5, label = spTopE[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[2,5], y = -24.5, label = spTopE[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[3,5], y = -24.5, label = spTopE[3,2], cex = 3, angle = 78) +
#   annotate("text", x =  spTopE[4,5], y = -24.5, label = spTopE[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopE[5,5], y = -24.5, label = spTopE[5,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("cyan4","cyan4"))
#  
# pdf("figures/4violinEastern2575PT.pdf", width = 10, height =24)
# plot_grid(bbSpaceE, aSpaceE,chillSpaceE,forceSpaceE, photoSpaceE, nrow = 5, align = "v")
# dev.off()


#####################################################################
## Western spp only #################################################

west <- subset(spInfo, transect != "east")
westSp <- unique(west$species.name)

dataWest <- data[data$species.name %in% westSp, ]


overlappingW <- c("spialb","betpap","popbal","rhoalb","alninc")
spMiniW <- west[!west$species %in% overlappingW,]
spTopW <- west[west$species %in% overlappingW,]

# bbSpaceW <-  ggplot() + 
#   geom_violin(dataWest, mapping=aes(x = meanBB, y = value, group = species.name, width =1, fill = type)) +
#   geom_point(west, mapping= aes(x = meanBB, y = meanBB)) +
#   geom_linerange(data = west,aes(x = meanBB, y = meanBB, ymin=bb2.5, ymax = bb97.5)) +
#   theme_classic() +  
#   ylim(-5,70) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
#   labs( x = "Species", y = "Estimated budburst", main = NA) +
#   annotate("text", x = 40, y = 70, label = "Western transect", cex =5) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 70, label = "b)", cex =5) +
#   annotate("text", x = spTopW[1,5], y = -1.5, label = spTopW[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[2,5], y =-1.5, label = spTopW[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[3,5], y = -1.5, label = spTopW[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[4,5], y = -1.5, label = spTopW[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[5,5], y = -1.5, label = spTopW[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[6,5], y = -1.5, label = spTopW[6,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("mediumpurple2","mediumpurple2"))
# 
# # intercept
# aSpaceW <-  ggplot() + 
#   geom_violin(dataWest, mapping=aes(x = meanBB, y = intercept, group = species.name, width =1, fill = type)) +
#   geom_point(west, mapping= aes(x = meanBB, y = Int)) +
#   geom_linerange(data = west,aes(x = meanBB, y = Int, ymin=Int2.5, ymax = Int97.5)) +
#   theme_classic() + 
#   ylim(-5,60)+
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
#   labs( x = "Species", y = "Species intercept", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 60, label = "d)", cex =5) +
#   annotate("text", x = spTopW[1,5], y = -1.5, label = spTopW[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[2,5], y = -1.5, label = spTopW[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[3,5], y = -1.5, label = spTopW[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[4,5], y = -1.5, label = spTopW[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[5,5], y = -1.5, label = spTopW[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[6,5], y = -1.5, label = spTopW[6,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("darkolivegreen4","darkolivegreen4"))
# 
# 
# chillSpaceW <-  ggplot() + 
#   geom_violin(dataWest, mapping=aes(x = meanBB, y = chill, group = species.name, width =1, fill = type)) +
#   geom_point(west, mapping= aes(x = meanBB, y = chill)) +
#   geom_linerange(data = west, aes(x = meanBB, y = chill, ymin=chill2.5, ymax = chill97.5)) +
#   theme_classic() +  
#   ylim(-45,10) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
#   labs( x = "Species", y = "Chilling response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 10, label = "b)", cex =5) +
#   annotate("text", x = 38, y = 10, label = "Western transect", cex =5) +
#   annotate("text", x = spTopW[1,5], y = -44.5, label = spTopW[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[2,5], y = -44.5, label = spTopW[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[3,5], y = -44.5, label = spTopW[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[4,5], y = -44.5, label = spTopW[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[5,5], y = -44.5, label = spTopW[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[6,5], y = -44.5, label = spTopW[6,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("#cc6a70ff","#cc6a70ff"))
# 
# forceSpaceW <- ggplot() + 
#   geom_violin(dataWest, mapping=aes(x = meanBB, y = force, group = species.name, width =1, fill = type)) +
#   geom_point(west, mapping= aes(x = meanBB, y = force)) +
#   geom_linerange(data = west,aes(x = meanBB, y = force, ymin=force2.5, ymax = force97.5)) +
#   theme_classic() +  
#   ylim(-25,10) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
#   labs( x = "Species", y = "Forcing response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 10, label = "d)", cex =5) +
#   annotate("text", x = spTopW[1,5], y = -24.5, label = spTopW[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[2,5], y = -24.5, label = spTopW[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[3,5], y = -24.5, label = spTopW[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[4,5], y = -24.5, label = spTopW[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[5,5], y = -24.5, label = spTopW[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[6,5], y = -24.5, label = spTopW[6,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("#f9b641ff","#f9b641ff"))
# 
# photoSpaceW <- ggplot() + 
#   geom_violin(dataWest, mapping=aes(x = meanBB, y = photo, group = species.name, width =1, fill = type)) +
#   geom_point(west, mapping= aes(x = meanBB, y = photo)) +
#   geom_linerange(data = west,aes(x = meanBB, y = photo, ymin=photo2.5, ymax = photo97.5)) +
#   theme_classic() +  
#   ylim(-25,10) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15),
#         legend.position = "none") + 
#   scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,62)) +
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 10, label = "f)", cex =5) +
#   annotate("text", x = spTopW[1,5], y = -24.5, label = spTopW[1,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[2,5], y = -24.5, label = spTopW[2,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[3,5], y = -24.5, label = spTopW[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[4,5], y = -24.5, label = spTopW[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[5,5], y = -24.5, label = spTopW[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTopW[6,5], y = -24.5, label = spTopW[6,2], cex = 3, angle = 78) +
#   scale_fill_manual(values = c("cyan4","cyan4"))

# pdf("figures/4violinWestern9725.pdf", width = 10, height =24)
# plot_grid(bbSpaceW, aSpaceW,chillSpaceW,forceSpaceW, photoSpaceW, nrow = 5, align = "v")
# dev.off()

#Making better figures for the MS:

# pdf("figures/violinBetaAlphaEW.pdf", width = 16, height = 8)
# plot_grid(bbSpaceE, bbSpaceW, aSpaceE, aSpaceW, nrow = 2, ncol = 2, align = "v")
# dev.off()
# 
# pdf("figures/violinCFPEW.pdf", width = 16, height = 16)
# plot_grid(chillSpaceE, chillSpaceW, forceSpaceE, forceSpaceW, photoSpaceE, photoSpaceW, nrow = 3, ncol = 2, align = "v")
# dev.off()

# Making lizzie's other suggested plot - 2 connected dots:

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
  labs( x = "Species ordered by predicted budburst date", y = "Estimated parameter (days/standardized units)", main = NA) +
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
  labs( x = "Species ordered by predicted budburst date under high cues", y = "Day of budburst (days/standardized units)", main = NA) +
  theme(legend.position = "none") +  annotate("text", x = 21.5, y = 60, label = "b)      Western transect", cex =8) +
  annotate("text", x = spTopW[1,5], y = -1, label = spTopW[1,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -1, label = spTopW[2,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -1, label = spTopW[3,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[4,5], y = -1, label = spTopW[4,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[5,5], y = -1, label = spTopW[5,2], cex = 5, angle = 78) +
  scale_shape_manual(values = c(1,2,16)) +
  annotate("text", x = 66, y = 60.5, label = "High cues", cex = 7) +  
  annotate("text", x = 66, y = 44.5, label = "Intercept", cex = 7) +
  annotate("text", x = 66, y = 27.5, label = "Low cues", cex = 7) + 
  geom_segment(aes(x = 61, y = 60.5, xend = 62.5 , yend = 60.5)) +
  geom_segment(aes(x = 61, y = 44.5, xend = 62.5 , yend = 44.5)) +
  geom_segment(aes(x = 61, y = 27.5, xend = 62.5 , yend = 27.5))

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
  labs( x = "Species ordered by predicted budburst date under high cues", y = "Day of budburst (days/standardized units)", main = NA) +
  theme(legend.position = "none") +  annotate("text", x = 18, y = 60, label = "b)      Western transect", cex =8) +
  annotate("text", x = spTopW[1,7], y = -1, label = spTopW[1,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[2,7], y = -1, label = spTopW[2,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[3,7], y = -1, label = spTopW[3,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[4,7], y = -1, label = spTopW[4,2], cex = 5, angle = 78) +
  annotate("text", x = spTopW[5,7], y = -1, label = spTopW[5,2], cex = 5, angle = 78) +
  scale_shape_manual(values = c(1,2,16)) 

+
  annotate("text", x = 66, y = 60.5, label = "High cues", cex = 7) +  
  annotate("text", x = 66, y = 44.5, label = "Intercept", cex = 7) +
  annotate("text", x = 66, y = 27.5, label = "Low cues", cex = 7) + 
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
  labs( x = "Species ordered by predicted budburst date", y = "Estimated parameter (days/standardized units)", main = NA) +
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
  theme(legend.position = "none") +  annotate("text", x = 29, y = 60, label = "a)      Eastern transect", cex =8) +
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
  annotate("text", x = 68, y = 63, label = "High cues", cex = 7) +  
  annotate("text", x = 68, y = 48, label = "Intercept", cex = 7) +
  annotate("text", x = 68, y = 36.5, label = "Low cues", cex = 7) + 
  geom_segment(aes(x = 64, y = 63, xend = 65 , yend = 63)) +
  geom_segment(aes(x = 64, y = 48, xend = 65.5 , yend = 48)) +
  geom_segment(aes(x = 64, y = 36.5, xend = 65 , yend = 36.5))
dotEBw


# pdf("figures/dotBetaAlphaLongColor.pdf", width = 12, height = 16)
# plot_grid(dotE, dotW, nrow = 2, ncol = 1, align = "v")
# dev.off()

pdf("figures/dotBetaAlphaLongBW_LH.pdf", width = 12, height = 16)
plot_grid(dotEBw, dotWBw, nrow = 2, ncol = 1, align = "v")
dev.off()
#############################################
# old plots
####### Old plot not spaced out ####################3

# bbSp <- ggplot() + 
#   stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = value), .width = c(.90, .5), cex = 0.75, fill = "mediumpurple2") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Estimated budburst", main = NA)+ 
#   scale_color_identity(name = "Model post",
#                        breaks = c("black"),
#                        labels = c("Model Posterior"),
#                        guide = guide_legend(override.aes = list(
#                          linetype = c(NA),
#                          shape = c(8)))) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 80, label = "a)", cex =5) 
# 
# 
# chillSp <- ggplot() + 
#   stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = chill), .width = c(.90, .5), cex = 0.75, fill = "#cc6a70ff") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Chilling Response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 
# # 
# # pdf("figures/testChill2.pdf", width = 10, height =5)
# # chillSp
# # dev.off()
# 
# # Forcing 
# 
# forceSp <- ggplot() + 
#   stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = force), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Forcing response", main = NA)+
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "c)", cex =5) 
# 
# # pdf("figures/testforce.pdf", width = 10, height =5)
# # forceSp
# # dev.off()
# 
# # Photoperiod
# 
# photoSp <- ggplot() + 
#   stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = photo), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 
# 
# pdf("figures/4panel.pdf", width = 10, height =20)
# plot_grid(bbSp, chillSp,forceSp, photoSp, nrow = 4, align = "v", rel_heights = c(1/4, 1/4, 1/4,1.2/3))
# dev.off()

################################################################################
#### Now can we somehow make this work with different spacing? 

# bbSpace <- ggplot() + 
#   stat_eye(data = data, aes(x = meanBB, y = value), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   #geom_errorbar(aes(xmin= bb_df$bb25, xmax = bb_df$bb75), width= 0) +
#  # xlim (15,80) +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Estimated budburst", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 80, label = "a)", cex =5) + scale_x_continuous( breaks = data$meanBB, labels = data$species.name,limits = c(15,80) )


# bbSpace <- ggplot() + 
#   geom_point(data = spInfo, aes(x = meanBB, y = value)) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   geom_pointrange(data = spInfo,aes(x = meanBB, y = value, ymin=bb25, ymax = bb75)) +
#   # xlim (15,80) +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size=10,
#                                 angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Estimated budburst", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 80, label = "a)", cex =5) + 
#   scale_x_continuous( breaks = spInfo$meanBB, labels = spInfo$species,limits = c(17.2,68)) +
#   annotate("text", x = 1, y = 80, label = "a)", cex =5)

overlapping <- c("aromel","prupen","vacmyr","spialb","ilemuc", "vibcas","betpap","betlen","symalb","acerub","lyolig","rhopri")
spMini <- spInfo[!spInfo$species %in% overlapping,]
spTop <- spInfo[spInfo$species %in% overlapping,]

# bbSpace <- ggplot() + 
#   geom_point(data = spInfo, aes(x = meanBB, y = value)) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   geom_pointrange(data = spInfo,aes(x = meanBB, y = value, ymin=bb5, ymax = bb90)) +
#   # xlim (15,80) +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size= 8.9,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Estimated budburst", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 17.5, y = 80, label = "a)", cex =5) + 
#   scale_x_continuous( breaks = spMini$meanBB, labels = spMini$species,limits = c(17.2,68)) +
#   annotate("text", x = spTop[1,5], y = 15, label = spTop[1,2], cex = 3, angle = 78) +
#   annotate("text", x = 24.5, y = 15, label = spTop[2,2], cex = 3, angle = 78) +
#   annotate("text", x = 40.3, y = 15, label = spTop[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[4,5], y = 15, label = spTop[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[5,5], y = 15, label = spTop[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[6,5], y = 15, label = spTop[6,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[7,5], y = 15, label = spTop[7,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[8,5], y = 15, label = spTop[8,2], cex = 3, angle = 78) +
#   annotate("text", x = 29.68, y = 15, label = spTop[9,2], cex = 3, angle = 78) +
#   annotate("text", x = 43.1, y = 15, label = spTop[10,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[11,5], y = 15, label = spTop[11,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[12,5], y = 15, label = spTop[12,2], cex = 3, angle = 78) 

# chillSpace <- ggplot() + 
#   geom_point(data = spInfo, aes(x = meanBB, y = chill)) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   geom_pointrange(data = spInfo,aes(x = meanBB, y = chill, ymin=chill5, ymax = chill95)) +
#   # xlim (15,80) +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size= 8.9,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Chilling response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 17.5, y = 0.1, label = "b)", cex =5) + 
#   scale_x_continuous( breaks = spMini$meanBB, labels = spMini$species,limits = c(17.2,68)) +
#   annotate("text", x = spTop[1,5], y = -36.5, label = spTop[1,2], cex = 3, angle = 78) +
#   annotate("text", x = 24.5, y = -36.5, label = spTop[2,2], cex = 3, angle = 78) +
#   annotate("text", x = 40.3, y = -36.5, label = spTop[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[4,5], y = -36.5, label = spTop[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[5,5], y = -36.5, label = spTop[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[6,5], y = -36.5, label = spTop[6,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[7,5], y = -36.5, label = spTop[7,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[8,5], y = -36.5, label = spTop[8,2], cex = 3, angle = 78) +
#   annotate("text", x = 29.68, y = -36.5, label = spTop[9,2], cex = 3, angle = 78) +
#   annotate("text", x = 43.1, y = -36.5, label = spTop[10,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[11,5], y = -36.5, label = spTop[11,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[12,5], y = -36.5, label = spTop[12,2], cex = 3, angle = 78) 
# 
# # stat_eye(data = data, aes(x = meanBB, y = chill), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
# # #geom_hline(yintercept = sum[4,1], linetype="dashed") +
# # theme_classic() +  
# # xlim (15,80) +
# # theme(axis.text.x = element_text( size=10,
# #                                   angle = 78, 
# #                                   hjust=1),
# #       axis.title.y=element_text(size = 12),
# #       axis.title=element_text(size=15) ) + # angle of 55 also works
# # #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
# # labs( x = "Species", y = "Estimated budburst", main = NA) +
# # theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 
# 
# #
# forceSpace <- ggplot() + 
#   geom_point(data = spInfo, aes(x = meanBB, y = force)) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   geom_pointrange(data = spInfo,aes(x = meanBB, y = force, ymin=force2.5, ymax = force97.5)) +
#   # xlim (15,80) +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size= 8.9,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Forcing response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 17.5, y = 1.75, label = "c)", cex =5) + 
#   scale_x_continuous( breaks = spMini$meanBB, labels = spMini$species,limits = c(17.2,68)) +
#   annotate("text", x = spTop[1,5], y = -20.2, label = spTop[1,2], cex = 3, angle = 78) +
#   annotate("text", x = 24.5, y = -20.2, label = spTop[2,2], cex = 3, angle = 78) +
#   annotate("text", x = 40.3, y = -20.2, label = spTop[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[4,5], y = -20.2, label = spTop[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[5,5], y = -20.2, label = spTop[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[6,5], y = -20.2, label = spTop[6,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[7,5], y = -20.2, label = spTop[7,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[8,5], y = -20.2, label = spTop[8,2], cex = 3, angle = 78) +
#   annotate("text", x = 29.68, y = -20.2, label = spTop[9,2], cex = 3, angle = 78) +
#   annotate("text", x = 43.1, y = -20.2, label = spTop[10,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[11,5], y = -20.2, label = spTop[11,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[12,5], y = -20.2, label = spTop[12,2], cex = 3, angle = 78) 
# 
# #   stat_eye(data = data, aes(x = meanBB, y = force), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
# #   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
# #   theme_classic() +  
# #   xlim (15,80) +
# #   theme(axis.text.x = element_blank(),
# #         axis.title.y=element_text(size = 12),
# #         axis.title=element_text(size=15) ) + # angle of 55 also works
# #   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
# #   labs( x = "Species", y = "Forcing response", main = NA)+
# #   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "c)", cex =5) 
# 
# # pdf("figures/testforce.pdf", width = 10, height =5)
# # forceSp
# # dev.off()
# 
# # Photoperiod
# 
# photoSpace <- ggplot() + 
#   geom_point(data = spInfo, aes(x = meanBB, y = photo)) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   geom_pointrange(data = spInfo,aes(x = meanBB, y = photo, ymin=photo2.5, ymax = photo97.5)) +
#   # xlim (15,80) +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size= 8.9,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 17.5, y = 1.75, label = "d)", cex =5) + 
#   scale_x_continuous( breaks = spMini$meanBB, labels = spMini$species,limits = c(17.2,68)) +
#   annotate("text", x = spTop[1,5], y = -10, label = spTop[1,2], cex = 3, angle = 78) +
#   annotate("text", x = 24.5, y = -10, label = spTop[2,2], cex = 3, angle = 78) +
#   annotate("text", x = 40.3, y = -10, label = spTop[3,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[4,5], y = -10, label = spTop[4,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[5,5], y = -10, label = spTop[5,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[6,5], y = -10, label = spTop[6,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[7,5], y = -10, label = spTop[7,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[8,5], y = -10, label = spTop[8,2], cex = 3, angle = 78) +
#   annotate("text", x = 29.68, y = -10, label = spTop[9,2], cex = 3, angle = 78) +
#   annotate("text", x = 43.1, y = -10, label = spTop[10,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[11,5], y = -10, label = spTop[11,2], cex = 3, angle = 78) +
#   annotate("text", x = spTop[12,5], y = -10, label = spTop[12,2], cex = 3, angle = 78) 
# #   stat_eye(data = data, aes(x = meanBB, y = photo), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
# #   theme_classic() +  xlim (15,80) +
# #   theme(axis.text.x = element_text( size=10,
# #                                     angle = 78, 
# #                                     hjust=1),
# #         axis.title.y=element_text(size = 12),
# #         axis.title=element_text(size=15) ) + # angle of 55 also works
# #   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
# #   labs( x = "Species", y = "Photoperiod response", main = NA) +
# #   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 
# 
# pdf("figures/4panelSpace.pdf", width = 10, height =20)
# plot_grid(bbSpace, chillSpace,forceSpace, photoSpace, nrow = 4, align = "v")
# dev.off()


# bbSpE <- ggplot() + 
#   stat_eye(data = dataEast, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = value), .width = c(.90, .5), cex = 0.75, fill = "mediumpurple2") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Estimated budburst", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 80, label = "a)", cex =5) 
# 
# 
# chillSpE <- ggplot() + 
#   stat_eye(data = dataEast, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = chill), .width = c(.90, .5), cex = 0.75, fill = "#cc6a70ff") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Chilling Response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 
# 
# #forcing
# forceSpE <- ggplot() + 
#   stat_eye(data = dataEast, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = force), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Forcing response", main = NA)+
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "c)", cex =5) 
# 
# # Photoperiod
# photoSpE <- ggplot() + 
#   stat_eye(data = dataEast, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = photo), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 
# 
# pdf("figures/4panelEastern.pdf", width = 10, height =20)
# plot_grid(bbSpE, chillSpE,forceSpE, photoSpE, nrow = 4, align = "v", rel_heights = c(1/4, 1/4, 1/4,1.2/3))
# dev.off()
# 
# bbSpW <- ggplot() + 
#   stat_eye(data = dataWest, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = value), .width = c(.90, .5), cex = 0.75, fill = "mediumpurple2") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Estimated budburst", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 80, label = "a)", cex =5) 
# 
# 
# chillSpW <- ggplot() + 
#   stat_eye(data = dataWest, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = chill), .width = c(.90, .5), cex = 0.75, fill = "#cc6a70ff") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Chilling Response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 
# 
# #forcing
# forceSpW <- ggplot() + 
#   stat_eye(data = dataWest, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = force), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Forcing response", main = NA)+
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "c)", cex =5) 
# 
# # Photoperiod
# photoSpW <- ggplot() + 
#   stat_eye(data = dataWest, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = photo), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 
# 
# pdf("figures/4panelWestern.pdf", width = 10, height =20)
# plot_grid(bbSpW, chillSpW,forceSpW, photoSpW, nrow = 4, align = "v", rel_heights = c(1/4, 1/4, 1/4,1.2/3))
# dev.off()

#########################################################
###### Remake fig but with dots and bars ################

chillPtW <- aggregate(dataWest[c("meanBB", "Int","chill","force","photo","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")], dataWest[c("species.name","type","transect")], FUN = mean)
# names(chillPtW) <- c("species.name","type","transect","Budburst","Intercept","Int5", "Int95", "Int25", "Int75","chill5","chill95","chill25","chill75","force5","force95","force25","force75","photo5","photo95","photo25","photo75")


chilldotW <- ggplot(chillPtW,aes(y= chill, x = meanBB, colour = "#cc6a70ff"), size = 7) + geom_point(size =7, colour ="#cc6a70ff", shape = 15) +
  geom_errorbar(aes(ymin= chill5, ymax = chill95,xmin= meanBB, xmax = meanBB), width= 0, size = 0.5, colour = "#cc6a70ff")+
  geom_errorbar(aes(ymin= chill25, ymax = chill75,xmin= meanBB, xmax = meanBB), width= 0, size = 1.5, colour = "#cc6a70ff") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-40,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size= 17),
        legend.position = "none") + 
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
  scale_fill_manual(values = c("maroon","maroon")) + 
  geom_segment(aes(x = 62, y = 0, xend = 62 , yend = -7),
  arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 62, y = -9, label = "Earlier", cex =5) 
chilldotW

forcedotW <- ggplot(chillPtW,aes(y= force, x = meanBB), size = 7) +
  geom_point(size = 7,color = "#f9b641ff", shape = 15) +
  geom_errorbar(aes(ymin= force5, ymax = force95,xmin= meanBB, xmax = meanBB),color = "#f9b641ff", width= 0, size =0.5) +
  geom_errorbar(aes(ymin= force25, ymax = force75,xmin= meanBB, xmax = meanBB),color = "#f9b641ff", width= 0, size =1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-25,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none") + 
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
  geom_point(size = 7, color = "cyan4", shape =15) +
  geom_errorbar(aes(ymin= photo5, ymax = photo95,xmin= meanBB, xmax = meanBB), width= 0, size = 0.5, color = "cyan4") +
  geom_errorbar(aes(ymin= photo25, ymax = photo75,xmin= meanBB, xmax = meanBB),color = "cyan4", width= 0, size =1.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-10,0) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,65)) +
  labs( x = "Species ordered by predicted budburst date", y = "Photoperiod response (days/standardized unit)",#"Days per standardized photoperiod",
        main = NA) +
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
  geom_point(size = 7, colour = "#cc6a70ff", shape  = 17) +
  geom_errorbar(aes(ymin= chill5, ymax = chill95,xmin= meanBB, xmax = meanBB), width= 0, size = 0.5, colour = "#cc6a70ff") + geom_errorbar(aes(ymin= chill25, ymax = chill75,xmin= meanBB, xmax = meanBB), width= 0, size = 1.5, colour = "#cc6a70ff") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + ylim(-45,5) +
  theme(axis.text.x = element_text( size = 17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none") + 
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
  geom_point(size = 7,color = "#f9b641ff", shape  = 17) +
  geom_errorbar(aes(ymin= force5, ymax = force95,xmin= meanBB, xmax = meanBB),color = "#f9b641ff", width= 0, size = 0.5) +
  geom_errorbar(aes(ymin= force25, ymax = force75,xmin= meanBB, xmax = meanBB),color = "#f9b641ff", width= 0, size = 1.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +  ylim(-25,5) +
  theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=17),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,65)) +
  labs( x = "", y = "Forcing response (days/standardized unit)"#"Days per standardized forcing")+#expression("Forcing response (days/"*~degree*C*")")
  , main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 2, label = "c)", cex =10) +
  annotate("text", x = spTopE[1,5], y = -23.5, label = spTopE[1,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -23.5, label = spTopE[2,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -23.5, label = spTopE[3,2], cex = 6, angle = 78) +
  annotate("text", x =  spTopE[4,5], y = -23.5, label = spTopE[4,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -23.5, label = spTopE[5,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[6,5], y = -23.5, label = spTopE[6,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[7,5], y = -23.5, label = spTopE[7,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[8,5], y = -23.5, label = spTopE[8,2], cex = 6, angle = 78) +
  annotate("text", x = spTopE[9,5], y = -23.5, label = spTopE[9,2], cex = 6, angle = 78) +
  scale_fill_manual(values = c("#f9b641ff","#f9b641ff")) +
  geom_segment(aes(x = 65, y = 2, xend = 65 , yend = -3),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
  annotate("text", x = 65, y = -5, label = "Earlier", cex =5) 
forcedotE


photodotE <- ggplot(chillPtE,aes(y= photo, x = meanBB, colour = "Photoperiod"), size = 7) +
  geom_point(size = 7, color = "cyan4", shape  = 17) +
  geom_errorbar(aes(ymin= photo5, ymax = photo95,xmin= meanBB, xmax = meanBB), width= 0, size = 0.5, color = "cyan4") +
  geom_errorbar(aes(ymin= photo25, ymax = photo75,xmin= meanBB, xmax = meanBB), width= 0, size = 1.5, color = "cyan4") +
  geom_segment(aes(x = 65, y = -0.2, xend = 65 , yend = -2),
               arrow = arrow(length = unit(0.5, "cm")), col = "black") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +  ylim(-10,0) +
  theme(axis.text.x = element_text( size=17,angle = 78,  hjust=1),
        axis.text.y=element_text(size = 15),
        axis.title=element_text(size=  17),
        legend.position = "none") + 
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

pdf("figures/dotCFPEWSizeColorUnitsShape.pdf", width = 20, height = 16)
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
     col = "maroon", 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "High cue rank",
     main = "Eastern transect",
     cex.lab =1.5,
     cex =1.5, cex.main = 2)
abline(0,1)

plot(rankE$rankLowC~rankE$rankInt, 
     col = "maroon", 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "Low cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)

plot(rankE$rankHighC~rankE$rankLowC, 
     col = "maroon", 
     pch = 19,
     xlab = "Low cue rank",
     ylab = "High cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)


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
     cex =1.5, cex.main = 2
     )
abline(0,1)

plot(rankW$rankLowC~rankW$rankInt, 
     col = "forestgreen", 
     pch = 19,
     xlab = "Intercept rank",
     ylab = "Low cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)

plot(rankW$rankHighC~rankW$rankLowC, 
     col = "forestgreen", 
     pch = 19,
     xlab = "Low cue rank",
     ylab = "High cue rank",
     cex.lab =1.5,
     cex =1.5)
abline(0,1)
dev.off()
