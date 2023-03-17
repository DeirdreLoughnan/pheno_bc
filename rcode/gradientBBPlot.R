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

load("output/bb_4sites_phylo_mini.Rda")
sum <- summary(mdl.4phyloMini)$summary 

post <- rstan::extract(mdl.4phyloMini)

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
a_sp97.5 = (sum[grep("a_sp", rownames(sum)), "97.5%"])
a_sp2.5 = (sum[grep("a_sp", rownames(sum)), "2.5%"])

b_photo = sum[grep("b_photo", rownames(sum)), 1]
b_chill = sum[grep("b_chill", rownames(sum)), 1]
b_force = sum[grep("b_warm", rownames(sum)), 1]

b_photo25 = sum[grep("b_photo", rownames(sum)), "25%"]
b_chill25 = sum[grep("b_chill", rownames(sum)), "25%"]
b_force25 = sum[grep("b_warm", rownames(sum)), "25%"]

b_photo75 = sum[grep("b_photo", rownames(sum)), "75%"]
b_chill75 = sum[grep("b_chill", rownames(sum)), "75%"]
b_force75 = sum[grep("b_warm", rownames(sum)), "75%"]

b_photo2.5 = sum[grep("b_photo", rownames(sum)), "2.5%"]
b_chill2.5 = sum[grep("b_chill", rownames(sum)), "2.5%"]
b_force2.5 = sum[grep("b_warm", rownames(sum)), "2.5%"]

b_photo97.5 = sum[grep("b_photo", rownames(sum)), "97.5%"]
b_chill97.5 = sum[grep("b_chill", rownames(sum)), "97.5%"]
b_force97.5 = sum[grep("b_warm", rownames(sum)), "97.5%"]



# par(mfrow = c(1,3))
# plot(spInfo$meanBB ~ spInfo$chill, xlab = "Chilling response", ylab = "Estimated budburst"); abline(lm(meanBB~chill, spInfo))
# plot(spInfo$meanBB ~ spInfo$force, xlab = "Forcing response", ylab = "Estimated budburst"); abline(lm(meanBB~force, spInfo))
# plot(spInfo$meanBB ~ spInfo$photo, xlab = "Photoperiod response", ylab = "Estimated budburst"); abline(lm(meanBB~photo, spInfo))

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# I think what we want is a loop that goes through each iteration of the posteriors and calculates the bb, but using 20 for forcing, 12 for photoperiod, 75 (75/10 when rescaled), and smithers to start
# 

#If we are using the old model, we will use the z-scored values for the parameters
photo <- -0.5041133 #8 h photo
siteSM <- 0
force <- -0.5077191 #15 C trt
chill <- -0.4023109 # 5/15

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



# now get the order of diff spp bb that I can use to order the figure
spInfo <- read.csv("input/species_list.csv")

spInfo <- spInfo[order(spInfo$species),]
head(spInfo)
spInfo$meanBB <- colMeans(m)
colnames(m) <- spInfo$species.name

spInfo$Int <- a_sp
spInfo$Int2.5 <- a_sp2.5
spInfo$Int97.5 <- a_sp97.5

spInfo$force <- b_force[2:48]
spInfo$chill <- b_chill[2:48]
spInfo$photo <- b_photo[50:96]

spInfo$force2.5 <- b_force2.5[2:48]
spInfo$chill2.5 <- b_chill2.5[2:48]
spInfo$photo2.5 <- b_photo2.5[50:96]

spInfo$force97.5 <- b_force97.5[2:48]
spInfo$chill97.5 <- b_chill97.5[2:48]
spInfo$photo97.5 <- b_photo97.5[50:96]

spInfo$force25 <- b_force25[2:48]
spInfo$chill25 <- b_chill25[2:48]
spInfo$photo25 <- b_photo25[50:96]

spInfo$force75 <- b_force75[2:48]
spInfo$chill75 <- b_chill75[2:48]
spInfo$photo75 <- b_photo75[50:96]

quantile2575 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.75, 0.25))
  return(returnQuanilte)
}

quantile97.52.5 <- function(x){
  returnQuanilte <- quantile(x, prob = c(0.975, 0.025))
  return(returnQuanilte)
}

bb_quan <- apply(m, 2, quantile2575)
bb_t <- t(bb_quan)
bb_df <- data.frame(bb_t)
colnames(bb_df)[colnames(bb_df) == "X75."] <- "bb25"
colnames(bb_df)[colnames(bb_df) == "X25."] <- "bb75"

bb_quan9725 <- apply(m, 2, quantile97.52.5)
bb_t9725 <- t(bb_quan9725)
bb_df9725 <- data.frame(bb_t9725)
colnames(bb_df9725)[colnames(bb_df9725) == "X97.5."] <- "bb97.5"
colnames(bb_df9725)[colnames(bb_df9725) == "X2.5."] <- "bb2.5"
# colnames(bb_df)[colnames(bb_df) == "X97.5."] <- "bb25"
# colnames(bb_df)[colnames(bb_df) == "X2.75."] <- "bb75"


spInfo <- cbind(spInfo, bb_df)
spInfo <- cbind(spInfo, bb_df9725)
spInfo$value <- spInfo$meanBB


m <- data.frame(m)

long <- melt(m)
names(long) <- c("species.name", "value")

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

names(data) <- c("species.name","value","species","type","transect","meanBB", "Int","Int2.5","Int97.5","spMeanForce", "spMeanChill", "spMeanPhoto","force2.5","chill2.5","photo2.5","force97.5", "chill97.5","photo97.5","force25","chill25","photo25","force75", "chill75","photo75","bb25","bb75","bb2.5","bb97.5", "spacing","chill", "force","photo","intercept")

### Adding intercept
####### Old plot not spaced out ####################3

bbSp <- ggplot() + 
  stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = value), .width = c(.90, .5), cex = 0.75, fill = "mediumpurple2") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Estimated budburst", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 80, label = "a)", cex =5) 


chillSp <- ggplot() + 
  stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = chill), .width = c(.90, .5), cex = 0.75, fill = "#cc6a70ff") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Chilling Response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 
# 
# pdf("figures/testChill2.pdf", width = 10, height =5)
# chillSp
# dev.off()

# Forcing 

forceSp <- ggplot() + 
  stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = force), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Forcing response", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "c)", cex =5) 

# pdf("figures/testforce.pdf", width = 10, height =5)
# forceSp
# dev.off()

# Photoperiod

photoSp <- ggplot() + 
  stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = photo), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Photoperiod response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 

pdf("figures/4panel.pdf", width = 10, height =20)
plot_grid(bbSp, chillSp,forceSp, photoSp, nrow = 4, align = "v", rel_heights = c(1/4, 1/4, 1/4,1.2/3))
dev.off()

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

bbSpace <- ggplot() + 
  geom_point(data = spInfo, aes(x = meanBB, y = value)) +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  geom_pointrange(data = spInfo,aes(x = meanBB, y = value, ymin=bb25, ymax = bb75)) +
  # xlim (15,80) +
  theme_classic() +  
  theme(axis.text.x = element_text( size= 8.9,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Estimated budburst", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 17.5, y = 80, label = "a)", cex =5) + 
  scale_x_continuous( breaks = spMini$meanBB, labels = spMini$species,limits = c(17.2,68)) +
  annotate("text", x = spTop[1,5], y = 15, label = spTop[1,2], cex = 3, angle = 78) +
  annotate("text", x = 24.5, y = 15, label = spTop[2,2], cex = 3, angle = 78) +
  annotate("text", x = 40.3, y = 15, label = spTop[3,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[4,5], y = 15, label = spTop[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[5,5], y = 15, label = spTop[5,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[6,5], y = 15, label = spTop[6,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[7,5], y = 15, label = spTop[7,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[8,5], y = 15, label = spTop[8,2], cex = 3, angle = 78) +
  annotate("text", x = 29.68, y = 15, label = spTop[9,2], cex = 3, angle = 78) +
  annotate("text", x = 43.1, y = 15, label = spTop[10,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[11,5], y = 15, label = spTop[11,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[12,5], y = 15, label = spTop[12,2], cex = 3, angle = 78) 

chillSpace <- ggplot() + 
  geom_point(data = spInfo, aes(x = meanBB, y = chill)) +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  geom_pointrange(data = spInfo,aes(x = meanBB, y = chill, ymin=chill2.5, ymax = chill97.5)) +
  # xlim (15,80) +
  theme_classic() +  
  theme(axis.text.x = element_text( size= 8.9,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Chilling response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 17.5, y = 0.1, label = "b)", cex =5) + 
  scale_x_continuous( breaks = spMini$meanBB, labels = spMini$species,limits = c(17.2,68)) +
  annotate("text", x = spTop[1,5], y = -36.5, label = spTop[1,2], cex = 3, angle = 78) +
  annotate("text", x = 24.5, y = -36.5, label = spTop[2,2], cex = 3, angle = 78) +
  annotate("text", x = 40.3, y = -36.5, label = spTop[3,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[4,5], y = -36.5, label = spTop[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[5,5], y = -36.5, label = spTop[5,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[6,5], y = -36.5, label = spTop[6,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[7,5], y = -36.5, label = spTop[7,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[8,5], y = -36.5, label = spTop[8,2], cex = 3, angle = 78) +
  annotate("text", x = 29.68, y = -36.5, label = spTop[9,2], cex = 3, angle = 78) +
  annotate("text", x = 43.1, y = -36.5, label = spTop[10,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[11,5], y = -36.5, label = spTop[11,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[12,5], y = -36.5, label = spTop[12,2], cex = 3, angle = 78) 

  # stat_eye(data = data, aes(x = meanBB, y = chill), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
  # #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  # theme_classic() +  
  # xlim (15,80) +
  # theme(axis.text.x = element_text( size=10,
  #                                   angle = 78, 
  #                                   hjust=1),
  #       axis.title.y=element_text(size = 12),
  #       axis.title=element_text(size=15) ) + # angle of 55 also works
  # #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  # labs( x = "Species", y = "Estimated budburst", main = NA) +
  # theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 

#
forceSpace <- ggplot() + 
  geom_point(data = spInfo, aes(x = meanBB, y = force)) +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  geom_pointrange(data = spInfo,aes(x = meanBB, y = force, ymin=force2.5, ymax = force97.5)) +
  # xlim (15,80) +
  theme_classic() +  
  theme(axis.text.x = element_text( size= 8.9,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Forcing response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 17.5, y = 1.75, label = "c)", cex =5) + 
  scale_x_continuous( breaks = spMini$meanBB, labels = spMini$species,limits = c(17.2,68)) +
  annotate("text", x = spTop[1,5], y = -20.2, label = spTop[1,2], cex = 3, angle = 78) +
  annotate("text", x = 24.5, y = -20.2, label = spTop[2,2], cex = 3, angle = 78) +
  annotate("text", x = 40.3, y = -20.2, label = spTop[3,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[4,5], y = -20.2, label = spTop[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[5,5], y = -20.2, label = spTop[5,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[6,5], y = -20.2, label = spTop[6,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[7,5], y = -20.2, label = spTop[7,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[8,5], y = -20.2, label = spTop[8,2], cex = 3, angle = 78) +
  annotate("text", x = 29.68, y = -20.2, label = spTop[9,2], cex = 3, angle = 78) +
  annotate("text", x = 43.1, y = -20.2, label = spTop[10,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[11,5], y = -20.2, label = spTop[11,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[12,5], y = -20.2, label = spTop[12,2], cex = 3, angle = 78) 

#   stat_eye(data = data, aes(x = meanBB, y = force), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   xlim (15,80) +
#   theme(axis.text.x = element_blank(),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Forcing response", main = NA)+
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "c)", cex =5) 

# pdf("figures/testforce.pdf", width = 10, height =5)
# forceSp
# dev.off()

# Photoperiod

 photoSpace <- ggplot() + 
geom_point(data = spInfo, aes(x = meanBB, y = photo)) +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  geom_pointrange(data = spInfo,aes(x = meanBB, y = photo, ymin=photo2.5, ymax = photo97.5)) +
  # xlim (15,80) +
  theme_classic() +  
  theme(axis.text.x = element_text( size= 8.9,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Photoperiod response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 17.5, y = 1.75, label = "d)", cex =5) + 
  scale_x_continuous( breaks = spMini$meanBB, labels = spMini$species,limits = c(17.2,68)) +
  annotate("text", x = spTop[1,5], y = -10, label = spTop[1,2], cex = 3, angle = 78) +
  annotate("text", x = 24.5, y = -10, label = spTop[2,2], cex = 3, angle = 78) +
  annotate("text", x = 40.3, y = -10, label = spTop[3,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[4,5], y = -10, label = spTop[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[5,5], y = -10, label = spTop[5,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[6,5], y = -10, label = spTop[6,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[7,5], y = -10, label = spTop[7,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[8,5], y = -10, label = spTop[8,2], cex = 3, angle = 78) +
  annotate("text", x = 29.68, y = -10, label = spTop[9,2], cex = 3, angle = 78) +
  annotate("text", x = 43.1, y = -10, label = spTop[10,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[11,5], y = -10, label = spTop[11,2], cex = 3, angle = 78) +
  annotate("text", x = spTop[12,5], y = -10, label = spTop[12,2], cex = 3, angle = 78) 
#   stat_eye(data = data, aes(x = meanBB, y = photo), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
#   theme_classic() +  xlim (15,80) +
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + # angle of 55 also works
#   #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 

pdf("figures/4panelSpace.pdf", width = 10, height =20)
plot_grid(bbSpace, chillSpace,forceSpace, photoSpace, nrow = 4, align = "v")
dev.off()

#####################################################################
## Eastern spp only #################################################

east <- subset(spInfo, transect != "west")
eastSp <- unique(east$species.name)

dataEast <- data[data$species.name %in% eastSp, ]

bbSpE <- ggplot() + 
  stat_eye(data = dataEast, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = value), .width = c(.90, .5), cex = 0.75, fill = "mediumpurple2") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Estimated budburst", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 80, label = "a)", cex =5) 


chillSpE <- ggplot() + 
  stat_eye(data = dataEast, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = chill), .width = c(.90, .5), cex = 0.75, fill = "#cc6a70ff") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Chilling Response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 

#forcing
forceSpE <- ggplot() + 
  stat_eye(data = dataEast, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = force), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Forcing response", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "c)", cex =5) 

# Photoperiod
photoSpE <- ggplot() + 
  stat_eye(data = dataEast, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = photo), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Photoperiod response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 

pdf("figures/4panelEastern.pdf", width = 10, height =20)
plot_grid(bbSpE, chillSpE,forceSpE, photoSpE, nrow = 4, align = "v", rel_heights = c(1/4, 1/4, 1/4,1.2/3))
dev.off()

#### Get the violin plots to be spaced #########################

# this seems promising
# ggplot() + 
#   geom_violin(data = dataEast, aes(x = meanBB, y = photo, group = species.name, width =1)) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + 
#   scale_x_continuous( breaks = east$meanBB, labels = east$species.name,limits = c(17.2,68)) +
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 
overlappingE <- c("aromel","vibcas","betlen","lyolig","rhopri")
spMiniE <- east[!east$species %in% overlappingE,]
spTopE <- east[east$species %in% overlappingE,]

bbSpaceE <- ggplot() + 
  geom_violin(dataEast, mapping=aes(x = meanBB, y = value, group = species.name, width =1, fill = type)) +
  geom_point(east, mapping= aes(x = meanBB, y = meanBB)) +
  geom_linerange(data = east,aes(x = meanBB, y = meanBB, ymin=bb2.5, ymax = bb97.5) ) +
  theme_classic() +  
  ylim(-1,80) +
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,70)) +
  labs( x = "Species", y = "Estimated budburst", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 80, label = "a)", cex =5) +
  annotate("text", x =40, y = 80, label = "Eastern transect", cex =5) +
  annotate("text", x = spTopE[1,5], y = -0.2, label = spTopE[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -0.2, label = spTopE[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -0.2, label = spTopE[3,2], cex = 3, angle = 78) +
  annotate("text", x = 40, y = -0.2, label = spTopE[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -0.2, label = spTopE[5,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("mediumpurple2","mediumpurple2"))

# intercept
aSpaceE<-  ggplot() + 
  geom_violin(dataEast, mapping=aes(x = meanBB, y = intercept, group = species.name, width =1, fill = type)) +
  geom_point(east, mapping= aes(x = meanBB, y = Int)) +
  geom_linerange(data = east,aes(x = meanBB, y = Int, ymin=Int2.5, ymax = Int97.5)) +
  theme_classic() +  
  ylim(-5,60)+
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,70)) +
  labs( x = "Species", y = "Species intercept", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 60, label = "c)", cex =5) +
  annotate("text", x = spTopE[1,5], y = -4, label = spTopE[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -4, label = spTopE[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -4, label = spTopE[3,2], cex = 3, angle = 78) +
  annotate("text", x = 40, y = -4, label = spTopE[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -4, label = spTopE[5,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("darkolivegreen4","darkolivegreen4"))

chillSpaceE <- ggplot() + 
  geom_violin(dataEast, mapping=aes(x = meanBB, y = chill, group = species.name, width =1, fill = type)) +
  geom_point(east, mapping= aes(x = meanBB, y = chill)) +
  geom_linerange(data = east, aes(x = meanBB, y = chill, ymin=chill2.5, ymax = chill97.5)) +
  theme_classic() +  
  ylim(-45,10) +
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,70)) +
  labs( x = "Species", y = "Chilling response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 10, label = "a)", cex =5) +
  annotate("text", x = 45, y = 10, label = "Eastern transect", cex =5) +
  annotate("text", x = spTopE[1,5], y = -43.7, label = spTopE[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -43.7, label = spTopE[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -43.7, label = spTopE[3,2], cex = 3, angle = 78) +
  annotate("text", x = 40, y = -43.7, label = spTopE[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -43.7, label = spTopE[5,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("#cc6a70ff","#cc6a70ff"))

forceSpaceE <-ggplot() + 
  geom_violin(dataEast, mapping=aes(x = meanBB, y = force, group = species.name, width =1, fill = type)) +
  geom_point(east, mapping= aes(x = meanBB, y = force)) +
  geom_linerange(data = east,aes(x = meanBB, y = force, ymin=force2.5, ymax = force97.5)) +
  theme_classic() +  
  ylim(-30,10) +
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,70)) +
  labs( x = "Species", y = "Forcing response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 10, label = "c)", cex =5) +
  annotate("text", x = spTopE[1,5], y = -29.5, label = spTopE[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -29.5, label = spTopE[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -29.5, label = spTopE[3,2], cex = 3, angle = 78) +
  annotate("text", x = 40, y = -29.5, label = spTopE[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -29.5, label = spTopE[5,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("#f9b641ff","#f9b641ff"))


photoSpaceE <- ggplot() + 
  geom_violin(dataEast, mapping=aes(x = meanBB, y = photo, group = species.name, width =1, fill = type)) +
  geom_point(east, mapping= aes(x = meanBB, y = photo)) +
  geom_linerange(data = east,aes(x = meanBB, y = photo, ymin=photo2.5, ymax = photo97.5)) +
  theme_classic() +  
  ylim(-15,10) +
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniE$meanBB, labels = spMiniE$species,limits = c(23,70)) +
  labs( x = "Species", y = "Photoperiod response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 23, y = 10, label = "e)", cex =5) +
  annotate("text", x = spTopE[1,5], y = -14.5, label = spTopE[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[2,5], y = -14.5, label = spTopE[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[3,5], y = -14.5, label = spTopE[3,2], cex = 3, angle = 78) +
  annotate("text", x = 40, y = -14.5, label = spTopE[4,2], cex = 3, angle = 78) +
  annotate("text", x = spTopE[5,5], y = -14.5, label = spTopE[5,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("cyan4","cyan4"))
 
pdf("figures/4violinEastern2575.pdf", width = 10, height =24)
plot_grid(bbSpaceE, aSpaceE,chillSpaceE,forceSpaceE, photoSpaceE, nrow = 5, align = "v")
dev.off()


#####################################################################
## Western spp only #################################################

west <- subset(spInfo, transect != "east")
westSp <- unique(west$species.name)

dataWest <- data[data$species.name %in% westSp, ]

bbSpW <- ggplot() + 
  stat_eye(data = dataWest, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = value), .width = c(.90, .5), cex = 0.75, fill = "mediumpurple2") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Estimated budburst", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 80, label = "a)", cex =5) 


chillSpW <- ggplot() + 
  stat_eye(data = dataWest, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = chill), .width = c(.90, .5), cex = 0.75, fill = "#cc6a70ff") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Chilling Response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "b)", cex =5) 

#forcing
forceSpW <- ggplot() + 
  stat_eye(data = dataWest, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = force), .width = c(.90, .5), cex = 0.75, fill = "#f9b641ff") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_blank(),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Forcing response", main = NA)+
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 5, label = "c)", cex =5) 

# Photoperiod
photoSpW <- ggplot() + 
  stat_eye(data = dataWest, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = photo), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15) ) + # angle of 55 also works
  #geom_point(data = meanTrophic, aes(x = meanSlope,y = trophic.level, colour = "purple"), shape = 8, size = 3)+
  labs( x = "Species", y = "Photoperiod response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 

pdf("figures/4panelWestern.pdf", width = 10, height =20)
plot_grid(bbSpW, chillSpW,forceSpW, photoSpW, nrow = 4, align = "v", rel_heights = c(1/4, 1/4, 1/4,1.2/3))
dev.off()

#### Get the violin plots to be spaced #########################

# this seems promising
# ggplot() + 
#   geom_violin(data = dataEast, aes(x = meanBB, y = photo, group = species.name, width =1)) +
#   #geom_hline(yintercept = sum[4,1], linetype="dashed") +
#   theme_classic() +  
#   theme(axis.text.x = element_text( size=10,
#                                     angle = 78, 
#                                     hjust=1),
#         axis.title.y=element_text(size = 12),
#         axis.title=element_text(size=15) ) + 
#   scale_x_continuous( breaks = east$meanBB, labels = east$species.name,limits = c(17.2,68)) +
#   labs( x = "Species", y = "Photoperiod response", main = NA) +
#   theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 4, label = "d)", cex =5) 
overlappingW <- c("spialb","betpap","popbal")
spMiniW <- west[!west$species %in% overlappingW,]
spTopW <- west[west$species %in% overlappingW,]

bbSpaceW <-  ggplot() + 
  geom_violin(dataWest, mapping=aes(x = meanBB, y = value, group = species.name, width =1, fill = type)) +
  geom_point(west, mapping= aes(x = meanBB, y = meanBB)) +
  geom_linerange(data = west,aes(x = meanBB, y = meanBB, ymin=bb2.5, ymax = bb97.5)) +
  theme_classic() +  
  ylim(-1,80) +
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,68)) +
  labs( x = "Species", y = "Estimated budburst", main = NA) +
  annotate("text", x =40, y = 80, label = "Western transect", cex =5) +
  theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 80, label = "b)", cex =5) +
  annotate("text", x = spTopW[1,5], y = -0.2, label = spTopW[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[2,5], y =-0.2, label = spTopW[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -0.2, label = spTopW[3,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("mediumpurple2","mediumpurple2"))

# intercept
aSpaceW <-  ggplot() + 
  geom_violin(dataWest, mapping=aes(x = meanBB, y = intercept, group = species.name, width =1, fill = type)) +
  geom_point(west, mapping= aes(x = meanBB, y = Int)) +
  geom_linerange(data = west,aes(x = meanBB, y = Int, ymin=Int2.5, ymax = Int97.5)) +
  theme_classic() + 
  ylim(-5,60)+
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,68)) +
  labs( x = "Species", y = "Species intercept", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 60, label = "d)", cex =5) +
  annotate("text", x = spTopW[1,5], y = -4, label = spTopW[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -4, label = spTopW[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -4, label = spTopW[3,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("darkolivegreen4","darkolivegreen4"))


chillSpaceW <-  ggplot() + 
  geom_violin(dataWest, mapping=aes(x = meanBB, y = chill, group = species.name, width =1, fill = type)) +
  geom_point(west, mapping= aes(x = meanBB, y = chill)) +
  geom_linerange(data = west, aes(x = meanBB, y = chill, ymin=chill2.5, ymax = chill97.5)) +
  theme_classic() +  
  ylim(-45,10) +
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,68)) +
  labs( x = "Species", y = "Chilling response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 10, label = "b)", cex =5) +
  annotate("text", x = 38, y = 10, label = "Western transect", cex =5) +
  annotate("text", x = spTopW[1,5], y =-43.7, label = spTopW[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -43.7, label = spTopW[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -43.7, label = spTopW[3,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("#cc6a70ff","#cc6a70ff"))

forceSpaceW <- ggplot() + 
  geom_violin(dataWest, mapping=aes(x = meanBB, y = force, group = species.name, width =1, fill = type)) +
  geom_point(west, mapping= aes(x = meanBB, y = force)) +
  geom_linerange(data = west,aes(x = meanBB, y = force, ymin=force2.5, ymax = force97.5)) +
  theme_classic() +  
  ylim(-30,10) +
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,68)) +
  labs( x = "Species", y = "Forcing response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 10, label = "d)", cex =5) +
  annotate("text", x = spTopW[1,5], y = -29.5, label = spTopW[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -29.5, label = spTopW[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -29.5, label = spTopW[3,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("#f9b641ff","#f9b641ff"))

photoSpaceW <- ggplot() + 
  geom_violin(dataWest, mapping=aes(x = meanBB, y = photo, group = species.name, width =1, fill = type)) +
  geom_point(west, mapping= aes(x = meanBB, y = photo)) +
  geom_linerange(data = west,aes(x = meanBB, y = photo, ymin=photo2.5, ymax = photo97.5)) +
  theme_classic() +  
  ylim(-15,10) +
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
        axis.title.y=element_text(size = 12),
        axis.title=element_text(size=15),
        legend.position = "none") + 
  scale_x_continuous( breaks = spMiniW$meanBB, labels = spMiniW$species,limits = c(15,70)) +
  labs( x = "Species", y = "Photoperiod response", main = NA) +
  theme(legend.title = element_blank()) +  annotate("text", x = 15, y = 10, label = "f)", cex =5) +
  annotate("text", x = spTopW[1,5], y = -14.5, label = spTopW[1,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[2,5], y = -14.5, label = spTopW[2,2], cex = 3, angle = 78) +
  annotate("text", x = spTopW[3,5], y = -14.5, label = spTopW[3,2], cex = 3, angle = 78) +
  scale_fill_manual(values = c("cyan4","cyan4"))

# pdf("figures/4violinWestern9725.pdf", width = 10, height =24)
# plot_grid(bbSpaceW, aSpaceW,chillSpaceW,forceSpaceW, photoSpaceW, nrow = 5, align = "v")
# dev.off()

#Making better figures for the MS:

pdf("figures/violinBetaAlphaEW.pdf", width = 12, height = 8)
plot_grid(bbSpaceE, bbSpaceW, aSpaceE, aSpaceW, nrow = 2, ncol = 2, align = "v")
dev.off()

pdf("figures/violinCFPEW.pdf", width = 12, height = 16)
plot_grid(chillSpaceE, chillSpaceW, forceSpaceE, forceSpaceW, photoSpaceE, photoSpaceW, nrow = 3, ncol = 2, align = "v")
dev.off()
