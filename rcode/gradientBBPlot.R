# started Feb 15, 2023 by Deirdre 

# Aim of this code is to make a cool figure that estimates BB from set conditions and them ranks the early to late bb individuals; do the cues responses similarly vary?

# the set cues will be: 12h photoperiod, 20C, high chill - 75/10
require(rstan)

if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  else{
  setwd("/home/deirdre/pheno_bc") # for midge
}

load("output/bb_4sites_phylo_mini.Rda")
sum <- summary(mdl.4phyloMini)$summary 

post <- rstan::extract(mdl.4phyloCont)

a_sp = mean(sum[grep("a_sp", rownames(sum)), 1])
mu_b_warm = sum[grep("b_warm", rownames(sum)), 1]
mu_b_photo = sum[grep("mu_b_photo", rownames(sum)), 1]
mu_b_chill1 = sum[grep("mu_b_chill1", rownames(sum)), 1]
mu_b_inter_pc1 = sum[grep("mu_b_inter_pc1", rownames(sum)), 1]
mu_b_inter_wp = sum[grep("mu_b_inter_wp", rownames(sum)), 1]
mu_b_inter_wc1 = sum[grep("mu_b_inter_wc1", rownames(sum)), 1]
mu_b_inter_ws2 = sum[grep("mu_b_inter_ws2", rownames(sum)), 1]
mu_b_inter_s2c1 = sum[grep("mu_b_inter_s2c1", rownames(sum)), 1]
mu_b_inter_ps2 = sum[grep("mu_b_inter_ps2", rownames(sum)), 1]
mu_b_inter_ws3 = sum[grep("mu_b_inter_ws3", rownames(sum)), 1]
mu_b_inter_s3c1 = sum[grep("mu_b_inter_s3c1", rownames(sum)), 1]
mu_b_inter_ps3 = sum[grep("mu_b_inter_ps3", rownames(sum)), 1]
mu_b_inter_ws4 = sum[grep("mu_b_inter_ws4", rownames(sum)), 1]
mu_b_inter_s4c1 = sum[grep("mu_b_inter_s4c1", rownames(sum)), 1]
mu_b_inter_ps4 = sum[grep("mu_b_inter_ps4", rownames(sum)), 1]
b_site2 = sum[grep("b_site2", rownames(sum)), 1]
b_site3 = sum[grep("b_site3", rownames(sum)), 1]
b_site4 = sum[grep("b_site4", rownames(sum)), 1]

# #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
# I think what we want is a loop that goes through each iteration of the posteriors and calculates the bb, but using 20 for forcing, 12 for photoperiod, 75 (75/10 when rescaled), and smithers to start
# 

photo <- 8
siteSM <- 0
force <- 5
chill <- 5.5

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
head(spInfo)
spInfo$meanBB <- colMeans(m)

colnames(m) <- sort(spInfo$species.name)

m <- data.frame(m)
long <- melt(m)
names(long) <- c("species.name", "value")

long <- merge(long,spInfo, by = "species.name")
spOrderData <- spInfo[order(spInfo$meanBB),]
spOrder <- as.factor(spOrderData$species.name)
# longPhotoInfo$mean <- rowMeans(longPhotoInfo[,c("Site1","Site2","Site3","Site4")], na.rm=TRUE)

data <- long[order(long$meanBB),]

data$species.name <- factor(data$species.name, levels=unique(data$species.name) )
data <- transform(data, variable=reorder(species.name, -meanBB) ) 

bbSp <- ggplot() + 
  stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = value), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
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
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "c)", cex =5) 



pdf("figures/test.pdf", width = 10, height =5)
bbSp
dev.off()

bChill <- data.frame(post$b_chill1[1:1000,])
colnames(bChill) <- sort(spInfo$species.name)
longChill <- melt(bChill)
names(longChill) <- c("species.name", "chill")


data <- cbind(data, longChill$chill)
head(longest)
names(data) <- c("species.name","value","species","type","transect","meanBB", "spp2","chill")

chillSp <- ggplot() + 
  stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = chill), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
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
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "c)", cex =5) 

pdf("figures/testChill.pdf", width = 10, height =5)
chillSp
dev.off()

# Forcing 
bForce <- data.frame(post$b_warm[1:1000,])
colnames(bForce) <- sort(spInfo$species.name)
longForce <- melt(bForce)
names(longForce) <- c("species.name", "force")


data <- cbind(data, longForce$force)

names(data) <- c("species.name","value","species","type","transect","meanBB", "spp2","chill","force")

forceSp <- ggplot() + 
  stat_eye(data = data, aes(x = fct_reorder(.f = species.name, .x = value, .fun = mean, na.rm = T), y = force), .width = c(.90, .5), cex = 0.75, fill = "cyan4") +
  #geom_hline(yintercept = sum[4,1], linetype="dashed") +
  theme_classic() +  
  theme(axis.text.x = element_text( size=10,
                                    angle = 78, 
                                    hjust=1),
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
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "c)", cex =5) 

pdf("figures/testforce.pdf", width = 10, height =5)
forceSp
dev.off()

# Photoperiod

bPhoto <- data.frame(post$b_photo[1:1000,])
colnames(bPhoto) <- sort(spInfo$species.name)
longPhoto <- melt(bPhoto)
names(longPhoto) <- c("species.name", "photo")


data <- cbind(data, longPhoto$photo)
head(longest)
names(data) <- c("species.name","value","species","type","transect","meanBB", "spp2","chill","force","photo")

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
  labs( x = "Species", y = "Estimated budburst", main = NA)+ 
  scale_color_identity(name = "Model fit",
                       breaks = c("black"),
                       labels = c("Model Posterior"),
                       guide = guide_legend(override.aes = list(
                         linetype = c(NA),
                         shape = c(8)))) +
  theme(legend.title = element_blank()) +  annotate("text", x = 1, y = 10, label = "c)", cex =5) 

pdf("figures/testphoto.pdf", width = 10, height =5)
photoSp
dev.off()