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


if(length(grep("deirdreloughnan", getwd()) > 0)) { 
  setwd("~/Documents/github/pheno_bc") 
}  

#load("output/final/ew_phylo_output_newpriors.Rda")
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

head(pheno.t)
# change the df to lowercase:

pheno.t <- merge(pheno.t, spInfo, by = "species")


b.force.both.ew <- sumew[grep("^b_warm\\[", rownames(sumew))]
b.photo.both.ew <- sumew[grep("^b_photo\\[", rownames(sumew))]
b.chill.both.ew <- sumew[grep("^b_chill1", rownames(sumew))]

species.both <- sort(unique(tolower(pheno.t$species)))
species.fact.both <-as.numeric( as.factor(unique(tolower(pheno.t$species))))

type.both <- c("tree", "tree", "tree","tree", "shrub", "shrub", "shrub", "tree", "tree", "tree", "tree", "tree",
               "tree", "shrub","shrub","tree","tree", "shrub", "shrub","shrub",  "shrub", "shrub", "shrub", "shrub", "shrub",
               "tree","tree","tree","tree","tree","tree","tree","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub")
transect.both <- c("western","eastern","eastern","eastern","both","western","western","eastern","eastern","eastern","both","eastern","western","eastern","eastern","eastern","eastern","eastern","eastern","western","eastern","western","eastern","both","eastern","western","eastern","eastern","eastern","eastern","eastern","western","eastern","western","western","western","western","western","western","western","western","western","western","eastern","eastern","western","eastern")

both.ew <- data.frame(species.both, transect.both, species.fact.both, b.force.both.ew, b.photo.both.ew, b.chill.both.ew, type.both)

color_scheme_set("viridis")
ggplot(both.ew, aes(x= b.chill.both.ew, y = b.force.both.ew, col = type.both, shape = transect.both)) +
  geom_point() +
  ylim (-15, -2) +
  xlim (-30, 0) +
  labs( y = "High forcing", x = "High chilling") +
  geom_text(aes(label=species.both),hjust= 0.5, vjust= 1.5, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_blank()) +
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff"))


ggplot(both.ew, aes(x= b.chill.both.ew, y = b.photo.both.ew, col = type.both)) +
  geom_point() +
  ylim (-6, -2) +
  xlim (-28, 0) +
  labs (x = "High chilling", y = "Long photoperiod") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1.5, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_blank())+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff"))


ggplot(both.ew, aes(x= b.force.both, y = b.photo.both, col = type.both)) +
  geom_point() +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1.5, show.legend = F)+
  ylim (-8, -1.5) +
  xlim (-20, 0) +
  labs(x = "High forcing", y = "Long photoperiod") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_blank())+
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff"))


############## BB plots ##################################

term.bb.both.ew <- ddply(pheno, c("species"), summarize, mean = mean(bb, na.rm = TRUE))
names(term.bb.both.ew) <- c("species.both", "mean")
term.bb.both.ew$species.both <- tolower(term.bb.both.ew$species.both)

term.both.ew <- merge(term.bb.both.ew, both.ew, by = "species.both", all =TRUE)
term.both.ew <- term.both.ew[,c("species.both","mean","b.force.both.ew","b.chill.both.ew","b.photo.both.ew","transect.both", "type.both")]
term.both.ew <- term.both.ew[complete.cases(term.both.ew), ] 

tf.both <-  ggplot(term.both.ew, aes(y = b.force.both.ew, x= mean,col = type.both, shape = transect.both)) +
  geom_point() +
  ylim (-15, -1.5) +
  xlim (0, 60) +
  labs(x = "Mean day of budburst", y = "High chilling") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_blank()) + scale_fill_manual(values=c( "#593d9cff","#cc6a70ff")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff"))



tc.both <- ggplot(term.both.ew, aes(y = b.chill.both.ew, x= mean,col = type.both, shape = transect.both)) +
  geom_point() +
  ylim (-25, -1.5) +
  xlim (0, 60) +
  labs(x = "Mean day of budburst", y = "High chilling") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_blank()) +
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff"))


tp.both <- ggplot(term.both.ew, aes(y = b.photo.both.ew, x= mean,col = type.both, shape = transect.both)) +
  geom_point() +
  ylim (-5.5, -2) +
  xlim (0, 60) +
  labs(x = "Mean day of budburst", y = "Long photoperiod")+
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1, show.legend = F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_blank()) +
  scale_fill_manual(values=c( "#593d9cff","#cc6a70ff")) + scale_color_manual(values=c( "#593d9cff","#cc6a70ff"))



plotblank()
text(5,5, "Leafout \n Change (days) due to 5° warming", font = 2, srt = 90)

plotlet("b_photo", "b_warm", 
        #    ylab = "Advance due to 5° warming", 
        #     xlab = "Advance due to 4 hr longer photoperiod", 
        ylim = c(-27, 0.5),
        xlim = c(-16, 0.5),
        group = treeshrub,
        data = sumerl)
legend("topleft", bty = "n", inset = 0.035, legend = "C.", text.font=2)
plotlet("b_chill1", "b_warm", 
        #   ylab = "Advance due to 5° warming", 
        #   xlab = "Advance due to 30d 4° chilling", 
        ylim = c(-27, 0.5),
        xlim = c(-28, -8),
        yaxt="n",
        group = treeshrub,
        data = sumerl)
axis(2, seq(0, -25, by = -5), labels = FALSE)
legend("topleft", bty = "n", inset = 0.035, legend = "D.", text.font=2)
plotblank()
