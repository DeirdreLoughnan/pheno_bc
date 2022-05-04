# Started May 3, 2022 by deirde

# the aim of this code is to generate the model output for my phenology ms
load("output/final/ew_phylo_output_newpriors.Rda")

sumt <- summary(mdl.ewphylo)$summary 
col4fig <- c("mean","sd","25%","50%","75%","Rhat")
col4table <- c("mean","sd","2.5%","50%","97.5%","Rhat")

# manually to get right order
mu_params <- c("a_z",
               "lam_interceptsa",
               "mu_b_warm",
               "mu_b_photo",
               "mu_b_chill1",
               "b_site",
               "mu_b_inter_wp",
               "mu_b_inter_wc1",
               "mu_b_inter_pc1",
               "mu_b_inter_ws",
               "mu_b_inter_ps",
               "mu_b_inter_sc1"
          
)
meanzew <- sumt[mu_params, col4fig]

rownames(meanzew) = c("a_z",
                      "lam_interceptsa",
                      "mu_b_warm",
                      "mu_b_photo",
                      "mu_b_chill1",
                      "b_site",
                      "mu_b_inter_wp",
                      "mu_b_inter_wc1",
                      "mu_b_inter_pc1",
                      "mu_b_inter_ws",
                      "mu_b_inter_ps",
                      "mu_b_inter_sc1"

)

meanzew.table <- sumt[mu_params, col4table]
row.names(meanzew.table) <- row.names(meanzew)
head(meanzew.table)
#write.table(meanzew.table , "output/term.mdl.esti.dldf.csv", sep = ",", row.names = FALSE)

# Begin by checking to see what cue is most important and whether there are strong correlations between cues:
df.mean.t <- data.frame(bb.force = sumt[grep("b_warm_ncp", rownames(sumt)), 1],
                        bb.photo = sumt[grep("b_photo_ncp", rownames(sumt)), 1],
                        bb.chill = sumt[grep("^b_chill1", rownames(sumt)), 1])

df.mean.t[which(df.mean.t$bb.force > df.mean.t$bb.photo), ] # species 25- rho alb
df.mean.t[which(df.mean.t$bb.chill > df.mean.t$bb.force), ] # 13


# all correlated
summary(lm(bb.force~bb.photo, data=df.mean.t))
summary(lm(bb.force~bb.chill, data=df.mean.t))
summary(lm(bb.force~bb.photo, data=df.mean.t))

pdf(file.path( "figures/changes.pheno.standardized.pdf"), width = 7, height = 8)
par(mfrow = c(1,1), mar = c(5, 10, 2, 1))
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
       col = "midnightblue")
arrows(meanzew[, "75%"], nrow(meanzew):1, meanzew[, "25%"], nrow(meanzew):1,
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
dev.off()

# add advance/delay arrows
# par(xpd=NA)
# arrows(1, 15.5, 6, 15.5, len=0.1, col = "black")
# legend(5, 16.5, legend="delay", bty="n", text.font = 1, cex=0.75)
# arrows(-1, 15.5, -6, 15.5, len=0.1, col = "black")
# legend(-12, 16.5, legend="advance", bty="n", text.font = 1, cex=0.75)
# legend(-2, 16.5, legend="0", bty="n", text.font = 1, cex=0.75)
# par(xpd=FALSE)


## Replicating Flynn Figure 2:

b.force.both <- sumt[grep("b_warm", rownames(sumt))]; b.force.both <- b.force.both[48:94]
b.photo.both <- sumt[grep("b_photo", rownames(sumt))]; b.photo.both <- b.photo.both[48:94]
b.chill.both <- sumt[grep("^b_chill1", rownames(sumt))]


shrubs.both = c("VIBLAN","RHAFRA","RHOPRI","SPIALB","VACMYR","VIBCAS", "AROMEL","ILEMUC", "KALANG", "LONCAN", "LYOLIG", "alninc","alnvir","amelan", "corsto","loninv", "menfer","rhoalb", "riblac","rubpar","samrac","shecan","sorsco","spibet","spipyr","symalb","vacmem","vibedu")
trees.both = c("ACEPEN", "ACERUB", "ACESAC", "BETALL", "BETLEN", "CORCOR", "FAGGRA", "FRANIG", "HAMVIR", "NYSSYL", "POPGRA", "PRUPEN", "QUEALB" , "QUERUB", "QUEVEL", "acegla","betpap", "poptre", "popbal")

pheno.term <- pheno[,c("bb", "chillport.z2", "force.z2", "photo.z2", "species", "lab2","transect")]
pheno.t <- pheno.term[complete.cases(pheno.term), ] # 1780 rows data 


species.both <- sort(unique(tolower(pheno.t$species)))
species.fact.both <-as.numeric( as.factor(unique(tolower(pheno.t$species))))
type.both <- c("tree", "tree", "tree","tree", "shrub", "shrub", "shrub", "tree", "tree", "tree", "tree", "tree",
               "tree", "shrub","shrub","tree","tree", "shrub", "shrub","shrub",  "shrub", "shrub", "shrub", "shrub", "shrub",
               "tree","tree","tree","tree","tree","tree","tree","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub", "shrub","shrub", "shrub", "shrub")
both <- data.frame(species.both, species.fact.both, b.force.both, b.photo.both, b.chill.both, type.both)


#pdf(file.path( "figures/chill_vs_force_dldf.pdf"), width = 7, height = 8)
cf.both <- ggplot(both, aes(x= b.chill.both, y = b.force.both, col = type.both)) +
  geom_point() +
  ylim (-25, 1) +
  xlim (-30, 0) +
  labs( y = "High forcing", x = "High chilling") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()

#pdf(file.path( "figures/chill_vs_photo_dldf.pdf"), width = 7, height = 8)
cp.both <- ggplot(both, aes(x= b.chill.both, y = b.photo.both, col = type.both)) +
  geom_point() +
  ylim (-5.5, 1) +
  xlim (-30, 0) +
  labs (x = "High chilling", y = "Long photoperiod") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#dev.off()
#legend.position = "none"

#pdf(file.path( "figures/force_vs_photo_dldf.pdf"), width = 7, height = 8)
fp.both <- ggplot(both, aes(x= b.force.both, y = b.photo.both, col = type.both)) +
  geom_point() +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1)+
  ylim (-5.5, 0.5) +
  xlim (-22, 0.5) +
  labs(x = "High forcing", y = "Long photoperiod") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#dev.off()

## Plotting the day to bb with the cues on the y-axis 
pheno$species <- tolower(pheno$species)
term.bb.both <- ddply(pheno, c("species"), summarize, mean = mean(bb, na.rm = TRUE))
names(term.bb.both) <- c("species.both", "mean")

term.both <- merge(term.bb.both, both, by = "species.both", all =TRUE)
term.both <- term.both[,c("species.both","mean","b.force.both","b.chill.both","b.photo.both")]
term.both <- term.both[complete.cases(term.both), ] 

tf.both <-  f

tc.both <- ggplot(term.both, aes(y = b.chill.both, x= mean,col = type.both)) +
  geom_point() +
  labs(x = "Mean day of budburst", y = "High chilling") +
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

tp.both <- ggplot(term.both, aes(y = b.photo.both, x= mean,col = type.both)) +
  geom_point() +
  labs(x = "Mean day of budburst", y = "Long photoperiod")+
  geom_text(aes(label=species.both),hjust=0.5, vjust= 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##### General boxplots across treatments:
# As per Lizzie's July 7 post: I should look at how linear these relationships are
# Plot raw data (bb~ chill) with the chill effect from the model plotted on top

plot(pheno.t$tbb ~ pheno.t$chillport.z2, pch = 19, col = "darkgreen")
abline(a = sumt[grep("mu_a", rownames(sumt)), "mean"], b = sumt[grep("mu_chill", rownames(sumt)), "mean"])

plot(sumt[grep("sigma_chill", rownames(sumt)), "mean"])
plot(sumt[grep("log", rownames(sumt)), "mean"])

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


