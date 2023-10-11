library(maptools)
library(mapdata)
library(geosphere)
library(mapproj)
library(igraph)
library(ggmap)
library(ggplot2)
library(rworldmap)
library(maps)
library(mapdata)
library(marmap) #
library(dplyr)
library(plyr)


library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

(west <- data.frame(Longitude = c(-127.1686, -120.7816,-122.1417, - 119.8891), Latitude = c(54.7824, 49.0646,52.1417,50.8837)))
(east <- data.frame(Longitude = c(-72.1900,-74.0248,-71.09543597, -71.37388296), Latitude = c(42.5315, 45.9310,44.92466697,43.99837498)))
                                                                       
# (72.1900,74.0248) (42.5315, 45.9310)
pdf("figures/phenoTraitMap.pdf", height = 6, width = 6)
ggplot(data = world) +
  geom_sf(fill = "beige", col = "black") +
  #geom_sf(data = sites, size = 4, shape = 23, fill = "darkred")+
  geom_point(data = west, aes(x = Longitude, y = Latitude), size = 4, 
             shape = 23, fill = "chartreuse4") +
  geom_point(data = east, aes(x = Longitude, y = Latitude), size = 4, 
             shape = 23, fill = "chartreuse4") +
  coord_sf(xlim = c(-140.15, -54.12), ylim = c(38, 71), expand = T)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()

ggplot(data = world) +
  geom_sf() +
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 30, 
             shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-140.15, -54.12), ylim = c(35, 80), expand = FALSE)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
#nacho says that we can use a different background map that looks a little better 
#(i.e. without having entire countries cut off)
mappath<-"figures"
mapname<-"photoshift_map.pdf"
pdf(file.path(mappath,mapname), width = 9, height = 6)

#quartz()
#mapDevice() #create world map shaped window
sites <- st_as_sf(data.frame(longitude = c(-150), latitude = c(50)), coords = c("longitude", "latitude"), crs = 4326, 
                  agr = "constant")

map("world", fill=TRUE
    ,col="grey65"
    ,boundary=F,interior=F
    ,ylim=c(38, 80), xlim=c(-140.15, -54.12),mar=c(0,0,0,0)
    ,projection='albers',par=c(0,0),wrap=T
    ,resolution=1,border="lightgrey",myborder=0)

points(c(photop_all$long), c(photop_all$lat), pch=21, bg="salmon4", cex=abs(.03*photop_all$time2))

points(c(timeER$long), c(timeER$lat), pch=8,cex=1.8, col="salmon4", lwd=2.5)

for (i in c(1:nrow(photop_all))){
  inter2 <- gcIntermediate(c(photop_all$long[i], photop_all$lat[i]),
                           c(photop_all$long[i], photop_all$lat[i]+photop_all$space2[i]), n=50, addStartEnd=TRUE)
  lines(inter2, col="darkblue",lwd=2)
}
#Make stars for points that  are "ER"
for (i in c(1:nrow(spacerER))){
  inter3 <- gcIntermediate(c(spacerER$long[i], spacerER$lat[i]),
                           c(spacerER$long[i], spacerER$lat[i]+spacerER$space2[i]), n=50, addStartEnd=TRUE)
  arrows(inter3[1,1],min(inter3[,2]),inter3[1,1],max(inter3[,2]),length=.10,code=2, angle=45,col="darkblue",lwd=2)
}
#how do i add a legend to the map?
legend(-140,15, # position
       legend = c("5 days earlier","55 days earlier","105 days earlier","Exceeds range possible with temporal shift","Equivalent spatial shift"), 
       title = "Equivalent shift with climate change",
       pch = c(21,21,21,8),
       pt.cex=c(abs(.03*c(-5,-55,-105)),1.1,NA),
       pt.bg="salmon4",
       col=c("black","black","black","salmon4","darkblue"),
       lty=c(NA,NA,NA,NA,1),
       lwd=c(1,1,1,2.5,2),
       cex = .9,
       bty = "n") # border
mtext("A)", side=3, adj=0)

dev.off()

library(tmap)
library(tmaptools)
library(sf)
library(raster)
library(viridis)

bavaria <- read_sf("northamerica.shp")
lat_map_talk <- tm_shape(North_America) +
  tm_polygons(col="white")+
  tm_shape(coord.neg10)+
  #tm_dots(col= "#593d9cff", size = 0.2, shape = 21)+
  tm_dots(col= "#20A387FF", size = 0.5, shape = 21)+
  tm_shape(coord.pos01) +
  tm_dots(col= "#287D8EFF", size = 0.5, shape = 21)+
  tm_shape(coord.neg21)+
  tm_dots(col= "#55C667FF", size = 0.5, shape = 21)+
  tm_shape(coord.neg32) +
  tm_dots(col= "#95D840FF", size = 0.5, shape = 21)+
  tm_shape(coord.pos12)+
  tm_dots(col= "#482677FF", size = 0.5, shape = 21)+
  tm_shape(coord.neg43) +
  tm_dots(col= "#DCE319FF", size = 0.5, shape = 21)+
  tm_add_legend(type ='symbol', col = c("#482677FF","#287D8EFF","#20A387FF","#55C667FF","#95D840FF","#DCE319FF"), labels = c("10 to 20", "0 to 10", "0 to -10", "-10 to -20", "-20 to -30","-30 to -40"), title = "Shift in phenology") 

