# Started Nov 25, 2021 by Deirdre 

# aim of this code is to get a phylogenetic tree for my pheno_bc species
# code adapted from Get.ospree.phylo.R
#https://github.com/FePhyFoFum/big_seed_plant_trees/releases

rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(tidyverse)
library(stringr)
library(ape)
library(phytools)
library(geiger)
library(pez)
library(caper)
library(phangorn)

## set your wd here:
setwd("~/Documents/github/pheno_bc")

bb <- read.csv("notes/species_list.csv", header=TRUE)
bb$temp <- bb$species.name
temp <- str_split_fixed(bb$temp, " ", 2)
bb$phylo.name <- paste(temp[,1], temp[,2], sep="_")
bb$genus <- temp[,1]
bb$species <- temp[,2]

sps.list <- sort(unique(bb$phylo.name))
genus.list=sort(unique(bb$genus))

## load phylo (from Smith and Brown 2019)
phy.plants<-read.tree("data/ALLMB.tre")

## getting a list of genera in S&B's phylo
phy.genera<-unlist(
  lapply(strsplit(phy.plants$tip.label, "_"),function(x){return(x[1])})
)
phy.genera.uniq<-sort(unique(phy.genera))

#temp <- sort(unique(phy.plants$tip.label)) #356305 species

## how many phenobc species are in the phylogeny?
phenosp.genus.inphylo<-genus.list[which(genus.list%in%phy.genera.uniq)]


## first prune the phylogeny to include only these genera
phy.genera.phenobc<-drop.tip(phy.plants,
                            which(!phy.genera %in% phenosp.genus.inphylo)) #8814 tips
rm(phy.plants)
View(sort(phy.genera.phenobc$tip.label))
# now prune just the species I want
phy.plants.phenobc<- drop.tip(phy.genera.phenobc,
                            which(!phy.genera.phenobc$tip.label %in% sps.list))

length(phy.plants.phenobc$tip.label)
sort(phy.plants.phenobc$tip.label)
# only 172 species are in the phylogeny


plot(phy.plants.phenobc,cex=.5)

# save phylogeny
#write.tree(phy.plants.phenobc,"data/SBphylo_phenobc.tre")

# Smith and Brown tree is missing Alnus viridis and Rhamnus frangula, maybe they are in the Zanne tree -- nope

# phy.zanne<-read.tree("data/Vascular_Plants_rooted.dated.tre")
# 
# missing <- c("Alnus_viridis", "Rhamnus_frangula")
# missing.genera <- c("Alnus", "Rhamnus")
# 
# phy.genera.zanne<-unlist(
#   lapply(strsplit(phy.zanne$tip.label, "_"),function(x){return(x[1])})
# )
# phy.genera.uniq<-sort(unique(phy.genera.zanne))
# 
# ## first prune the phylogeny to include only these genera
# phy.genera.zanne.phenobc <- drop.tip(phy.plants,
#                     which(!phy.genera.zanne %in% phenosp.genus.inphylo)) #8814 tips
# View(sort(phy.genera.zanne.phenobc$tip.label))
