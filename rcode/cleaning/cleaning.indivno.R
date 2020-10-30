# Cleaning the names of different sites and the labels for individuals

if(length(grep("deirdreloughnan", getwd())>0)) { 
  setwd("~/Documents/github/pheno_bc") 
} else {
  setwd("~/Documents/github/pheno_bc")
}

require(plyr)
require(tidyr)
require(stringr)
# rm(list=ls()) 
# options(stringsAsFactors = FALSE)

mp<-read.csv("pheno_data/MP_indivno.csv",na.strings= "")

sm<-read.csv("pheno_data/SM_indivno.csv", na.strings = "")

mp1<-mp %>% 
  separate(flask_id, c("species","chill","photo","force","flask"), "_")

head(mp1)

unique(mp1$indiv)

#sf = strawberry flats
#eg = east gate
#ll= lightening lake
#wjt= windy joe trail
#tml = twenty minute lake
# ss = ski slope
# ski = nordic ski site

mp1$indiv[mp1$indiv=="Ase"] <- "a_sf"
mp1$indiv[mp1$indiv=="bll"] <- "b_ll"
mp1$indiv[mp1$indiv=="bsf"] <- "b_sf"
mp1$indiv[mp1$indiv=="A"] <- "a"
mp1$indiv[mp1$indiv=="Tml"] <- "tml"
mp1$indiv[mp1$indiv=="sg"] <- "sf"
mp1$indiv[mp1$indiv=="G"] <- "g_eg"
mp1$indiv[mp1$indiv=="csf"] <- "c_sf"
mp1$indiv[mp1$indiv=="ALL"] <- "a_ll"
mp1$indiv[mp1$indiv=="Bski"] <- "b_ski"
mp1$indiv[mp1$indiv=="Atml"] <- "a_tml"
mp1$indiv[mp1$indiv=="Aeg"] <- "a_eg"
mp1$indiv[mp1$indiv=="Geg"] <- "g_eg"
mp1$indiv[mp1$indiv=="G eg"] <- "g_eg"
mp1$indiv[mp1$indiv=="all"] <- "a_ll"
mp1$indiv[mp1$indiv=="Asslope"] <- "a_ski"
mp1$indiv[mp1$indiv=="wjgt"] <- "wjt"
mp1$indiv[mp1$indiv=="Aski"] <- "a_ski"
mp1$indiv[mp1$indiv=="asf"] <- "a_sf"
mp1$indiv[mp1$indiv=="Asfd"] <- "a_sf"
mp1$indiv[mp1$indiv=="feg"] <- "f_eg"
mp1$indiv[mp1$indiv=="Feg"] <- "f_eg"
mp1$indiv[mp1$indiv=="Heg"] <- "h_eg"
mp1$indiv[mp1$indiv=="tmj"] <- "tml"
mp1$indiv[mp1$indiv=="C sf"] <- "c_sf"
mp1$indiv[mp1$indiv=="Feg'"] <- "f_eg"
mp1$indiv[mp1$indiv=="Bsf"] <- "b_sf"
mp1$indiv[mp1$indiv=="BLL"] <- "b_ll"
mp1$indiv[mp1$indiv=="Asf"] <- "a_sf"
mp1$indiv[mp1$indiv=="bLL"] <- "b_ll"
mp1$indiv[mp1$indiv=="sfA"] <- "a_sf"
mp1$indiv[mp1$indiv=="A eg"] <- "a_eg"
mp1$indiv[mp1$indiv=="Dsf"] <- "d_sf"
mp1$indiv[mp1$indiv=="dsf"] <- "d_sf"
mp1$indiv[mp1$indiv=="Beg-juv"] <- "b_eg_juv"
mp1$indiv[mp1$indiv=="atml"] <- "a_tml"
mp1$indiv[mp1$indiv=="WJT"] <- "wjt"
mp1$indiv[mp1$indiv=="A -ss"] <- "a_ss"
mp1$indiv[mp1$indiv=="beg'"] <- "b_eg"
mp1$indiv[mp1$indiv=="ALL\\"] <- "a_ll"
mp1$indiv[mp1$indiv=="SG"] <- "sf"
mp1$indiv[mp1$indiv=="120"] <- "12"
mp1$indiv[mp1$indiv=="sfc"] <- "c_sf"
mp1$indiv[mp1$indiv=="beg"] <- "b_eg"
mp1$indiv[mp1$indiv=="bski"] <- "b_ski"
mp1$indiv[mp1$indiv==""] <- "NA"
mp1$indiv[mp1$indiv=="x"] <- "NA"
mp1$indiv[mp1$indiv=="LL"] <- "ll"
mp1$indiv[mp1$indiv=="aLl"] <- "a_ll"

#########################################################################
#########################################################################
sm1<-sm %>% 
  separate(trt, c("chill","photo","force"), "_")

unique(sm1$indiv)

sm1$indiv[sm1$indiv==" A"] <- "a"
sm1$indiv[sm1$indiv=="A"] <- "a"
sm1$indiv[sm1$indiv==" a"] <- "a"
sm1$indiv[sm1$indiv==" A38"] <- "a"
sm1$indiv[sm1$indiv==" A28"] <- "a"
sm1$indiv[sm1$indiv==" A39"] <- "a"
sm1$indiv[sm1$indiv==" A1"] <- "a"
sm1$indiv[sm1$indiv==" A2"] <- "a"
sm1$indiv[sm1$indiv==" A3"] <- "a"
sm1$indiv[sm1$indiv==" A4"] <- "a"
sm1$indiv[sm1$indiv==" A5"] <- "a"
sm1$indiv[sm1$indiv==" A6"] <- "a"
sm1$indiv[sm1$indiv==" A7"] <- "a"
sm1$indiv[sm1$indiv==" A8"] <- "a"
sm1$indiv[sm1$indiv==" A9"] <- "a"
sm1$indiv[sm1$indiv==" A10"] <- "a"
sm1$indiv[sm1$indiv==" a11"] <- "a"
sm1$indiv[sm1$indiv==" A12"] <- "a"
sm1$indiv[sm1$indiv==" A13"] <- "a"
sm1$indiv[sm1$indiv==" A14"] <- "a"
sm1$indiv[sm1$indiv==" A15"] <- "a"
sm1$indiv[sm1$indiv==" A16"] <- "a"
sm1$indiv[sm1$indiv==" A17"] <- "a"
sm1$indiv[sm1$indiv==" A18"] <- "a"
sm1$indiv[sm1$indiv==" A19"] <- "a"
sm1$indiv[sm1$indiv==" A20"] <- "a"
sm1$indiv[sm1$indiv==" A21"] <- "a"
sm1$indiv[sm1$indiv==" A22"] <- "a"
sm1$indiv[sm1$indiv==" A23"] <- "a"
sm1$indiv[sm1$indiv==" A24"] <- "a"
sm1$indiv[sm1$indiv==" A25"] <- "a"
sm1$indiv[sm1$indiv==" A26"] <- "a"
sm1$indiv[sm1$indiv==" A27"] <- "a"
sm1$indiv[sm1$indiv==" A28"] <- "a"
sm1$indiv[sm1$indiv==" A29"] <- "a"
sm1$indiv[sm1$indiv==" A30"] <- "a"
sm1$indiv[sm1$indiv==" A31"] <- "a"
sm1$indiv[sm1$indiv==" A32"] <- "a"
sm1$indiv[sm1$indiv==" A33"] <- "a"
sm1$indiv[sm1$indiv==" A34"] <- "a"
sm1$indiv[sm1$indiv==" A35"] <- "a"
sm1$indiv[sm1$indiv==" A36"] <- "a"
sm1$indiv[sm1$indiv==" A37"] <- "a"
sm1$indiv[sm1$indiv==" A38"] <- "a"
sm1$indiv[sm1$indiv==" A39"] <- "a"

sm1$indiv[sm1$indiv==" B"] <- "b"
sm1$indiv[sm1$indiv=="B"] <- "b"
sm1$indiv[sm1$indiv==" b"] <- "b"
sm1$indiv[sm1$indiv==" BB1"] <- "b"
sm1$indiv[sm1$indiv==" b1"] <- "b"
sm1$indiv[sm1$indiv==" b2"] <- "b"
sm1$indiv[sm1$indiv==" b3"] <- "b"
sm1$indiv[sm1$indiv==" B4"] <- "b"
sm1$indiv[sm1$indiv==" b5"] <- "b"
sm1$indiv[sm1$indiv==" B6"] <- "b"
sm1$indiv[sm1$indiv==" b7"] <- "b"
sm1$indiv[sm1$indiv==" b8"] <- "b"
sm1$indiv[sm1$indiv==" b9"] <- "b"
sm1$indiv[sm1$indiv==" b10"] <- "b"
sm1$indiv[sm1$indiv==" b11"] <- "b"
sm1$indiv[sm1$indiv==" b12"] <- "b"
sm1$indiv[sm1$indiv==" b13"] <- "b"
sm1$indiv[sm1$indiv==" B14"] <- "b"
sm1$indiv[sm1$indiv==" B15"] <- "b"
sm1$indiv[sm1$indiv==" B16"] <- "b"
sm1$indiv[sm1$indiv==" b17"] <- "b"
sm1$indiv[sm1$indiv==" B18"] <- "b"
sm1$indiv[sm1$indiv==" B19"] <- "b"
sm1$indiv[sm1$indiv==" B20"] <- "b"
sm1$indiv[sm1$indiv==" b21"] <- "b"
sm1$indiv[sm1$indiv==" B22"] <- "b"
sm1$indiv[sm1$indiv==" B23"] <- "b"
sm1$indiv[sm1$indiv==" b24"] <- "b"
sm1$indiv[sm1$indiv==" b25"] <- "b"
sm1$indiv[sm1$indiv==" B26"] <- "b"
sm1$indiv[sm1$indiv==" B27"] <- "b"
sm1$indiv[sm1$indiv==" b28"] <- "b"
sm1$indiv[sm1$indiv==" B29"] <- "b"
sm1$indiv[sm1$indiv==" b30"] <- "b"
sm1$indiv[sm1$indiv==" B31"] <- "b"
sm1$indiv[sm1$indiv==" B32"] <- "b"
sm1$indiv[sm1$indiv==" B33"] <- "b"
sm1$indiv[sm1$indiv==" B34"] <- "b"
sm1$indiv[sm1$indiv==" b35"] <- "b"
sm1$indiv[sm1$indiv==" B36"] <- "b"
sm1$indiv[sm1$indiv==" B37"] <- "b"
sm1$indiv[sm1$indiv==" b38"] <- "b"
sm1$indiv[sm1$indiv==" B39"] <- "b"
sm1$indiv[sm1$indiv==" b40"] <- "b"
sm1$indiv[sm1$indiv==" b41"] <- "b"
sm1$indiv[sm1$indiv==" 3B"] <- "3b"
sm1$indiv[sm1$indiv==" 3b"] <- "3b"

sm1$indiv[sm1$indiv==" onion b"] <- "onion_b"
sm1$indiv[sm1$indiv==" onion B"] <- "onion_b"
sm1$indiv[sm1$indiv==" onionB"] <- "onion_b"
sm1$indiv[sm1$indiv==" onionE"] <- "onion_e"
sm1$indiv[sm1$indiv==" onion E"] <- "onion_b"

sm1$indiv[sm1$indiv==" oc 23"] <- "oct_23"
sm1$indiv[sm1$indiv==" oc 25"] <- "oct_25"
sm1$indiv[sm1$indiv==" oc25"] <- "oct_25"
sm1$indiv[sm1$indiv==" oc 21"] <- "oct_21"
sm1$indiv[sm1$indiv==" oc21"] <- "oct_21"
sm1$indiv[sm1$indiv==" oc21"] <- "oct_21"
sm1$indiv[sm1$indiv==" oc21`"] <- "oct_21"
sm1$indiv[sm1$indiv==" Oc21"] <- "oct_21"
sm1$indiv[sm1$indiv==" 0c21"] <- "oct_21"
sm1$indiv[sm1$indiv==" oc21/8"] <- "oct_21"

sm1$indiv[sm1$indiv==" spibe7_01"] <- "NA"
sm1$indiv[sm1$indiv==" 12(skinny"] <- "12"
sm1$indiv[sm1$indiv==" xx"] <- "NA"
sm1$indiv[sm1$indiv==" XX"] <- "NA"
sm1$indiv[sm1$indiv==" ?"] <- "NA"
sm1$indiv[sm1$indiv==" "] <- "NA"

sm1$indiv[sm1$indiv==" w/ 10"] <- "w/10"
sm1$indiv[sm1$indiv==" w10"] <- "w/10"
sm1$indiv[sm1$indiv==" wih 10"] <- "w/10"
sm1$indiv[sm1$indiv==" w/10"] <- "w/10"
sm1$indiv[sm1$indiv==" area w/9"] <- "w/9"
sm1$indiv[sm1$indiv==" area w/ 9"] <- "w/9"

sm1$indiv[sm1$indiv==" D"] <- "d"
sm1$indiv[sm1$indiv==" b"] <- "b"
sm1$indiv[sm1$indiv==" bC"] <- "b"
sm1$indiv[sm1$indiv==" c"] <- "c"
sm1$indiv[sm1$indiv==" cx"] <- "c"

sm1$indiv[sm1$indiv==" C"] <- "c"
sm1$indiv[sm1$indiv==" E"] <- "e"

sm1$indiv[sm1$indiv==" clw"] <- "cl"
sm1$indiv[sm1$indiv==" Cl"] <- "cl"
sm1$indiv[sm1$indiv==" CL"] <- "cl"
sm1$indiv[sm1$indiv==" cl"] <- "cl"
sm1$indiv[sm1$indiv==" CLA"] <- "cl"

sm1$indiv[sm1$indiv==" #2"] <- "2"
sm1$indiv[sm1$indiv==" #3"] <- "3"
sm1$indiv[sm1$indiv==" #1"] <- "1"

unique(sm1$indiv)


head(sm1)
head(mp1)

indiv<-rbind(sm1,mp1)

indiv$species[indiv$species=="rupar"] <- "rubpar"
indiv$species[indiv$species=="spipyr-missing"] <- "spipyr"
indiv$species[indiv$species=="poptre-missing"] <- "poptre"
indiv$species[indiv$species=="spibe"] <- "spibet"
indiv$species[indiv$species=="popre"] <- "poptre"
indiv$species[indiv$species=="bepap"] <- "betpap"
indiv$species[indiv$species=="shecan8"] <- "shecan"
indiv$species[indiv$species=="Shecan"] <- "shecan"
indiv$species[indiv$species=="corso"] <- "corsto"

head(indiv)

unique(indiv$species)

#creating a unique label for each row, one that deals with the issue of having multiples of some species in each flask
indiv$lab<-paste(indiv$site,indiv$chill, indiv$photo,indiv$force,indiv$flask, indiv$species, sep="_")

#start by identifying the samples that are duplicates, demarcated with T or F
indiv$dup<-duplicated(indiv[,"lab"])

indiv<-indiv %>% 
  group_by(dup) %>% 
  mutate(ref=ifelse(dup, "b", "a"))

indiv$lab2<-paste(indiv$lab, indiv$ref, sep="_")

head(indiv)
#write.csv(indiv, "indiv.no.cleaned.csv")

## GOOO ####
# indiv %>%
#   group_by_at(vars(8)) %>%
#   mutate(ref = ifelse(dup, paste0("id", first(id)), NA_character_))
# 
# test<- indiv %>% 
#   group_by(dup) %>% 
#   mutate(ref=ifelse(dup, paste("b",1:n()), NA_character_))
