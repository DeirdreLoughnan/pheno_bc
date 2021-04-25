# I want to merge the DFlynn data with my data, so first I have to rename to columns etc. 

if(length(grep("deirdreloughnan", getwd()) > 0)) {
  setwd("~/Documents/github/pheno_bc/")
} else {
  setwd("~/Documents/github/pheno_bc")
}

require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
# require(chillR)

rm(list = ls()) 
options(stringsAsFactors = FALSE)

# read in the cleaning phenology data:
df <- read.csv("input/dflynn.data.csv", header=TRUE, na.strings=c("","NA"))
head(df)

dl <- read.csv("input/bc_phenology_Feb52021.csv", header=TRUE, na.strings=c("","NA"))
head(dl)

# In the DF data they used their own scale, where bbch 7 is called 3, 

df <- df[, c("id", "sp", "rep","site","ind","treatcode","warm","photo","chill", "tleaf","lleaf", "day")]

# changing column names 
names(dl)
names(df)

colnames(df)[colnames(df) == "id"] <- "lab2"
colnames(df)[colnames(df) == "site"] <- "population"
colnames(df)[colnames(df) == "sp"] <- "species"
colnames(df)[colnames(df) == "ind"] <- "indiv"
colnames(df)[colnames(df) == "treatcode"] <- "treatment"
colnames(df)[colnames(df) == "warm"] <- "force"
colnames(df)[colnames(df) == "tleaf"] <- "bbch.t"
colnames(df)[colnames(df) == "lleaf"] <- "bbch.l"

# removing the zero chill treatment:
#df <- subset(df, chill != "chill0")

# Changing the labels of the different treatments:
unique(df$force) # cool warm
df$force[df$force == "cool"] <- "LF"
df$force[df$force == "warm"] <- "HF"

# unique(df$chill) # cool warm
# df$chill[df$chill == "chill2"] <- "HC"
# df$chill[df$chill == "chill1"] <- "LC"

unique(df$photo) # cool warm
df$photo[df$photo == "long"] <- "HP"
df$photo[df$photo == "short"] <- "LP"

df$lab3 <- paste(df$lab2, df$chill, df$photo, df$force, sep = "_")
head(df)

# Now I am going to calculate the day of bb for ther terminal bud and the first lateral bud, using my own method, but recall the difference in the DF scale, 3= bbch stage 7 
bursted  <- subset(df, bbch.t >= 3)

# for every unique, label, pop, trt, flask, sp, it is getting the minimum (ie first) day that an obs of 7 or greater was observed
terminalbb <- aggregate(bursted["day"],
                        bursted[c("lab3","lab2", "population", "treatment","chill","force","photo", "species")], 
                        FUN = min)
names(terminalbb)[names(terminalbb) == "day"] <- "tbb"
head(terminalbb)

# Because df data is only for the day the first lat bud bb, I can use the same method and don't have to use the complex methods I used for my data 

lat.bursted  <- subset(df, bbch.l >= 3)

# for every unique, label, pop, trt, flask, sp, it is getting the minimum (ie first) day that an obs of 7 or greater was observed
lateralbb <- aggregate(lat.bursted["day"],
                       lat.bursted[c("lab3","lab2", "population", "treatment","chill","force","photo", "species")], 
                        FUN = min)
names(lateralbb)[names(lateralbb) == "day"] <- "latbb1"
head(terminalbb)

pheno <- merge(terminalbb, lateralbb, by = c("lab2"), all.x = TRUE, all.y = TRUE) 

# There are some the only have lateral bb and others that only have terminal bb and I don't know how to merge these without loose 2 rows of data, here is my round about attempt
pheno.both <- subset(pheno, tbb > 0 & latbb1 > 0)
pheno.both <- pheno.both[ , c("lab2","lab3.x","population.x", "treatment.x", "chill.x", "force.x", "photo.x", "species.x","tbb","latbb1")]
names(pheno.both) <- c("lab2","lab3","population", "treatment", "chill", "force", "photo", "species","tbb","latbb1")

pheno.lat <- subset(pheno, is.na(tbb))
pheno.lat <- pheno.lat[, c("lab2", "lab3.y", "population.y", "treatment.y", "chill.y", "force.y", "photo.y", "species.y","tbb", "latbb1")]
names(pheno.lat) <- c("lab2", "lab3", "population", "treatment", "chill", "force", "photo", "species","tbb", "latbb1")

pheno.term <- subset(pheno, is.na(latbb1))
pheno.term <- pheno.term[ , c("lab2","lab3.x","population.x", "treatment.x", "chill.x", "force.x", "photo.x", "species.x","tbb","latbb1")]
names(pheno.term) <- c("lab2","lab3","population", "treatment", "chill", "force", "photo", "species","tbb","latbb1")

pheno.df <- rbind(pheno.both, pheno.term, pheno.lat)


write.csv(pheno.df, "input/day.of.bb.DFlynn.chill0.csv", row.names = FALSE)
