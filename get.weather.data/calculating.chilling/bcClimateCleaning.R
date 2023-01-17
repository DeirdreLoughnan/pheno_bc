# started Jan 16, 2023 by DL

# Clean climate data and calculate average field chilling for BC climate data:

setwd("~/Documents/github/pheno_bc")

smithers <- read.csv("input/smithersDaily_1943_2018.csv")
smithers$pop <- "smithers"

hope <- read.csv("input/hopeDaily_1975_2018.csv")
hope$pop <- "hope"

clim <- rbind(smithers, hope)

# subset to 1975, convert blanks to NA etc

clim75 <- subset(clim, year > 1974)
range(clim75$year)

sort(unique(clim75$tempMax))

clim75$tempMax <- gsub("LegendEE","", clim75$tempMax)
clim75$tempMax[clim75$tempMax == ""] <- "NA"
clim75$tempMax[clim75$tempMax == "LegendMM"] <- "NA"

sort(unique(clim75$tempMin))

clim75$tempMin <- gsub("LegendEE","", clim75$tempMin)
clim75$tempMin[clim75$tempMin == ""] <- "NA"
clim75$tempMin[clim75$tempMin == "LegendMM"] <- "NA"

sort(unique(clim75$tempMean))

clim75$tempMean <- gsub("LegendEE","", clim75$tempMean)
clim75$tempMean[clim75$tempMean == ""] <- "NA"
clim75$tempMean[clim75$tempMean == "LegendMM"] <- "NA"

write.csv(clim75, "winterTemp.csv", row.names = F)

## Now calculate the chilling:
sm <- subset(clim75, pop == "smithers")
mp <- subset(clim75, pop == "hope")

calibration = list(
  Average = data.frame(time_min = rep(5, 12),
                       time_max = rep(14, 12),
                       time_suns = rep(17, 12),
                       C_m = rep(0.35, 12))
)

year = as.numeric(sm$year)
month = as.numeric(sm$month)
day = sm$day

Tmin = data.frame(year, month, day, T = as.numeric(sm$tempMin))
Tmax = data.frame(year, month, day, T = as.numeric(sm$tempMax))

hrly_SM = vector()

for(j in 1:nrow(sm)){
  
  xy_SM <- Th_interp(Tmin, Tmax,
                     day = j,
                     tab_calibr = calibration$Average)
  
  hrly_SM = rbind(hrly_SM,
                  data.frame(
                    DAY = clim75[j,'day'],
                    Temp = xy_SM$Th,
                    Year = Tmin$year[j], 
                    JDay = as.numeric(clim75[j,'day'], "day"),
                    month = Tmin$month[j],
                    day = Tmin$day[j],
                    Hour = 1:24
                  )
  )
  
}

# now manning park
year = as.numeric(mp$year)
month = as.numeric(mp$month)
day = mp$day

Tmin = data.frame(year, month, day, T = as.numeric(mp$tempMin))
Tmax = data.frame(year, month, day, T = as.numeric(mp$tempMax))

hrly_MP = vector()

for(j in 1:nrow(mp)){
  
  xy_MP <- Th_interp(Tmin, Tmax,
                     day = j,
                     tab_calibr = calibration$Average)
  
  hrly_MP = rbind(hrly_MP,
                  data.frame(
                    DAY = clim75[j,'day'],
                    Temp = xy_MP$Th,
                    Year = Tmin$year[j], 
                    JDay = as.numeric(clim75[j,'day'], "day"),
                    month = Tmin$month[j],
                    day = Tmin$day[j],
                    Hour = 1:24
                  )
  )
 
}

write.csv(hrly_MP, "manningparkHourly.csv", row.names = F)


hrly_MP <- read.csv("input/manningparkHourly.csv")
hrly_SM <- read.csv("input/smithersHourly.csv")

head(hrly_MP)
hrly_MP$JDay <- do.call(paste, list(hrly_MP$month, hrly_MP$day, hrly_MP$Year))
hrly_MP$dateCode <- as.Date(hrly_MP$JDay, format=c("%m %d %Y"))
hrly_MP$julian <- as.numeric(format(hrly_MP$dateCode, "%j"))

if(apply(hrly_MP, 2, function(x) all(!is.na(x)))["Temp"]) {
  
  chillcalcsHope <- chilling(hrly_MP, hrly_MP$julian[1], hrly_MP$day[nrow(hrly_MP)]) # 
} else { chillcalcsHopeope <- data.frame("Season"=NA,"End_year"=NA,"Chilling_Hours"=NA, "Utah_Model"=NA, "Chill_portions"=NA) }

chillingHope <- chillcalcsHope[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")]

# Smithers:
head(hrly_SM)
hrly_SM$JDay <- do.call(paste, list(hrly_SM$month, hrly_SM$day, hrly_SM$Year))
hrly_SM$dateCode <- as.Date(hrly_SM$JDay, format=c("%m %d %Y"))
hrly_SM$julian <- as.numeric(format(hrly_SM$dateCode, "%j"))

if(apply(hrly_SM, 2, function(x) all(!is.na(x)))["Temp"]) {
  
  chillcalcsSmithers <- chilling(hrly_SM, hrly_SM$julian[1], hrly_SM$day[nrow(hrly_SM)]) # 
} else { chillcalcsSmithers <- data.frame("Season"=NA,"End_year"=NA,"Chilling_Hours"=NA, "Utah_Model"=NA, "Chill_portions"=NA) }

chillingSmithers <- chillcalcsSmithers[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")]
