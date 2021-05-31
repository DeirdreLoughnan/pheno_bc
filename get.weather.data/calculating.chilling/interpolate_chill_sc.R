rm(list=ls())
options(stringsAsFactors = FALSE)

library(chillR)
library (Interpol.T)

###load fullTable from weatherdataextraction.R first
if(length(grep("deirdreloughnan", getwd())>0)) {
  setwd("~/Documents/github/pheno_bc")
} else{
  setwd("~/deirdre/pheno") #
}

###load fullTable from weatherdataextraction.R first
fullTable <- read.csv("get.weather.data/fullweatherdata.csv")

## Hope calcs

### Hope field

calibration_hope = list(
  Average = data.frame(time_min = rep(5, 12),
                       time_max = rep(14, 12),
                       time_suns = rep(17, 12),
                       C_m = rep(0.35, 12))
)
#chillcalcsHope <- vector()

hope <- subset(fullTable, site== "Hope" & X <132)
hope$DAY<-strptime(hope$DAY,"%Y-%m-%d", tz="GMT")


year = as.numeric(format(hope$DAY, "%Y"))
month = as.numeric(format(hope$DAY, "%m"))
day = as.numeric(format(hope$DAY, "%d"))

Min_Temp_degC_MP= data.frame(year, month, day, T = hope$Min_Temp_degC)
Max_Temp_degC_MP = data.frame(year, month, day, T = hope$Max_Temp_degC)

hrly_MP = vector()

for(j in 1:nrow(hope)){
  
  xy_MP <- Th_interp(Min_Temp_degC_MP, Max_Temp_degC_MP, #function that creates 24 values of hourly temperature from minimum and maximum daily values.
                  day = j,
                  tab_calibr = calibration_hope$Average)
  
  hrly_MP = rbind(hrly_MP,
               data.frame(
                 DAY = hope[j,'DAY'],
                 Temp = xy_MP$Th,
                 Year = Min_Temp_degC_MP$year[j], 
                 JDay = as.numeric(format(hope[j,'DAY'], "%j")),
                 month = Min_Temp_degC_MP$month[j],
                 day = Min_Temp_degC_MP$day[j],
                 Hour = 1:24
               )
  )
  
}
# Skip interpolation if NA for temperature data
if(apply(hrly_MP, 2, function(x) all(!is.na(x)))["Temp"]) {
  
  chillcalcsHope <- chilling(hrly_MP, hrly_MP$JDay[1], hrly_MP$JDay[nrow(hrly_MP)]) # 
} else { chillcalcsHopeope <- data.frame("Season"=NA,"End_year"=NA,"Chilling_Hours"=NA, "Utah_Model"=NA, "Chill_portions"=NA) }

chillingHope <- chillcalcsHope[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")]


### Hope experimental #################################################################

# Chilling in experimental setting. 
#For high chill 55 days at 4 deg C (Hope)
dayseq = seq(as.numeric(format(as.POSIXlt("2019-10-29", "%Y-%m-%d"), "%j")),
             as.numeric(format(as.POSIXlt("2019-12-22", "%Y-%m-%d"), "%j")))
length(dayseq)
chill.hc.hope <- data.frame(
  Year = as.numeric(rep("2019", 55*24)),
  JDay = as.numeric(rep(dayseq, each = 24)),
  Hour = rep(0:23,55),
  Temp = 4
)

chill.hc.hope.calc <- chilling(chill.hc.hope,chill.hc.hope$JDay[1] 
                       ,chill.hc.hope$JDay[nrow(chill.hc.hope)])

# For 28 low chill days at 4 deg C (Hope)
dayseq = seq(as.numeric(format(as.POSIXlt("2019-10-29", "%Y-%m-%d"), "%j")),
             as.numeric(format(as.POSIXlt("2019-11-25", "%Y-%m-%d"), "%j")))
length(dayseq)
chill.lc.hope <- data.frame(
  Year = as.numeric(rep("2019", 28*24)),
  JDay = as.numeric(rep(dayseq, each = 24)),
  Hour = rep(0:23,28),
  Temp = 4
)

chill.lc.hope.calc <- chilling(chill.lc.hope,chill.lc.hope$JDay[1] 
                               ,chill.lc.hope$JDay[nrow(chill.lc.hope)])

colz = c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")

hopeCalc <- rbind(chill.hc.hope.calc[colz], chill.lc.hope.calc[colz], chillingHope)

allHope <- data.frame(Treatment = c("HC", "LC","Field"), hopeCalc)
allHope$site <- "MP"

allHope
### Smithers calcs   ###############################################################

### Smithers field

calibration_l = list(
  Average = data.frame(time_min = rep(5, 12),
                       time_max = rep(14, 12),
                       time_suns = rep(17, 12),
                       C_m = rep(0.35, 12))
)
#chillcalcsSmithers <- vector()

smithers <- subset(fullTable, site== "Smithers" & X < 59)
smithers$DAY<-strptime(smithers$DAY,"%Y-%m-%d", tz="GMT")

year = as.numeric(format(smithers$DAY, "%Y"))
month = as.numeric(format(smithers$DAY, "%m"))
day = as.numeric(format(smithers$DAY, "%d"))

Min_Temp_degC_SM= data.frame(year, month, day, T = smithers$Min_Temp_degC)
Max_Temp_degC_SM = data.frame(year, month, day, T = smithers$Max_Temp_degC)

hrly_SM = vector()

for(j in 1:nrow(smithers)){
  
  xy_SM <- Th_interp(Min_Temp_degC_SM, Max_Temp_degC_SM, #function that creates 24 values of hourly temperature from minimum and maximum daily values.
                  day = j,
                  tab_calibr = calibration_l$Average)
  
  hrly_SM = rbind(hrly_SM,
               data.frame(
                 DAY = smithers[j,'DAY'],
                 Temp = xy_SM$Th,
                 Year = Min_Temp_degC_SM$year[j], 
                 JDay = as.numeric(format(smithers[j,'DAY'], "%j")),
                 month = Min_Temp_degC_SM$month[j],
                 day = Min_Temp_degC_SM$day[j],
                 Hour = 1:24
               )
  )
  
}
# Skip interpolation if NA for temperature data
if(apply(hrly_SM, 2, function(x) all(!is.na(x)))["Temp"]) {
  
  chillcalcS <- chilling(hrly_SM, hrly_SM$JDay[1], hrly_SM$JDay[nrow(hrly_SM)]) # 
} else { chillcalcS <- data.frame("Season"=NA,"End_year"=NA,"Chilling_Hours"=NA, "Utah_Model"=NA, "Chill_portions"=NA) }

chillingSmithers <- chillcalcS[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")]

### Smithers experimental
#For high chill 59 days at 4 deg C (Smithers)
dayseq = seq(as.numeric(format(as.POSIXlt("2019-10-25", "%Y-%m-%d"), "%j")),
             as.numeric(format(as.POSIXlt("2019-12-22", "%Y-%m-%d"), "%j")))
length(dayseq)
chill.hc.smithers <- data.frame(
  Year = as.numeric(rep("2019", 59*24)),
  JDay = as.numeric(rep(dayseq, each = 24)),
  Hour = rep(0:23,59),
  Temp = 4
)

chill.hc.smithers.calc <- chilling(chill.hc.smithers,chill.hc.smithers$JDay[1] 
                               ,chill.hc.smithers$JDay[nrow(chill.hc.smithers)])

# For high chill 32 days at 4 deg C (Hope)
dayseq = seq(as.numeric(format(as.POSIXlt("2019-10-25", "%Y-%m-%d"), "%j")),
             as.numeric(format(as.POSIXlt("2019-11-25", "%Y-%m-%d"), "%j")))
length(dayseq)
chill.lc.smithers <- data.frame(
  Year = as.numeric(rep("2019", 32*24)),
  JDay = as.numeric(rep(dayseq, each = 24)),
  Hour = rep(0:23,32),
  Temp = 4
)

chill.lc.smithers.calc <- chilling(chill.lc.smithers,chill.lc.smithers$JDay[1] 
                               ,chill.lc.smithers$JDay[nrow(chill.lc.smithers)])

colz = c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")

smithersCalc <- rbind(chill.hc.smithers.calc[colz], chill.lc.smithers.calc[colz], chillingSmithers)

allsmithers <- data.frame(Treatment = c("HC", "LC","Field"), smithersCalc)
allsmithers$site <- "SM"

allsmithers
#biind smithers and hope together  

allChill <- rbind(allsmithers, allHope)
allChill
write.csv(allChill, "chilling_values_Hope_Smithers.csv", row.names=FALSE, eol="\r\n")
