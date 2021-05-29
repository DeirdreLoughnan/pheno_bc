
library(chillR)
library (Interpol.T)
###load fullTable from weatherdataextraction.R first


## Hope calcs

### Hope field

calibration_l = list(
  Average = data.frame(time_min = rep(5, 12),
                       time_max = rep(14, 12),
                       time_suns = rep(17, 12),
                       C_m = rep(0.35, 12))
)
chillcalcsH <- vector()

xx <- filter(fullTable, site== "Hope")
xx$DAY<-strptime(xx$DAY,"%Y-%m-%d", tz="GMT")


year = as.numeric(format(xx$DAY, "%Y"))
month = as.numeric(format(xx$DAY, "%m"))
day = as.numeric(format(xx$DAY, "%d"))

Min_Temp_degC= data.frame(year, month, day, T = xx$Min_Temp_degC)
Max_Temp_degC = data.frame(year, month, day, T = xx$Max_Temp_degC)

hrly = vector()

for(j in 1:nrow(xx)){
  
  xy <- Th_interp(Min_Temp_degC, Max_Temp_degC, #function that creates 24 values of hourly temperature from minimum and maximum daily values.
                  day = j,
                  tab_calibr = calibration_l$Average)
  
  hrly = rbind(hrly,
               data.frame(
                 DAY = xx[j,'DAY'],
                 Temp = xy$Th,
                 Year = Min_Temp_degC$year[j], 
                 JDay = as.numeric(format(xx[j,'DAY'], "%j")),
                 month = Min_Temp_degC$month[j],
                 day = Min_Temp_degC$day[j],
                 Hour = 1:24
               )
  )
  
}
# Skip interpolation if NA for temperature data
if(apply(hrly, 2, function(x) all(!is.na(x)))["Temp"]) {
  
  chillcalcH <- chilling(hrly, hrly$JDay[1], hrly$JDay[nrow(hrly)]) # 
} else { chillcalcH <- data.frame("Season"=NA,"End_year"=NA,"Chilling_Hours"=NA, "Utah_Model"=NA, "Chill_portions"=NA) }

chillcalcsH <- rbind(chillcalcsH, chillcalcH[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")])


### Hope experimental 

# Chilling in experimental setting. For 53 days at 4 deg C (Hope)

dayseq = seq(as.numeric(format(as.POSIXlt("2019-01-02", "%Y-%m-%d"), "%j")),
             as.numeric(format(as.POSIXlt("2019-02-23", "%Y-%m-%d"), "%j")))

chill2 <- data.frame(
  Year = as.numeric(rep("2019", 53*24)),
  JDay = as.numeric(rep(dayseq, each = 24)),
  Hour = rep(0:23,53),
  Temp = 4
)


chill2calc <- chilling(chill2,chill2$JDay[1] 
                       ,chill2$JDay[nrow(chill2)])

colz = c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")

hopeCalc <- rbind(chill2calc[colz], chillcalcsH)

allHope <- data.frame(Treatment = c("Hope Experimental", "Hope Field"), hopeCalc)



### Smithers calcs  


### Smithers field

calibration_l = list(
  Average = data.frame(time_min = rep(5, 12),
                       time_max = rep(14, 12),
                       time_suns = rep(17, 12),
                       C_m = rep(0.35, 12))
)
chillcalcsS <- vector()

xx <- filter(fullTable, site== "Smithers")
xx$DAY<-strptime(xx$DAY,"%Y-%m-%d", tz="GMT")


year = as.numeric(format(xx$DAY, "%Y"))
month = as.numeric(format(xx$DAY, "%m"))
day = as.numeric(format(xx$DAY, "%d"))

Min_Temp_degC= data.frame(year, month, day, T = xx$Min_Temp_degC)
Max_Temp_degC = data.frame(year, month, day, T = xx$Max_Temp_degC)

hrly = vector()

for(j in 1:nrow(xx)){
  
  xy <- Th_interp(Min_Temp_degC, Max_Temp_degC, #function that creates 24 values of hourly temperature from minimum and maximum daily values.
                  day = j,
                  tab_calibr = calibration_l$Average)
  
  hrly = rbind(hrly,
               data.frame(
                 DAY = xx[j,'DAY'],
                 Temp = xy$Th,
                 Year = Min_Temp_degC$year[j], 
                 JDay = as.numeric(format(xx[j,'DAY'], "%j")),
                 month = Min_Temp_degC$month[j],
                 day = Min_Temp_degC$day[j],
                 Hour = 1:24
               )
  )
  
}
# Skip interpolation if NA for temperature data
if(apply(hrly, 2, function(x) all(!is.na(x)))["Temp"]) {
  
  chillcalcS <- chilling(hrly, hrly$JDay[1], hrly$JDay[nrow(hrly)]) # 
} else { chillcalcS <- data.frame("Season"=NA,"End_year"=NA,"Chilling_Hours"=NA, "Utah_Model"=NA, "Chill_portions"=NA) }

chillcalcsS <- rbind(chillcalcsS, chillcalcS[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")])

### Smithers experimental
# Chilling in experimental setting. For 28 days at 4 deg C (manning)

dayseq = seq(as.numeric(format(as.POSIXlt("2019-01-02", "%Y-%m-%d"), "%j")),
             as.numeric(format(as.POSIXlt("2019-01-29", "%Y-%m-%d"), "%j")))

chill1 <- data.frame(
  Year = as.numeric(rep("2019", 28*24)),
  JDay = as.numeric(rep(dayseq, each = 24)),
  Hour = rep(0:23,28),
  Temp = 4
)


chill1calc <- chilling(chill1,chill1$JDay[1] 
                       ,chill1$JDay[nrow(chill1)])

colz = c("Season", "End_year", "Chilling_Hours","Utah_Model","Chill_portions")

smithCalc <- rbind(chill1calc[colz], chillcalcsS)

#bind smith field and experimental
allSmith <- data.frame(Treatment = c("Smithers Experimental", "Smithers Field"), smithCalc)



#biind smithers and hope together  

allChill <- rbind(allSmith, allHope)

write.csv(allChill, "chilling_values_Hope_Smithers.csv", row.names=FALSE, eol="\r\n")
