library(chillR)
library (Interpol.T)
###load fullTable from weatherdataextraction.R first


## Hope calcs
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
    
    chillcalc <- chilling(hrly, hrly$JDay[1], hrly$JDay[nrow(hrly)]) # 
  } else { chillcalc <- data.frame("Season"=NA,"End_year"=NA,"Chilling_Hours"=NA, "Utah_Model"=NA, "Chill_portions"=NA) }
  
  chillcalcsH <- rbind(chillcalcsH, chillcalc[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")])
  
  chillcalcsH <-cbind(Site= "Hope", chillcalcsH)  

  
  
### Smithers calcs  
  
  
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
    
    chillcalc <- chilling(hrly, hrly$JDay[1], hrly$JDay[nrow(hrly)]) # 
  } else { chillcalc <- data.frame("Season"=NA,"End_year"=NA,"Chilling_Hours"=NA, "Utah_Model"=NA, "Chill_portions"=NA) }
  
  chillcalcsS <- rbind(chillcalcsS, chillcalc[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")])
  
  chillcalcsS <-cbind(Site= "Smithers", chillcalcsS)
  
  
#biind smithers and hope together  
  
allChill <- rbind(chillcalcsH, chillcalcsS)
  
write.csv(allChill, "chilling_values_Hope_Smithers.csv", row.names=FALSE, eol="\r\n")
  
  
  
  
  
