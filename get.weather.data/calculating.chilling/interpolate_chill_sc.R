library(chillR)

###load fullTable from weatherdataextraction.R first

calibration_l = list(
  Average = data.frame(time_min = rep(5, 12),
                       time_max = rep(14, 12),
                       time_suns = rep(17, 12),
                       C_m = rep(0.35, 12))
)
chillcalcs <- vector()

  xx <- filter(fullTable, site== "Hope")
  xx$DAY<-strptime(xx$DAY,"%Y-%m-%d", tz="GMT")
  
  #this step is unnecessary?
  
  #add interpolated climate data for studies with warming treatments (ambient plus 0.76, ambient plus 4 degrees)
 #if(length(grep("ambplus0.76",i))==1){xx$Min_Temp_degC<-xx$Min_Temp_degC+0.76;xx$Max_Temp_degC<-xx$Max_Temp_degC+0.76}# pagter15
# if(length(grep("ambplus4",i))==1){xx$Min_Temp_degC<-xx$Min_Temp_degC+4;xx$Max_Temp_degC<-xx$Max_Temp_degC+4}#skre08
  #if(length(grep("ambplus2.25",i))==1){xx$Min_Temp_degC<-xx$Min_Temp_degC+2.25;xx$Max_Temp_degC<-xx$Max_Temp_degC+2.25}
 # if(length(grep("ambplus4.5",i))==1){xx$Min_Temp_degC<-xx$Min_Temp_degC+4.5;xx$Max_Temp_degC<-xx$Max_Temp_degC+4.5}
 # if(length(grep("ambplus6.75",i))==1){xx$Min_Temp_degC<-xx$Min_Temp_degC+6.75;xx$Max_Temp_degC<-xx$Max_Temp_degC+6.75}
# if(length(grep("ambplus9",i))==1){xx$Min_Temp_degC<-xx$Min_Temp_degC+9;xx$Max_Temp_degC<-xx$Max_Temp_degC+9}
  
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
  
  chillcalcs <- rbind(chillcalcs, chillcalc[c("Season","End_year","Chilling_Hours","Utah_Model","Chill_portions")])
  
 write.csv(chillcalcs, "testchill.csv", row.names=FALSE, eol="\r\n")
  
# what is ht2 ??????
#ht2$Year <- as.numeric(substr(ht2$DAYtime, 1, 4))
#ht2$JDay <- as.numeric(format(strptime(substr(ht2$DAYtime, 1, 10), "%Y-%m-%d"), "%j"))
#ht2$Hour <- as.numeric(substr(ht2$DAYtime, 12, 13))
#names(ht2)[2] = "Temp"

#chill0calc <- chilling(ht2, 273, 21) # 56 chill portions by Jan 21 last year.
