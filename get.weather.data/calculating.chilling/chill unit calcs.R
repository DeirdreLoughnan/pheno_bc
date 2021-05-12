# Started in 2016 by Dan Flynn #
# Calculating chill units for 2015 budburst experiment #
# Updated 2017 by Lizzie to my filepath # 

setwd("/Users/Lizzie/Documents/git/projects/treegarden/budexperiments/analyses")
toload <- sort(dir("./input")[grep("Budburst Data", dir('./input'))], T)[1]
load(file.path("input", toload))

loadHFdata <- TRUE

# Dynamic Model

library(chillR)
# Luedeling recommends Dynamic Model (chill portions)

if(!loadHFdata){ # Dan had this commented out but it is where he builds the HF data ... so I adjusted it
# HF chilling from Oct 1 2014
htemp <- read.csv("..//data/hf001-10-15min-m.csv")

htemp$datetime <- as.POSIXlt(htemp$datetime, format = "%Y-%m-%dT%H:%M")

# select years. There are 44 na's
htemp1 <- htemp[htemp$datetime > "2014-09-30" & htemp$datetime < "2015-05-31" & !is.na(htemp$datetime),]
htemp2 <- htemp[htemp$datetime > "2015-11-01" & !is.na(htemp$datetime),]

# aggregate to hourly averages 
ht <- aggregate(airt ~ format(htemp1$datetime, "%Y-%m-%d %H"), mean, data = htemp1)
names(ht)[1] = 'datetime'
save(file = "HF Temp Data to May.Rdata", list=c("ht", "htemp1", "htemp2"))
}

if(FALSE){ # This is Dan's original code, mine above extends to end of May which seems to have no effect on answers and lets me look at chilling later in year
htemp <- read.csv("..//data/hf001-10-15min-m.csv")

htemp$datetime <- as.POSIXlt(htemp$datetime, format = "%Y-%m-%dT%H:%M")

# select years. There are 44 na's
htemp1 <- htemp[htemp$datetime > "2014-09-30" & htemp$datetime < "2015-01-22" & !is.na(htemp$datetime),]
htemp2 <- htemp[htemp$datetime > "2015-11-01" & !is.na(htemp$datetime),]

# aggregate to hourly averages 
ht <- aggregate(airt ~ format(htemp1$datetime, "%Y-%m-%d %H"), mean, data = htemp1)
names(ht)[1] = 'datetime'
save(file = "HF Temp Data.Rdata", list=c("ht", "htemp1", "htemp2")) # to save a lot of time in Sweaving
    }

if(loadHFdata){
# load("HF Temp Data.Rdata")
load("HF Temp Data to May.Rdata")
}

ht$Year <- as.numeric(substr(ht$datetime, 1, 4))
ht$JDay <- as.numeric(format(strptime(substr(ht$datetime, 1, 10), "%Y-%m-%d"), "%j"))
ht$Hour <- as.numeric(substr(ht$datetime, 12, 13))
names(ht)[2] = "Temp"

# Chilling starting in Sept 2014
ht2 <- aggregate(airt ~ format(htemp1$datetime, "%Y-%m-%d %H"), mean, data = htemp1)
names(ht2)[1] = 'datetime'

ht2$Year <- as.numeric(substr(ht2$datetime, 1, 4))
ht2$JDay <- as.numeric(format(strptime(substr(ht2$datetime, 1, 10), "%Y-%m-%d"), "%j"))
ht2$Hour <- as.numeric(substr(ht2$datetime, 12, 13))
names(ht2)[2] = "Temp"

chill0calc <- chilling(ht2, 273, 21) # 56 chill portions by Jan 21 last year.

# What about by Mar 1, April 1 or May 10?
chill60 <- chilling(ht2, 273, 60) # 57 chill portions by 1 March.
chill91 <- chilling(ht2, 273, 91) # 79 chill portions by 1 April.
chill121 <- chilling(ht2, 273, 121) # 100 chill portions by 1 May (10 May is another option, see Richardonson et al. 2009, Tree Phys)

# Helpful code to convert dates ... 
jtomod <- as.Date("130", format = "%j")
format(jtomod, "%d-%b")
modtoj <- as.Date("Apr-01", format = "%b-%d")
format(modtoj, "%j")

# Chilling in experimental setting. From Jan 22 - Feb 12 at 4, then to either 4 or 1.5 for additional 30 d. total: 22 d at 4, then 30 d extra.
# diff(as.POSIXlt(c("2015-01-22", "2015-02-12", "2015-03-14"), "%Y-%m-%d"))

dayseq = seq(as.numeric(format(as.POSIXlt("2015-01-22", "%Y-%m-%d"), "%j")),
             as.numeric(format(as.POSIXlt("2015-03-14", "%Y-%m-%d"), "%j")))

chill1 <- data.frame(
  Year = as.numeric(rep("2015", 52*24)),
  JDay = as.numeric(rep(dayseq, each = 24)),
  Hour = rep(0:23, 52),
  Temp = 4
  )

# bug: Jday vs JDay, returns error for "chillout not found"

chill1calc <- chilling(chill1, 0, 80)

chill2 <- chill1
chill2$Temp = c(rep(4, 22), rep(1.5, 30))

chill2calc <- chilling(chill2, 0, 80) # FEWER chill portions and fewer Utah Model chill hours

colz = c("Chilling_Hours","Utah_Model","Chill_portions")

allcalc <- rbind(chill0calc[colz], chill91[colz], chill121[colz],
    chill0calc[colz]+chill1calc[colz], chill0calc[colz]+chill2calc[colz])

hfallcalc <- data.frame(Treatment = c("Field chilling", "Field chilling to 1 Apr", "Field chilling to 1 May", "Field chill + 4.0° x 30 d", "Field chill + 1.5° x 30 d"), allcalc)
rownames(hfallcalc)=NULL

## St. Hipp data 

stemp <- read.csv("..//data/St. Hip Weather Data.csv")

stemp$datetime <- as.POSIXlt(stemp$Date, format = "%Y-%m-%d")

# interporlate to hourly, based on max min 
# install.packages("Interpol.T")
library(Interpol.T)

# Build a calibration table, here we don't actually have hourly data, use best guess

calibration_l = list(
  Average = data.frame(time_min = rep(5, 12),
                       time_max = rep(14, 12),
                       time_suns = rep(17, 12),
                       C_m = rep(0.35, 12))
                      )

year = as.numeric(substr(stemp$Month, 1, 4))
month = as.numeric(substr(stemp$Month, 6, 7))
day = stemp$Day

Tmin = data.frame(year, month, day, T = stemp$Tmin)
Tmax = data.frame(year, month, day, T = stemp$Tmax)

hrly = vector()

for(i in 1:nrow(stemp)){
  
  xx <- Th_interp(Tmin, Tmax, 
          day = i,
          tab_calibr = calibration_l$Average)
  
  hrly = rbind(hrly,
            data.frame(
              date = strptime(stemp$Date[i], "%Y-%m-%d"),
              Temp = xx$Th,
              Year = Tmin$year[i], 
              JDay = as.numeric(format(strptime(stemp$Date[i], "%Y-%m-%d"), "%j")),
              month = Tmin$month[i],
              day = Tmin$day[i],
              Hour = 1:24
              )
            )
  
  }

st <- hrly
st <- st[st$date > "2014-08-30" & st$date < "2015-03-01",]

stmay <- hrly[hrly$date > "2014-08-30" & hrly$date < "2015-05-31",]

chill0calc.SH <- chilling(st, 273, 21) # 39 chill portions by Jan 20 last year.
calc.SH60 <- chilling(stmay, 273, 60) # 44
calc.SH91 <- chilling(stmay, 273, 91) # 64
calc.SH121 <- chilling(stmay, 273, 121) # 84 chill portions

colz = c("Chilling_Hours","Utah_Model","Chill_portions")


allcalc <- rbind(chill0calc[colz], chill0calc[colz]+chill1calc[colz], chill0calc[colz]+chill2calc[colz],
            chill0calc.SH[colz], chill0calc.SH[colz]+chill1calc[colz], chill0calc.SH[colz]+chill2calc[colz])


allcalc <- data.frame(Site = c("Harvard Forest", "","","St. Hippolyte","",""), Treatment = rep(c("Field chilling","4.0° x 30 d", "1.5° x 30 d"), 2), allcalc)
rownames(allcalc)=NULL

allcalc <- rbind(chill0calc[colz], chill91[colz], chill121[colz],
    chill0calc[colz]+chill1calc[colz], chill0calc[colz]+chill2calc[colz],
    chill0calc.SH[colz], calc.SH91[colz], calc.SH121[colz],
    chill0calc.SH[colz]+chill1calc[colz], chill0calc.SH[colz]+chill2calc[colz])

allcalc <- data.frame(Site = c("Harvard Forest", "","","", "", "St. Hippolyte","","", "", ""),
    Treatment = rep(c("Field chilling until collection","Field chilling to 1 Apr",
       "Field chilling to 1 May", "Field chilling + 4.0° x 30 d",
       "Field chilling + 1.5° x 30 d"), 2), allcalc)

rownames(allcalc)=NULL

# Three models: Chilling Hours, Utah, Dynamic. Chilling hours: all temperatures equally effective (0-7.2 deg C). 
# Utah: Richardson EA, Seeley SD, Walker DR (1974) A model for estimating the completion of rest for Redhaven and Elberta peach trees. Hortscience 9(4):331–332
# No chill units below 1.4 deg C, half a unit from 1.4 to 2.4, full unit from 2.4 to 9.1, HALF from 9.1 to 12.4, FULL from 12.4 to 15.9??!
# Dynamic: Erez 1990, Fishman 1987.


# Sequential, parallel, and unified models: need multiple years of data?
# Ecodormancy, from external factors, Endodormancy from internal factors (broken by chilling)

