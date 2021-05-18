### LOAD SMITHERS DATA
library(rvest)

#load sept data
url <- "https://climate.weather.gc.ca/climate_data/daily_data_e.html?timeframe=2&Year=2019&Month=9&Day=12&hlyRange=2010-02-02%7C2021-05-12&dlyRange=2010-04-08%7C2021-05-12&mlyRange=%7C&StationID=48628&Prov=BC&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2021&selRowPerPage=25&Line=0&searchMethod=contains&txtStationName=Smithers"
septTable <- url %>%
  xml2::read_html() %>%
  html_nodes(xpath='/html/body/main/div[5]/table') %>%
  html_table()
septTable <- septTable[[1]]
head(septTable)

#change to date format
date<- paste(septTable$DAY,"Sep", "2019", sep = "-")
head(date)
septTable$DAY <- as.Date(date,format="%d-%b-%Y")

#load oct data
url2 <- "https://climate.weather.gc.ca/climate_data/daily_data_e.html?timeframe=2&Year=2019&Month=10&Day=12&hlyRange=2010-02-02%7C2021-05-12&dlyRange=2010-04-08%7C2021-05-12&mlyRange=%7C&StationID=48628&Prov=BC&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2021&selRowPerPage=25&Line=0&searchMethod=contains&txtStationName=Smithers"
octTable <- url2 %>%
  xml2::read_html() %>%
  html_nodes(xpath='/html/body/main/div[5]/table') %>%
  html_table()
octTable <- octTable[[1]]
head(octTable)

#change to date format
date<- paste(octTable$DAY,"Oct", "2019", sep = "-")
head(date)
octTable$DAY <- as.Date(date,format="%d-%b-%Y")
octTable$DAY

#bind sept and oct tables to create full data table
smithTable <- rbind(septTable,octTable)
smithTable


#clean up names
names(smithTable)<- gsub(" ", "_", names(smithTable))
names(smithTable)<- gsub("°C", "degC", names(smithTable))
names(smithTable)<- gsub("Definition", "", names(smithTable))
names(smithTable)

#subset for date, max, min, mean temp
smithTable <-subset(smithTable, select=c("DAY", "Max_Temp_degC","Min_Temp_degC", "Mean_Temp_degC"  ))

#convert to numeric data from character
smithTable$Max_Temp_degC <- as.numeric(smithTable$Max_Temp_degC)
smithTable$Min_Temp_degC <- as.numeric(smithTable$Min_Temp_degC)
smithTable$Mean_Temp_degC <- as.numeric(smithTable$Mean_Temp_degC)

#add site column
smithTable$site <- "Smithers"

smithTable



##LOAD HOPE SLIDE DATA

#load sept data 
url3 <- "https://climate.weather.gc.ca/climate_data/daily_data_e.html?hlyRange=1969-12-31%7C2014-03-20&dlyRange=1975-05-01%7C2021-05-12&mlyRange=1975-01-01%7C2007-02-01&StationID=951&Prov=BC&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2021&selRowPerPage=25&Line=19&searchMethod=contains&txtStationName=hope&timeframe=2&Day=13&Year=2019&Month=9#"
septTableH <- url3 %>%
  xml2::read_html() %>%
  html_nodes(xpath='/html/body/main/div[5]/table') %>%
  html_table()
  
septTableH <- septTableH[[1]]
head(septTableH)

#change to date format
septTableH$DAY <- gsub("†", "", septTableH$DAY)
date <- gsub('\\s+', '', septTableH$DAY)
date <- paste(date,"Sep", "2019", sep = "-")
head(date)
septTableH$DAY <- as.Date(date,format="%d-%b-%Y")
septTableH$DAY


#load oct data
url4 <- "https://climate.weather.gc.ca/climate_data/daily_data_e.html?hlyRange=1969-12-31%7C2014-03-20&dlyRange=1975-05-01%7C2021-05-12&mlyRange=1975-01-01%7C2007-02-01&StationID=951&Prov=BC&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2021&selRowPerPage=25&Line=19&searchMethod=contains&txtStationName=hope&timeframe=2&Day=13&Year=2019&Month=10#"
octTableH <- url4 %>%
  xml2::read_html() %>%
  html_nodes(xpath='/html/body/main/div[5]/table') %>%
  html_table()
octTableH <- octTableH[[1]]
head(octTableH)

#change to date format
octTableH$DAY <- gsub("†", "", octTableH$DAY)
date <- gsub('\\s+', '', octTableH$DAY)
date <- paste(date,"Oct", "2019", sep = "-")
head(date)
octTableH$DAY <- as.Date(date,format="%d-%b-%Y")
octTableH$DAY

#bind sept and oct tables to create full data table
hopeTable <- rbind(septTableH,octTableH)
hopeTable


#clean up names
names(hopeTable)<- gsub(" ", "_", names(hopeTable))
names(hopeTable)<- gsub("°C", "degC", names(hopeTable))
names(hopeTable)<- gsub("Definition", "", names(hopeTable))
names(hopeTable)

#subset for date, max, min, mean temp
hopeTable <-subset(hopeTable, select=c("DAY", "Max_Temp_degC","Min_Temp_degC", "Mean_Temp_degC"  ))

#convert to numeric data from character
hopeTable$Max_Temp_degC <- as.numeric(hopeTable$Max_Temp_degC)
hopeTable$Min_Temp_degC <- as.numeric(hopeTable$Min_Temp_degC)
hopeTable$Mean_Temp_degC <- as.numeric(hopeTable$Mean_Temp_degC)

#add site column
hopeTable$site <- "Hope"

hopeTable


#combine both tables
fullTable <- rbind(smithTable,hopeTable)

#get rid of empty rows (where year=NA)
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
fullTable <- completeFun(fullTable, "DAY")


write.csv(fullTable, 'fullweatherdata.csv')

#Create plots

#Max temp
fullTable %>%
  ggplot( aes(x=DAY, y=Max_Temp_degC, group=site, color=site)) +
  geom_line() +
  ggtitle("Daily Max Temperature") +
  ylab("Temperature °C") +
  xlab("Date")+
  theme_bw()

#Min temp
fullTable %>%
  ggplot( aes(x=DAY, y=Min_Temp_degC, group=site, color=site)) +
  geom_line() +
  ggtitle("Daily Min Temperature") +
  ylab("Temperature °C") +
  xlab("Date")+
  theme_bw()


#Mean temp
fullTable %>%
  ggplot( aes(x=DAY, y=Mean_Temp_degC, group=site, color=site)) +
  geom_line() +
  ggtitle("Daily Mean Temperature") +
  ylab("Temperature °C") +
  xlab("Date")+
  theme_bw()

