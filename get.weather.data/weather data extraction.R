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
fullTable <- rbind(septTable,octTable)
fullTable


#clean up names
names(fullTable)<- gsub(" ", "_", names(fullTable))
names(fullTable)<- gsub("°C", "degC", names(fullTable))
names(fullTable)<- gsub("Definition", "", names(fullTable))
names(fullTable)

#subset for date, max, min, mean temp
fullTable <-subset(fullTable, select=c("DAY", "Max_Temp_degC","Min_Temp_degC", "Mean_Temp_degC"  ))

#convert to numeric data from character
fullTable$Max_Temp_degC <- as.numeric(fullTable$Max_Temp_degC)
fullTable$Min_Temp_degC <- as.numeric(fullTable$Min_Temp_degC)
fullTable$Mean_Temp_degC <- as.numeric(fullTable$Mean_Temp_degC)

fullTable


