#Scraping weather data for bc weather stations
#started by faith jones Nov 9 2020
rm(list = ls())



#
library(rvest)#read_html


#Scrape extreme miunimum temp data for PENTICTON monthly 1953-2012 from website - takes a few minutes to run  
#-----------------------------------------------


#1 year
	webpage <- read_html("https://climate.weather.gc.ca/climate_data/monthly_data_e.html?hlyRange=1953-01-01%7C2012-05-10&dlyRange=1941-04-01%7C2012-05-10&mlyRange=1941-01-01%7C2012-05-01&StationID=1053&Prov=BC&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2020&selRowPerPage=25&Line=1&searchMethod=contains&Month=12&Day=31&txtStationName=penticton&timeframe=3&Year=1953")
    nodes <- html_nodes(webpage, "table")
	dataTableAll <- data.frame(html_table(nodes))
	names(dataTableAll)
	dataTableAll$Extr.Min.Temp.Definition.C[1:12]
 	#ExtreamMinTempsYear  <- dataTableAll$Extr.Min.Temp.Definition.C[!dataTableAll$Extr.Min.Temp.Definition.C == "Summary, average and extreme values are based on the data above."]

#Loop through all years 

firstYear <- 1953  	
lastYear <- 2012
nyear <- lastYear - firstYear
years <- as.character(c(1953: 2012))

websiteShort <- "https://climate.weather.gc.ca/climate_data/monthly_data_e.html?hlyRange=1953-01-01%7C2012-05-10&dlyRange=1941-04-01%7C2012-05-10&mlyRange=1941-01-01%7C2012-05-01&StationID=1053&Prov=BC&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2020&selRowPerPage=25&Line=1&searchMethod=contains&Month=12&Day=31&txtStationName=penticton&timeframe=3&Year="
str(websiteShort)

#make an empty data frame for loop data
ExtreamMinTempsPenticton <- data.frame(matrix(NA,(nyear+1)*12,3))
names(ExtreamMinTempsPenticton) <- c("Year", "Month", "temp1")
ExtreamMinTempsPenticton$Year <- rep(years, each = 12)
ExtreamMinTempsPenticton$Month <- rep(c(1:12), times = (nyear+1))


for (i in 1:(nyear+1)){
	websiteYear <- paste0(websiteShort, years[i], sep = "") # cobine website urkl with year
	webpagei <- read_html(websiteYear)
	nodesi <- html_nodes(webpagei, "table") # select the table of data from teh web page
	dataTablei <- data.frame(html_table(nodesi))#extract table 
	minTempsi <- dataTablei$Extr.Min.Temp.Definition.C[1:12]#select only monthly extream min temps
	ExtreamMinTempsPenticton$temp1[ExtreamMinTempsPenticton$Year == years[i]] <- minTempsi # put data into a big dataframe

}

#Clean the data of text*
ExtreamMinTempsPenticton$temp  <-as.numeric(gsub("Legend.*","",ExtreamMinTempsPenticton$temp1))
ExtreamMinTempsPenticton$temp1 <- NULL

#Get min temp per year
minTempYearPen <- aggregate(ExtreamMinTempsPenticton$temp, by = list(ExtreamMinTempsPenticton$Year), FUN = min,  na.rm = TRUE)
names(minTempYearPen) <- c("Year", "minTemp")
plot(as.numeric(minTempYearPen$Year), minTempYearPen$minTemp, type = "l")

#number of days below -23
ExtreamMinTempsPenticton$tempBelow20 <- ExtreamMinTempsPenticton$temp
ExtreamMinTempsPenticton$tempBelow20 [ExtreamMinTempsPenticton$temp <= -22 ] <- 1
ExtreamMinTempsPenticton$tempBelow20 [ExtreamMinTempsPenticton$temp > -22 ] <- 0
penDaysBelow20 <- aggregate(ExtreamMinTempsPenticton$tempBelow20, by = list(ExtreamMinTempsPenticton$Year), FUN = sum, na.rm = TRUE)
names(penDaysBelow20) <- c("Year", "nDayBelow20")
plot(penDaysBelow20$nDayBelow20 ~ as.numeric(penDaysBelow20$Year), type = "l")

#Scrape  miunimum daily temp data for OSOYOOS datily 1953-2012 from website - takes a few minutes to run  
#-----------------------------------------------


osURL <- "https://climate.weather.gc.ca/climate_data/daily_data_e.html?hlyRange=%7C&dlyRange=1967-04-01%7C2009-09-15&mlyRange=1967-01-01%7C
2007-02-01&StationID=1043&Prov=BC&urlExtension=_
e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2020&selRowPerPage=25&Line=3&searchMethod=contains&Month=9&Day=8&txtStationName=Osoyoos&timeframe=2&Year=2009"

#Split the url so i can specify month and year 
osURLshort1 <- "https://climate.weather.gc.ca/climate_data/daily_data_e.html?hlyRange=%7C&dlyRange=1967-04-01%7C2009-09-15&mlyRange=1967-01-01%7C2007-02-01&StationID=1043&Prov=BC&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2020&selRowPerPage=25&Line=3&searchMethod=contains&Month="
isURLshort2 <- "&Day=8&txtStationName=Osoyoos&timeframe=2&Year="

#1 month (select min temp each month)
	webpage <- read_html(osURL)
    nodes <- html_nodes(webpage, "table")
	dataTableAll <- data.frame(html_table(nodes))
	names(dataTableAll)
	minTemp <- min(as.numeric(dataTableAll$Min.Temp.Definition.C[1:31]), na.rm = TRUE)
 	#ExtreamMinTempsYear  <- dataTableAll$Extr.Min.Temp.Definition.C[!dataTableAll$Extr.Min.Temp.Definition.C == "Summary, average and extreme values are based on the data above."]


#Loop through all data (trial 1 year)

yearsos <- c(1967:2009)
nyearos <- length(yearsos+1)
yearsosChar <- as.character(yearsos)


#Dataframe to inpout data
#make an empty data frame for loop data
MinTempsOs<- data.frame(matrix(NA,nyearos*12,3))
names(MinTempsOs) <- c("Year", "Month", "temp")
MinTempsOs$Year <- rep(yearsos, each = 12)
MinTempsOs$Month <- rep(c(1:12), times = nyearos)


for(yeari in 1:nyearos){


	for(monthi in 1:12){

		#Make up teh url with correct day and month 
		osURLMonth <- paste(osURLshort1, monthi, sep = "" )
		osURLMonth2 <- paste(osURLMonth, isURLshort2, sep = "")
		osURLi<- paste(osURLMonth2, yearsos[yeari], sep = "")

		#read in data table
		webpage <- read_html(osURLi) #set html
	    nodes <- html_nodes(webpage, "table") #specify i want to see the table on teh website
		dataTableAll <- data.frame(html_table(nodes)) # read in table
		minTemp <- min(as.numeric(dataTableAll$Min.Temp.Definition.C[1:31]), na.rm = TRUE) # get minimum temp that month 
	 	MinTempsOs$temp[MinTempsOs$Year == yearsos[yeari] & MinTempsOs$Month == monthi] <- minTemp # put data into a big dataframe

	}

}


minTempYearOs <- aggregate(MinTempsOs$temp, by = list(MinTempsOs$Year), FUN = min,  na.rm = TRUE)
names(minTempYearOs) <- c("Year", "minTemp")
plot(minTempYearOs)
plot(as.numeric(minTempYearOs$Year), minTempYearOs$minTemp, type = "l")
