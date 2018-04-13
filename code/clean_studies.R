## Functions to clean each study individually

## TODO: Implemented the Bjonerass cleaning criteria for each study.

library(data.table)
library(ggplot2)
source("pigfxns.R")



clean_tejon = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	tcdat = dat[study == "tejon"]
	rm(dat) # Free up some space

	tplot = ggplot(tcdat) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)

	if(plotit) tplot;

	# There are a few pigs where we seem to have errant movements
	badpigs = c("tejonF02_F12", "tejonM314")

	badout = array(NA, dim=2)
	fullind = list()

	for(i in 1:length(badpigs)){
	  
	  bp = badpigs[i]
	  bpdat = tcdat[pigID == bp, list(longitude, latitude)]
	  kmeans.result = kmeans(bpdat, 1)
	  centers = kmeans.result$centers[kmeans.result$cluster, ]
	  distances <- sqrt(rowSums((bpdat - centers)^2))
	  
	  outliers <- order(distances, decreasing=T)[1:6]

	  # print(outliers) 
	  badout[i] = outliers[1]
	  fullind[[i]] = sapply(1:length(outliers), function(j) which((tcdat$pigID == bp)
	  												 & (tcdat$longitude == bpdat$longitude[outliers[j]]) & 
	  												 (tcdat$latitude == bpdat$latitude[outliers[j]])))
	  
	  if(plotit){
		  plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1)
		  points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
		  points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
		}
	  
	}

	# Remove 6 outliers for "tejonF02_F12" and 1 outlier for "tejonM314"
	rminds = c(fullind[[1]], fullind[[2]][1])

	trimdat = tcdat[-rminds, ]
	tp = ggplot(trimdat) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)

	if(plotit) tp;

	return(trimdat)

}


clean_txcamp = function(plotit=FALSE){

	# For each pig, remove unrealistic observations
	dat = fread("../data/formatted/full_pig_data.csv")
	tcdat = dat[study == "txcamp"]
	rm(dat) # Free up some space

	# At least three pigs have some errant movements
	tplot = ggplot(tcdat) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# There are a few pigs where we seem to have errant movements: txcamp20141, txcamp20150, txcamp20169.  
	# For these pigs, let's identify potential outliers in movement

	badpigs =c("txcamp20141", "txcamp20150", "txcamp20169")

	badout = array(NA, dim=3)
	fullind = array(NA, dim=3)

	for(i in 1:length(badpigs)){
	  
	  bp = badpigs[i]
	  bpdat = tcdat[pigID == bp, list(longitude, latitude)]
	  kmeans.result = kmeans(bpdat, 1)
	  centers = kmeans.result$centers[kmeans.result$cluster, ]
	  distances <- sqrt(rowSums((bpdat - centers)^2))
	  
	  outliers <- order(distances, decreasing=T)[1:5]

	  # print(outliers) 
	  badout[i] = outliers[1]
	  fullind[i] = which((tcdat$pigID == bp) & (tcdat$longitude == bpdat$longitude[outliers[1]]) & (tcdat$latitude == bpdat$latitude[outliers[1]]))
	  
	  if(plotit){
		  plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1)
		  points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
		  points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
		}
	  
	}

	# Remove the substantial outliers.
	trimdat = tcdat[-fullind, ]
	tp = ggplot(trimdat) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tp;
	return(trimdat)

}

clean_fl_raoul = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "fl_raoul"]
	rm(dat) # Free up some space

	# None of the pigs visually look to have errant points
	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Look at fix time distribution...looks reasonable
	fl$datetime = as.POSIXct(fl$datetime)
	fixtimes = fl[, list(diffs=median(diffunits(datetime))), by=pigID]

	# There are 12 of 18 pigs that meet the criteria
	trimdat = fl
	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=200)), by=pigID]

	return(trimdat)
}


clean_tx_tyler_w2 = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "tx_tyler_w2"]
	rm(dat) # Free up some space

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Look at fix time distribution...looks reasonable
	fl$datetime = as.POSIXct(fl$datetime)
	fixtimes = fl[, list(diffs=median(diffunits(datetime))), by=pigID]

	# There are 12 of 18 pigs that meet the criteria
	trimdat = fl
	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=200)), by=pigID]

	return(trimdat)
}


clean_srel_contact = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "srel_contact"]
	rm(dat) # Free up some space

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# # There are a few pigs where we seem to have errant movements
	# badpigs = c("srel_contactP702", "srel_contactP705", "srel_contactP708")

	# badout = array(NA, dim=2)
	# fullind = list()

	# for(i in 1:length(badpigs)){
	  
	#   bp = badpigs[i]
	#   bpdat = fl[pigID == bp, list(longitude, latitude)]
	#   kmeans.result = kmeans(bpdat, 1)
	#   centers = kmeans.result$centers[kmeans.result$cluster, ]
	#   distances <- sqrt(rowSums((bpdat - centers)^2))
	  
	#   outliers <- order(distances, decreasing=T)[1:8]

	#   # print(outliers) 
	#   badout[i] = outliers[1]
	#   fullind[[i]] = sapply(1:length(outliers), function(j) which((fl$pigID == bp)
	#   												 & (fl$longitude == bpdat$longitude[outliers[j]]) & 
	#   												 (fl$latitude == bpdat$latitude[outliers[j]])))
	  
	#   if(plotit){
	#   	dev.new()
	# 	  plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1, main=bp)
	# 	  points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
	# 	  points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
	# 	}
	  
	# }

	# Look at fix time distribution...looks reasonable
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))
	fixtimes = fl[, list(diffs=median(diffunits(datetime))), by=pigID]

	# There are 12 of 18 pigs that meet the criteria
	trimdat = fl
	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=200)), by=pigID]

	return(trimdat)
}

clean_tx_tyler_w1 = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "tx_tyler_w1"]
	rm(dat) # Free up some space
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	# Trim the final 23 fixes as many are post-capture fixes.
	numtrim = 23
	fldrop = fl[order(pigID, datetime)][, lapply(.SD, function(x) x[-((length(longitude) - numtrim):(length(longitude)))]), by=pigID]

	tplot = ggplot(fldrop) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Look at fix time distribution...looks reasonable
	fixtimes = fldrop[, list(diffs=median(diffunits(datetime))), by=pigID]

	trimdat = fldrop
	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=150)), by=pigID]

	return(trimdat)
}

clean_florida = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "florida"]
	rm(dat) # Free up some space
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	# Trim the final 23 fixes as many are post-capture fixes.
	numtrim = 23
	fldrop = fl[order(pigID, datetime)][, lapply(.SD, function(x) x[-((length(longitude) - numtrim):(length(longitude)))]), by=pigID]

	tplot = ggplot(fldrop) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;


	# There are a few pigs where we seem to have errant movements
	badpigs = c("florida11_06", "florida33_326", "florida44_319")

	badout = array(NA, dim=2)
	fullind = list()

	for(i in 1:length(badpigs)){
	  
	  bp = badpigs[i]
	  bpdat = fldrop[pigID == bp, list(longitude, latitude)]
	  kmeans.result = kmeans(bpdat, 1)
	  centers = kmeans.result$centers[kmeans.result$cluster, ]
	  distances <- sqrt(rowSums((bpdat - centers)^2))
	  
	  outliers <- order(distances, decreasing=T)[1:6]

	  # print(outliers) 
	  badout[i] = outliers[1]
	  fullind[[i]] = sapply(1:length(outliers), function(j) which((fldrop$pigID == bp)
	  												 & (fldrop$longitude == bpdat$longitude[outliers[j]]) & 
	  												 (fldrop$latitude == bpdat$latitude[outliers[j]])))
	  
	  if(plotit){
	  	dev.new()
		  plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1, main=bp)
		  points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
		  points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
		}
	  
	}

	rminds = c(fullind[[1]][1], fullind[[2]][1], fullind[[3]][c(1, 6)])
	trimdat = fldrop[-rminds, ]

	tplot = ggplot(trimdat) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Look at fix time distribution...looks reasonable
	fixtimes = trimdat[, list(diffs=median(diffunits(datetime))), by=pigID]

	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=150)), by=pigID]

	return(trimdat)
}


clean_cali2 = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "cali2"]
	rm(dat) # Free up some space
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Look at fix time distribution...looks reasonable
	trimdat = fl
	fixtimes = trimdat[, list(diffs=median(diffunits(datetime))), by=pigID]

	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=150)), by=pigID]

	return(trimdat)
}

clean_cali0 = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "cali0"]
	rm(dat) # Free up some space
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Look at fix time distribution...looks reasonable
	trimdat = fl
	fixtimes = trimdat[, list(diffs=median(diffunits(datetime))), by=pigID]

	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=150)), by=pigID]

	return(trimdat)
}

clean_mo_kurt0 = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "mo_kurt0"]
	rm(dat) # Free up some space
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Look at fix time distribution...looks reasonable
	trimdat = fl
	fixtimes = trimdat[, list(diffs=median(diffunits(datetime))), by=pigID]

	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=150)), by=pigID]

	return(trimdat)
}


clean_tx_susan = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "tx_susan"]
	rm(dat) # Free up some space
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Trim the first 10 and last 10 fixes as they include errant movements
	numtrim = 10
	fldrop = fl[order(pigID, datetime)][, lapply(.SD, function(x) x[-c((1:10), ((length(longitude) - numtrim):(length(longitude))))]), by=pigID]
	ggplot(fldrop) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)

	# 	# There are a few pigs where we seem to have errant movements
	# badpigs = c("florida11_06", "florida33_326", "florida44_319")

	# badout = array(NA, dim=2)
	# fullind = list()

	# for(i in 1:length(badpigs)){
	  
	#   bp = badpigs[i]
	#   bpdat = fldrop[pigID == bp, list(longitude, latitude)]
	#   kmeans.result = kmeans(bpdat, 1)
	#   centers = kmeans.result$centers[kmeans.result$cluster, ]
	#   distances <- sqrt(rowSums((bpdat - centers)^2))
	  
	#   outliers <- order(distances, decreasing=T)[1:6]

	#   # print(outliers) 
	#   badout[i] = outliers[1]
	#   fullind[[i]] = sapply(1:length(outliers), function(j) which((fldrop$pigID == bp)
	#   												 & (fldrop$longitude == bpdat$longitude[outliers[j]]) & 
	#   												 (fldrop$latitude == bpdat$latitude[outliers[j]])))
	  
	#   if(plotit){
	#   	dev.new()
	# 	  plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1, main=bp)
	# 	  points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
	# 	  points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
	# 	}
	  
	# }


	# Look at fix time distribution...looks reasonable
	trimdat = fldrop



	fixtimes = trimdat[, list(diffs=median(diffunits(datetime))), by=pigID]

	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=150)), by=pigID]



	return(trimdat)
}

clean_tx_tyler_k1 = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "tx_tyler_k1"]
	rm(dat) # Free up some space
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# Look at fix time distribution...looks reasonable
	trimdat = fl

	fixtimes = trimdat[, list(diffs=median(diffunits(datetime))), by=pigID]

	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=150)), by=pigID]

	return(trimdat)
}


clean_michigan = function(plotit=FALSE){

	dat = fread("../data/formatted/full_pig_data.csv")
	fl = dat[study == "michigan"]
	rm(dat) # Free up some space
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	# 
	fldrop = fl[latitude > 43.75] #fl[order(pigID, datetime)][, lapply(.SD, function(x) x[-c((1:10), ((length(longitude) - numtrim):(length(longitude))))]), by=pigID]
	ggplot(fldrop) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)

		# There are a few pigs where we seem to have errant movements
	badpigs = c("michigan36635", "michigan38313")

	badout = array(NA, dim=2)
	fullind = list()

	for(i in 1:length(badpigs)){
	  
	  bp = badpigs[i]
	  bpdat = fldrop[pigID == bp, list(longitude, latitude)]
	  kmeans.result = kmeans(bpdat, 1)
	  centers = kmeans.result$centers[kmeans.result$cluster, ]
	  distances <- sqrt(rowSums((bpdat - centers)^2))
	  
	  outliers <- order(distances, decreasing=T)[1:20]

	  # print(outliers) 
	  badout[i] = outliers[1]
	  fullind[[i]] = sapply(1:length(outliers), function(j) which((fldrop$pigID == bp)
	  												 & (fldrop$longitude == bpdat$longitude[outliers[j]]) & 
	  												 (fldrop$latitude == bpdat$latitude[outliers[j]])))
	  
	  if(plotit){
	  	dev.new()
		  plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1, main=bp)
		  points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
		  points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
		}
	  
	}


	# Look at fix time distribution...looks reasonable
	trimdat = fldrop

	fixtimes = trimdat[, list(diffs=median(diffunits(datetime))), by=pigID]

	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=150)), by=pigID]



	return(trimdat)
}

	
	# dat = fread("../data/formatted/full_pig_data.csv")
	# dat$datetime = as.POSIXct(strptime(dat$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	# # For each study, look at the distributions of speeds
	# get_speed = function(lon, lat, datetime){

	# 	mindiffs = as.numeric(diffunits(datetime)) / 60 # to hours
	# 	dists = as.vector(get_distance_vect(as.data.table(data.frame(longitude=lon, latitude=lat)))) / 1000 # to miles

	# 	return(dists / mindiffs)
	# }

	# temp = dat[pigID %in% c("florida11_06", "florida33_326", "florida44_319")]
	# speed = temp[order(datetime, pigID), list(speed=get_speed(longitude, latitude, datetime), 
	# 																					dists=as.vector(get_distance_vect(as.data.table(data.frame(longitude=longitude, latitude=latitude)))) / 1000), by=pigID]


	# for(studynm in unique(dat$study)){
	# 	fl = dat[study == studynm]

	# 	# Look at fix time distribution...looks reasonable
	# 	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))
	# 	fixtimes = fl[, list(diffs=median(diffunits(datetime))), by=pigID]

	# 	# There are 12 of 18 pigs that meet the criteria
	# 	trimdat = fl
	# 	goodpigs = trimdat[, list(goodpig=runs(datetime, ctime=130, clength=100)), by=pigID]

	# 	print(studynm)
	# 	print(sum(goodpigs$goodpig))
	# 	print(nrow(goodpigs))
	# }





