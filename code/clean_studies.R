## Functions to clean each study individually

## TODO: Implemented the Bjonerass cleaning criteria for each study.

library(data.table)
library(ggplot2)
library(yaml)
source("pigfxns.R")

anal_params = yaml.load_file("analysis_parameters.yml")
maxtime = anal_params$maxtime
minlength = anal_params$minlength

clean_pigs = function(studynm, plotit=FALSE){
	# Clean the pig movement data for a particular study
	#
	# Parameters
	# ----------
	# studynm : str
	# 	Name of study to be cleaned
	#
	# Notes
	# -----
	# All cleaning proceeds with the following steps
	# 1. Visualize the study data
	# 2. Perform study-specific cleaning
	# 3. Remove all points where pigs move faster than 40 km per hour (Mayer et al. book)

	dat = fread("../data/formatted/full_pig_data.csv")

	fl = dat[study == studynm]
	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	tplot = ggplot(fl) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
	if(plotit) tplot;

	################################
	### Study specific cleaning ####
	################################

	if(studynm == 'tejon'){

		trimdat = clean_tejon(fl)

	} else if(studynm == 'txcamp'){

		trimdat = clean_txcamp(fl)

	} else if(studynm == "fl_raoul"){

		trimdat = clean_fl_raoul(fl)

	} else if(studynm == "tx_tyler_w2"){

		trimdat = clean_tx_tyler_w2(fl)

	} else if(studynm == "srel_contact"){

		trimdat = clean_srel_contact(fl)

	} else if(studynm == "tx_tyler_w1"){

		trimdat = clean_tx_tyler_w1(fl)

	} else if(studynm == "florida"){

		trimdat = clean_florida(fl)

	} else if(studynm == 'cali0'){

		trimdat = clean_cali0(fl)

	}  else if(studynm == 'cali1'){

		trimdat = clean_cali1(fl)

	} else if(studynm == 'cali2'){

		trimdat = clean_cali2(fl)

	} else if(studynm == 'cali3'){

		trimdat = clean_cali3(fl)

	} else if(studynm == 'cali4'){

		trimdat = clean_cali4(fl)

	} else if(studynm %like% 'mo_kurt'){

		trimdat = clean_mo_kurt0(fl)

	} else if(studynm == "tx_susan"){

		trimdat = clean_tx_susan(fl)

	} else if(studynm == "tx_tyler_k1"){

		trimdat = clean_tx_tyler_k1(fl)

	} else if(studynm == "michigan"){

		trimdat = clean_michigan(fl)	
	} else{
		trimdat = fl
	}

	################################
	### General cleaning 				####
	################################

	# Remove pigs that are moving too fast
	cat("Removing unrealistic pigs speeds", "\n")
	trimdat = trim_speed(trimdat, maxspeed=40)

	# Examine how many pigs will actually be useful in the analysis
	goodpigs = trimdat[, list(goodpig=runs(datetime, 
																				ctime=maxtime, 
																				clength=minlength)), by=pigID]

	return(trimdat)

}

#############################################
##### Study specific cleaning functions #####
#############################################

clean_tejon = function(fl){
	# Cleaning for tejon data

	# There are a few pigs where we seem to have errant movements
	badpigs = c("tejonF02_F12", "tejonM314")

	badout = array(NA, dim=2)
	fullind = list()

	for(i in 1:length(badpigs)){
	  
	  bp = badpigs[i]
	  bpdat = fl[pigID == bp, list(longitude, latitude)]
	  kmeans.result = kmeans(bpdat, 1)
	  centers = kmeans.result$centers[kmeans.result$cluster, ]
	  distances <- sqrt(rowSums((bpdat - centers)^2))
	  
	  outliers <- order(distances, decreasing=T)[1:6]

	  # print(outliers) 
	  badout[i] = outliers[1]
	  fullind[[i]] = sapply(1:length(outliers), function(j) which((fl$pigID == bp)
	  												 & (fl$longitude == bpdat$longitude[outliers[j]]) & 
	  												 (fl$latitude == bpdat$latitude[outliers[j]])))
	  
	 #  if(plotit){
		#   plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1)
		#   points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
		#   points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
		# }
	  
	}

	# Remove 6 outliers for "tejonF02_F12" and 1 outlier for "tejonM314"
	rminds = c(fullind[[1]], fullind[[2]][1])

	trimdat = fl[-rminds, ]

	return(trimdat )

	trimdat = trim_speed(trimdat, maxspeed=40)
	tp = ggplot(trimdat) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)

	#if(plotit) tp;

	return(trimdat)

}


clean_txcamp = function(fl){

	# There are a few pigs where we seem to have errant movements: txcamp20141, txcamp20150, txcamp20169.  
	# For these pigs, let's identify potential outliers in movement

	badpigs =c("txcamp20141", "txcamp20150", "txcamp20169")

	badout = array(NA, dim=3)
	fullind = array(NA, dim=3)

	for(i in 1:length(badpigs)){
	  
	  bp = badpigs[i]
	  bpdat = fl[pigID == bp, list(longitude, latitude)]
	  kmeans.result = kmeans(bpdat, 1)
	  centers = kmeans.result$centers[kmeans.result$cluster, ]
	  distances <- sqrt(rowSums((bpdat - centers)^2))
	  
	  outliers <- order(distances, decreasing=T)[1:5]

	  # print(outliers) 
	  badout[i] = outliers[1]
	  fullind[i] = which((fl$pigID == bp) & (fl$longitude == bpdat$longitude[outliers[1]]) & (fl$latitude == bpdat$latitude[outliers[1]]))
	  
	 #  if(plotit){
		#   plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1)
		#   points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
		#   points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
		# }
	  
	}

	# Remove the substantial outliers.
	trimdat = fl[-fullind, ]
	return(trimdat)

}

clean_fl_raoul = function(fl){
	return(fl)
}


clean_tx_tyler_w2 = function(fl){
	return(fl)
}


clean_srel_contact = function(fl){
	return(fl)
}

clean_tx_tyler_w1 = function(fl){

	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	# Trim the final 23 fixes as many are post-capture fixes.
	numtrim = 23
	fldrop = fl[order(pigID, datetime)][, lapply(.SD, function(x) x[-((length(longitude) - numtrim):(length(longitude)))]), by=pigID]

	# Look at fix time distribution...looks reasonable
	fixtimes = fldrop[, list(diffs=median(diffunits(datetime))), by=pigID]
	return(fldrop)
}

clean_florida = function(fl){

	fl$datetime = as.POSIXct(strptime(fl$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

	# Trim the final 23 fixes as many are post-capture fixes.
	numtrim = 23
	fldrop = fl[order(pigID, datetime)][, lapply(.SD, function(x) x[-((length(longitude) - numtrim):(length(longitude)))]), by=pigID]

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
	  
	 #  if(plotit){
	 #  	dev.new()
		#   plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1, main=bp)
		#   points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
		#   points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
		# }
	  
	}

	rminds = c(fullind[[1]][1], fullind[[2]][1], fullind[[3]][c(1, 6)])
	trimdat = fldrop[-rminds, ]

	return(trimdat)
}

clean_cali0 = function(fl){
	return(fl)
}

clean_cali1 = function(fl){
	return(fl)
}


clean_cali2 = function(fl){
	return(fl)
}

clean_cali3 = function(fl){

	# Only include pig that has 3 months of data
	fldrop = fl[pigID == "caliCC1a"]
	return(fldrop)
}

clean_cali4 = function(fl){
	return(fl)
}

clean_mo_kurt0 = function(fl){
	return(fl)
}

clean_tx_susan = function(fl){

	# Trim the first 10 and last 10 fixes as they include errant movements
	numtrim = 10
	fldrop = fl[order(pigID, datetime)][, lapply(.SD, function(x) x[-c((1:10), ((length(longitude) - numtrim):(length(longitude))))]), by=pigID]

	return(fldrop)
}

clean_tx_tyler_k1 = function(fl){
	return(fl)
}

clean_michigan = function(fl){

	fldrop = fl[latitude > 43.75] #fl[order(pigID, datetime)][, lapply(.SD, function(x) x[-c((1:10), ((length(longitude) - numtrim):(length(longitude))))]), by=pigID]

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
	  
	 #  if(plotit){
	 #  	dev.new()
		#   plot(bpdat[,list(longitude, latitude)], pch=19, col=kmeans.result$cluster, cex=1, main=bp)
		#   points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
		#   points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)
		# }
	  
	}

	return(fldrop)
}


