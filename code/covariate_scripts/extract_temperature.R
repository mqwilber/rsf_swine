## Extract monthly temperature data for the studies used in the analyses 
## 
## Description
## -----------
## This script will extract monthly covariates of temperature for the region 
## under study

library(data.table)
library(lubridate)
library(ncdf4)
library(raster)
library(parallel)
source("covariate_fxns.R")

# Formatted summary of the studies
study_sum = fread("../../data/formatted/study_summary.csv")
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)

# Load in netcdf temperature file
# nc = nc_open("../../data/covariate_data/temperature/air.mon.mean.nc")
bftemp = unstack(brick("../../data/covariate_data/temperature/air.mon.mean.nc", 
											varname="air"))
		
startdate = as.POSIXct("1948-01-01 GMT", tz="GMT")
enddate = as.POSIXct("2017-12-01 GMT", tz="GMT")
nmdates = seq(startdate, enddate, by="months")


for(studynm in unique(study_sum$study)){
	#studynm = "ga_steve"

	print(paste("Extracting temperature data for", studynm))

	# Range of indices for study
	minmax = study_sum[study == studynm, list(minmonth=month(datetime_mindate), 
									 		  minyear=year(datetime_mindate),
									 	      maxmonth=month(datetime_maxdate),
									 	      maxyear=year(datetime_maxdate))]

	mindate = strptime(paste(minmax$minyear, minmax$minmonth, "01", sep="-"), format="%Y-%m-%d", tz="GMT")
	maxdate = strptime(paste(minmax$maxyear, minmax$maxmonth, "01", sep="-"), format="%Y-%m-%d", tz="GMT")
	dates = seq(mindate, maxdate, by="months")

	tempras = bftemp[nmdates %in% dates]

	ind = study_sum$study == studynm

	cellsize = 0.07
	extobj = extent(c(xmin=study_sum$longitude_min[ind] - cellsize, 
										xmax=study_sum$longitude_max[ind] + cellsize,
					 					ymin=study_sum$latitude_min[ind] - cellsize, 
					 					ymax=study_sum$latitude_max[ind] + cellsize))

	for(i in 1:length(dates)){

		tras = tempras[[i]]
		tras_agg = disaggregate(tras, fact=c(8, 8)) # Dissaggregate data to help cropping
		cras = crop(rotate(tras_agg), extobj)

		tfp = file.path("../../data/covariate_data/temperature", studynm)
		dir.create(tfp, showWarnings = FALSE)

		fnm = paste(studynm, "_temperature_", month(dates)[i], "_", year(dates)[i], 
												".tif", sep="")
		writeRaster(cras, file.path(tfp, fnm), format="GTiff", overwrite=TRUE)
	}

}	

#mclapply(unique(study_sum$study), raster_convert, nc, study_sum, timedf, mc.cores=4)




