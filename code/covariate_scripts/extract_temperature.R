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
nc = nc_open("../../data/covariate_data/temperature/air.mon.mean.nc")

# Time units: hours since hours since 1800-01-01 00:00:0.0
hourvals = ncvar_get(nc, "time")
secondvals = hourvals * 60 * 60 # Convert to seconds
datevals = as.POSIXct(secondvals, origin="1800-01-01 00:00:00", tz="GMT")
timedf = data.table(date=datevals, month=month(datevals), year=year(datevals))

#raster_convert = function(studynm, nc, study_sum, timedf){

for(studynm in unique(study_sum$study)){
	#studynm = "ga_steve"


	print(paste("Extracting temperature data for", studynm))

	# Range of indices for study
	minmax = study_sum[study == studynm, list(minmonth=month(datetime_mindate), 
									 		  minyear=year(datetime_mindate),
									 	      maxmonth=month(datetime_maxdate),
									 	      maxyear=year(datetime_maxdate))]

	#  Find the time-dependent indexes
	timeinds = which(((timedf$year >= minmax$minyear) & (timedf$year <= minmax$maxyear)) & 
					  !((timedf$year == minmax$minyear) & (timedf$month < minmax$minmonth)) &
					  !((timedf$year == minmax$maxyear) & (timedf$month > minmax$maxmonth)))

	print(nrow(timedf[timeinds, ]))

	ind = study_sum$study == studynm

	cellsize = 0
	extobj = extent(study_sum$longitude_min[ind] - cellsize, study_sum$longitude_max[ind] + cellsize,
					 study_sum$latitude_min[ind] - cellsize, study_sum$latitude_max[ind] + cellsize)

	res = ncdf_to_raster(nc, timeinds, "air", extobj)

	# Write temperature rasters
	tfp = file.path("../../data/covariate_data/temperature", studynm)
	dir.create(tfp, showWarnings = F)


	for(i in 1:length(res)){

		monthyear = paste(timedf[timeinds[i], month], timedf[timeinds[i], year], sep="_")
		fname = paste(studynm, "_temperature_", monthyear, ".grd", sep="")
		writeRaster(res[[i]], file.path(tfp, fname), format="raster", overwrite=TRUE)

	}
}	

#mclapply(unique(study_sum$study), raster_convert, nc, study_sum, timedf, mc.cores=4)




