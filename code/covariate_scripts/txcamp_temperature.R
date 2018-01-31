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
source("covariate_scripts/covariate_fxns.R")

# Formatted summary of the studies
study_sum = fread("../data/formatted/study_summary.csv")
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)


# Load in netcdf temperature file
nc = nc_open("../data/covariate_data/temperature/air.mon.mean.nc")

# Time units: hours since hours since 1800-01-01 00:00:0.0
hourvals = ncvar_get(nc, "time")
secondvals = hourvals * 60 * 60 # Convert to seconds
datevals = as.POSIXct(secondvals, origin="1800-01-01 00:00:00", tz="GMT")
timedf = data.table(date=datevals, month=month(datevals), year=year(datevals))

studynm = "txcamp"

# Range of indices for study
minmax = study_sum[study == studynm, list(minmonth=month(datetime_mindate), 
								 		  minyear=year(datetime_mindate),
								 	      maxmonth=month(datetime_maxdate),
								 	      maxyear=year(datetime_maxdate))]

#  Find the time-dependent indexes
timeinds = which(((timedf$month >= minmax$minmonth) & 
				  (timedf$month <= minmax$maxmonth)) &
				 ((timedf$year >= minmax$minyear) & 
				  (timedf$year <= minmax$maxyear)))

ind = study_sum$study == studynm

extobj = extent(study_sum$longitude_min[ind], study_sum$longitude_max[ind],
				 study_sum$latitude_min[ind], study_sum$latitude_max[ind])

res = ncdf_to_raster(nc, timeinds, "air", extobj)

# TODO: write to files



