## Extract and format the NDVI data for each study

library(data.table)
library(lubridate)
library(ncdf4)
library(raster)
library(parallel)
source("covariate_fxns.R")

study_sum = fread("../../data/formatted/study_summary.csv")
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)

# Perform this for each study. Just txcamp now

studynm = "txcamp"
sind = study_sum$study == studynm

# Step 1: Parse the filenames to extract dates
flnames = Sys.glob(file.path("../../data/covariate_data/ndvi/downloaded", studynm, "*.tif"))

month1 = sapply(strsplit(basename(flnames), "_"), function(x) month(strptime(x[3], format="%Y%m%d")))
month2 = sapply(strsplit(basename(flnames), "_"), function(x) month(strptime(x[4], format="%Y%m%d")))

# Step 2: Compile the biweekly rasters into monthly rasters. To do this, take the
# mean NDVI for all rasters that overlap with a given month. 
minmonth = month(study_sum[sind]$datetime_mindate)
minyear = year(study_sum[sind]$datetime_mindate)
maxmonth = month(study_sum[sind]$datetime_maxdate)
maxyear = year(study_sum[sind]$datetime_maxdate)

# monthndvi = list()

for(yr in minyear:maxyear){

	if(minyear == maxyear)
		range = minmonth:maxmonth
	else if(yr == minyear)
		range = minmonth:12
	else if(yr == maxyear)
		range = 1:maxmont
	else
		range = 1:12

	for(mn in range){
		ind = (mn == month1) | (mn == month2)

		if(sum(ind) != 0){

			# Monthly average of NDVI for biweekly NDVI measures
			rasstack = stack(flnames[ind])
			raslay = mean(rasstack)
			rascrop = crop(raslay, extent(study_sum[sind]$longitude_min, study_sum[sind]$longitude_max, 
																		study_sum[sind]$latitude_min, study_sum[sind]$latitude_max))
			# monthndvi[[mn]] = rascrop

			fname = file.path("../../data/covariate_data/ndvi/", studynm, paste(studynm, "_ndvi", "_", mn, "_", yr, ".grd", sep=""))
			writeRaster(rascrop, fname, format="raster", overwrite=TRUE)

		}
	}
}

