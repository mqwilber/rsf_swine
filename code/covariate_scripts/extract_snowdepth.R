## Extract monthly snow depth data for the studies used in the analyses 
## 
## Description
## -----------
## This script will extract monthly covariates of snowdepth for the region 
## under study

library(data.table)
library(lubridate)
library(ncdf4)
library(raster)
library(parallel)
source("covariate_fxns.R")
snodasdir = "../../data/covariate_data/snowdepth"

# Formatted summary of the studies
study_sum = fread("../../data/formatted/study_summary.csv")
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)


for(studynm in paste0("la_steve", 0:5)){ #unique(study_sum$study)) {

	if(studynm != "canada"){
	#studynm = "ga_steve"

		print(paste("Extracting snowdepth data for", studynm))

		# Range of indices for study
		minmax = study_sum[study == studynm, list(minmonth=month(datetime_mindate), 
										 		  										minyear=year(datetime_mindate),
										 	      									maxmonth=month(datetime_maxdate),
										 	      									maxyear=year(datetime_maxdate))]

		mindate = strptime(paste(minmax$minyear, minmax$minmonth, "01", sep="-"), 
																										format="%Y-%m-%d", tz="GMT")
		maxdate = strptime(paste(minmax$maxyear, minmax$maxmonth, "01", sep="-"), 
																										format="%Y-%m-%d", tz="GMT")

		dates = seq(mindate, maxdate, by="months")
		ind = study_sum$study == studynm

		cellsize = 0.02
		extobj = extent(c(xmin=study_sum$longitude_min[ind] - cellsize, 
											xmax=study_sum$longitude_max[ind] + cellsize,
						 					ymin=study_sum$latitude_min[ind] - cellsize, 
						 					ymax=study_sum$latitude_max[ind] + cellsize))

		for(i in 1:length(dates)){
			# canada_temperature_3_2015.tif

			# Load in the appropriate raster
			tyear = year(dates)[i]
			tmonth = month(dates)[i]
			mnthdir = ifelse(tmonth < 10, paste("0", tmonth, sep=""), tmonth)
			studydir = file.path(snodasdir, "downloaded", tyear, mnthdir)
			tras = raster(Sys.glob(file.path(studydir, "*.tif")))

			# Crop data
			cras = crop(tras, extobj)

			# Save cropped raster
			tfp = file.path(snodasdir, studynm)
			dir.create(tfp, showWarnings = FALSE)
			fnm = paste(studynm, "_snowdepth_", tmonth, "_", tyear, ".tif", sep="")
			writeRaster(cras, file.path(tfp, fnm), format="GTiff", overwrite=TRUE)
		}
	}

}	




