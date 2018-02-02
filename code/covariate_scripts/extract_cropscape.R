## Script to extract and format the CropScape data

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

for(studynm in study_sum$study){

	if(studynm == "txcamp"){ # Just process txcamp for now

		tras = raster(file.path("../../data/covariate_data/croplayer/", studynm, 
											paste(studynm, "_croplayer_raw.tif", sep="")))

		ind = study_sum$study == studynm
		extobj = extent(study_sum$longitude_min[ind], study_sum$longitude_max[ind], 
										study_sum$latitude_min[ind], study_sum$latitude_max[ind])

		cras = crop(tras, extobj)
		forestind = (values(cras) == 142) # Forest habitat
		values(cras)[forestind] = 1 # Forest habitat
		values(cras)[!forestind] = 0 # Not forest

		writeRaster(cras, filename=file.path("../../data/covariate_data/croplayer/", 
										studynm, paste(studynm, "_croplayer.grd", sep="")), 
										format="raster", overwrite=TRUE)


	}

}

