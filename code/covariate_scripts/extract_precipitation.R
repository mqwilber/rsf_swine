## Extract monthly precipitation data for the studies used in the analyses 
## 
## Description
## -----------
## This script will extract monthly covariates of precipitation for the region 
## under study

library(data.table)
library(lubridate)
library(pbdNCDF4) # To parallelize the ncdf file reading
library(raster)
library(parallel)
source("covariate_fxns.R")

base = "/Users/mqwilber/Repos/rsf_swine/data/covariate_data/precipitation"
rasfolder = "raster_precip_data"

# Formatted summary of the studies
study_sum = fread("../../data/formatted/study_summary.csv")
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)

# If month_year rasters are not yet built, build them. Otherwise, use what is there.
#sink(file.path(base, "log_ras.txt"))
if(!file.exists(file.path(base, rasfolder))){

	tfp = file.path(base, rasfolder)
	dir.create(tfp, showWarnings = F)

	precipfiles = Sys.glob(file.path(base, "raw_precip", "precip.*.nc"))

	for(i in 1:length(precipfiles)){

		flname = precipfiles[i]
		yr = strsplit(basename(flname), ".", fixed=TRUE)[[1]][2]
		cat("Working on file", basename(flname), "\n")

		tnc = nc_open(flname)

		# Extract specific date values
		hourvals = ncvar_get(tnc, "time")
		secondvals = hourvals * 60 * 60 # Convert to seconds
		datevals = as.POSIXct(secondvals, origin="1900-01-01 00:00:00", tz="GMT")
		timedf = data.table(date=datevals, month=month(datevals), year=year(datevals))

		rasbricks = unstack(brick(flname, varname="precip"))

		# Get mean precip raster for each month
		unqmonths = unique(timedf$month)

		for(month in unqmonths){

			cat(basename(flname), ":", "month", month, "\n")

			tis = which(timedf$month == month)
			raslist = rasbricks[tis]

			# Get total rainfall over the month.
			rasMeans = sum(stack(raslist))

			fname = paste("precip_", month, "_", yr, ".grd", sep="")
			writeRaster(rasMeans, file.path(tfp, fname), format="raster", overwrite=TRUE)

		}
	}

}
# sink()

####################################
### Format precip for each study ###
####################################

for(studynm in unique(study_sum$study)){
	#studynm = "ga_steve"

	print(paste("Extracting precipitation data for", studynm))

	# Range of indices for study
	minmax = study_sum[study == studynm, list(minmonth=month(datetime_mindate), 
									 		  minyear=year(datetime_mindate),
									 	      maxmonth=month(datetime_maxdate),
									 	      maxyear=year(datetime_maxdate))]

	mindate = strptime(paste(minmax$minyear, minmax$minmonth, "01", sep="-"), format="%Y-%m-%d", tz="GMT")
	maxdate = strptime(paste(minmax$maxyear, minmax$maxmonth, "01", sep="-"), format="%Y-%m-%d", tz="GMT")
	dates = seq(mindate, maxdate, by="months")


	ind = study_sum$study == studynm
	cellsize = 0.07
	extobj = extent(c(xmin=study_sum$longitude_min[ind] - cellsize, 
										xman=study_sum$longitude_max[ind] + cellsize,
					 					ymin=study_sum$latitude_min[ind] - cellsize, 
					 					ymax=study_sum$latitude_max[ind] + cellsize))

	pfiles = paste0("precip_", month(dates), "_", year(dates), ".grd", sep="")

	for(j in 1:length(pfiles)){

		fl = pfiles[j]
		tras = raster(file.path(base, rasfolder, fl))
		tras_agg = disaggregate(tras, fact=c(8, 8)) # Dissaggregate data to help cropping
		cras = crop(rotate(tras_agg), extobj)

		tfp = file.path(base, studynm)
		dir.create(tfp, showWarnings = FALSE)

		fnm = paste(studynm, "_precipitation_", month(dates)[j], "_", year(dates)[j], 
												".tif", sep="")
		writeRaster(cras, file.path(tfp, fnm), format="GTiff", overwrite=TRUE)

	}

}	




