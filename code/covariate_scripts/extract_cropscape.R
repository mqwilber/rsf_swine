## Script to extract and format the CropScape data
## 
## Script reads in Cropscape data and breaks the data into 10-11 indicator variables
## based on the crop groups defined in cropgroupings.csv
##
## Each of the attributes is then saved as a separate covariate file for each study.
## Note that they are all saved under the covariate "croplayer/" with the format
## study_habitattype_year.grd
##
## Author: Mark Wilber

library(data.table)
library(lubridate)
library(ncdf4)
library(raster)
library(parallel)
source("covariate_fxns.R") # File with useful functions for formatting covariates
base = "/Users/mqwilber/Repos/rsf_swine"

# Formatted summary of the studies
study_sum = fread(file.path(base, "data/formatted/study_summary.csv"))
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)

# Extract 4 crop identifiers from metadata
#		crops, grain_hay_seeds, forest, open_water
csmeta = fread(file.path(base, "data/covariate_data/croplayer/cropgroupings.csv"))

unqgroups = unique(csmeta$group_name)


cropvals = lapply(1:length(unqgroups), function(x) csmeta[group_name == unqgroups[x], value])


# Loop through different studies to format croplayer covariates
for(studynm in study_sum$study){

	if(studynm == "txcamp" | studynm == "michigan"){ # Just process txcamp and michigan for now

		cat("Processing", studynm, "\n")

		ind = study_sum$study == studynm
		minyear = year(study_sum[ind]$datetime_mindate)
		maxyear = year(study_sum[ind]$datetime_maxdate)

		# Loop through different croplayer years
		for(yr in minyear:maxyear){

			tras = raster(file.path(base, "data/covariate_data/croplayer", studynm, 
												paste(studynm, "_croplayer_", yr, "_raw.tif", sep="")))

			extobj = extent(study_sum$longitude_min[ind], study_sum$longitude_max[ind], 
											study_sum$latitude_min[ind], study_sum$latitude_max[ind])
			cras = crop(tras, extobj)

			# Make separate rasters for each crop type
			for(j in 1:length(unqgroups)){

				tempras = cras
				ctype = unqgroups[j]

				# Nothing in non-agricultural land. Don't make raster for this type
				if(ctype != "nothing"){

					ctypeind = (values(cras) %in% cropvals[[j]]) # Indicator TRUE if values match ctype
					values(tempras)[ctypeind] = 1 # ctype habitat
					values(tempras)[!ctypeind] = 0 # Not ctype forest

					# Save formated raster file
					writeRaster(tempras, filename=file.path(base, "data/covariate_data/croplayer", 
												studynm, paste(studynm, "_", ctype, "_", yr, ".grd", sep="")), 
												format="raster", overwrite=TRUE)
				}

			} # End croptype
		} # End if
	} # End year
} # End study

