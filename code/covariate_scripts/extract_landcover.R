#!/usr/bin/env Rscript

# Extract command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

## Script to extract and format NLCD landcover data
## -----------------------------------------------
## 
## The script reads in the National landcover data and and crops and projects the
## the data for each pig study (excluding `canada`).
##
## Generates a the following types of forest cover
## 	1. Deciduous forest
##	2. Evergreen forest
##  3. Mixed forest
##  4. Other trees
##
## Author: Mark Wilber

library(data.table)
library(lubridate)
library(raster)
library(parallel)
source("covariate_fxns.R") # File with useful functions for formatting covariates
base = "/Users/mqwilber/Repos/rsf_swine"
crsbase = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Formatted summary of the studies
study_sum = fread(file.path(base, "data/formatted/study_summary.csv"))
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)

# Nationwide NLCD data from 2011
nlcd = raster(file.path(base, "data/covariate_data/landcover/raw_nlcd", 
													"nlcd_2006_landcover_2011_edition_2014_10_10.img"))

nlcdTrans = projectExtent(nlcd, crsbase)

# Nationwide NLCD percent canopy layer from 2011.
cover = raster(file.path(base, "data/covariate_data/landcover/raw_canopycover", 
													"nlcd2011_usfs_conus_canopy_analytical.img"))
coverTrans = projectExtent(cover, crsbase)


# NLCD metadata specifies groupings
nlcdmeta = fread(file.path(base, "data/covariate_data/landcover/raw_nlcd/", 
														"nlcd_landcover_groupings.csv"))

# Fine-tuning NLCD groupings
#		- Considering all wetlands and developed types as identical
nlcdmeta$grouping[nlcdmeta$grouping %like% "wetland"] = "wetland"
nlcdmeta$grouping[nlcdmeta$grouping %like% "developed"] = "developed"


unqgroups = unique(nlcdmeta$grouping)

# Loop through different studies to format landcover and canopy cover covariates
for(studynm in study_sum$study){

	# Select a few studies to process
	if(any(studynm %in% args)) { #"txcamp", "tejon", "tx_tyler_w2","fl_raoul", "srel_contact"))) {

		cat("Processing", studynm, "\n")

		ftp = file.path(base, "data/covariate_data/landcover", studynm)
		dir.create(ftp, showWarnings = FALSE)

		ind = study_sum$study == studynm

		buffer = 0.005
		extobj = extent(c(xmin=study_sum$longitude_min[ind] - buffer, 
											xmax=study_sum$longitude_max[ind] + buffer, 
											ymin=study_sum$latitude_min[ind] - buffer, 
											ymax=study_sum$latitude_max[ind] + buffer))

		# Crop and project the landcover groupings
		cras_extent = projectExtent(crop(nlcdTrans, extobj), crs(nlcd))
		cras = projectRaster(crop(nlcd, extent(cras_extent)), crs=crsbase, method="ngb")

		# Make separate rasters for each landcover type
		for(j in 1:length(unqgroups)){

			tempras = cras
			ctype = unqgroups[j]
			landvals = nlcdmeta[grouping == ctype]$values

			ctypeind = (values(cras) %in% landvals) # Indicator TRUE if values match ctype
			values(tempras)[ctypeind] = 1 # ctype habitat
			values(tempras)[!ctypeind] = 0 # Not ctype forest

			# Aggregate raster for some studies to reduce computational time
			if(studynm == "cali0" & ctype == "developed"){
				tempras = raster::aggregate(tempras, fact=4)
			}

			# Save formated raster file
			fname = file.path(ftp, paste(studynm, "_", ctype, ".tif", sep=""))
			writeRaster(tempras, filename=fname, format="GTiff", overwrite=TRUE)

			if(ctype %in% c("developed", "barren_land")){ #"wetland", "developed", 

				fnameshp = paste(strsplit(fname, ".", fixed=TRUE)[[1]][1], ".shp", sep="")
				system(paste("gdal_polygonize.py ", fname, " -f 'ESRI Shapefile' ", fnameshp, sep=""))

				# Dummy raster
				buffer = 0.005
				dumras = raster()
				extent(dumras) = extent(cras)
				projection(dumras) = projection(cras)
				res(dumras) = res(cras)

				shp = shapefile(fnameshp)
				distras = dist_nearest_neighbor(dumras, shp)
				writeRaster(distras, filename=file.path(ftp, 
									paste(studynm, "_", ctype, "_nndistance.tif", sep="")), 
									format="GTiff", overwrite=TRUE)
			}


		}

		# Crop and project the tree cover density
		lcras_extent = projectExtent(crop(coverTrans, extobj), crs(cover))
		lcras = projectRaster(crop(cover, extent(lcras_extent)), crs=crsbase, method="ngb")

		writeRaster(lcras, filename=file.path(ftp, paste(studynm, "_canopycover.tif", sep="")), 
										format="GTiff", overwrite=TRUE)

	} # End if
} # End study

