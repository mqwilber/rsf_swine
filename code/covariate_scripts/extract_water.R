#!/usr/bin/env Rscript

# Extract command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

## Creates distance to nearest water rasters for pig studies
##
## Description
## -----------
## Given a shapefile from the NWI:
## 1. Clip the shapefile using gdal and save the minimal shapefile
## 2. Compute the distance 2 nearest water over a raster

library(data.table)
library(sp)
library(raster)
source("covariate_fxns.R")


summarydat = fread("/Users/mqwilber/Repos/rsf_swine/data/formatted/study_summary.csv")
base = "/Users/mqwilber/Repos/rsf_swine/data/covariate_data/water"

## STEP 1: Crop the state-level wetland shape files

# List of NWI shape files referenced by study
nwi_shpfiles = list("tejon"=file.path(base, "downloaded/CA_Wetlands_South.shp"),
										"txcamp"=file.path(base, "downloaded/TX_Wetlands_Central.shp"),
										"fl_raoul"=file.path(base, "downloaded/FL_Wetlands.shp"),
										"tx_tyler_w2"=file.path(base, "downloaded/TX_Wetlands_Central.shp"),
										"srel_contact"=file.path(base, "downloaded/SC_Wetlands.shp"),
										"tx_tyler_w1"=file.path(base, "downloaded/TX_Wetlands_Central.shp"),
										"cali2"=file.path(base, "downloaded/CA_Wetlands_North.shp"),
										"florida"=file.path(base, "downloaded/FL_Wetlands.shp"),
										"cali0"=file.path(base, "downloaded/CA_Wetlands_NorthCentral.shp"),
										"mo_kurt0"=file.path(base, "downloaded/MO_Wetlands.shp"))

for(studynm in summarydat$study){

	if(studynm %in% args) {

		## STEP 1: Crop the state-level wetland shape files

		# NWI data comes as an AEA projection. Re-project the Lat Lon study bounds
		# to crop appropriately
		tdat = summarydat[study == studynm]
		longpoints = rep(c(tdat$longitude_min, tdat$longitude_max), c(2, 2))
		latpoints = rep(c(tdat$latitude_min, tdat$latitude_max), 2)

		ll = data.frame(x=longpoints, y=latpoints)
		spdat = SpatialPoints(ll)
		proj4string(spdat) = "+proj=longlat +datum=WGS84"

		aea_dat = spTransform(spdat, 
								CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))

		buffer = 500
		coords = as.data.frame(aea_dat)
		xmin = min(coords$x) - buffer
		xmax = max(coords$x) + buffer
		ymin = min(coords$y) - buffer
		ymax = max(coords$y) + buffer

		# MUST HAVE GDAL INSTALLED AND ON PATH
		dirpath = file.path(base, studynm)
		dir.create(dirpath, showWarnings = F)
		dir.create(file.path(dirpath, "raw"), showWarnings = F)
		savefl = file.path(dirpath, "raw", paste(studynm, "_water_nwi_raw.shp", sep=""))

		# Only do this if the file doesn't yet exist
		if(!file.exists(savefl)){
			gdal_call = paste("ogr2ogr -overwrite -clipsrc", xmin, ymin, xmax, ymax, 
														savefl, nwi_shpfiles[[studynm]], sep=" ")

			cat("Cropping shapefile for", studynm, "...", "\n")
			system(gdal_call)
		}

		## Step 2: Compute distance to nearest water polygon using GDAL cropped data

		# Load cropped shape file and transform to Lat Long.
		cshp = shapefile(savefl)
		cshp_trans = spTransform(cshp, CRS("+proj=longlat +datum=WGS84"))

		# Only consider permanent and semi-permanent water sources: "H" and "F" in NWI lingo: https://www.fws.gov/wetlands/Data/Wetland-Codes.html
		peren = cshp_trans[grepl("H", cshp_trans$ATTRIBUTE) | grepl("F", cshp_trans$ATTRIBUTE), ]

		# Create a dummy raster
		buffer = 0.005
		ras = raster()
		extent(ras) = extent(c(xmin=min(ll$x) - buffer,
													 xmax=max(ll$x) + buffer,
													 ymin=min(ll$y) - buffer,
													 ymax=max(ll$y) + buffer))

		projection(ras) = "+proj=longlat +datum=WGS84"
		res(ras) = c(0.0003258969, 0.0003258508)

		cat("Computing distance to nearest water for", studynm, "\n")

		tif_path = file.path(dirpath, paste(studynm, "_water_nndistance.tif", sep=""))

		#if(!file.exists(tif_path)){
			distras = dist_nearest_neighbor(ras, peren)
			writeRaster(distras, tif_path, format="GTiff", overwrite=TRUE)
		#}
	}
}
