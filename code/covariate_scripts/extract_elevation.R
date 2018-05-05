#!/usr/bin/env Rscript

# Extract command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

## Extracts elevation rasters from the USGS NED data downloaded for each
## particular study.  If necessary, merges rasters.

library(raster)
library(data.table)
library(lubridate)
library(xml2)
library(plyr)

# Formatted summary of the studies
study_sum = fread("../../data/formatted/study_summary.csv")
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)

ned_crs_string = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"
crs_string = "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"

for(studynm in unique(study_sum$study)){

	if(studynm %in% args) {

		cat("Processing elevation for", studynm, "\n")

		# Check for image file
		imgfls = Sys.glob(file.path("../../data/covariate_data/elevation", 
																										studynm, "raw", "*.img"))
		ras_list = list()
		if(length(imgfls) != 0){

			# Load all image rasters
			for(img in imgfls){
				tras = raster(img)

				# Raster is by default in ned_crs_string
				cras = projectRaster(tras, crs=crs_string)
				ras_list[[img]] = cras
			}

		} else{

			# If no image file, there should be jpgs with metadata
			jpgfls = Sys.glob(file.path("../../data/covariate_data/elevation", 
																										studynm, "raw", "*.jpg"))

			ras_list = list()
			for(jpg in jpgfls){

				# Extract associated metadata from XML files
				metaroot = substring(basename(jpg), 4, 12)
				fp = file.path("../../data/covariate_data/elevation", studynm, "raw",
												paste(metaroot, "_meta.xml", sep=""))

				xmldoc = read_xml(fp)
				spdom = xml_children(xml_children(xml_children(xmldoc)[[1]])[[5]])[[1]]

				bbox = list()
				for(child in xml_children(spdom)){
					bbox[[xml_name(child)]] = as.numeric(xml_text(child))
				}

				# Rename extent names
				newname = c("westbc"="xmin",  "eastbc"="xmax", 
										"southbc"="ymin", "northbc"="ymax")
				names(bbox) = revalue(names(bbox), newname)
				bbox = unlist(bbox)[c("xmin", "xmax", "ymin", "ymax")]

				# Load raster and set spatial extent
				tras = raster(jpg)
				crs(tras) = ned_crs_string
				extent(tras) = extent(bbox)
				pras = projectRaster(tras, crs=crs_string)
				ras_list[[metaroot]] = pras
			}
		} 

		# If necessary merge the elevation rasters.  Only can handle two rasters!
		if(length(ras_list) > 1){
			ras_list = lapply(ras_list, function(x) {
																		origin(x) = origin(ras_list[[1]])
																		return(x)})
			fullras = merge(ras_list[[1]], ras_list[[2]])

			# Special case of 4 rasters
			if(length(ras_list) == 4){

				fullras2 = merge(ras_list[[3]], ras_list[[4]])
				fullras = merge(fullras, fullras2)

			}
		} else{
			fullras = ras_list[[1]]
		}

		# Save the formatted elevation raster
		raspath = file.path("../../data/covariate_data/elevation", studynm, 
															paste(studynm, "_elevation.tif", sep=""))
		writeRaster(fullras, raspath, format="GTiff", overwrite=TRUE)
	}
}
