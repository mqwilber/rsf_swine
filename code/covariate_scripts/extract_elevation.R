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

crs_string = "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"

for(studynm in unique(study_sum$study)){

	if(studynm %in% c('tejon', 'txcamp')) {

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
			crs(tras) = crs_string
			extent(tras) = extent(bbox)
			ras_list[[metaroot]] = tras
		}

		# If necessary merge the elevation rasters.  Only can handle two rasters!
		if(length(ras_list) > 1){
			ras_list = lapply(ras_list, function(x) {
																		origin(x) = origin(ras_list[[1]])
																		return(x)})
			fullras = merge(ras_list[[1]], ras_list[[2]])
		} else{
			fullras = ras_list[[1]]
		}

		# Save the formatted elevation raster
		raspath = file.path("../../data/covariate_data/elevation", studynm, 
															paste(studynm, "_elevation.tif", sep=""))
		writeRaster(fullras, raspath, formt="GTiff", overwrite=TRUE)
	}
}
