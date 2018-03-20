# Make and save an extent shapefiles for each study.
#
# Shape files are saved in the ndvi/downloaded covariate folder
#
# Author: Mark Wilber

library(raster)
library(data.table)

base = "../data/covariate_data/ndvi/downloaded"

studysum = fread("../data/formatted/study_summary.csv")

for(studynm in studysum$study){

	# Make extent object
	tstud = studysum[study == studynm]
	buffer = 0.007
	extobj = extent(c(xmin=tstud$longitude_min - buffer, xmax=tstud$longitude_max + buffer,
										ymin=tstud$latitude_min - buffer, ymax=tstud$latitude_max + buffer))

	sp = as(extobj, "SpatialPolygons")
	crs(sp) = "+proj=longlat +datum=WGS84"

	# Create directory and save shape file
	fp = file.path(base, studynm)
	dir.create(fp)
	shapefile(sp, file.path(fp, paste(studynm, "_extent.shp", sep="")), overwrite=TRUE)

}