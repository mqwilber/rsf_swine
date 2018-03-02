## Script to format and extract the masting layers used in analysis

library(raster)
library(data.table)
library(lubridate)


# Formatted summary of the studies
study_sum = fread("../../data/formatted/study_summary.csv")
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)

# Format the masting layer
dens = raster("../../data/covariate_data/masting/NA_density_raster.txt")
crs(dens) = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
densproj = projectRaster(dens, crs="+proj=longlat +datum=WGS84 +ellps=WGS84")

# The projections are weird for the species richness...
# spp = raster("../../data/covariate_data/masting/NA_spp_rich_raster.txt")
# spp = crop(spp, extent(dens))
# crs(spp) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
# sppproj = projectRaster(spp, crs="+proj=longlat +datum=WGS84 +ellps=WGS84")

# Loop through studies
for(studynm in study_sum$study){

	cat("Working on", studynm, "\n")
	ind = study == studynm
	minlon = study_sum$longitude_min[ind]
	maxlon = study_sum$longitude_max[ind]
	minlat = study_sum$latitude_min[ind]
	maxlat = study_sum$latitude_max[ind]

	buffer = 0.02
	extobj = extent(c(xmin=minlon - buffer, xmax=maxlon + buffer, 
										ymin=minlat - buffer, ymax=maxlat + buffer))

	tras = crop(densproj, extobj)
	tfp = file.path("../../data/covariate_data/masting", studynm) 
	dir.create(tfp, showWarnings=FALSE)

	rasname = paste(studynm, "_masting.grd", sep="")
	writeRaster(tras, file.path(tfp, rasname), format="raster", overwrite=TRUE)

}