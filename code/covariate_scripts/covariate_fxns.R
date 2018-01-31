## Functions for manipulating covariates
library(ncdf4)
library(rgdal)
library(sp)
library(plotKML)
library(raster)

ncdf_to_raster = function(nc, timeindex, varname, extobj, latname="lat", 
							lonname="lon", convertLon=TRUE){
	# Given a ncdf object, extract the values for a given attribute
	# at particular times. 
	#
	# Parameters
	# ----------
	# nc : ncdf4 object
	# timeindex : int or vector
	#	An index (or indexes) referring the times to extract for the attribute data
	# varname : str
	# 	The name of the variable to extract
	# extobj: an extent object from the raster package
	#	Used the crop the generated raster objects
	# latname : str
	#	The name of the y/latitude dimension
	# lonname : str
	# 	The name of the x/longitude dimension
	# convertLon : bool
	# 	If TRUE, assumes longitude is measured from 0 - 360 degrees and converts
	#   to -180 to 180.
	#
	# Notes
	# -----
	# Assumes all times are sequential

	v1 = nc$var[[varname]]
	varsize = v1$varsize
	ndims = v1$ndims
	nt = varsize[ndims] # Time steps

	start = rep(1, ndims)
	start[ndims] = timeindex[1] # Start reading at a particular time step
	count = varsize # begin w/count=(nx,ny,nz,...,nt), reads entire var
	count[ndims] = length(timeindex) # Read this many time steps sequentially

	# Extract the appropriate variable data
	vardata = ncvar_get(nc, varname, start=start, count=count)
	# vectorize the data for each year
	varlist = lapply(1:length(timeindex), function(x) as.vector(vardata[, , x]))
	varvect = do.call(rbind, varlist)

	# Extract the matching dimension data
	lat = ncvar_get(nc, latname)
	lon = ncvar_get(nc, lonname)

	if(convertLon)
		lon = ifelse(lon > 180, lon - 360, lon)
	
	lonlat = as.data.frame(expand.grid(lon=lon, lat=lat))
	
	rasterlist = list()

	# Loop through the different time indexes
	for(i in 1:nrow(varvect)){

		out_dat = data.frame(cbind(varvect[i, ], lonlat))
		colnames(out_dat) = c(paste(varname, timeindex[i], sep=""), lonname, 
																		latname)

		# Convert to a raster object
		pnt_dat = SpatialPointsDataFrame(coords = out_dat[,c(lonname,latname)], 
							data=out_dat,
		                    proj4string = CRS(paste("+proj=longlat +datum=WGS84", 
		                    					"+ellps=WGS84 +towgs84=0,0,0")))
		out_rast = vect2rast(pnt_dat, cell.size=.5)
		out_rast = raster(out_rast)
		rasterlist[[i]] = crop(out_rast, extobj)
	}

	return(rasterlist)

}

