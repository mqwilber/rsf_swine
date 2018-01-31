## Functions for manipulating covariates
library(ncdf4)
library(rgdal)
library(sp)
library(plotKML)
library(raster)

raster_from_ncdf = function(nc, timeindex, varname, latname="lat", lonname="lon"){
	# Given a ncdf object, extract the values for a given attribute
	# at a particular time
	#
	# Parameters
	# ----------
	# nc : ncdf4 object
	# timeindex : int
	#	An index specifying the time dimension to extract
	# varname : str
	# 	The name of the variable to extract

	v1 = nc$var[[varname]]
	varsize = v1$varsize
	ndims = v1$ndims
	nt = varsize[ndims] # Time steps

	start = rep(1, ndims)
	start[ndims] = timeindex # Read only the a particular time step
	count = varsize # begin w/count=(nx,ny,nz,...,nt), reads entire var
	count[ndims] = 1 # Change to only read the a single time step

	# Extract the appropriate variable data
	vardata = as.vector(ncvar_get(nc, varname, start=start, count=count))

	# Extract the matching dimension data
	lat = ncvar_get(nc, latname)
	lon = ncvar_get(nc,lonname)
	lonlat = as.matrix(expand.grid(lon=lon, lat=lat))
	out_dat = data.frame(cbind(vardata, lonlat))
	colnames(out_dat) = c(varname, lonname, latname)

	# Convert to a raster object
	pnt_dat = SpatialPointsDataFrame(coords = out_dat[,c(lonname,latname)], data=out_dat,
	                           proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
	out_rast = vect2rast(pnt_dat, cell.size=.5)
	out_rast = raster(out_rast)

	return(out_rast)

}

