## Functions for manipulating covariates
library(pbdNCDF4)
library(rgdal)
library(sp)
library(plotKML)
library(raster)
library(rgeos)

ncdf_to_raster = function(nc, timeindex, varname, extobj, latname="lat", 
							lonname="lon", convertLon=TRUE, cellsize=0.5){
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
	# cellsize : float
	# 	When converting to a raster, specify in the cell size to use (e.g. in degrees)
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
	count[ndims] = length(timeindex) # Read many time steps sequentially

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
		out_rast = vect2rast(pnt_dat, cell.size=cellsize)
		out_rast = raster(out_rast)
		rasterlist[[i]] = crop(out_rast, extobj)
	}

	return(rasterlist)

}

dist_nearest_neighbor = function(ras, shp, center_dist=TRUE){
	# For each cell in the given raster, computes distance to the nearest polygon
	# center in the shapefile
	#
	# Parameters
	# ----------
	# ras : RasterLayer
	# shp : SpatialPolygonsDataFrame 
	# center_dist : bool
	#		If TRUE, computes distance to nearest polygon centers (much faster but less accurate).
	#		If FALSE, computes distance to nearest perimeter point

	polys = shp
	tras = ras
	polys = crop(polys, ras)

	# Polygons are present
	if(dim(polys)[1] != 0){

		# For each cell in raster, compute distance to nearest

		if(center_dist)
			objectpts = gCentroid(polys, byid=TRUE)
		else
			objectpts = as(as(polys, "SpatialLinesDataFrame"), "SpatialPoints")

		pts = rasterToPoints(tras, spatial=TRUE)

		numcalcs = as.numeric(length(pts)) * as.numeric(length(objectpts))
		cat("Distance to nearest calculation will require", 
								numcalcs, "calculations", "\n")

		# If there are too many points for nearest neighbor, subsample 
		if(numcalcs > 7e8){

			numneeded = floor(7e8 / length(pts))
			numsamp = length(objectpts) - numneeded
			cat("Sub-sampling points to 7e8 calculations", "\n")

			set.seed(123)
			objectpts = objectpts[-sample(1:length(objectpts), numsamp, replace=F)]
		}

		# Distance is in kilometers, convert to meters. 
		cat("Calculating distances...", "\n")
		mindists = apply(spDists(pts, objectpts, longlat=TRUE), 1, min)
		values(tras) = mindists*1000

	} else{
		values(tras) = NA
	}

	return(tras)
}

min_distance_buffered = function(pt, sl, buffers=c(0.5, 1, 2, 3, 4, 5),
                                 crs_str="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"){

  # Tries to be smart about computing nearest distance.  
  # Computes a range of rectangular buffers and only computes distance for points
  # in the buffer. For small buffers, this is mush faster than computing the 
  # distance to all points and taking the minimum.
  #
  # If none of the buffers hold any points, return the distance of maximum 
  # buffer for that raster.
  #
  # Parameters
  # ----------
  # pt : 1 X 2 matrix
  #		A point for which to find the minimum distance
  # sl : SpatialPointsDataFrame
  #		All the points against which to compute nearest distance
  # buffers: vector
  #		Vector of sizes buffers to draw around points. In km
  # crs_str : string
  # 	Projection of the raster
  #
  #	Returns
  # -------
  # : distance to nearest point or max buffer.

  xy = pt
  buffersrng = 0.0003258508*(buffers*100 / 3) # Range of buffers
  
  inpoly = array(NA,  dim=length(buffers))
  inds = list()
  
  for(i in 1:length(buffersrng)){
    
    buffer = buffersrng[i]
    xs = c(xmin=xy[, 1] - buffer, xmax=xy[, 1] + buffer)
    ys = c(ymin=xy[, 2] - buffer, ymax=xy[, 2] + buffer)
    ind = point.in.polygon(sl@coords[, 1], sl@coords[, 2], rep(xs, 2), rep(ys, c(2, 2)))
    
    inpoly[i] = any(ind > 0)
    inds[[i]] = ind
  }
  
  if(!any(inpoly)){
    dist = spDistsN1(pt, c(xs[1], ys[1]), longlat = T)
  } else{
    ptsin = sl[inds[[which(inpoly)[1]]] > 0, ]
    dist = min(spDists(pt, ptsin, longlat=T))
  }
  return(dist)
}





