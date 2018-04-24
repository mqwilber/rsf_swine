## Compute home range estimates for all pigs

library(raster)
library(sp)
library(data.table)
library(adehabitatHR)
library(lubridate)
library(ggmap)

mcp_calc = function(longitude, latitude){
	# Given lat lon, compute MCP home range

	crsalb = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
	crslatlon = "+proj=longlat +datum=WGS84 +ellps=WGS84"
	spdat = SpatialPoints(data.frame(longitude, latitude), 
												proj4string=CRS(crslatlon))
	sptrans = spTransform(spdat, CRS(crsalb))
	mcpres = mcp(sptrans, percent=95, unout="km2")
	return(mcpres$area)
}

# Load in the  full data
allpigdt = fread("../data/formatted/full_pig_data.csv")
pigattributes = fread("../data/formatted/pig_attributes.csv")

# This conversion is a bit slow...
allpigdt[, month:=month(datetime)]

# Drop canada and judas pigs
allpigdt_nc = allpigdt[!(study %in% c("canada"))]

# Drop all pigs with less 30 fixes
shortpigs = allpigdt_nc[, list(short=length(datetime) < 30), by=pigID]
goodpigs = shortpigs[short == FALSE, pigID]
pigdt = allpigdt_nc[pigID %in% goodpigs]

# For all pigs drop the first and last 50 observations
pigdt = pigdt[order(pigID, datetime)]
numtrim = 20
trimpig = pigdt[, lapply(.SD, function(x) x[-c(1:numtrim, (length(latitude) - numtrim):length(latitude))]), by=pigID]

shortpigs = trimpig[, list(short=length(datetime) < 10), by=pigID]
goodpigs = shortpigs[short == FALSE, pigID]
trimpig = trimpig[pigID %in% goodpigs]


# Calculate MCP for all pigs across all studies
mcpdt = trimpig[, list(MCPareakm2=mcp_calc(longitude, latitude), 
										 total_months=length(unique(month)),
										 mean_lat=round(mean(latitude), 6),
										 mean_lon=round(mean(longitude), 6)), by=list(study, pigID)]
# Join sex data
mcpdt = merge(mcpdt, pigattributes[, list(pigID, sex)], by="pigID")
fwrite(mcpdt, "../results/mcp_est.csv", row.names = F)

# FAST PYTHON GEO CODER
# import pandas as pd
# import reverse_geocoder as rg

# dt = pd.read_csv("../results/mcp_est.csv")
# coords = [(lat, lon) for lat, lon in zip(dt.mean_lat, dt.mean_lon)]
# geocodes = rg.search(coords)  
# states = [g['admin1'] for g in geocodes]    
# dt = dt.assign(state=states)

# # Drop some of the judas pig data
# dtdrop = dt[~(dt.state.isin(["Oregon", "South Carolina", 
# 													 "Kansas", "Louisiana", 
# 													 "New Mexico"]) & (dt.study == "judas_pig"))]

# dtdrop.to_csv("../results/mcp_est.csv")
