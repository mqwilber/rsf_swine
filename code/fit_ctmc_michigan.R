library(lubridate)
library(data.table)
library(ggplot2)
library(ggmap)
library(parallel)
source("pigfxns.R")

# Load in the full data
dat = fread("../data/formatted/full_pig_data.csv")
mdat = dat[study == "michigan"]
mdat$datetime = as.POSIXct(mdat$datetime, format="%Y-%m-%d %H:%M:%S", tz="GMT")


# Plot the Michigan data
us = c(left=min(mdat$longitude), bottom=min(mdat$latitude), 
			 right=max(mdat$longitude), top = max(mdat$latitude))
map = get_map(location=geocode("Michigan"), zoom=7)

front = 200
back = 300
tdat = mdat[, list(longitude=longitude[front:(length(longitude) - back)],
									 latitude=latitude[front:(length(latitude) - back)],
									 datetime=datetime[front:(length(datetime) - back)]), by=pigID]
plt = ggmap(map) + geom_point(data=tdat, aes(x=longitude, y=latitude), color="red")
plt

# Check the datetime ranged
diffunit = function(x){
	dt = diff(x)
	units(dt) = "mins"
	return(dt)
}

tdat[, list(mindate=min(datetime), maxdate=max(datetime), 
								meddelta=median(diffunit(datetime)),
								mindelta=min(diffunit(datetime)),
								maxdelta=max(diffunit(datetime)),
								numfixes=length(datetime)), by=pigID]

# A few pigs with some longer stretches...but the median is all at 30 mins.

# From the map, there seem to be two or three outliers
badpigs = unique(tdat$pigID)

badout = array(NA, dim=3)
fullind = array(NA, dim=3)

for(i in 1:length(badpigs)){
  
  bp = badpigs[i]
  bpdat = tdat[pigID == bp, list(longitude, latitude)]
  kmeans.result = kmeans(bpdat, 1)
  centers = kmeans.result$centers[kmeans.result$cluster, ]
  distances <- sqrt(rowSums((bpdat - centers)^2))
  
  outliers <- order(distances, decreasing=T)[1:5]

  print(outliers) 
  badout[i] = outliers[1]
  fullind[i] = which((tdat$pigID == bp) & (tdat$longitude == bpdat$longitude[outliers[1]]) & (tdat$latitude == bpdat$latitude[outliers[1]]))
  
  dev.new()
  plot(bpdat[,list(longitude, latitude)], type="b", pch=19, col=kmeans.result$cluster, cex=1, main=bp)
  points(kmeans.result$centers[, c("longitude", "latitude")], col=1:3, pch=15, cex=2)
  points(bpdat[outliers, c("longitude", "latitude")], pch="+", col=4, cex=3)

  
}


# Clearly, point 13280 has some problems...drop it
cleandat = tdat[-fullind[3]][!(pigID == "michigan36635" & latitude < 44.2 & longitude < -84.4)]
plt = ggmap(map) + geom_point(data=cleandat[pigID == "michigan36635"][-c((length(datetime) - 10):length(datetime))], aes(x=longitude, y=latitude), color="red")
plt

with(cleandat[order(datetime)][pigID == "michigan36635"], plot(longitude, latitude, type="b"))


# Looks pretty good.  TODO: Need to consistently remove capture effects.
# vectorField



ppigs = function(i, allpigs, numpoints, buffer, timestep, studynm){
  
  unqpig = unique(allpigs$pigID)
  
  cat("Working on pig", i, "of", length(unqpig), ":", unqpig[i], "\n")
        
  pigdata = allpigs[pigID == unqpig[i], ]
  
  extobj = extent(min(pigdata$longitude) - buffer, max(pigdata$longitude) + buffer, 
                  min(pigdata$latitude) - buffer, max(pigdata$latitude) + buffer)
  
  loc_stack = process_covariates(c("temperature", "croplayer"), studynm, extobj, 
                                 min(pigdata$datetime), max(pigdata$datetime), ext="loc")
  grad_stack = process_covariates(c("temperature"), studynm, extobj, 
                                 min(pigdata$datetime), max(pigdata$datetime), ext="grad")
  gradxy_stack = process_covariates(c("croplayer"), studynm, extobj, 
                                 min(pigdata$datetime), max(pigdata$datetime), ext="grad",
                                 distgrad=TRUE)
  
  gradx_stack = lapply(gradxy_stack, function(m) m$xgrad)
  grady_stack = lapply(gradxy_stack, function(m) m$ygrad)
  
  # Ensure projections are the same
  locnames = names(loc_stack)
  yr = year(min(pigdata$datetime))
  frname = paste('forest_', yr, "_loc", sep="")
  loc_stackproj = stack(lapply(loc_stack, function(x) projectRaster(x, loc_stack[[frname]],
                                                                                method="ngb")))
  grad_stackproj = stack(lapply(grad_stack, function(x) projectRaster(x, loc_stack[[frname]],
                                                                                method="bilinear")))
  gradx_stackproj = stack(lapply(gradx_stack, function(x) projectRaster(x, loc_stack[[frname]],
                                                                                method="ngb")))
  grady_stackproj = stack(lapply(grady_stack, function(x) projectRaster(x, loc_stack[[frname]],
                                                                                method="ngb")))
  
  pigdata = pigdata[, list(x=longitude, y=latitude, datetime=datetime)]
  tglmdat = fit_ctmcmodel(pigdata[1:numpoints, ], unqpig[i], loc_stackproj, grad_stackproj, 
                          method="interp", impute=1, timestep=timestep,
                          buffer=0.002, mc.cores=1, path2ctmcMethod="ShortestPath",
                          xygrad=TRUE, xgrad.stack = gradx_stackproj, ygrad.stack = grady_stackproj)
  
  cat("Done with pig", i, "of", length(unqpig), ":", unqpig[i], "\n")
  return(tglmdat)
}

unqpig = unique(cleandat$pigID)
#sink("log_ctmc.txt")
allglmdat = lapply(1:1, ppigs, cleandat, 100, 0.002, "15 mins", "michigan")
#sink()


















