## This script is performing the following preliminary movement analysis
## on the Texas pig data. 
##
##  Steps
## -------
## 1. Extract 5 pigs with equally spaced (or roughly so) time steps
## 2. Format the datetime so that it will work the packages
## 3. Extract relevant covariates from raster and shp files
##     a. Cultivated land
##     b. Forest cover
##     c. Distance to water
##     d. Distance to nearest road
##     e. Distance to neighboring pigs (intraspecific interactions)
## 4. Run a discrete time movement model analysis and a continuous time movement
##    model analysis to understand how environmental characteristics affect pig
##    movement.


# Load in necessary packages
library(data.table)
library(raster)
library(momentuHMM)
library(ggplot2)
library(lubridate)
library(crawl)
library(sp)
library(rgdal)

## Useful functions
source('/Users/mqwilber/Repos/rsf_swine/code/pigfxns.R')


##### Step 1 #######

dt = fread("/Users/mqwilber/Repos/rsf_swine/data/formatted/full_pig_data_cleaned.csv")

# Extract texas pigs
txpigs = dt[study == "txcamp"]
rm(dt)

# Convert datetime
txpigs$datetimef = as.POSIXct(strptime(txpigs$datetime, format='%Y-%m-%d %H:%M:%S'), tz="GMT")

# Extract pigs with the most observations
cntpigs = txpigs[, list(counts=length(datetime)), by=list(pigID)][order(counts, decreasing = T)]
head(cntpigs)


# Look for overlap
runsdt = txpigs[, list(runs=runs(datetimef, clength=500, ctime=30, ctimemin=0)), by="pigID"]
runsIDs = runsdt$pigID[runsdt$runs]

# extract the specific IDs
ind = array(FALSE, dim=nrow(txpigs))
for(ID in runsIDs){
	ind = ind | txpigs$pigID == ID
}

pigdt = txpigs[ind]
dim(pigdt)
unique(pigdt$pigID)

# Look for where the pigs are overlapping in time
minmax = pigdt[, list(min=min(datetimef), max=max(datetimef)), by="pigID"]
minmax

# Let's just focus on 03-03 for the following analysis
redpigs = minmax[month(min) == 3 & day(min) == 3]$pigID
redpigs

ind = array(FALSE, dim=nrow(pigdt))
for(ID in redpigs){
	ind = ind | pigdt$pigID == ID
}

# Get the test pigs
testpigsdt = pigdt[ind]
dim(testpigsdt)
unique(testpigsdt$pigID)
sum(table(testpigsdt$pigID))

ind = (testpigsdt$pigID == testpigsdt$pigID[1])
pig1 = testpigsdt[ind]



# length(diff(testpigsdt$datetimef[ind]))
# sum(ind)

# testpigsdt[, list(diffs=diff(datetimef)), by=list(pigID)]
# testpigsdt[, list(min=min(datetimef), max=max(datetimef)), by="pigID"]
# testpigsdt[, list(maxdiff=max(diff(datetimef)), 
# 				  meandiff=mean(diff(datetimef)),
# 				  mindiff=min(diff(datetimef))), by="pigID"]

# for(ID in unique(testpigsdt$pigID)){

# 	ind = testpigsdt$pigID == ID
# 	print(ID)
# 	print(rle(diff(testpigsdt[ind]$datetimef) < 100)) # In minutes

# }

############### Try fitting splines to Lat and Long positions ###################

pig1$latsc = scale(pig1$latitude)
pig1$longsc = scale(pig1$longitude)
plot(pig1$longsc, pig1$latsc, type="l")

library(splines)

pig1t = pig1[1:1000, ]
fitlong = lm(longsc ~ bs(datetimef, df=300), data=pig1t)
fitlongsp = smooth.spline(pig1t$datetimef, pig1t$longsc)

plot(pig1t$datetimef, pig1t$longsc, type="l")
lines(pig1t$datetimef, predict(fitlong), col="green")
lines(fitlongsp, col="red")
summary(fitlong)

fitlat = lm(latsc ~ bs(datetimef, df=300), data=pig1t)
fitlatsp = smooth.spline(pig1t$datetimef, pig1t$latsc)
plot(pig1t$datetimef, pig1t$latsc, type="l")
lines(pig1t$datetimef, predict(fitlat), col="green")
lines(fitlatsp, col="red")
summary(fitlat)


plot(predict(fitlongsp)$y, predict(fitlatsp)$y, type="l", col="red")
lines(predict(fitlong), predict(fitlat), type="l", col="blue")
lines(pig1t$longsc, pig1t$latsc, type="l")

pred = predict(fitlong, newdata=data.frame(datetimef=seq(min(pig1t$datetimef), max(pig1t$datetimef), len=2000)))
plot(seq(min(pig1t$datetimef), max(pig1t$datetimef), len=2000), pred, type="l")
lines(pig1t$datetimef, pig1t$longsc, type="l", col="red")


# A very heavy tailed distribution...but if we aren't making inference than
# that is ok. If we are, we could let this 
shapiro.test(residuals(fitlong))

## Load in raster data for and assign the latlongs of the texas pigs to a raster
## value

## 1. Load in cultivated land raster

r = raster("~/Downloads/2016_Cultivated_Layer/2016_Cultivated_Layer.img")
inMemory(r)

# Project pig lat longs to the appropriate projection
coordinates(pig1) <- ~longitude+latitude
proj4string(pig1) <- CRS("+proj=longlat")
pig1 <- spTransform(pig1, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

# Match the pig lat longs to the cultivated land values
vals = extract(r, pig1)
pig1$cult_land = vals

## 2. Load in cropscape data.
crop = raster("~/Downloads/CDL_2016_clip_20180118163303_1712778437/CDL_2016_clip_20180118163303_1712778437.tif")
pig1 <- spTransform(pig1, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

cvals = extract(crop, pig1)
pig1$crop_vals = cvals


## 3. Load in water body data
# The input file geodatabase
fgdb <- "/Users/mqwilber/Downloads/usgs-rivers_tx/Rivers_Streams_Waterbodies.gdb"

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)

# Read the feature class
library(rgeos)
fc <- readOGR(dsn=fgdb,layer="Waterbodies")
pig1 <- spTransform(pig1, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))

# Only look at the area of interest
fcred = crop(fc, extent(-98.8802, -98.07772, 29.10497, 29.78581))
min_dists = apply(gDistance(pig1, fcred, byid=TRUE), 2, min)
pig1$dist_water = min_dists


##### Try fitting the CTMC model to the pig data ######

pig1df = as.data.table(pig1)






