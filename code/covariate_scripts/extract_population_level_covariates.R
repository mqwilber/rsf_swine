## Extract the population level covariates and append them to the summary file
##
## Population-level covariates include
##
## 1. Level 2 ecoregion
## 2. Pig density
## 3. Genetic similarity

library(raster)
library(sp)
library(rgeos)
library(data.table)

studysum = fread("../../data/formatted/study_summary.csv")
minsum = as.data.frame(studysum[, list(study, longitude=longitude_mean, 
																							latitude=latitude_mean)])

spdat = SpatialPointsDataFrame(minsum[, c("longitude", "latitude")],  
											 minsum[, "study", drop=F],
											 proj4string=CRS("+proj=longlat"))

# 1. Link level 2 ecoregions to studies

ecoshp = shapefile("../../data/covariate_data/ecoregion/downloads/NA_CEC_Eco_Level2.shp")
spdat = spTransform(spdat, crs(ecoshp))

pointatts = over(spdat, ecoshp)
minsum$l2name = pointatts$NA_L2NAME

# One Louisiana population doesn'tt have a value. Set it to TEXAS-LOUISIANA COASTAL PLAIN
minsum$l2name[minsum$study == "la_steve4"] = "TEXAS-LOUISIANA COASTAL PLAIN"

# Rename Mississippi alluvial and southeast usa coastal plains
ms = "MISSISSIPPI ALLUVIAL AND SOUTHEAST USA COASTAL PLAINS"
minsum$l2name[minsum$l2name == ms] = "SOUTHEAST USA COASTAL PLAINS"

studysum$l2ecoregion = minsum$l2name
fwrite(studysum, "../../data/formatted/study_summary.csv")
