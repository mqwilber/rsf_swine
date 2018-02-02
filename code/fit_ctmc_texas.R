## Script for fitting the CTMC pig movement models to the Texas data


library(data.table)
library(raster)
library(ggplot2)
library(lubridate)
library(crawl)
library(sp)
library(rgdal)
library(rgeos)
library(ctmcmove)
library(splines)
source('/Users/mqwilber/Repos/rsf_swine/code/pigfxns.R')

## Step 0: Load in and format the pig data.
## a. For each pig, remove unrealistic observations

dat = fread("../data/formatted/full_pig_data.csv")
tcdat = dat[study == "txcamp"]
rm(dat) # Free up some space

# At least three pigs have some errant movements
tplot = ggplot(tcdat) + geom_path(aes(x=longitude, y=latitude)) + facet_wrap(~pigID)
ggsave("temppigID.pdf", width=10, height=10)

# Look at the distance between steps for each pig...slow
distdat = tcdat[, list(dist=get_distance_vect(data.frame(longitude=longitude, 
                                latitude=latitude))), by=pigID]
ggplot(distdat) + geom_boxplot(aes(x=pigID, y=dist))

## Step 1: Load in and format the necessary covariate data

## Step 2: Fit a continuous time basis function and get the continuous path

## Step 3: Fit the CTMC model using the ctmcmove package

## Step 4: Look at how resources effect movement
##      - How the the utilization distribution changes with time
##      - How introduction to an area leads affects probability of spread.

