## Script for analyzing movement trajectories for the pigs in Camp Bullis

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
library(ncdf4)

## Step 1: Load in data and extract Camp Bullis
dt = fread("../data/formatted/full_pig_data.csv")
pig_attr = fread("../data/formatted/pig_attributes.csv")
campBull = dt[dt$study == "txcamp"]
rm(dt) # To free up space

## Step 2: Format the covariates to prepare for analysis



## Step 3: Run the movement analysis

## Step 4: Visualize the results