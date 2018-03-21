## Compiling Covariates

There are 36 "studies" for which we need covariates. In general, the goal is to make a .tif raster file for each covariate for each study.

### Covariates that need extracting for each study

1. NDVI/EVI (Plant productivity)
2. Distance to water
3. Landcover (Forest canopy cover)
4. Croptypes (Agriculture)
5. Drought (Drought index)
6. Elevation
7. Distance to roads

The file `covariate_table.xlsx` contains a description of how to extract each of these covariates under the column "Extraction Methodology".

### Extraction overview

1. All of the covariates live in the directory `rsf_swine/data/covariate_data/`.
2. There is a separate folder for each covariate.
3. The structure of this folder is important in that when the raw covariate files are downloaded they need to be put in a specific place in a specific way.
4. After the covariate data is downloaded and placed in the correct directory in the correct format (see the excel file), the appropriate `extract_*.R` script should be run with the current working directory as `rsf_swine/code/covariate_scripts`. 

*Running an `extract_*.R` script for a study*

1. Open the `extract_*.R` script
2. Look for a line similar to `if(studynm %in% c('srel_contact'))`.
3. Replace the `srel_contact`, in this example, with the study name you are considering (e.g. `tejon`, `txcamp`, etc.).  
4. Covariates will be formatted for all studies that are named and this can take a long time for some studies.  So I have been just running one, or a few, at a time.

### Studies that need all covariates extracted

1. cali0, ..., cali4
2. canada
3. florida
4. ga_bill
5. ga_steve
6. judas_pig
	- Ignore until we figure out the capturing times
7. la_hartley
8. la_steve
9. michigan
10. mo_kurt0, ..., mo_kurt8
11. sc_jim2
12. scjim1
13. srel_vacuum
14. srs_kilgo
15. tx_christy_14r
16. tx_chirsty_15r
17. tx_susan
18. tx_tyler_k1
19. tx_tyler_w1
	- Need NDVI for this study still

**Priorities for extraction**

1. ga_bill
2. ga_steve
3. la_hartley
4. scjim1
5. sc_jim2
6. canada
	- This is one could be tricky as a lot of the crop data isn't consistent into Canada. We might have to exclude this one...
7. srs_kilgo
8. tx_christy_14r
9. tx_chirsty_15r
10. Others


## EVI extraction: Top priority

Currently we are using NDVI for a proxy for plant productivity.  However, EVI seems to be a much better proxy, but is harder to get a hold of in a formatted way.  All the EVI data is freely available on MODIS, but to access it you need to query MODIS, stitch relevant MODIS tiles together to get the area that you want, and reproject the MODIS data.  Ideally, we want to do this for each study site.

There are a few R packages that seem useful

- MODIStsp
- MODISTools

They both require GDAL or the MODIS tool

There is also a Python package, `pyMODIS` that would work.  This will take some fiddling around to figure out the best way to extract this from MODIS.

**Goals**

Extract EVI for all studies into 250m by 250m rasters that span the extent of the study. 

## Distance to road metric: Second priority

1. Come up with a recipe for "distance to nearest road" and ensure consistent road data across studies. Ideally, this would be the same data source for all studies.
1a. I think the recipe should be something like
	- Download road shp file (from where, I am not sure)
	- Crop road shp file to the extent of each study
	- Compute distance to nearest road point for each cell in raster.
2. Generate the distance to nearest road rasters for all studies.

## Update "distance to nearest" metric: Second priority

Currently, this metric is being computed as distance to a centroid of a polygon. The function is in `code/covariate_scripts/covariate_fxns.R` and is named `dist_nearest_neighbor`. This really should be being computed as distance to the nearest perimeter point.  Because there are many more points than polygons, this could be quite slow for some of the data.  Would be useful if we could think about ways to speed up this calculation.




