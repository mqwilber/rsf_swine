
## Compiling Covariates

There are 24 "studies" for which we need covariates.  There will actually be more as some of them are pretty dispersed and we will need to group these into smaller areas so we can compile and extract covariates in a reasonable amount of time and feel confident making inference on the population-level as opposed to the meta-population level.

In general, the goal is to make a .tif raster file for each covariate for each study.

### Covariates that need extracting for each study

1. NDVI (Plant productivity)
2. Distance to water
3. Landcover (Forest canopy cover)
4. Croptypes (Agriculture)
5. Drought (Drought index)
6. Elevation
7. Distance to roads

The file `covariate_table_Feb2018.xlsx` contains a description of how to extract each of these covariates under the column "Extraction Methodology".

### Extraction overview

1. All of the covariates live in the directory `rsf_swine/data/covariate_data/`.
2. There is a separate folder for each covariate.
3. The structure of this folder is important in that when the raw covariate files are downloaded they need to be put in a specific place in a specific way.
4. After the covariate data is downloaded and placed in the correct directory in the correct format (see the exfel , the appropriate `extract_*.R` script should be run with the current working directory as `rsf_swine/code/covariate_scripts`. 

*Running an `extract_*.R` script for a study*

1. Open the `extract_*.R` script
2. Look for a line similar to `if(studynm %in% c('srel_contact'))`.
3. Replace the `srel_contact`, in this example, with the study name you are considering (e.g. `tejon`, `txcamp`, etc.).  
4. Covariates will be formatted for all studies that are named and this can take a long time for some studies.  So I have been just running one, or a few, at a time.

## Studies that need covariate extraction

1. cali
	- Will likely need to be broken into subpopulations for analyses
2. canada
	- Hold off on this one as getting comparable data for this study could be tough.
3. florida
4. ga_bill
5. ga_steve
6. judas_pig
	- Ignore until we figure out the capturing times
7. la_hartley
8. la_steve
9. michigan
10. mo_kurt
	- Will likely need to be broken into sub-populations
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


**Priorities**

1. cali
2. florida
3. michigan
4. mo_kurt
5. tx_tyler_k1
6. tx_susan
7. la_steve












