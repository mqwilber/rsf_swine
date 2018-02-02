## Data notes as of 1-24-2018

1. Get Florida pigs attribute data [Done]
2. Judas attribute data [Done]
3. Cleaning up Canada data [Done]
4. Data descriptions 
	- Close but still missing Judas and California
5. Checking for and identifying "errant" points in the GPS data.

## Data notes as of 1-25-2018

1. Data descriptions for Judas and California [For Sarah]
3. Get GPS projections for most studies [For Sarah]
	- Still need Michigan, SRS, Judas, Movement, and Tejon
4. Metadata for completed pig attribute data, GPS data, and study data [For Mark]
5. Why is SRS_Kilgo weird? [For Sarah]
	- Looks like only some of the points are weird...Might be a transient collar thing again [DONE...the weird points were there, so chopped them out]

6. Missing attribute data for 4 Canada pigs.
	- set(['canada690089A', 'canada679060A', 'canada679070A', 'canada690091A'])
7. Move pig study summary into the formatted folder
8. Get the most recent data file from Sarah for storage
9. Talk about covariates and where we can get them from.

## Date notes as of 1-26-2018

1. Need data description for California, Judas, and Canada.
2. Still need GPS projections for Canada and Movement.
3. Need to write metadata/README for files and folders in archival/ [Done]
4. Contact Ryan about missing attribute data for 4 Canada pigs: set(['canada690089A', 'canada679060A', 'canada679070A', 'canada690091A'])
5. Get the most recent data file from Sarah for storage. [Done]


## Data notes as of 2-1-2018

The raster files will all be pre-processed and will have the following format

**Organizing raster data**

covariate_data/
	elevation/
		txcamp/
			txcamp_elevation.grd
			txcamp_elevation.gri
		judas_pig/
			judas_pig_elevation.grd
			judas_pig_elevation.gri

	temperature/
		txcamp/
			txcamp_temperature_1_2016.grd # Where 1 is the month and 2016 is the year
			txcamp_temperature_1_2016.gri
			txcamp_temperature_2_2016.grd
			txcamp_temperature_2_2016.gri

The pre-processing scripts will be stored in code/covariate_scripts/
