#!/bin/bash

for studynm in "$@"
do
	printf "Extracting NDVI...\n"
	Rscript extract_ndvi.R $studynm
	#printf "Extracting water...\n"
	#Rscript extract_water.R $studynm
	printf "Extracting landcover...\n"
	Rscript extract_landcover.R $studynm
	printf "Extracting drought...\n"
	Rscript extract_drought.R $studynm
	printf "Extracting croplayer...\n"
	Rscript extract_croplayer.R $studynm
	printf "Extracting elevation..\n"
	Rscript extract_elevation.R $studynm
done