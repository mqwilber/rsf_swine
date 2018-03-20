#!/bin/bash

YEAR=$1
MONTH=$2

cd ~/Repos/rsf_swine/data/covariate_data/snowdepth/downloaded/$YEAR/$MONTH

# Untar all files
for filename in *.tar
do
	tar -zxvf $filename
done

# Unzip all files
gunzip *.gz

# Remove all .Hdr files
rm *.Hdr

# Remove everything but the snowpack data 1036
find . -type f ! -name "*1036*.dat" | xargs rm

# For each file make the appropriate .hdr file
for filename in *.dat
do
	printf "ENVI\nsamples = 6935\nlines   = 3351\nbands   = 1\nheader offset = 0\nfile type = ENVI Standard\ndata type = 2\ninterleave = bsq\nbyte order = 1" > "${filename%.*}".hdr
done


# Use gdal to convert all files to GTiffs
for filename in *.dat
do
	gdal_translate -of GTiff -a_srs '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' -a_nodata -9999 -a_ullr  -124.73333333 52.87500000 -66.94166667 24.95000000 $filename "${filename%.*}".tif
done

# Clean up extra files
rm *.dat
rm *.hdr
