## Script builds monthly rasters of Snow Depth data over the United States
## using the SNODAS database.
##
## The resulting file is a GeoTIFF file with the monthly average snowpack for 
## each year specified in the script.
##
## Data is at a 1 km by 1 km scale
##
## Author: Mark Wilber
library(parallel)
library(raster)

snodas_url = "ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked"
base = "/Users/mqwilber/Repos/rsf_swine/data/covariate_data"

download_snodas_files = function(yr, month, months, base, snodas_url, yrdir){

		minday = TRUE

		for(day in 1:31){

				# Download file the file exists
				strday = ifelse(day >= 10, as.character(day), paste("0", day, sep=""))

		    minday = tryCatch({

		    	fpath = file.path(snodas_url, yr, paste(month, months[[month]], sep="_"), 
		    										paste("SNODAS_", yr, month, strday, ".tar", sep=""))
		    	dir.create(file.path(yrdir, month), showWarnings=FALSE)
		    	download.file(fpath, file.path(yrdir, month, paste("SNODAS_", yr, month, strday, ".tar", sep="")))
		      minday = TRUE
		      minday

      	}, error = function(err){

	      	minday = FALSE
	      	return(minday)

	      })

		}
}


years = 2004:2005
months = list("01"="Jan", "02"="Feb", "03"="Mar", "04"="Apr", "05"="May", 
							"06"="Jun", "07"="Jul", "08"="Aug", "09"="Sep", 
							"10"="Oct", "11"="Nov", "12"="Dec")

month_fxn = function(month, yr, months, base, snodas_url, yrdir){

		# Get all the SNODAS files for a particular month
	res = download_snodas_files(yr, month, months, base, snodas_url, yrdir)

	# Untar all of these files and convert them to tifs
	system(paste("./convert_snodas.sh", yr, month))

	# Load all .tifs, take mean, and save ras
	alltifs = Sys.glob(file.path(yrdir, month, "*.tif"))
	rasStack = do.call(stack, lapply(alltifs, raster))
	meanras = mean(rasStack)
	writeRaster(meanras, file.path(yrdir, month, paste("snowdepth_", month, "_", yr, ".tif", sep="")), 
													format="GTiff", overwrite=TRUE)
	rm(meanras)

	system(paste("rm", file.path(yrdir, month, "us_ssmv*.tif")))

}


for(yr in years){

	yrdir = file.path(base, "snowdepth", "downloaded", yr)
	dir.create(yrdir, showWarnings = FALSE)

	res = mclapply(names(months), month_fxn, yr, months, base, snodas_url, yrdir, mc.cores=4)

}