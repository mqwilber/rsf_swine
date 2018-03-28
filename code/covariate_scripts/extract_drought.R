#!/usr/bin/env Rscript

# Extract command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

## Extract and convert drought data

library(data.table)
library(lubridate)

# Formatted summary of the studies
study_sum = fread("../../data/formatted/study_summary.csv")
study_sum$datetime_mindate = as.POSIXct(study_sum$datetime_mindate)
study_sum$datetime_maxdate = as.POSIXct(study_sum$datetime_maxdate)

for(studynm in study_sum$study){

	if(studynm %in% args){
		# Assuming only one file per study...TODO
		files = Sys.glob(file.path("../../data/covariate_data/drought/", studynm, 
																			paste(studynm, "_drought*.csv", sep="")))

		drdat = fread(files[1])

		# Compute DSCI as given by http://droughtmonitor.unl.edu/AboutUSDM/DroughtClassification.aspx
		drdat[, DSCI:=(1*D0 + 2*D1 + 3*D2 + 4*D3 + 5*D4)]
		fwrite(drdat, files[1])
	}


}
