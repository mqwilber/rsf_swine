# Load in the pig data

library(data.table)
dt = fread("full_pig_data.csv")

# Remove positions that don't exist
ind = is.na(dt$latitude) | is.na(dt$longitude)
dt_trun = dt[!ind, ]

# Make a unique pigID
dt_trun[, pigID:=paste(collarID, study, sep="")]

# Date are not in the same format
lmt_date = tstrsplit(dt_trun$lmt_date, " ")[[1]]
dt_trun$lmt_date_strp = lmt_date


datetime = as.character(paste(dt_trun$lmt_date_strp, dt_trun$lmt_time))

# studies michigan and srel_jim don't have ant seconds on the datetimes..add 00
for(study in c("michigan", "srel_jim")){
  datetime[dt_trun$study == study] = paste(datetime[dt_trun$study == study], ":00", sep="")
}

dt_trun$datetime = as.POSIXct(strptime(datetime, format="%m/%d/%y %H:%M:%S", tz="GMT"))

# Save the cleaned data
write.csv(dt_trun, "full_pig_data_cleaned.csv")