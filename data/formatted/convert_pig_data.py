import pandas as pd
import numpy as np
import glob
import os
import sys

"""
Description
-----------
Python script for concatenating the various GPS data

"""

file_names = glob.glob("../archival/*.csv")
file_names.remove("../archival/judas_pig.csv")

# Read and convert the Judas pig files because the current file has some issues
judas_files = glob.glob("/Users/mqwilber/Desktop/Data_for_Mark/GPS_data/JudasIN/GPS_Collar*.csv")

# This file has some WONKY dates...drop it
judas_files.remove("/Users/mqwilber/Desktop/Data_for_Mark/GPS_data/JudasIN/GPS_Collar20129_20161018112203.csv")

judas_dat = pd.concat([pd.read_csv(nm) for nm in judas_files])

# Make columns a bit easier to read
judas_dat.columns = [u'No', u'collarID', u'UTC_Date', u'UTC_Time', 
	   u'lmt_date', u'lmt_time',
       u'Origin', u'SCTS_Date', u'SCTS_Time', u'latitude',
       u'longitude', u'Heightm', u'dop', u'fixtype', u'3D_Error',
       u'Mort. Status', u'Activity', u'Main', u'Beacon', u'Temp',
       u'Easting', u'Northing', u'AnimalID', u'GroupID']

judas_dat.loc[:, "study"] = "judas_pig"

min_cols = ['collarID', 'study', 'lmt_date', 'lmt_time', 'dop', 'fixtype', 'latitude', 'longitude']
jd_trun = judas_dat[min_cols]

# Convert datetime...a bit slow but it works
# 
date_time = pd.to_datetime(jd_trun.lmt_date + " " +  jd_trun.lmt_time)
jd_trun.loc[:, "lmt_time"] = [str(a.time()) for a in date_time]
jd_trun.loc[:, "lmt_date"] = [b.replace("2016", "16") for b in jd_trun.lmt_date]



# Load in the non Judas studies
dfs = [pd.read_csv(fl) for fl in file_names]

for i, df in enumerate(dfs):
	df.loc[:, 'Study'] = os.path.split(file_names[i])[-1].split(".")[0]

full_dat = pd.concat(dfs)

full_dat.columns = [u'collarID', u'study', u'utc_date', u'utc_time', 
	   u'lmt_date', u'lmt_time', u'ecef_x_m', u'ecef_y_m', u'ecef_Z_m',
       u'latitude', u'longitude', u'dop', u'fixtype', u'sats',
       u'height_altitude', u'temp_C', u'easting', u'northing', u'easting1',
       u'northing1', u'habitat', u'sex']

full_dat_trun = full_dat[min_cols]

comb_dat = pd.concat([full_dat_trun, jd_trun])

# Save the resulting pig data
comb_dat.to_csv("full_pig_data.csv", index=False)



