import pandas as pd
import numpy as np
import glob
import os
import sys

"""
Description
-----------
Python script converts, and formats all GPD data for the pig RSF study.
All the raw datasets are stored in ../archival and are converted in this script.

"""

##############################################################
### Look at the Tejon and California datasets specifically ###
##############################################################

print("Working on Tejon and California data...")

tejon = pd.read_csv("../archival/tejon_data.csv")
cali = pd.read_csv("../archival/california_lindsey_data.csv")

del tejon['Unnamed: 0']
tejon.columns = ['collarID', 'datetime', 'latitude',  'longitude', 'altitude', 
                 'fixtype', 'dop', 'temp_C', 'sex']

tejon.loc[:, "study"] = "tejon"

cali.columns = [u'PigIDforms', u'collarID', u'UID', u'Date', u'Time', u'latitude',
       u'longitude', u'Ecoregion', u'Ecoregion_code', u'Site', u'Site_Code',
       u'County', u'County_code', u'Sex', u'Sex_code', u'Age', u'Age_code',
       u'Sounder', u'ClusterID', u'GroupID', u'Year', u'Month', u'Week',
       u'Day', u'Night', u'Offset_hr', u'Offset_min', u'Temp', u'Elev',
       u'Slope', u'NDVI_notscaled', u'NDVI', u'STEPLENGTH', u'Dist_hr',
       u'Near_locrd', u'Near_majrd', u'Near_NHD', u'Near_hylin', u'Near_fevwa',
       u'Near_Ag', u'Near_Barre', u'Near_Conif', u'Near_Hardw', u'Near_Herb',
       u'Near_urban', u'Near_shrub']

cali.loc[:, 'datetime'] = cali.Date + " " + cali.Time
cali.loc[:, "study"] = "caliLind"
cali.loc[:, 'dop'] = np.nan
cali.loc[:, 'fixtype'] = np.nan

# Extract the necessary columns from the data
tejon_red = tejon[['collarID', 'study', 'latitude', 'longitude', 'dop', 'fixtype', 'datetime']]
cali_red = cali[['collarID', 'study', 'latitude', 'longitude', 'dop', 'fixtype', 'datetime']]

tejon_red.loc[:, "pigID"] = tejon_red.study + tejon_red.collarID.astype("S21")
cali_red.loc[:, "pigID"] = cali_red.study + cali_red.collarID.astype("S21")

# Get the cali datatimes in the right format...slow because all different formats
cali_red.loc[:, "datetime"] = pd.to_datetime(cali.datetime)
print("Done")

#########################################################
### Judas pig files need to be loaded in individually ### 
#########################################################

print("Working on Judas pig data...")

# Read and convert the Judas pig files because the current file has some issues
judas_files = glob.glob("../archival/judas_pigs/GPS_Collar*.csv")

# This file has some WONKY dates...drop it
judas_files.remove("../archival/judas_pigs/GPS_Collar20129_20161018112203.csv")

judas_dat = pd.concat([pd.read_csv(nm) for nm in judas_files])

# Make columns a bit easier to read
judas_dat.columns = [u'No', u'collarID', u'UTC_Date', u'UTC_Time', 
	   u'lmt_date', u'lmt_time',
       u'Origin', u'SCTS_Date', u'SCTS_Time', u'latitude',
       u'longitude', u'Heightm', u'dop', u'fixtype', u'3D_Error',
       u'Mort. Status', u'Activity', u'Main', u'Beacon', u'Temp',
       u'Easting', u'Northing', u'AnimalID', u'GroupID']

judas_dat.loc[:, "study"] = "judas_pig"
judas_dat.loc[:, "datetime"] = pd.to_datetime(judas_dat.lmt_date + " " + judas_dat.lmt_time)
judas_red = judas_dat[['collarID', 'study', 'latitude', 'longitude', 'dop', 'fixtype', 'datetime']]
judas_red.loc[:, "pigID"] = judas_red.study + judas_red.collarID.astype("S21")

print("Done")

###########################################################################
###  Load in the non-Judas studies and non-tejon/non-california studies ###
###########################################################################

print("Working on combined data...")

file_names = glob.glob("../archival/*.csv")
file_names.remove("../archival/california_lindsey_data.csv") # Not yet formatted
file_names.remove("../archival/tejon_data.csv") # Not yet formatted
file_names.remove("../archival/contact_study_all_pigs.csv") # Not yet formatted
file_names.remove("../archival/srel_jim.csv")

dfs = [pd.read_csv(fl) for fl in file_names]

# Loop through each of the files
for i, df in enumerate(dfs):

    study = os.path.split(file_names[i])[-1].split(".")[0]
    df.loc[:, 'Study'] = study

    df.columns = [u'collarID', u'study', u'utc_date', u'utc_time', 
       u'lmt_date', u'lmt_time', u'ecef_x_m', u'ecef_y_m', u'ecef_Z_m',
    u'latitude', u'longitude', u'dop', u'fixtype', u'sats',
    u'height_altitude', u'temp_C', u'easting', u'northing', u'easting1',
    u'northing1', u'habitat', u'sex']

    df.loc[:, "pigID"] = df.study + df.collarID.astype("S21")

    print(study)
    if study == "michigan":
           df.loc[:, "lmt_time"] = df.lmt_time + ":00"

    # # Convert datetime
    # if study != "srel_jim":
    #     df.loc[:, "datetime"] = pd.to_datetime(df.lmt_date + " " + df.lmt_time, 
                                format="%m/%d/%y %H:%M:%S")
    df.loc[:, "datetime"] = pd.to_datetime(df.lmt_date + " " + df.lmt_time, 
                                format="%d/%m/%y %H:%M:%S")


full_dat = pd.concat(dfs)

print("Done")

##############################################
### Load and format the Hartley Hogs files ###
##############################################

# print("Working on the Hartley Hogs...")

# hart = pd.read_csv("../archival/hartley_hogs.csv")
# hart.columns = [u'objectid', u'collarID', u'year', u'datetime', u'latitude',
#        u'longitude', u'temp_f', u'hour24', u'sex', u'distance']

# hart.loc[:, "study"] = "hartley"
# hart.loc[:, "pigID"] = hart.study + hart.collarID.astype("S21")
# hart.loc[:, 'dop'] = np.nan
# hart.loc[:, "fixtype"] = np.nan

# print("Working on hartley datetime...")
# hart.loc[:, "datetime"] = pd.to_datetime(hart.datetime)

# print("Done")


##########################################
### Load and format the SREL Jim study ###
##########################################

# Each file needs to be loaded individually

# print("Working on SREL files...")
# srel_files = glob.glob("../archival/srel_jim/*.xlsx")
# sreldfs = [pd.read_excel(sf, sheetname="cleaned") for sf in srel_files]

# # 171M
# for i, df in enumerate(sreldfs):

#     cid = os.path.split(srel_files[i])[-1].split("_")[0]

#     # Regex

#     if cid in ["P208M", "P163F", "P164F", "P166M", "P167F", "P169F"]:

#         df.columns = [u'fix_num', u'Unnamed: 1',  u'Unnamed: 2', u'datetime',
#                       u'Unnamed: 4', u'Unnamed: 5', u'Unnamed: 6', u'Unnamed: 7', 
#                       u'latitude', u'longitude',   u'altitude ', u'Time',
#                       u'Temp ',  u'fixtype',  u'sats', u'dop']

#         df.loc[:, 'collarID'] = os.path.split(srel_files[i])[-1].split("_")[0]
#         df.loc[:, "study"] = "srel_jim"
#         df.loc[:, "pigID"] = df.study + df.collarID

#     elif cid in ["P156F", "P161M", "P168M", "P180M", "P194F", "P197F", ]:

#     else cid in ["P172M", "P202M"]


#########################################
### Load and format the contact study ###
#########################################

print("Working on contact study...")

contact = pd.read_csv("../archival/contact_study_all_pigs.csv")

contact.columns = [u'collarID', u'LocNum', u'Acquisition Time', u'Acquisition Start Time',
       u'GPS Fix Time', u'GPS Date', u'GPS Time', u'fixtype',
       u'latitude', u'longitude', u'GPS UTM Zone', u'GPS UTM Northing',
       u'GPS UTM Easting', u'GPS Altitude', u'GPS Speed', u'GPS Heading',
       u'GPS Horizontal Error', u'dop',
       u'GPS Satellite Bitmap', u'GPS Satellite Count', u'GPS Navigation Time',
       u'Activity Count', u'Temperature', u'Predeployment Data', u'Error']

contact.loc[:, "datetime"] = pd.to_datetime(contact['GPS Fix Time'], format="%Y.%m.%d %H:%M:%S")
contact.loc[:, "study"] = "contact"
contact.loc[:, "pigID"] = contact.study + contact.collarID

# Drop non-pig columns
contact_red = contact[~((contact.pigID == "contacthard drive") | (contact.pigID == "contact**time") | 
                     (contact.pigID == "contactdata on ") | (contact.pigID == "contacton original"))]


print("Done")

##########################################
### Format and combine all of the data ###
##########################################

cols = ['pigID', 'collarID', 'study', 'latitude', 'longitude', 'dop', 'fixtype', 'datetime']
full_dat_trun = full_dat[cols]
judas_red = judas_red[cols]
tejon_red = tejon_red[cols]
cali_red = cali_red[cols]
contact_red = contact_red[cols]
# hart_red = hart[cols]

comb_dat = pd.concat([full_dat_trun, judas_red, tejon_red, cali_red, contact_red])

# Save the resulting pig data
comb_dat.to_csv("full_pig_data.csv", index=False)



