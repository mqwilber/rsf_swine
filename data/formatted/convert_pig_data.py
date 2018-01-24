import pandas as pd
import numpy as np
import glob
import os
import sys

"""
Description
-----------
Python script converts, and formats all GPs data for the pig RSF study.
All the raw datasets are stored in ../archival and are converted in this script.

Pig attribute data are converted to match the the format of the GPS data.

Both converted files are saved in the formatted data folder

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

cali.columns = [u'collarID', u'PID', u'UID', u'Date', u'Time', u'latitude',
       u'longitude', u'Ecoregion', u'Ecoregion_code', u'Site', u'Site_Code',
       u'County', u'County_code', u'Sex', u'Sex_code', u'Age', u'Age_code',
       u'Sounder', u'ClusterID', u'GroupID', u'Year', u'Month', u'Week',
       u'Day', u'Night', u'Offset_hr', u'Offset_min', u'Temp', u'Elev',
       u'Slope', u'NDVI_notscaled', u'NDVI', u'STEPLENGTH', u'Dist_hr',
       u'Near_locrd', u'Near_majrd', u'Near_NHD', u'Near_hylin', u'Near_fevwa',
       u'Near_Ag', u'Near_Barre', u'Near_Conif', u'Near_Hardw', u'Near_Herb',
       u'Near_urban', u'Near_shrub']

cali.loc[:, 'datetime'] = cali.Date + " " + cali.Time
cali.loc[:, "study"] = "cali"
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


judas_dat = pd.read_csv("../archival/judas_pigs.csv")
judas_dat.columns = [u'collarID', u'study', u'UTC_Date', u'UTC_Time', u'LMT_Date',
       u'LMT_Time', u'x', u'y', u'z', u'latitude', u'longitude',
       u'dop', u'fixtype', u'Height [m]', u'Temp', u'Easting',
       u'Northing', u'Unnamed: 17', u'Unnamed: 18', u'Unnamed: 19',
       u'Unnamed: 20', u'Unnamed: 21', u'Unnamed: 22', u'Unnamed: 23',
       u'Unnamed: 24', u'Unnamed: 25', u'Unnamed: 26', u'Unnamed: 27']
judas_dat.loc[:, "study"] = "judas_pig"
judas_dat.loc[:, "datetime"] = pd.to_datetime(judas_dat.LMT_Date + " " + judas_dat.LMT_Time,
                                      format="%m/%d/%y %H:%M:%S")
judas_dat.loc[:, "pigID"] = judas_dat.study + pd.Series([str(cID).split("_")[-1] for cID in judas_dat.collarID])

print("Done")

###########################################################################
###  Load in the non-Judas studies and non-tejon/non-california studies ###
###########################################################################

print("Working on combined data...")

file_names = ["txcamp.csv", "florida.csv", "movepig.csv", "michigan.csv", 
                    "srel_vacuum.csv", "kilgo.csv"]

dfs = [pd.read_csv("../archival/" + fl) for fl in file_names]

# Loop through each of the files
for i, df in enumerate(dfs):

    study = os.path.split(file_names[i])[-1].split(".")[0]

    if study == "kilgo":
      df.loc[:, 'Study'] = "srs_" + study
    else:
      df.loc[:, 'Study'] = study

    df.columns = [u'collarID', u'study', u'utc_date', u'utc_time', 
       u'lmt_date', u'lmt_time', u'ecef_x_m', u'ecef_y_m', u'ecef_Z_m',
    u'latitude', u'longitude', u'dop', u'fixtype', u'sats',
    u'height_altitude', u'temp_C', u'easting', u'northing', u'easting1',
    u'northing1', u'habitat', u'sex']

    if study == "srel_vacuum":
      df.loc[:, "pigID"] = [study + str(ID).strip("M").strip("F") for ID in df.collarID]
    elif study == "movepig":

      # Clean up ID to match attributes
      df.loc[:, "collarID"] = df.collarID.str.replace("LA_Hart_", "LASt_")
      df.loc[:, "collarID"] = [ID.strip("n").replace("n2", "") for ID in df.collarID]
      df.loc[:, "pigID"] = df.study + df.collarID.astype("S21")
    else:
      df.loc[:, "pigID"] = df.study + df.collarID.astype("S21")

    print(study)
    if study == "michigan":
           df.loc[:, "lmt_time"] = df.lmt_time + ":00"

    df.loc[:, "datetime"] = pd.to_datetime(df.lmt_date + " " + df.lmt_time, 
                                format="%m/%d/%y %H:%M:%S")


full_dat = pd.concat(dfs)

print("Done")

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
contact.loc[:, "study"] = "srel_contact"
contact.loc[:, "pigID"] = contact.study + "P" + contact.collarID.str.strip("M").str.strip("F")

# Drop non-pig columns
contact_red = contact[~((contact.pigID == "srel_contactPhard drive") | (contact.pigID == "srel_contactP**time") | 
                     (contact.pigID == "srel_contactPdata on ") | (contact.pigID == "srel_contactPon original"))]


print("Done")

##########################################
### Format and combine all of the data ###
##########################################

cols = ['pigID', 'collarID', 'study', 'latitude', 'longitude', 'dop', 'fixtype', 'datetime']
full_dat_trun = full_dat[cols]
judas_red = judas_dat[cols]
tejon_red = tejon_red[cols]
cali_red = cali_red[cols]
contact_red = contact_red[cols]
# hart_red = hart[cols]

comb_dat = pd.concat([full_dat_trun, judas_red, tejon_red, cali_red, contact_red])

# Drop all Null lats and longs
comb_dat = comb_dat[~(comb_dat.longitude.isnull() | comb_dat.latitude.isnull())]

# Drop non int lats and longs
inds = []
for i in comb_dat.latitude:
  try:
    np.float(i)
    inds.append(True)
  except Exception:
    inds.append(False)

comb_dat = comb_dat[inds]

# Save the resulting pig data
comb_dat.to_csv("full_pig_data.csv", index=False)

### Link the attributes pigIDs and study IDs ###

# TODO: Add in FLORIDA data
attrib = pd.read_csv("../archival/all_pig_attributes.csv")
attrib.replace({'Study' : {'SREL-Contact' : 'srel_contact', 
                           'Judas' : 'judas_pig',
                           'Camp Bullis' : 'txcamp',
                           'SREL-vacuum' : 'srel_vacuum',
                           'SRS' : 'srs_kilgo',
                           'Michigan': 'michigan',
                           'Movement': "movepig",
                           'Tejon': 'tejon',
                           'CA': 'cali'}}, inplace=True)

attrib.columns = [u'collarID', u'sex', u'age_class', u'weight_lb', u'start_date',
       u'end_date', u'study', u'notes']

#comb_dat.groupby('study')['pigID'].unique()['cali']
attrib.loc[:, "pigID"] = attrib.study + attrib.collarID

# Clean up weights
attrib.loc[:, "weight_lb"] = attrib.weight_lb.str.replace(".*\(estimated\)", "").str.strip("+").str.strip("-").str.replace(" ", "")
attrib = attrib.replace('', np.nan, regex=True)
attrib.loc[:, "weight_lb"] = attrib.weight_lb.astype(np.float)

attrib.to_csv("pig_attributes.csv", index=False)

# Check that all IDs match...all IDs are matching!
gb = comb_dat.groupby('study') 

for study in attrib.study.unique():
  cvals = gb['pigID'].unique()[study]
  atvals = attrib.loc[attrib.study == study, "pigID"]

  print(study)
  print(set(cvals) - set(atvals))

  if study == "judas_pig":
    print(set(cvals) - set(atvals))

  if study == "movepig":

    print(len(cvals))
    print(len(atvals))

    # Pigs with no attribute data
    diff = set(cvals) - set(atvals)

    # Attributes with no pig data
    diff2 = set(atvals) - set(cvals)
    print("Unique: {0}".format(len(set(atvals) - set(cvals)) + len(set(cvals) - set(atvals))))
    print("Intersect: {0}".format(len(np.intersect1d(atvals, cvals))))


comb_dat[comb_dat.study == "movepig"].pigID.value_counts().loc[diff]











