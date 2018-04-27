import pandas as pd
import numpy as np
import glob
import os
import sys
import string
from sklearn.cluster import KMeans
import pytz 

"""
Description
-----------
Python script converts, and formats all GPs data for the pig RSF study.
All the raw datasets are stored in ../archival and are converted in this script.

Pig attribute data are converted to match the the format of the GPS data.

Both converted files are saved in the formatted data folder.

Throughout this script various adjustments are made to the data. Some of the
important adjustments include

1. Points from SRS Kilgo are dropped that are clearly erroneous based on outlying
GPS locations.

2. Michigan pigs that were sampled in the year 2000 are removed as there are no
H:M:S fixes on the GPS dates.

2.0. All michigan points below 43.75 latitude are remove as these are pre-release
and post-capture data.

2.1. Errant txcamp points are removed.

2.2. Errant srel_contact points are removed

3. Data are cleaned and IDs are changed to match with pig attributes file

NOTE
----
Must be run with Python 2.7

"""

# Time-zone information of the raw data

# tzs = {'cali0': "America/Los_Angeles", # Check
#  'cali1': "America/Los_Angeles", # Check
#  'cali2': "America/Los_Angeles", # Check
#  'cali3': "America/Los_Angeles",  # Check
#  'cali4': "America/Los_Angeles", # Check
#  'canada': "America/Regina", # Don't know...need to ask Ryan
#  'fl_raoul': "America/New_York", # Seems to be UTC time
#  'florida': "America/New_York", 
#  'ga_bill': "America/New_York",
#  'ga_steve': "America/New_York",
#  'judas_pig': "LMT", # Check
#  'la_hartley': "America/Chicago",
#  'la_steve': "America/Chicago", 
#  'michigan': "America/New_York", # Check
#  'mo_kurt0': "America/Chicago",  
#  'mo_kurt1': "America/Chicago",
#  'mo_kurt2' : "America/Chicago",
#  'mo_kurt3' : "America/Chicago", 
#  'mo_kurt4' : "America/Chicago",
#  'mo_kurt5': "America/Chicago", 
#  'mo_kurt6': "America/Chicago", 
#  'mo_kurt7': "America/Chicago",
#  'mo_kurt8': "America/Chicago", 
#  'sc_jim2': "America/New_York", 
#  'scjim1' : "America/New_York", 
#  'srel_contact': "America/New_York", # Check
#  'srel_vacuum': "America/New_York", # Looks like LMT time...
#  'srs_kilgo': "America/New_York", # Check
#  'tejon': "America/Los_Angeles", # Check
#  'tx_christy_14r': "America/Chicago", 
#  'tx_christy_15r': "America/Chicago",
#  'tx_susan': "America/Chicago", 
#  'tx_tyler_k1': "America/Chicago", 
#  'tx_tyler_w1': "America/Chicago", 
#  'tx_tyler_w2': "America/Chicago", 
#  'txcamp': "America/Chicago" #check}

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

# Split the California populations into separate populations that are spatially
# proximate
kobj = KMeans(n_clusters=5)
X = cali_red[['longitude', 'latitude']].values
kfit = kobj.fit(X)
labels = kfit.labels_

cali_red.loc[:, "study"] = np.array("cali" + pd.Series(labels).astype(np.str))

# Relabel studies to be consistent upon regeneration of data
# Because I already created rasters for some of these files, 
# can't do a systematic reordering.
cali_latlon = cali_red.groupby('study').agg({'longitude':np.mean, 
                                             'latitude': np.mean})

cali0ind = (cali_latlon.latitude.round(6) == 37.888436) & \
        (cali_latlon.longitude.round(6) == -121.865549)

cali1ind = (cali_latlon.latitude.round(6) == 34.986554) & \
        (cali_latlon.longitude.round(6) == -118.690084)

cali2ind = (cali_latlon.latitude.round(6) == 38.928151) & \
        (cali_latlon.longitude.round(6) == -123.315835)

cali3ind = (cali_latlon.latitude.round(6) == 38.357676) & \
        (cali_latlon.longitude.round(6) == -122.105563)

cali4ind = (cali_latlon.latitude.round(6) == 37.455387) & \
        (cali_latlon.longitude.round(6) == -121.806334)


nms = ["cali" + str(i) for i in range(5)]
nminds = [i[i].index[0] for i in  [cali0ind, cali1ind, cali2ind, cali3ind, cali4ind]]
repdict = {key : val for key, val in zip(nms, nminds)}
revalue = cali_red.study.map(repdict)
cali_red = cali_red.assign(study=revalue)

# Relabel the cali data for consistency
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

cIDs = pd.Series([str(cID).split("_")[-1] for cID in judas_dat.collarID])
judas_dat.loc[:, "pigID"] = judas_dat.study + judas_dat.collarID.astype("S21") #pd.Series([str(cID).split("_")[-1] for cID in judas_dat.collarID])

print("Done")

##############################
### Format the Kilgo Study ###
##############################

print("Working on SRS_Kilgo...")

kilgo = pd.read_csv("../archival/kilgo.csv")

kilgo.columns = [u'study', u'collarID', u'sex', u'age', u'eartag', u'gmt_date', u'gmt_time',
     u'lmt_date', u'lmt_time', u'ecef_x', u'ecef_y', u'ecef_z', u'latitude',
     u'longitude', u'height', u'dop', u'fixtype', u'validated', u'sats_used',
     u'ch1_satid', u'ch1_c/n', u'ch2_satid', u'ch2_c/n', u'ch3_satid',
     u'ch3_c/n', u'ch4_satid', u'ch4_c/n', u'ch5_satid', u'ch5_c/n',
     u'ch6_satid', u'ch6_c/n', u'ch7_satid', u'ch7_c/n', u'ch8_satid',
     u'ch8_c/n', u'ch9_satid', u'ch9_c/n', u'ch10_satid', u'ch10_c/n',
     u'ch11_satid', u'ch11_c/n', u'ch12_satid', u'ch12_c/n', u'main_vol',
     u'bu_vol', u'temp', u'remarks']

kilgo.loc[:, 'study'] = "srs_kilgo"
kilgo.loc[:, 'pigID'] = kilgo.study + kilgo.collarID
kilgo.loc[:, "datetime"] = pd.to_datetime(kilgo.lmt_date + " "  + kilgo.lmt_time, 
                              format="%m/%d/%y %H:%M:%S")

# Clean errant kilgo points
kilgo_trun = kilgo[(kilgo.latitude > 32) & (kilgo.latitude < 35) & 
               (kilgo.longitude < -79) & (kilgo.latitude > -82.5)]

# Based on John Kilgo's suggestion remove all dates prior to 12/8/2015 for ID srs_kilgo16-4-2
kilgo_trun = kilgo_trun[~((kilgo_trun.datetime < pd.to_datetime("12/8/2015")) & 
														(kilgo_trun.pigID == "srs_kilgo16-4-2"))]

print("Done")

###########################################################################
###  Load in the non-Judas studies and non-tejon/non-california studies ###
###########################################################################

print("Working on combined data...")

file_names = ["txcamp.csv", "florida.csv", "michigan.csv", 
                    "srel_vacuum.csv"]

dfs = [pd.read_csv("../archival/" + fl) for fl in file_names]

# Loop through each of the files
for i, df in enumerate(dfs):

    study = os.path.split(file_names[i])[-1].split(".")[0]

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

# For the michigan study, exclude all the points from year 2000 as they don't have
# daily times. Also exclude all points below 43.75 latitude
full_dat.loc[:, "year"] = [a.year for a in full_dat.datetime] 
full_dat = full_dat[~((full_dat.year == 2000) & (full_dat.study == "michigan"))]

# For the txcamp study, remove all points < -98.8 longitude and > 29.80 lat.
full_dat = full_dat[~((full_dat.study == "txcamp") & 
                    ((full_dat.longitude < -98.8) | (full_dat.latitude > 29.80)))]

print("Done")

##########################################
### Load and format the movement study ###
##########################################

print("Working on the movement study...")

move = pd.read_csv("../archival/movepig.csv")

move.columns = [u'collarID', u'date', u'time', u'study', u'latitude', u'longitude', 
					u'dop', u'temp', u'yr', u'mo', u'hab2']

move.loc[:, "study"] = move.study.str.lower()
move.loc[:, "collarID"] = move.collarID.str.lower()

# It is not clear which studies are UTC datetime, though fl raoul is definitely UTC
move.loc[:, "datetime"] = pd.to_datetime(move.date + " " + move.time, format="%m/%d/%y %H:%M:%S")

# Convert tx_tyler to LMT
central = pytz.timezone("America/Chicago")
baserl = move[move.study == 'tx_tyler_w2'].datetime
updts = pd.Index(baserl).tz_localize(pytz.utc).tz_convert(central).tz_localize(None)
move.loc[move.study == "tx_tyler_w2", "datetime"] = pd.Series(updts).values

move.loc[:, "fixtype"] = np.nan
move.loc[move.study == "sc_jim", "study"] = "scjim1"

# Format pig IDs
idmapper = {'la_hart': 'la_hartley', 'scji': 'scjim1', 'gabi' : 'ga_bill',
			'gast': 'ga_steve', 'flra': 'fl_raoul', 'last': 'la_steve', 
			'moku': 'mo_kurt', 'sc_jim': 'sc_jim2', 'txch14': 'tx_christy_14r',
			'txch15': 'tx_christy_15r', 'txsu': 'tx_susan', 'txtyk1': 'tx_tyler_k1',
			'txtyw1': 'tx_tyler_w1', 'txtyw2': 'tx_tyler_w2'}

for idm in idmapper.iterkeys():
	move.loc[:, "collarID"] = move.collarID.str.replace(idm, idmapper[idm])

move.loc[:, "collarID"] = move.collarID.str.replace("_ ", "_").str.replace("__", "_")
move.loc[:, "collarID"] = move.collarID.str.replace("f_ho", "fho").str.replace("g_ho", "gho").str.replace("h_ho", "hho").str.replace("d_ho", "dho")

# Find and replace...
move.loc[:, "collarID"] = move.collarID.str.replace(r"([0-9]+)_([a-d])", r"\1\2")

move.loc[:, "collarID"] = [val.split("_")[-1] for val in move.collarID]
move.loc[:, "pigID"] = move.study + move.collarID

# Break mo_kurt into multiple spatially proximate studies
mo_kurt = move[move.study == "mo_kurt"]
kobj = KMeans(n_clusters=9)
X = mo_kurt[['longitude', 'latitude']].values
kfit = kobj.fit(X)
labels = kfit.labels_
move.loc[move.study == "mo_kurt", "study"] = np.array("mo_kurt" + pd.Series(labels).astype(np.str))

mo_kurt = move[move.study.str.find("mo_kurt") != -1]
mo_kurtll = mo_kurt.groupby('study').agg({k : np.mean for k in ['latitude', 'longitude']})

# Fix mo_kurt0
mokurt0ind = (mo_kurtll.latitude.round(6) == 37.583448) & \
             (mo_kurtll.longitude.round(6) == -90.897654)

# Label in a consistent order
othermk = mo_kurtll[~mokurt0ind].sort_values(['longitude', 'latitude'])
othermk.loc[:, 'newstudy'] = ["mo_kurt" + str(i) for i in range(1, 9)]

repdict = {mokurt0ind[mokurt0ind].index[0] : 'mo_kurt0'}
repdict.update({key : val for key, val in zip(othermk.index, othermk.newstudy)})

revalue = mo_kurt.study.map(repdict)

# Rename mo_kurt studies consistently
move.loc[mo_kurt.index, "study"] = revalue

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

# This is in UTC time, convert to Local eastern time
contact.loc[:, "datetime"] = pd.to_datetime(contact['GPS Fix Time'], format="%Y.%m.%d %H:%M:%S")
east = pytz.timezone("America/New_York")
localdatetime = pd.Index(contact.datetime).tz_localize(pytz.utc).tz_convert(east).tz_localize(None)
contact.loc[:, "datetime"] = pd.Series(localdatetime)


contact.loc[:, "study"] = "srel_contact"
contact.loc[:, "pigID"] = contact.study + "P" + contact.collarID.str.strip("M").str.strip("F")

# Drop non-pig columns
contact_red = contact[~((contact.pigID == "srel_contactPhard drive") | (contact.pigID == "srel_contactP**time") | 
                     (contact.pigID == "srel_contactPdata on ") | (contact.pigID == "srel_contactPon original"))]


# Remove errant srel_contact points
ind = (((contact.pigID == "srel_contactP702") & (contact.latitude > 33.33)) | \
      ((contact.pigID == "srel_contactP705") & (contact.latitude > 33.33)) | \
      ((contact.pigID == "srel_contactP708") & (contact.latitude < 33.2)))

contact_red = contact_red[~ind]

print("Done")


#####################################################
### Load and format the various Canada data files ###
#####################################################

print("Working on Canada data...")

filenames = glob.glob("../archival/canada_pigs/6*.csv")
idlinker = pd.read_csv("../archival/canada_pigs/animalID_to_collarID.csv")
idlinker.set_index("collar_ID", inplace=True)

dfs = []
for i, fl in enumerate(filenames):

  # get the collar IDs
  collarID = os.path.split(fl)[-1].split("_")[0]

  # If they exist, extract the animal ID
  try:
    animalID = idlinker.loc[str(collarID)][0]
  except:
    animalID = ""

  print("Canada pigID {0}".format(collarID + animalID))
  df = pd.read_csv(fl, skiprows=23) # Ignore the headers
  df.columns = [u'Acquisition Time', u'Acquisition Start Time', u'Iridium CEP Radius',
       u'Iridium Latitude', u'Iridium Longitude', u'datetime',
       u'fixtype', u'latitude', u'longitude', u'GPS UTM Zone',
       u'GPS UTM Northing', u'GPS UTM Easting', u'GPS Horizontal Error',
       u'GPS Navigation Time', u'Activity Count', u'Temperature',
       u'Satellite Uplink', u'Receive Time', u'Repetition Count',
       u'Low Voltage', u'Mortality', u'Iridium Command', u'Predeployment Data',
       u'Error']

  # Just include successful fixes
  df.loc[:, "collarID"] = collarID + animalID
  df.loc[:, "study"] = "canada"
  df.loc[:, "pigID"] = "canada" + collarID + animalID
  df.loc[:, "dop"] = np.nan
  df.loc[:, "dtpigID"] = df.pigID + df.datetime
  df.loc[:, "datetime"] = pd.to_datetime(df.datetime, format="%Y.%m.%d %H:%M:%S")
  dfs.append(df)

# TODO: What to do with duplicate times
canada_dat = pd.concat(dfs)
canada_dat = canada_dat[canada_dat.fixtype == "Succeeded"]
canada_dat = canada_dat.drop_duplicates("dtpigID", keep='first')

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
canada_red = canada_dat[cols]
kilgo_red = kilgo_trun[cols]
move_red = move[cols]

comb_dat = pd.concat([full_dat_trun, judas_red, tejon_red, cali_red, contact_red, 
                        canada_red, kilgo_red, move_red])

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

# Make sure all date times are formatted correctly
comb_dat.loc[:, "datetime"] = pd.to_datetime(comb_dat.datetime)
comb_dat.loc[:, "longitude"] = comb_dat.longitude.astype(np.float)
comb_dat.loc[:, "latitude"] = comb_dat.latitude.astype(np.float)


# Clean up errant points in Michigan after converting lat long to floats
comb_dat = comb_dat[~(((comb_dat.latitude < 43.75) | (comb_dat.longitude < -84.6)) & (comb_dat.study == "michigan"))]


# Reformat all Fixes

# Set all NULL fixes to 3D. An assumption, but given these were cleaned studies
# not unreasonable
comb_dat.loc[comb_dat.fixtype.isnull(), "fixtype"] = "3D"
comb_dat.loc[:, "fixtype"] = comb_dat.fixtype.str.replace(".*3D.*", "3D").str.replace(".*2D.*", "2D").str.replace("Succeeded", "3D")

# Remove the 2D fixes
comb_dat = comb_dat[comb_dat.fixtype == "3D"]

# Save the resulting pig data
comb_dat.to_csv("full_pig_data.csv", index=False)

################################################
### Link the attributes pigIDs and study IDs ###
################################################

attrib = pd.read_csv("../archival/all_pig_attributes.csv")
#attrib = pd.read_excel("../archival/all_pig_attributes_1_25_18.xlsx", sheetname="Sheet1")
attrib.replace({'Study' : {'SREL-Contact' : 'srel_contact', 
                           'Judas' : 'judas_pig',
                           'Camp Bullis' : 'txcamp',
                           'SREL-vacuum' : 'srel_vacuum',
                           'SRS' : 'srs_kilgo',
                           'Michigan': 'michigan',
                           'Movement': "movepig",
                           'Tejon': 'tejon',
                           'CA': 'cali',
                           'Florida': 'florida', 
                           'Canada': 'canada'}}, inplace=True)

attrib.columns = [u'collarID', u'sex', u'age_class', u'weight_lb', u'start_date',
       u'end_date', u'study', u'notes', 'smc_notes', 'cleaned']

# Reformatting movepig study
ind = attrib.study == "movepig"
attrib.loc[ind, "collarID"] = attrib.collarID[ind].str.lower()

# Replace ID names in attributes
for idm in idmapper.iterkeys():
	attrib.loc[ind, "collarID"] = attrib.collarID[ind].str.replace(idm, idmapper[idm])

attrib.loc[ind, "collarID"] = attrib.collarID[ind].str.replace("_ ", "_").str.replace("__", "_")

for let in string.lowercase:
	attrib.loc[ind, "collarID"] = attrib.collarID[ind].str.replace("{0}_ho".format(let), "{0}ho".format(let))

attrib.loc[ind, "collarID"] = attrib.collarID[ind].str.replace(r"([0-9]+)_([a-d])", r"\1\2")

attrib.loc[ind, "study"] = ["_".join(v[:-1]) for v in attrib[ind].collarID.str.split("_")]
attrib.loc[ind, "collarID"] = [v[-1] for v in attrib[ind].collarID.str.split("_")]

#comb_dat.groupby('study')['pigID'].unique()['cali']
attrib.loc[:, "pigID"] = attrib.study + attrib.collarID.astype("S21")

# Clean up weights
attrib.loc[:, "weight_lb"] = attrib.weight_lb.str.replace(".*\(estimated\)", "").str.strip("+").str.strip("-").str.replace(" ", "")
attrib = attrib.replace('', np.nan, regex=True)
attrib.loc[:, "weight_lb"] = attrib.weight_lb.astype(np.float)

# Clean up sexes
attrib.sex[attrib.sex == "male"] = "M"
attrib.sex[attrib.sex == "female"] = "F"

attrib.to_csv("pig_attributes.csv", index=False)

# Check that all IDs match...all IDs are matching!
# gb = comb_dat.groupby('study') 

# for study in comb_dat.study.unique():

#   cvals = gb['pigID'].unique()[study]
#   atvals = attrib.loc[attrib.study == study, "pigID"]

#   print(study)
#   print(set(cvals) - set(atvals))


##########################################
### Make a summary file for each study ###
##########################################

study_sum = comb_dat.groupby("study").agg({'latitude' : {'min': np.min, 'max': np.max, 'mean': np.mean}, 
                               'longitude' : {'min': np.min, 'max': np.max, 'mean': np.mean},
                               'pigID': {'num_pigs': lambda x: len(np.unique(x)),
                                         'num_fixes': len},
                                'datetime': {'mindate': np.min, 
                                             'maxdate': np.max}})

study_sum.columns = ['_'.join(col).strip() for col in study_sum.columns.values]

# TODO: Update as necessary. Add time zone information for each study
tzs = {'cali0': "America/Los_Angeles", 
 'cali1': "America/Los_Angeles", 
 'cali2': "America/Los_Angeles", 
 'cali3': "America/Los_Angeles", 
 'cali4': "America/Los_Angeles", 
 'canada': "America/Regina", 
 'fl_raoul': "America/New_York",
 'florida': "America/New_York", 
 'ga_bill': "America/New_York",
 'ga_steve': "America/New_York",
 'judas_pig': "UTC",
 'la_hartley': "America/Chicago",
 'la_steve': "America/Chicago", 
 'michigan': "America/New_York", 
 'mo_kurt0': "America/Chicago",  
 'mo_kurt1': "America/Chicago",
 'mo_kurt2' : "America/Chicago",
 'mo_kurt3' : "America/Chicago", 
 'mo_kurt4' : "America/Chicago",
 'mo_kurt5': "America/Chicago", 
 'mo_kurt6': "America/Chicago", 
 'mo_kurt7': "America/Chicago",
 'mo_kurt8': "America/Chicago", 
 'sc_jim2': "America/New_York", 
 'scjim1' : "America/New_York", 
 'srel_contact': "America/New_York", 
 'srel_vacuum': "America/New_York",
 'srs_kilgo': "America/New_York",
 'tejon': "America/Los_Angeles", 
 'tx_christy_14r': "America/Chicago", 
 'tx_christy_15r': "America/Chicago",
 'tx_susan': "America/Chicago", 
 'tx_tyler_k1': "America/Chicago", 
 'tx_tyler_w1': "America/Chicago", 
 'tx_tyler_w2': "America/Chicago", 
 'txcamp': "America/Chicago"}

tzsdt = pd.DataFrame([(s, t) for s, t in tzs.items()])
tzsdt.columns = ['study', 'tz']
study_sum = study_sum.join(tzsdt.set_index("study"))

study_sum.to_csv("study_summary.csv")

