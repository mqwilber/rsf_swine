# Schedule and notes for 1/22 - 1/29 RSF meeting

Preliminary schedule
---

## Monday

**Objectives**: Data cleaning, research questions, and intro to R coding

- Meet at 9:00 am at CSU biology building
- Discuss potential research questions that we are interested in getting out of this project
- Look at the data
	- Ideally, I would like to have a script that takes in all of the raw data and converts it to the "master" file.  This is the file on which we will be performing the various analyses.
- Discuss challenges with the data.
	- Where do we want to store the data? FTP? Box?
	- We can put less sensitive material (docs and stuff) in a shared Dropbox
	folder or something.  That way we can both actively edit.
- Start formatting the data.
	- Compile a single file that has all of the data description in a single place.  This file should include any idiosyncracies associated with the files.  We can discuss these together and figure out who we need to email to figure out more about the data.
	- Another file that has questions that we will be asking various owners of the datasets.
	
*Data challenges*

1. Ensure locations and times are in the same format across datasets.
2. Make separate files that include pig characteristics and site characteristics. 
	- This could be most easily done with an SQL/SQLlite database. But that might be overkill.
3. Make sure each pig has a unique ID (study)

**Goal**: Have comparable data across each different study site

I think we would like to have three datatables in our database

1. The GPS datatable: Will contain all of the GPS tracking data 

| pigID | CollarNumber | Study | LMT_time | LMT_date | lat | long | DOP |


- All date and times should be in the same format
- All locations should be lats and longs


2. Pig attributes: Will contain attributes of each pig in the study

| pigID | Study | weight | sex | other |


3. Study site attributes: Will contain attributes of the specific study

| Study | Extent | Duration | meanLat| meanLong | 


These could be stored as separate csv's but ideally they would be SQLlite datatables, so we didn't always have to have all of the data in memory.


**R coding**

Go through the first part of some of the Software Carpentry R lessons to help with familiarity and data management.

### Notes from Monday

Spent the majority of the day getting the data prepped for analysis.  All of the data is ready EXCEPT for SREL Jim data which we will finish formatting today.

Also, worked on a bit of R coding for the last 30 minutes or so.


## Tuesday

**Objectives**: Finish data cleaning, discuss questions, and more R coding

## Wednesday

**Objectives**: Discuss modeling approaches, discuss relevant covariates, search for and identify potential covariates, more R coding

## Thursday

**Objectives**: Meet with larger group to discuss covariates, clearly define tasks

## Friday

**Objectives**: TBD




