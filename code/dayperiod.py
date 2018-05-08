from astral import Astral
import pytz
import pandas as pd
import numpy as np
import os


def suncalc(a, date, lat, lon):
    # Calculate relavent time peroids at date given location
    res = a.sun_utc(date, lat, lon)
    res.update({'nadir1': a.solar_midnight_utc(date, lon)})
    res.update({'nadir2': a.solar_midnight_utc(date + 
                                        pd.Timedelta(value=1, unit="D"), lon)})
    return(res)

def get_timeperiod(data, tz):
    """
    Compute the time of day period given lat, lon, and date
    
    Parameters
    ----------
    data : DataFrame, with columns
        datetime (LMT), date, latitude, longitude at a minimum
    tz : str
        A valid timezone to convert to
    
    Returns
    -------
    : Series
        Day period for each observation in data
    
    Notes
    -----
    Only calculates unique values in data to speed things up. 
    """
    tz = pytz.timezone(tz)
    a = Astral()
    
    # Drop duplicates
    no_dup = (data[~data.date.duplicated()][['date', 'latitude', 'longitude']]
                    .reset_index(drop=True))
    
    # Calculate sun periods for each non-duplicate item
    num = len(no_dup)
    timecalc = pd.DataFrame([suncalc(a, no_dup.date[i], 
                                        no_dup.latitude[i], 
                                        no_dup.longitude[i]) for i in range(num)])

    # Convert to appropriate time zone. 
    # Could be off by an hour given daylight savings
    timecalc = timecalc.apply(lambda x: pd.Index(x).tz_convert(tz)
                                                   .tz_localize(None))

    # Combine date info
    fullcalc = pd.concat((no_dup.date, timecalc), axis=1)
    dt_join = (data.set_index("date")[['datetime']]
                   .join(fullcalc.set_index("date")))
    
    suntimes = (dt_join[["nadir1", "dawn", "noon", "dusk", "nadir2"]]
                        .reset_index(drop=True))
    cumtime = dt_join.datetime.values[:, np.newaxis] > suntimes
    sumres = cumtime.sum(axis=1)

    # Categories of the day based on solar time
    sumdict = {0: 'dusk-nadir', 
               1: "nadir-dawn", 
               2: "dawn-solarNoon", 
               3: "solarNoon-dusk", 
               4: "dusk-nadir",
               5: "nadir-dawn"}
    
    # Map day periods back to data
    dayperiod = sumres.map(sumdict)
    
    return(dayperiod)

if __name__ == '__main__':

    """
    This script must be run after compute_ctmc_data.R is run, but before
    allstudy_analysis.R.  This script adds a day period column that allows
    the time of day to be comparable across populations.
    """
    
    basedir = "/Users/mqwilber/Repos/rsf_swine/results/glmdata_by_study/"
    studysum = pd.read_csv("../data/formatted/study_summary.csv")

    for studynm in ["tx_tyler_w2"]: #studysum.study:

        filepath = os.path.join(basedir, studynm + ".csv")

        if os.path.isfile(filepath):

            print("Loading {0}".format(studynm))
            tdat = pd.read_csv(filepath)

            print("Converting {0}".format(studynm))
            datetime = pd.to_datetime(tdat.t, unit="s")
            date = datetime.dt.date
            tdat = tdat.assign(datetime=datetime, date=date)
            form_tdat = tdat[['date', 'datetime', 'x.current', 'y.current']]
            form_tdat.columns = ['date', 'datetime', 'longitude', 'latitude']
            tz = studysum.loc[studysum.study == studynm, 'tz'].values[0]

            dayperiod = get_timeperiod(form_tdat, tz)
            tdat = tdat.assign(dayperiod=dayperiod)

            
            print("Saving {0}".format(studynm))
            tdat.to_csv(filepath, index=False)


        

