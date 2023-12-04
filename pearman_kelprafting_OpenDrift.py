#!/usr/bin/env python

#imports
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely
import multiprocessing

from datetime import timedelta, datetime
from joblib import Parallel, delayed
from kelp_raft.models.leeway_NZ_killzone import Leeway #import custom Leeway function which terminates rafts when they enter the inferred source of origin
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_netCDF_CF_generic

#define the raftTracker function, which takes a row number as an input, uses that to get raft collection location and time, and then runs drift model backwards from that location
def raftTracker(i):
    reader_landmask = reader_global_landmask.Reader()  #read in a landmask for collisions
    raft_data = pd.read_csv('kelp_raft/raft_details.csv') #read in raft data csv file
    o = Leeway(loglevel=50) #reporting level #0 for debug; 20 for info; 50 for none
    o.add_reader(reader_landmask) # add landmask to model
    #add ocean/wind data required to run model - multiple sources for each, uses highest res available but can fall back on coarser resolution data where needed
    #*** Need to set up credentials to access each of these servers - detailed in the OpenDrift documentation
    o.add_reader(reader_netCDF_CF_generic.Reader('https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1D-m'))
    o.add_reader(reader_netCDF_CF_generic.Reader('https://my.cmems-du.eu/thredds/dodsC/global-reanalysis-phy-001-031-grepv2-daily'))
    o.add_reader(reader_netCDF_CF_generic.Reader('https://my.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy_my_0.083_P1M-m'))
    o.add_reader(reader_netCDF_CF_generic.Reader('https://my.cmems-du.eu/thredds/dodsC/global-reanalysis-phy-001-031-grepv2-monthly'))
    o.add_reader(reader_netCDF_CF_generic.Reader('http://www.psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis2/pressure'))
    o.add_reader(reader_netCDF_CF_generic.Reader('http://www.psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis2/pressure'))
    o.add_reader(reader_netCDF_CF_generic.Reader('http://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis2/Monthlies/pressure/uwnd.mon.mean.nc'))
    o.add_reader(reader_netCDF_CF_generic.Reader('http://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis2/Monthlies/pressure/vwnd.mon.mean.nc'))
    raftID = raft_data["SampleName"][i] #set info in model - taken from raft_details.csv
    raft_lon = raft_data["Long"][i] #set info in model - taken from raft_details.csv
    raft_lat = raft_data["Lat"][i] #set info in model - taken from raft_details.csv
    raft_date = raft_data["DateCollected"][i] #set info in model - taken from raft_details.csv
    raft_shp = "".join(["kelp_raft/shapefiles_2km_buffer/", raft_data["SHP_File"][i], ".shp"]) #read in shapefill to be used as a 'target zone' - i.e. the region from which molecular data suggest the raft originated from
    outfile = "".join(["kelp_raft/outputs/", raftID, ".nc"]) #set output file location
    #seed elements/particles at surface of collection site ± 2 days
    o.seed_elements(lon=raft_lon, lat=raft_lat, number=33333, radius=1000, time=[datetime.strptime(raft_date, '%d/%m/%Y') - timedelta(days = 2), datetime.strptime(raft_date, '%d/%m/%Y') + timedelta(days = 2)], origin_marker = i+1, object_type = 1) #PIW-1
    o.seed_elements(lon=raft_lon, lat=raft_lat, number=33333, radius=1000, time=[datetime.strptime(raft_date, '%d/%m/%Y') - timedelta(days = 2), datetime.strptime(raft_date, '%d/%m/%Y') + timedelta(days = 2)], origin_marker = i+1, object_type = 5) #PIW-5
    o.seed_elements(lon=raft_lon, lat=raft_lat, number=33333, radius=1000, time=[datetime.strptime(raft_date, '%d/%m/%Y') - timedelta(days = 2), datetime.strptime(raft_date, '%d/%m/%Y') + timedelta(days = 2)], origin_marker = i+1, object_type = 6) #PIW-6
    #set coastline interaction - sit on coast until moved off again #change to 'stranding' to terminate element on contact w/any coast - 'previous' to keep things bouncing.
    o.set_config('general:coastline_action', 'previous')
    o.tzone = gpd.read_file(raft_shp) #add target zone to the model
    #run model
    o.run(duration=timedelta(days=-730), time_step = timedelta(hours=-1), time_step_output=timedelta(hours=12), outfile = outfile) #run backwards for 2 years, hourly time steps but just saving the output every 12hrs
    #once model has run extract time taken for succesful particles to reach the target zone
    res = []
    for j in range(0, len(o.history)):
        if np.where(o.history['status'][j] == 1)[0] and "target" in o.status_categories:
            res.append(o.history['age_seconds'][j][np.where(o.history['status'][j] == o.status_categories.index("target"))[0][0]])
        else :
            res.append(np.nan)
    with open('kelp_raft/outputs/target_times.txt', 'a') as f: #write time take as a txt file for simple summary/analysis after the fact
        f.write(raftID)
        f.write(", ")
        f.write(", ".join(map(str, res)))
        f.write("\n")
    return(i)


#run the raftTracker in parallel to process multiple rafts at a time
inputs = range(0, 35) #all rows of raft_details.csv
num_cores = multiprocessing.cpu_count() #get number of cours
results = Parallel(n_jobs=num_cores-2)(delayed(raftTracker)(i) for i in inputs) #run in parallell - keeping 2 cores free just in case

