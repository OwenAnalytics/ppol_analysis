import obspy
from obspy import read
from obspy.io.nordic.core import read_nordic
import numpy as np
import matplotlib.pyplot as plt
from obspy.core.event import read_events
import fnmatch
import os
import glob
import matplotlib.gridspec as gridspec
import sys
from obspy.geodetics import locations2degrees
from obspy.core.inventory.station import station
from obspy.core.inventory.network import network
from obspy.io.sac.sacpz import attach_paz
import seaborn as sns
from obspy.geodetics import gps2dist_azimuth
import ppol_functions

# Set plot style
sns.set_style("white")
my_cmap = sns.cubehelix_palette(3, rot = -.4, as_cmap=true)

# Particle motion analysis on single event or whole event catalog? 
evi = 'single'  
#evi = 'all'

# Set following parameters
ns = 7  # Station numer
fmi = 6  # Lower corner frequency of bandpass filter
fma = 12  # Upper corner frequency of bandpass filter
taille = 0.3  # Length of analysis window
d_con = 0.02  # Distance condition

# Response file
resp_file = '/path/to/response/file' 
# Event catalog
catalog = '/path/to/event/catalog'
# Directory containing waveforms
path_to_wav = '/path/to/wav_directory/wav' 
# Directory where figures should be saved
path_to_fig = '/path/to/figures'
# Directory containing file with all information from ppol analysis (event, particle motion, parameters,...)
path_to_log = '/path/to/directory'

# Create network with station coordinates
lsb1=station('lsb1',37.32167,-32.27771,-1859)
lsb2=station('lsb2',37.29850,-32.32375,-2099)
lsb3=station('lsb3',37.29178,-32.28796,-1776)
lsb4=station('lsb4',37.28270,-32.24102,-2097)
lsb7=station('lsb7',37.26059,-32.29746,-1939)
my_net=network('4g',[lsb1,lsb2,lsb3,lsb4,lsb7])

# Set station coordinates from network
if ns == 7:
    sta_lon = my_net.stations[4].longitude
    sta_lat = my_net.stations[4].latitude
else:
    sta_lon = my_net.stations[ns - 1].longitude
    sta_lat = my_net.stations[ns - 1].latitude


# read catalog with events bigger m = 0.8
cat = read_events(catalog)

# Create list of events with either manual input of event number or whole catalog
if evi == 'single':
	evis = [int(input('event: '))]
else:
	evis = np.arange(len(cat))


# For every chosen event (single or catalog)
for i in evis:  
    event= cat[i]  # Get event
    if event.magnitudes == []:
        mag = 1
    else:
        mag = event.magnitudes[0].mag
    event_lat = event.origins[0].latitude  # Get event location
    event_lon = event.origins[0].longitude
    depth = event.origins[0].depth / 1000 
    # Calculate distance event - station
    epi_dist, az, baz = gps2dist_azimuth(event_lat, event_lon, sta_lat, sta_lon)
    epi_dist = epi_dist / 1000
    dist = obspy.geodetics.kilometer2degrees(epi_dist, radius = 6371) 

    # Choose only events for the analyis which are farther than 0.02 degrees away
    if dist > d_con:         
        et = event.origins[0].time #get eventtime
        path = path_to_wav + '/%d/%02d/%d-%02d-%02d-*'\
                %(et.year, et.month, et.year, et.month, et.day) #read in corresponding waveform

        # Find file with corresponding file name that fits to the eventtime
        for fname in glob.glob(path):
            stream = obspy.read(fname)
            if stream[0].stats.starttime > et-40 and stream[0].stats.starttime < et:
                st1 = stream
                #print('filename: %s' %fname)

        # Remove instrument response
        for tr in st1:
        	attach_paz(tr, path_to_response_file)
        st1.simulate(paz_remove = tr.stats.paz)

        # Preprocessing
        st1.detrend('demean')
        st1.detrend('simple')
        st1.taper(max_percentage = 0.01, type = 'hann')
        st1.filter('bandpass', freqmin = fmi, freqmax = fma)

        # Select traces from specific station
        st = obspy.core.stream.stream()        
        for tr in st1.select(station = 'lsb%s' %ns):
        	st.append(tr)

        #print(st)

        # Assign traces
        tr_z = st.select(channel = 'shz')[0]
        tr_n = st.select(channel = 'shy')[0]
        tr_e = st.select(channel = 'shx')[0]

        # Read in pick on specific station
        for pick in event.picks: # Read picked times 
            pick_split = list(pick.waveform_id.channel_code)
            if pick.waveform_id.station_code == tr_z.stats.station and pick.phase_hint == 'p':
                pick_time = pick.time
        # check if length of trace and eventtime-picktime are reasonable
        if len(tr_z.data) < 1000:
            print('something is wrong with the trace for event ' + i)
        elif (et - pick_time) > 100:
            print('wrong picktime for event ' + i)
        else:
            try:
                pick_time
            except nameerror:
                print("no picks for this station and event " + i)
            else:
                # Plot waveforms and analysis window, change window if needed
                pt_change = plot_waveforms(tr_z, tr_n, tr_e, pick_time, my_cmap, path_to_fig)  # Plot data and change pick time

                # Prepare data stream which will be analysed
                trim_z = tr_z.trim(pick_time + pt_change, pick_time + pt_change + taille)
                trim_n = tr_n.trim(pick_time + pt_change, pick_time + pt_change + taille)
                trim_e = tr_e.trim(pick_time + pt_change, pick_time + pt_change + taille)
                
                # Particle motion analysis
                SNR, Theta, ThetaH, CpH, CpZ, err3d, err2d, baz = \
                        ppol_anaylsis(trim_e.data, trim_n.data, trim_z.data, i, taille, my_cmap)

                with open(path_to_log+'/LSb%s_ppol.txt' %(diri,ns), 'a') as myfile:
                    myfile.write('%s %s %03.1f %05.2f %s %02d %02d %04.2f %05.1f %06.2f %06.2f %04.2f %04.2f %05.2f %05.2f %06.2f\n' \
                            %(st[0].stats.station, et, mag, dist, pick_time,  fmi, fma, taille, \
                            SNR, Theta, ThetaH, CpH, CpZ, err3d, err2d, baz))
