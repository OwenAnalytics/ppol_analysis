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
from obspy.core.inventory.station import Station
from obspy.core.inventory.network import Network
from obspy.io.sac.sacpz import attach_paz
import seaborn as sns
from obspy.geodetics import gps2dist_azimuth

sns.set_style("white")


# Set following parameters (uncommand whatever needed)
evi = 'single' #Particle motion analysis on single station
#evi = 'all' #Particle motion analysis on whole event catalog
ns = 7 #Station numer
fmi = 6 #Lower corner frequency of bandpass filter
fma = 12 #Upper corner frequency of bandpass filter
taille = 0.3 #Length of analysis window
d_con = 0.2 #Distance condition

diri = 'quakes'



#Create network with station coordinates
LSb1=Station('LSb1',37.32167,-32.27771,-1859)
LSb2=Station('LSb2',37.29850,-32.32375,-2099)
LSb3=Station('LSb3',37.29178,-32.28796,-1776)
LSb4=Station('LSb4',37.28270,-32.24102,-2097)
LSb7=Station('LSb7',37.26059,-32.29746,-1939)
my_net=Network('4G',[LSb1,LSb2,LSb3,LSb4,LSb7])

#Read catalog with events bigger M = 0.8
if diri == 'quakes':
    cat = read_events('quake08.xml')
    wav_path = 'LSFIR_unrot'


#Creat list of events with either manual input of event number or whole catalog
if evi == 'single':
	evis = [int(input('Event: '))]
else:
	evis = np.arange(len(cat))

#Set station coordinates from network
if ns == 7:
    sta_lon = my_net.stations[4].longitude
    sta_lat = my_net.stations[4].latitude
else:
    sta_lon = my_net.stations[ns - 1].longitude
    sta_lat = my_net.stations[ns - 1].latitude

#For every chosen event (single or catalog)
for i in evis:
    ne = i    
    event= cat[ne] #Get event
    if event.magnitudes == []:
        mag = 1
    else:
        mag = event.magnitudes[0].mag
    event_lat = event.origins[0].latitude #Get event location
    event_lon = event.origins[0].longitude
    depth = event.origins[0].depth / 1000 
    epi_dist, az, baz = gps2dist_azimuth(event_lat, event_lon, sta_lat, sta_lon)
    epi_dist = epi_dist / 1000
    dist = obspy.geodetics.kilometer2degrees(epi_dist, radius = 6371) 
    #dist = locations2degrees(sta_lat, sta_lon, event_lat, event_lon) #Calculate distance
    #print(dist)
    if dist > d_con: #Choose only events which are farther than 0.02 degrees away

        '''
        READ AND PROCESS DATA STREAM
        '''

        et = event.origins[0].time #Get eventtime
        #print('Event time: %s' %et)
        path = '../WAV/%s/%d/%02d/%d-%02d-%02d-*'\
                %(wav_path, et.year, et.month, et.year, et.month, et.day) #Read in corresponding waveform

        #Find file with corresponding file name
        for fname in glob.glob(path):
            stream = obspy.read(fname)
            if stream[0].stats.starttime > et-40 and stream[0].stats.starttime < et:
                st1 = stream
                #print('Filename: %s' %fname)

        #Remove instrument response
        for tr in st1:
        	attach_paz(tr, '/Volumes/LS_M2/PpolPackage-master/response/SAC_PZs_XX_SPOBS2_SHZ_00_2001.001.00.00.00.0000_2021.001.24.60.60.99999')
        st1.simulate(paz_remove = tr.stats.paz)

        #Preprocessing
        st1.detrend('demean')
        st1.detrend('simple')
        st1.taper(max_percentage = 0.01, type = 'hann')
        st1.filter('bandpass', freqmin = fmi, freqmax = fma)

        #Select traces from specific station
        st = obspy.core.stream.Stream()        
        for tr in st1.select(station = 'LSb%s' %ns):
        	st.append(tr)

        #print(st)

        #Assign traces
        tr_z = st.select(channel = 'SHZ')[0]
        tr_n = st.select(channel = 'SHY')[0]
        tr_e = st.select(channel = 'SHX')[0]

        #Read in pick on specific station
        for pick in event.picks: #read picked times 
            pick_split = list(pick.waveform_id.channel_code)
            if pick.waveform_id.station_code == tr_z.stats.station and pick.phase_hint == 'P':
                tp = pick.time + 0.1 

        if len(tr_z.data) < 1000:
            print('Something is wrong with the trace...')
            print(ne)
        elif (et - tp) > 100:
            print('Wrong picktime')
            print(ne)
        else:
            try:
                tp
            except NameError:
                print("No picks for this event and station")
                print(ne)
            else:
                '''
                PLOT WAVEFORMS
                '''


                #Create colormap
                my_cmap = sns.cubehelix_palette(3, rot = -.4, as_cmap=True)
                xtime = np.linspace(0, tr_z.stats.npts * tr_z.stats.delta, tr_z.stats.npts)

                dt = np.abs(tp - tr_z.stats.starttime) #Time beteen stream starrtime and arrivaltime

                #PLOT FILTERED WAVEFORM
                f, axarr = plt.subplots(3, sharex = True)
                axarr[0].plot(xtime, tr_z.data, color = my_cmap.colors[250])
                axarr[0].set_ylabel('%s' %tr_z.stats.channel)
                axarr[0].axvline(dt, ymin = min(tr_z.data), ymax = max(tr_z.data), color = my_cmap.colors[20])
                axarr[0].set_title('Filtered waveform \n %s - %s' %(tr_z.stats.starttime, tr_z.stats.endtime))
                axarr[1].plot(xtime, tr_n.data, color = my_cmap.colors[250])
                axarr[1].set_ylabel('%s' %tr_n.stats.channel)
                axarr[1].axvline(dt, ymin = min(tr_z.data), ymax = max(tr_z.data), color = my_cmap.colors[20])
                axarr[2].plot(xtime, tr_e.data, color = my_cmap.colors[250])
                axarr[2].set_ylabel('%s' %tr_e.stats.channel)
                axarr[2].axvline(dt, ymin = min(tr_z.data), ymax = max(tr_z.data), color = my_cmap.colors[20])
                axarr[2].set_xlabel('time (s)')
                for i in range(len(axarr)):
                        axarr[i].set_xlim(0, tr_z.stats.npts * tr_z.stats.delta )

                #f.savefig('%s/LSb%s_figs/Event%03d_filt' %(diri,ns, ne), format = 'PS')
                #plt.show()

                #PLOT CHOSEN ANALYSIS WINDOW
                f, axarr = plt.subplots(3, sharex = True)

                for i in range(len(axarr)):
                        axarr[i].set_xlim(dt - taille, dt + 2*taille)
                        lims = axarr[i].get_xlim()
                        we = np.where( (xtime > lims[0]) & (xtime < lims[1]))
                        lim = np.max(np.abs(tr_z.data[we]))
                        axarr[i].set_ylim( -lim, lim)
                        #axarr[i].set_xlim(dt - 2*taille, dt + 5*taille)
                        #axarr[i].set_ylim( tr_z.data[we].min(), tr_z.data[i].max())

                axarr[0].plot(xtime, tr_z.data, color = my_cmap.colors[250])
                axarr[0].set_ylabel('%s' %tr_z.stats.channel)
                axarr[0].axvline(dt, ymin = 0, ymax = 1, color = my_cmap.colors[20])
                axarr[0].axvline(dt + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])
                axarr[0].set_title('Analysis window of filtered waveform')

                axarr[1].plot(xtime, tr_n.data, color = my_cmap.colors[250])
                axarr[1].set_ylabel('%s' %tr_n.stats.channel)
                axarr[1].axvline(dt, ymin = 0, ymax = 1, color = my_cmap.colors[20])
                axarr[1].axvline(dt + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])

                axarr[2].plot(xtime, tr_e.data, color = my_cmap.colors[250])
                axarr[2].set_ylabel('%s' %tr_e.stats.channel)
                axarr[2].axvline(dt, ymin = 0, ymax = 1, color = my_cmap.colors[20])
                axarr[2].axvline(dt + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])
                axarr[2].set_xlabel('time (s)')

                #plt.show()
                #f.savefig('%s/LSb%s_figs/Event%03d_window' %(diri,ns, ne))

                #ADJUST ANALYSIS WINDOW IF NECESSARY
                ch = 0
                #ch = float(input('Change in arrival time?'))

                f, axarr = plt.subplots(3, sharex = True)

                for i in range(len(axarr)):
                        axarr[i].set_xlim(dt - taille, dt + 2*taille)
                        lims = axarr[i].get_xlim()
                        we = np.where( (xtime > lims[0]) & (xtime < lims[1]))
                        lim = np.max(np.abs(tr_z.data[we]))
                        axarr[i].set_ylim( -lim, lim)
                        #axarr[i].set_xlim(dt - 2*taille, dt + 5*taille)
                        #axarr[i].set_ylim( tr_z.data[we].min(), tr_z.data[i].max())

                axarr[0].plot(xtime, tr_z.data, color = my_cmap.colors[250])
                axarr[0].set_ylabel('%s' %tr_z.stats.channel)
                axarr[0].axvline(dt + ch, ymin = 0, ymax = 1, color = my_cmap.colors[20])
                axarr[0].axvline(dt + ch + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])
                axarr[0].set_title('Analysis window of filtered waveform')

                axarr[1].plot(xtime, tr_n.data, color = my_cmap.colors[250])
                axarr[1].set_ylabel('%s' %tr_n.stats.channel)
                axarr[1].axvline(dt + ch, ymin = 0, ymax = 1, color = my_cmap.colors[20])
                axarr[1].axvline(dt + ch +taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])

                axarr[2].plot(xtime, tr_e.data, color = my_cmap.colors[250])
                axarr[2].set_ylabel('%s' %tr_e.stats.channel)
                axarr[2].axvline(dt + ch, ymin = 0, ymax = 1, color = my_cmap.colors[20])
                axarr[2].axvline(dt + ch + taille, ymin = 0, ymax = 1, color = my_cmap.colors[50])
                axarr[2].set_xlabel('time (s)')

                
                #f.savefig('%s/LSb%s_figs/Event%03d_window' %(diri,ns, ne), format = 'PS')
                #plt.show()

                #PREPARE DATA STREAM WHICH WILL BE ANALYSED
                trim_z = tr_z.trim(tp + ch, tp + ch + taille)
                trim_n = tr_n.trim(tp + ch, tp + ch + taille)
                trim_e = tr_e.trim(tp + ch, tp + ch + taille)

                #if ns == 7:
                #    continue
                #else:
                #    trim_n.data = -1 * trim_n.data
                #    trim_e.data = -1 * trim_e.data
                #trim_z.plot()
           

                #print(trim_e)

                '''
                BEGIN OF PARTICLE MOTION ANALYSIS
                '''

                #put traces in mean centered e,n,z matrix
                mat = np.column_stack((trim_e.data, trim_n.data, trim_z.data)) 

                for j in range(3):
                        mat[:,j] = mat[:,j] - np.mean(mat[:,j])

                #calculate covariance
                cov = np.dot((np.transpose(mat)), mat) * (1/trim_z.stats.npts)

                #Eigenvalues and vectors
                val, vec = np.linalg.eig(cov)


                args = np.argsort(val)[::-1]  #get sorted args
                val = np.sort(val)[::-1] #sort eigenvalues
                vecs = np.column_stack((vec[:,args[0]],vec[:,args[1]],\
                        vec[:,args[2]])) #sort eigenvectors according to eigenvalues

                #Normalize first eigenvector - should already be done by numpy
                vec1 = vecs[:,0]
                #nvec1 =  np.sqrt(vec1[0]**2 + vec1[1]**2 + vec1[2]**2) 
                #vec1 = vec1 / nvec1

                #Calculate angle of polarization

                if vec1[0] <= 0 and vec1[1] < 0:
                    Theta = (np.pi + np.arctan(vec1[0] / vec1[1])) * 180 / np.pi
                elif vec1[0] < 0  and vec1[1] > 0:
                    Theta = ( 2 * np.pi + np.arctan(vec1[0] / vec1[1] )) * 180 / np.pi
                elif vec1[0] >= 0 and vec1[1] > 0:
                    Theta = (np.arctan(vec1[0] / vec1[1] )) * 180 / np.pi
                elif vec1[0] >= 0 and vec1[1] <  0:
                    Theta = ( np.pi + np.arctan(vec1[0] / vec1[1] )) * 180 / np.pi

                #Degree of rectilinarity 3D

                rect = 1 - (val[1] + val[2]) / (2*val[0])

                #Horizontal covariance matrix
                cov_h = np.zeros((2,2))
                cov_h[0,0] = cov[0,0]
                cov_h[0,1] = cov[0,1]
                cov_h[1,0] = cov[1,0]
                cov_h[1,1] = cov[1,1]

                valh, vech = np.linalg.eig(cov_h) #eigenvalues and eigenvectors

                #Degree of horizontal rectilinarity
                CpH = 1 - np.min(valh)/np.max(valh)

                #2D azimuth
                argsh = np.argsort(valh)[::-1]  #get sorted args
                valh = np.sort(valh)[::-1] #sort eigenvalues
                vecsh = np.column_stack((vech[:,argsh[0]],vech[:,argsh[1]])) #sort eigenvectors according to eigenvalues

                #Normalize first eigenvector - should already be done by numpy

                vec1h = vecsh[:,0]
                
                if vec1h[0] >= 0 and vec1h[1] > 0:
                    ThetaH = (np.arctan(vec1h[0] / vec1h[1])) * 180 / np.pi
                elif vec1h[0] >= 0  and vec1h[1] < 0:
                    ThetaH = (np.pi + np.arctan(vec1h[0] / vec1h[1])) * 180 / np.pi
                elif vec1h[0] < 0 and vec1h[1] > 0:
                    ThetaH = (2 * np.pi + np.arctan(vec1h[0] / vec1h[1])) * 180 / np.pi
                elif vec1h[0] < 0 and vec1h[1] <  0:
                    ThetaH = ( np.pi + np.arctan(vec1[0] / vec1[1])) * 180 / np.pi
                

                #Resolution of the ambiguity between thetaH and theta
                ep1 = Theta - ThetaH
                ep2 = ThetaH - Theta

                if ep1 > 140:
                    Theta = Theta - 180
                if ep2 > 140:
                    Theta = Theta + 180

                #Signal to noise ration
                SNR = (np.max(valh) - np.min(valh)) / np.min(valh)

                #Angle of vertical polarization

                #Longitudinatl component
                dnorthr = (ThetaH + 180.) * np.pi / 180
                L = np.cos(dnorthr) * (trim_n.data) + np.sin(dnorthr) * (trim_e.data)
                #To define real value of the backazimuth????


                f, axarr = plt.subplots(2)
                #Horizontal particle motion
                axarr[0].plot(trim_e.data - np.mean(trim_e.data), trim_n.data - np.mean(trim_n.data), color = my_cmap.colors[250])
                axarr[0].set_title('Horizontal particle motion')
                axarr[0].set_ylabel('BHN')
                axarr[0].set_xlabel('BHE')
                limsx = axarr[0].get_xlim()
                limsy = axarr[0].get_ylim()
                limxy = np.max([np.max(np.abs(limsx)), np.max(np.abs(limsy))])
                axarr[0].set_ylim(-limxy, limxy)
                axarr[0].set_xlim(-limxy, limxy)
                x_lin = np.linspace(-limxy, limxy, 10)
                m = vec1h[1] / vec1h[0]
                axarr[0].plot(x_lin, m * x_lin, color = my_cmap.colors[50])
                

                #Vertical particle motion
                axarr[1].plot(L, trim_z.data, color = my_cmap.colors[250])
                axarr[1].set_title('Vertical particle motion')
                axarr[1].set_ylabel('BHZ')
                axarr[1].set_xlabel('L')
                
                limsx = axarr[1].get_xlim()
                limsy = axarr[1].get_ylim()
                limxy = np.max([np.max(np.abs(limsx)), np.max(np.abs(limsy))])
                axarr[1].set_ylim(-limxy, limxy)
                axarr[1].set_xlim(-limxy, limxy)
                #Set subfigures in quadratic form
                for ax in axarr:
                    ax.set(adjustable='box-forced', aspect='equal')
                plt.tight_layout()
                plt.show()
                f.savefig('%s/LSb%s_figs/Event%03d_motion' %(diri, ns, ne), format = 'PS')



                Vpol = np.arccos(vec1[2]) * 180 / np.pi

                if Vpol <= 90 and Vpol >= 0:
                    Vpol = Vpol
                elif vec1[0] < 0 and vec1[1] < 0  and vec1[2] > 0 and Vpol < 90:
                    Vpol = 180 - Vpol;
                elif vec1[0] > 0 and vec1[1] > 0 and vec1[2] > 0 and Vpol < 90:
                    Vpol = 180 - Vpol
                elif vec1[0] < 0 and vec1[1] > 0 and vec1[2] > 0 and Vpol < 90:
                    Vpol = 180 - Vpol


                if Vpol <= 90 and Vpol >= 0:
                    Vpol = -1 * ( Vpol - 180 )
                    Theta = Theta - 180
                    ThetaH = ThetaH - 180

                VV = np.column_stack((L - np.mean(L), trim_z.data - np.mean(trim_z.data)))
                cov_v = np.dot(np.transpose(VV), VV) * (1/trim_z.stats.npts)

                valv, vecv = np.linalg.eig(cov_v)

                CpZ = 1 - np.min(valv) / np.max(valv)

                err2d = np.rad2deg(np.arctan(np.sqrt(np.min(valh) / np.max(valh)))) #Calculate error
                errV = np.rad2deg(np.arctan(np.sqrt(np.min(valv) / np.max(valv))))
                err3d = np.rad2deg(np.arctan(np.sqrt(val[2] / (val[1] + val[0]))))

                if Theta < 0:
                    Theta += 360
                if ThetaH < 0:
                    ThetaH += 360

                #print('BAZmeas3D %s' %Theta)
                #print('BAZmeas2D %s' %ThetaH)
                #print('CpH %s' %CpH)
                #print('CpZ %s' %CpZ)
                #print('ErH2D %s' %err2d)
                #print('ErH3D %s' %err3d)
                plt.close('all')

                with open('%s/LSb%s_ppol.txt' %(diri,ns), 'a') as myfile:
                    myfile.write('%s %s %03.1f %05.2f %s %02d %02d %04.2f %05.1f %06.2f %06.2f %04.2f %04.2f %05.2f %05.2f %06.2f\n' %(st[0].stats.station, et, mag, dist, tp,  fmi, fma, taille, SNR, Theta, ThetaH, CpH, CpZ, err3d, err2d, baz))

                



