1. ppol_analysis_v0 and v1 are essentially the same code. v0 is the version I
used for my Masterthesis. v1 is a slightly cleaned up version which has not
been tested (I do not have the data anymore). 
Hence I am sure that there are a few errors, but they should be
easy to fix.

2. ppol_analysis_v1 uses two functions from file ppol_functions.py
The two functions plot the waveforms (plot_waveforms) to check wheter the
analysis window is OK, and perfom the particle motion analysis after Fontaine
et al. (2009) 

3. Paths need to be changed (response file, WAV directory, catalog, ...)

4. The event catalog quake08.xml was converted from a SEISAN event file and
contains all events with M > 0.8 

5. The response file was provided by Wayne Crawford (should be stored somewhere on
LS_M2/PpolPackage-master/response/SAC_PZs_XX_SPOBS2_SHZ_00_2001.001.00.00.00.0000_2021.001.24.60.60.99999
for those having access to the harddrive)

6. Make sure to only use events that are within the station network (all
events in quake08.xml should be within) as the location error is too big for
any events outside the array

7. Whales.xml is an event catalog with different whale songs that were located
with SEISAN 

8. The results from the particle motion analysis are writte in a file
containing all important information. This can then be used to plot the
difference between backzimuth from the particle motion and the backazimuth
from the event/station location for all events

9. The mean for all events can than be used rotate the seismograms and use the
amplitudes to improve the focal mechanisms. 

10. For a more detailed description, please read my Masterthesis (also in
directory

11. I am happy to answer any further questions! 
Contact me via mwenner@vaw.baug.ethz.ch

References:
Fontaine, F. R., Barruol, G., Kennett, B. L., Bokelmann, G. H. & Reymond, D. (2009), ‘Upper mantle anisotropy beneath australia and tahiti from p wave polarization: Implications for real-time earthquake location’, Journal of Geophysical Research: Solid Earth 114(B3).
Scholz, J.-R., Barruol, G., Fontaine, F. R., Sigloch, K., Crawford, W. C. & Deen, M. (2017), ‘Orient- ing ocean-bottom seismometers from p-wave and rayleigh wave polarizations’, Geophysical Journal International 208(3), 1277–1289.
