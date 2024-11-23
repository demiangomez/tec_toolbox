# TEC Toolbox
## A Matlab toolbox to process GNSS data and obtain total electron content (TEC) and perform FK analysis
### Author: Demián D. Gómez
This code is the realization of the theory explain in:  
*Gómez, Demián, Robert Smalley, Charles A. Langston, Terry Wilson, Michael Bevis, Ian Dalziel, Eric Kendrick, et al. “Virtual Array Beamforming of GPS TEC Observations of Co-Seismic Ionospheric Disturbances near the Geomagnetic South Pole Triggered by Teleseismic Megathrusts.” Journal of Geophysical Research: Space Physics, January 1, 2015, 2015JA021725. https://doi.org/10.1002/2015JA021725.*

Part of the RINEX reader and orbit determination are modified from PPPH by Berkay Bahadur

To use it, add the folder to your Matlab path. Create an ini file with the following content:

```
; variables that work in all fields are:  
; {ustation} : uppercase 4 or 9 character station name  
; {lstation} : lowercase 4 or 9 character station name  
; {year}     : 4 digit year  
; {doy}      : day of year  
; {month}    : month of the year  
; {day}      : day of the month  
; {session}  : session number as established in the rinex section  
; {year2d}   : 2 digit year e.g. 2019 = 19  

[general]  
; specify the date range to run the code  
; allowed formats: yyyy_doy, yyyy/mm/dd, gpswk-gpswkday, fyear  
; must be given in order  
; use a + sign to merge two days and process as a full arc  
; this is useful when an event occurs near the day boundary  
; e.g. 2020_347, 2020_348 + 2020_349, 2020_350  
date = 2020_349

; list of stations to process  
stations = CACU, CARD, HUEC, RMLS, PARR, PEHU, ALUA, SMAN, ZPLA;, PAGL, CTRO, LRNC, CNTR

; range of slowness for the FK filter  
slowness = -15:0.5:15

; Space vehicle to use  
fk_sv = G26

; time window for the analysis  
time_window = 16:00:00, 18:00:00

; interval of your GNSS data  
interval = 5

; polynomial for detrending  
polynomial = 7

; filter in mHz to bandpass the data  
filter = 3.10,3.25

; If you have a model to detrend the data, specify here  
model = /home/demian/Dropbox/OSU/Projects/Eclipses/2020/Models/sami3_eclipse_121420_1deg.nc  
```

The execute the code, use the example called "run_tec.m"
