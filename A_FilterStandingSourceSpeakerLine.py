#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code sets for all intervalls which are outside of a angle interval microphonedata to 0
As a result, the script always only calculates the sound volume maps for angles within the angle interval
The trigger data is for this script inside the same file with the microphonedata
"""
# ------------------------------- Imports -------------------------------
from os import path, getcwd, environ

environ["QT_API"] = "pyqt5"

import time as t
from math import pi, sqrt

from acoular import config, PowerSpectra, SteeringVector, BeamformerBase, MicGeom, RectGrid, __file__, RotatingFlow, \
    GeneralFlowEnvironment, MaskedTimeSamples
from acoular.tprocess import Trigger, AngleTracker, SpatialInterpolatorRotation

from numpy import save, array

from B_Acoular_SourceAndTprocess import SpatialInterpolatorRotationZeroing, ZeroedMaskedTimeSamples
from C_FilterFunctionality import getIntervallsByDegreeOverlapping, invertIntervalls, plotMaps, integrateSources,\
    averageAndMinimum

# ------------------------------- Initialize and set variables -------------------------------

config.global_caching = "none"

#------ Filenames and folder paths
h5Folder = path.join(getcwd(), 'micData', 'Messdaten')
h5savefileName = "2021-03-17_11-42-31_270037.h5"      # Name of the h5 file, in which the microph. data is stored
microdataFileName = path.join(h5Folder, h5savefileName)
micgeoName = 'tub_vogel64.xml'  # filename of the used microphone geometry
micgeofile = path.join(path.split(__file__)[0], 'xml', micgeoName)
mapFolder = path.join(getcwd(), 'Maps')     # Folder in which the generated maps should be stored


# ------ Geometry
inc = 0.01          # Size of the map pixels in cm?
nrIntervalls = 10   # Sets in how many pieces the Maps will be divided
angleRes = 2 * pi * (1 / nrIntervalls)  # Size of each angle interval (e.g. = pi/2 -> ventilator gets divided in
                                        # 4 segments) distance in m between microphone array and sound source
Z = 0.991               # Distance: microphone array <-> sound source   (in meter)
print(f"Z= {Z}m!!")
radius = 0.8                    # radius of the map in m
g = RectGrid(x_min=-radius, x_max=radius, y_min=-radius, y_max=radius, z=Z, increment=inc)
sectors = [[0.3, 0.3, 0.1]]     # [a,b,c], a=[x0,y0,r0], b=... areas over which the sound volume will get integrated

# ------ Sound source analysis
bSize = 1024        # block size for calculating the csm (should be a power of 2)
freq = 1500         # frequency of interest
bandwidth = 3       # bandwidth (0 = single frequency line, 3= third octave band)

temp = 16.9            # room temperature during measurement
c0 = sqrt(1.4 * 8.314462 * (273.15 + temp) / 0.02896)

startSample = 486
stopSample = 804690
trackerChannel = 64

useCorrectionFactor = False     # Multiply maps with the factor equivalent to nrIntervalls
if useCorrectionFactor:
    print("Maps are gonna be multiplied with the factor nrIntervalls! If you dont want that set 'useCorrectionFactor'"
          " to False")

# ------ Microphone, Trigger, Generators
mg = MicGeom(from_file=micgeofile)

# Micdata
microphonedata = ZeroedMaskedTimeSamples(name=microdataFileName)
microphonedata.start = startSample
microphonedata.stop = stopSample
microphonedata.invalid_channels = [trackerChannel]

# Angletracker
tr = MaskedTimeSamples(name=microdataFileName)      # AngleTracker data
invalid = list(range(trackerChannel))
tr.invalid_channels = invalid    # channels in which the microphone data is stored
tr.start = startSample
tr.stop = stopSample

trigger = Trigger(source=tr, threshold=35, trigger_type='dirac', multiple_peaks_in_hunk='first')
print(f"Trigger threshold is set to {trigger.threshold}")
angletracker = AngleTracker(source=microphonedata, trigger=trigger, rot_direction=-1, interp_points=5)

splittedMapName = h5savefileName[0:-3] + "freq" + str(int(freq)) + "_nrInts"+str(nrIntervalls)+"_"

mvirt = MicGeom()
mvirt.mpos_tot = mg.mpos_tot

rotdata = SpatialInterpolatorRotation(mics= mg, mics_virtual=mvirt, method='linear', array_dimension='2D', interp_at_zero=False,
                                      source=microphonedata, angle_source=angletracker)

# ----- Rotating air, steering vector
rotfield = RotatingFlow(rpm=int(angletracker.average_rpm), v0=0,
                        origin=array((0., 0., 0.)))
envRot = GeneralFlowEnvironment(c=c0, N=1000, ff=rotfield, Om=12)

steerRot = SteeringVector(grid=g, mics=mvirt, env=envRot)

# ------------------------------- Calculate -------------------------------
nIntervalls = 0

for nr, intervalls in enumerate(getIntervallsByDegreeOverlapping(angleRes, angletracker)):

    rotdataZero = SpatialInterpolatorRotationZeroing(mics= mg, mics_virtual=mvirt, source=microphonedata, angle_source=angletracker,
                                                     method='linear', array_dimension='2D', interp_at_zero=False)
    rotdataZero.intsToZero = invertIntervalls(intervalls, rotdata.numsamples)  # e.g. intervalls =[[0,100] ,[1000,1200]
                # describes which intervalls microphonedata will set to 0

    psRot = PowerSpectra(time_data=rotdataZero, window='Hanning', overlap='50%', block_size=bSize)
    bfRot = BeamformerBase(freq_data=psRot, steer=steerRot, r_diag=True)
    T = t.time()
    splittedRotMap = bfRot.synthetic(freq, bandwidth)   # calculates the sound maps

    print(f"[{nr+1}/{nrIntervalls}]\033[92m bfRot.synthetic time={format(t.time() - T, '.2f')}s \n\n \u001b[0m")

    save(path.join(mapFolder, splittedMapName + str(nIntervalls)), splittedRotMap)
    nIntervalls += 1

averageAndMinimum(mapFolder, splittedMapName, nrIntervalls-1, 10, nrIntervalls)
print(f"Maps are saved in {mapFolder}")

# ------------------------------- Plot results -------------------------------
plotMaps(mapFolder, splittedMapName, g, nrIntervalls, sectors)
integrateSources(sectors, mapFolder, splittedMapName, nrIntervalls - 1, g)
