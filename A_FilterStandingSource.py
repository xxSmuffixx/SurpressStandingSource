#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code sets for all intervalls which are outside of a angle interval microphonedata to 0
As a result, the script always only calculates the sound volume maps for angles within the angle interval
The trigger data is loaded from a diffrent file than the microphonedata
"""
# ------------------------------- Imports -------------------------------
from os import path, getcwd, environ

environ["QT_API"] = "pyqt5"

import time as t
from math import pi

from acoular import config, PowerSpectra, SteeringVector, BeamformerBase, MicGeom, RectGrid, __file__, RotatingFlow, \
    GeneralFlowEnvironment, MaskedTimeSamples
from acoular.tprocess import Trigger, AngleTracker

from numpy import save, array

from B_Acoular_SourceAndTprocess import SpatialInterpolatorRotationZeroing, ZeroedMaskedTimeSamples
from C_FilterFunctionality import getIntervallsByDegreeOverlapping, invertIntervalls, plotMaps, integrateSources,\
    averageAndMinimum

# ------------------------------- Initialize and set variables -------------------------------
config.global_caching = "none"

# ------ Filenames and folder paths
h5Folder = path.join(getcwd(), 'micData', 'Simuliert')
h5savefileName = "R_rpm500.h5"      # Name of the h5 file, in which the microph. data is stored
triggerFileName = 'trigger_rpm500.h5'
micgeoName = 'tub_vogel63.xml'  # filename of the used microphone geometry
micgeofile = path.join(path.split(__file__)[0], 'xml', micgeoName)
mapFolder = path.join(getcwd(), 'Maps')     # Folder in which the generated maps should be stored

# ------ Geometry
inc = 0.01          # Size of the map pixels in m
nrIntervalls = 10   # Sets in how many pieces the Maps will be divided
angleRes = 2 * pi * (1 / nrIntervalls)  # Size of each angle interval (e.g. = pi/2 -> ventilator gets divided in
                                        # 4 segments) distance in m between microphone array and sound source
Z = 0.991               # Distance: microphone array <-> sound source   (in meter)4
sideLenGrid = 1.6                    # length of each side of the scanning grid
g = RectGrid(x_min=-sideLenGrid/2, x_max=sideLenGrid/2, y_min=-sideLenGrid/2, y_max=sideLenGrid/2, z=Z, increment=inc)
sectors = [[0, 0.5, 0.1]]     # [a,b,c], a=[x0,y0,r0], b=... areas over which the sound volume will get integrated

# ------ Sound source analysis
bSize = 1024        # block size for calculating the csm (should be a power of 2 for FFT)
freq = 1500         # frequency of interest
bandwidth = 3       # bandwidth (0 = single frequency line, 3= third octave band)
c0 = 343            # Speed of sound For h5 files written with "WriteH5File.py"

# ------ Microphone, Trigger, Generators
mg = MicGeom(from_file=micgeofile)

tr = MaskedTimeSamples(name=path.join(h5Folder, triggerFileName))  # AngleTracker data
tr.load_timedata()
trigger = Trigger(source=tr, threshold=4)
print(f"Trigger threshold is set to {trigger.threshold}")

microdataFileName = path.join(h5Folder, h5savefileName)
microphonedata = ZeroedMaskedTimeSamples(name=microdataFileName)
microphonedata.load_timedata()

splittedMapName = h5savefileName[0:-3] +"freq"+str(int(freq)) + "_nrInts"+str(nrIntervalls)+"_bnd" + str(bandwidth)+"_"

rotdata = SpatialInterpolatorRotationZeroing(mics=mg, method='linear', array_dimension='2D', interp_at_zero=False,
                                      source=microphonedata)

angletracker = AngleTracker(trigger=trigger, source=rotdata)
rotdata.angle_source = angletracker

# ----- Rotating air, steering vector
rotfield = RotatingFlow(rpm=int(angletracker.average_rpm), v0=0,
                        origin=array((0., 0., 0.)))
envRot = GeneralFlowEnvironment(c=c0, N=1000, ff=rotfield, Om=12)
steerRot = SteeringVector(grid=g, mics=mg, env=envRot)

# ------------------------------- Calculate -------------------------------
nIntervalls = 0

for nr, intervalls in enumerate(getIntervallsByDegreeOverlapping(angleRes, angletracker)):

    rotdataZero = SpatialInterpolatorRotationZeroing(mics=mg, source=microphonedata, angle_source=angletracker,
                                                     method='linear', array_dimension='2D', interp_at_zero=False)
    rotdataZero.intsToZero = invertIntervalls(intervalls, rotdata.numsamples)  # e.g. intervalls =[[0,100] ,[1000,1200]
                # describes which intervalls microphonedata will set to 0

    psRot = PowerSpectra(time_data=rotdataZero, window='Hanning', overlap='50%', block_size=bSize)
    bfRot = BeamformerBase(freq_data=psRot, steer=steerRot, r_diag=True)
    T = t.time()
    splittedRotMap = bfRot.synthetic(freq, bandwidth)   # triggers the calculation

    print(f"[{nr+1}/{nrIntervalls}]\033[92m bfRot.synthetic time={format(t.time() - T, '.2f')}s \n\n \u001b[0m")

    save(path.join(mapFolder, splittedMapName + str(nIntervalls)), splittedRotMap)
    nIntervalls += 1

averageAndMinimum(mapFolder, splittedMapName, nrIntervalls-1, 10, nrIntervalls)
print(f"Maps are saved in {mapFolder}")

# ------------------------------- Plot results -------------------------------
plotMaps(mapFolder, splittedMapName, g, nrIntervalls, sectors)
integrateSources(sectors, mapFolder, splittedMapName, nrIntervalls - 1, g)
