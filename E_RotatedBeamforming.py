#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculates virtually rotated beamformer
"""
# ------------------------------- Imports -------------------------------

from os import path, getcwd, environ

environ["QT_API"] = "pyqt5"

import time as t
from math import pi

from acoular import config, PowerSpectra, SteeringVector, BeamformerBase, MicGeom, RectGrid, __file__, RotatingFlow, \
    GeneralFlowEnvironment, MaskedTimeSamples
from acoular.tprocess import Trigger, AngleTracker, SpatialInterpolatorRotation

from numpy import save, array



# ------------------------------- Initialize and set variables -------------------------------

config.global_caching = "none"

#------ Filenames and folder paths
h5Folder = path.join(getcwd(), 'micData', 'Simuliert')
h5savefileName = "S_rpm500.h5"      # Name of the h5 file, in which the microph. data is stored
triggerFileName = 'trigger_rpm500.h5'
micgeoName = 'tub_vogel63.xml'  # filename of the used microphone geometry
micgeofile = path.join(path.split(__file__)[0], 'xml', micgeoName)
mapFolder = path.join(getcwd(), 'Maps')     # Folder in which the generated maps should be stored

# ------ Geometry
inc = 0.01          # Size of the map pixels in cm?

Z = 1               # Distance: microphone array <-> sound source   (in meter)
print(f"Z= {Z}m!!")
radius = 0.8                    # radius of the map in m
g = RectGrid(x_min=-radius, x_max=radius, y_min=-radius, y_max=radius, z=Z, increment=inc)

# ------ Sound source analysis
bSize = 1024        # block size for calculating the csm (should be a power of 2)
freq = 1500         # frequency of interest
bandwidth = 3       # bandwidth (0 = single frequency line, 3= third octave band)
c0 = 343            # Speed of Sound,for h5 files written with "WriteH5File.py"
print(f"WARNING: c0 (speed of sound) set to {c0}. You may want to change that in 'E_RotatedBeamforming.py'\n")


# ------ Microphone, Trigger, Generators
mg = MicGeom(from_file=micgeofile)

tr = MaskedTimeSamples(name=path.join(h5Folder, triggerFileName))  # AngleTracker Daten
tr.load_timedata()
trigger = Trigger(source=tr, threshold=4)
print(f"Trigger threshold is set to {trigger.threshold}")

microdataFileName = path.join(h5Folder, h5savefileName)
microphonedata = MaskedTimeSamples(name=microdataFileName)

splittedMapName = h5savefileName[0:-3] + "freq" + str(int(freq)) + "_rotated"  #

rotdata = SpatialInterpolatorRotation(mics=mg, method='linear', array_dimension='2D', interp_at_zero=False,
                                      source=microphonedata)

angletracker = AngleTracker(trigger=trigger, source=rotdata)
rotdata.angle_source = angletracker


# ----- Rotating air, steering vector
rotfield = RotatingFlow(rpm=int(angletracker.average_rpm), v0=0,
                        origin=array((0., 0., 0.)))
envRot = GeneralFlowEnvironment(c=c0, N=1000, ff=rotfield, Om=12)
steerRot = SteeringVector(grid=g, mics=mg, env=envRot)

# ------------------------------- Calculate -------------------------------


psRot = PowerSpectra(time_data=rotdata, window='Hanning', overlap='50%', block_size=bSize)
bfRot = BeamformerBase(freq_data=psRot, steer=steerRot, r_diag=True)
T = t.time()
splittedRotMap = bfRot.synthetic(freq, bandwidth)   # calculates the sound maps

save(path.join(mapFolder, splittedMapName), splittedRotMap)

print(f"Maps are saved in {mapFolder}")
