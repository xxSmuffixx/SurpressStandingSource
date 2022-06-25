#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates conventional Beamformer
"""

# ------------------------------- Imports -------------------------------
from os import path, getcwd, environ
environ["QT_API"] = "pyqt5"
from acoular import config, PowerSpectra, SteeringVector, BeamformerBase, MicGeom, RectGrid, __file__, TimeSamples
from numpy import save
import time as t

# ------------------------------- Initialize and set variables -------------------------------
config.global_caching = "none"

#------ Filename and folder paths
h5Folder = path.join(getcwd(), 'micData', 'Simuliert')
h5savefileName = "S_rpm500.h5"      # Name of the h5 file, in which the microph. data is stored
micgeoName = 'tub_vogel63.xml'  # filename of the used microphone geometry
micgeofile = path.join(path.split(__file__)[0], 'xml', micgeoName)
mapFolder = path.join(getcwd(), 'Maps')     # Folder in which the generated maps should be stored

# ------ Geometry
inc = 0.01          # Size of the map pixels in cm
Z = 1               # Distance: microphone array <-> sound source   (in meter)
radius = 0.8        # radius of the map in m
g = RectGrid(x_min=-radius, x_max=radius, y_min=-radius, y_max=radius, z=Z, increment=inc)

# ------ Sound source analysis
bSize = 1024        # block size for calculating the csm (should be a power of 2)
freq = 1500         # frequency of interest
bandwidth = 3       # bandwidth (0 = single frequency line, 3= third octave band)
c0 = 343            # Speed of Sound,for h5 files written with "WriteH5File.py"

# ------ Microphone, Generator
mg = MicGeom(from_file=micgeofile)
microdataFileName = path.join(h5Folder, h5savefileName)
microphonedata = TimeSamples(name=microdataFileName)
splittedMapName = h5savefileName[0:-3] + "freq" + str(int(freq)) + "_notRotated"

# ----- steering vector
steer = SteeringVector(grid=g, mics=mg)
T = t.time()
# ------------------------------- Calculate -------------------------------
psRot = PowerSpectra(time_data=microphonedata, window='Hanning', overlap='50%', block_size=bSize)
bfRot = BeamformerBase(freq_data=psRot, steer=steer, r_diag=True)
map = bfRot.synthetic(freq, bandwidth)   # calculates the sound maps
save(path.join(mapFolder, splittedMapName), map)

D = t.time() -T
print(f"Time elapsed: {D}")
print(f"Maps are saved in {mapFolder}")
