#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use this script to plot the three map (standing, average, min) and the integrated map
"""
from os import path, getcwd

from acoular import RectGrid, MicGeom

from C_FilterFunctionality import plotMaps, integrateSources


# ------ Filenames and folder paths
splittedMapName = "R_rpm500freq1500_nrInts10_"      # Name of the h5 file, in which the microph. data is stored
#splittedMapName = "2021-03-17_11-42-31_270037freq1500_nrInts10_"

mapFolder = path.join(getcwd(), 'Maps')     # Folder in which the generated maps should be stored

# ------ Geometry
inc = 0.01          # Size of the map pixels in m
Z = 1               # Distance: microphone array <-> sound source   (in meter)
sideLenGrid = 1.6                    # length of each side of the scanning grid
g = RectGrid(x_min=-sideLenGrid/2, x_max=sideLenGrid/2, y_min=-sideLenGrid/2, y_max=sideLenGrid/2, z=Z, increment=inc)
#sectors = [[-0.3, 0.07, 0.1], [0,0.18,0.1]] # Zeile old
#sectors = [[-0.3, 0.07, 0.1], [0.23,0,0.05]] # Zeile new
#sectors = [[-0.165,-0.47,0.1],[0.165,0.425,0.1], [0.5,0,0.1]] #RRS
#sectors =[[0.5,0,0.1]] # S
#sectors = [[-0.165,-0.47,0.3]]#R
#sectors =[[0,0,2]] # R Fulls
sectors= [[-0.165,-0.47,0.3]]
#sectors =[]

plotMaps(mapFolder, splittedMapName, g, 10, sectors)
integrateSources(sectors, mapFolder, splittedMapName,9,g)
