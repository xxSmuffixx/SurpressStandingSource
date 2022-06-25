#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
from os import path, getcwd
from matplotlib.patches import Circle
from matplotlib.pyplot import xlabel, ylabel, gca
os.environ["QT_API"] = "pyqt5"
from acoular import L_p, RectGrid
from matplotlib.pylab import figure, imshow, colorbar, show
from numpy import load

npyMapFolder = path.join(getcwd(), 'Maps')
MapName = 'R_rpm500freq1500_notRotatedRot'
t = "R_rpm500"                     # title
sectors = [[-0.165,-0.47,0.3]]     # R
radius = 0.8
Z = 1
inc = 0.01
g = RectGrid(x_min=-radius, x_max=radius, y_min=-radius, y_max=radius, z=Z, increment=inc)
fig = figure("plots")
MapDB = L_p(load(path.join(npyMapFolder, MapName + ".npy")))
averageMapDBmx = MapDB.max()
print(f"averageMapDBmx = {averageMapDBmx} dB")
DBRange = 30
imshow(MapDB, vmax=averageMapDBmx + 2, vmin=averageMapDBmx-DBRange, interpolation='nearest',
       extent=g.extend(),
       origin='lower')
xlabel("x [m]")
ylabel("y [m]")
a = gca()
for s in sectors:
       circle = Circle((s[0], s[1]), radius=s[2], fill=False, linestyle='--', color='w', linewidth=2.5)
       a.add_patch(circle)

title = t

colorbar(label="Schalldruckpegel [dB]")
show()
