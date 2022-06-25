#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use this script to plot a single map inside the Maps folder
"""

from __future__ import print_function

from os import path, getcwd
import os

from matplotlib.patches import Circle
from matplotlib.pyplot import xlabel, ylabel, gca
from acoular import integrate

os.environ["QT_API"] = "pyqt5"



from acoular import L_p, RectGrid
from matplotlib.pylab import figure, imshow, colorbar, title, show
from numpy import load

npyMapFolder = path.join(getcwd(), 'Maps')
MapName = 'S_rpm500freq1500_rotated'
#sectors = [[0., 0.5, 2]]
sectors = []
radius = 0.8
Z = 1
inc = 0.01
g = RectGrid(x_min=-radius, x_max=radius, y_min=-radius, y_max=radius, z=Z, increment=inc)

fig = figure("plots")
Map = load(path.join(npyMapFolder, MapName + ".npy"))
MapDB =  L_p(Map)

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



colorbar(label="Schalldruckpegel [dB]")
for s in sectors:
       print({L_p(integrate(Map, g, s))})
show()




