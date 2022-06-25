#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from os import path
import h5py
from acoular import MicGeom, WNoiseGenerator, Trajectory, MovingPointSource, \
    SourceMixer, TimeCache, WriteH5, TimeSamples, PointSource, config
from acoular import __file__ as bpath
from numpy import pi, linspace, cos, sin, int32, arange, \
    zeros, less, string_, newaxis
from scipy.signal import argrelextrema

# ----- Initialize Variables -----
config.global_caching = 'none'

# position of sources
z1 = z2 = z3 = 1  # in [m]
r1 = r2 = 0.5
# starting angle -- NOT angle at trigger instant
phi1 = 2 * pi * (0 / 4)
phi2 = 2 * pi * (1 / 2)

rpm = -500          # rotational rate , anti - clockwise
rps = rpm / 60.
delta_t = 1. / abs(rps) / 16.0

tmax = 2            # measurement time in [s]
sfreq = 51200       # sampling rate in [Hz]
nsamples = int(tmax * sfreq)    # Number of samples

micgeofile2 = path.join(path.split(bpath)[0], 'xml', 'tub_vogel63.xml')
m3 = MicGeom(from_file=micgeofile2)

# source signals ( white noise )
n1 = WNoiseGenerator(sample_freq=sfreq, numsamples=nsamples, seed=1, rms=1)  # moving source
n2 = WNoiseGenerator(sample_freq=sfreq, numsamples=nsamples, seed=2, rms=1)  # moving source
n3 = WNoiseGenerator(sample_freq=sfreq, numsamples=nsamples, seed=3, rms=0)  # standing source

# ----- export trigger signal -----
def write_trigger(trigger, case="default"):
    h5f = h5py.File("% sTimeSeriesOpt.h5 " % case, "w")
    TachoData = h5f.create_group("TachoData ")
    tachoDataV = TachoData.create_dataset("tachoDataV",
                                          shape=(nsamples, 1),
                                          chunks=(nsamples, 1))
    tachoDataV[:, :] = trigger.reshape((nsamples, 1))
    tachoDataV.attrs.create("sampleCount", nsamples, dtype=int32)
    tachoDataV.attrs.create("sampleRateHz", sfreq)
    tachoDataV.attrs.create("units", string_("V"))
    h5f.close()

# jitter
d_rpm = 1
d_rps = d_rpm / 60.
f_jitter = 0.2      # Hz

# trajectory of source
tr1 = Trajectory()
tr2 = Trajectory()

for t in arange(0, tmax * 1.001, delta_t):
    jitter = d_rps * sin(2 * pi * f_jitter * t)
    phi = t * (rps + jitter) * 2 * pi  # angle
    # define points for trajectory spline
    tr1.points[t] = (r1 * cos(phi + phi1), r1 * sin(phi + phi1), z1)
    tr2.points[t] = (r2 * cos(phi + phi2), r2 * sin(phi + phi2), z2)

# generate trigger signal -> trigger when sin ( phi ) = 0 and cos ( phi ) < 0
t = linspace(0, tmax, nsamples, endpoint=False)
jitter = d_rps * sin(2 * pi * f_jitter * t)
phi = t * (rps + jitter) * 2 * pi  # angle
sinphi = sin(phi)
cosphi = cos(phi)
# find zeros in sin phi where cos phi is negative
indmin = argrelextrema(sinphi ** 2, less)  # sin ( phi )=0 where s i n ( phi ) has minima
trigger = zeros((nsamples,))
trigger[indmin] = 5.0
trigger *= cosphi < 0  # only keep minima where cos ( phi ) <0
trigger[0] = 0
trigger[-1] = 0

# point source
p1 = MovingPointSource(signal=n1, mpos=m3, trajectory=tr1)
p2 = MovingPointSource(signal=n2, mpos=m3, trajectory=tr2)
p3 = PointSource(signal=n3, mpos=m3, loc=(0, 0.5, z3))
p = SourceMixer(sources=[p1, p2, p3])

# time data here only consists of source 1
ta3 = TimeCache(source=p)
wh52 = WriteH5(source=ta3, name='RR_rpm500.h5')
wh52.save()
tr = TimeSamples(data=trigger[:, newaxis], sample_freq=sfreq, numsamples=nsamples, numchannels=1)
wh53 = WriteH5(source=tr, name='trigger_rpm500.h5')
wh53.save()
