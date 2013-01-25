#!/usr/bin/env python

from matplotlib.pylab import *
from scipy.special import erf

ndims = 3

def readdata(frame):
    tfile = 'zzz.temperature.%06d' % frame
    xfile = 'zzz.coord.%06d' % frame
    info = fromfile('zzz.info', sep=' ')
    info.shape = (-1, 8)
    time = info[frame, 2]
    t = fromfile(tfile)
    x = fromfile(xfile)
    x.shape = (-1, ndims)
    return x, t, time


def myplot(t, x, xx, frame, time):
    diffusivity = 3./3000/1e3
    plot(t, x[:,-1]*1e-3, 'ok',
         273-(1600-273)*erf(xx/sqrt(4*diffusivity*time)), xx*1e-3, markerfacecolor='w')
    return



clf()
xx = linspace(0, -250e3)
myr2sec = 1e6 * 86400 * 365.2422
for frame in (0, 5, 14):
    x, t, time = readdata(frame)
    myplot(t, x, xx, frame, time+myr2sec)


#axis([250, 1700, -50, 0])
xlabel('Temperature (Kelvin)')
ylabel('Depth (km)')
show()

