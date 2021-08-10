#!/usr/bin/python3
from __future__ import print_function
import sys, os
import numpy as np
#sys.path.append('/home/chaseshyu/Programs/DynEarthSol')
from Dynearthsol import Dynearthsol

def first_invariant(t):
    nstr = t.shape[1]
    ndims = 2 if (nstr == 3) else 3
    return np.sum(t[:,:ndims], axis=1) / ndims


def second_invariant(t):
    '''The second invariant of the deviatoric part of a symmetric tensor t,
    where t[:,0:ndims] are the diagonal components;
      and t[:,ndims:] are the off-diagonal components.'''
    nstr = t.shape[1]

    # second invariant: sqrt(0.5 * t_ij**2)
    if nstr == 3:  # 2D
        return np.sqrt(0.25 * (t[:,0] - t[:,1])**2 + t[:,2]**2)
    else:  # 3D
        a = (t[:,0] + t[:,1] + t[:,2]) / 3
        return np.sqrt( 0.5 * ((t[:,0] - a)**2 + (t[:,1] - a)**2 + (t[:,2] - a)**2) +
                        t[:,3]**2 + t[:,4]**2 + t[:,5]**2)


class Stuff():
    pass


def read_data(des, frame):
    stuff = Stuff()

    stuff.T = des.read_field(frame,'temperature')
    coordinate = des.read_field(frame, 'coordinate')
    stuff.x = np.array(coordinate[:,0])
    stuff.z = np.array(coordinate[:,-1])
    velocity = des.read_field(frame, 'velocity')
    stuff.vx = np.array(velocity[:,0])
    stuff.vz = np.array(velocity[:,-1])
    stuff.pls = des.read_field(frame,'plastic strain')

    stress = des.read_field(frame, 'stress')
    stuff.tI = first_invariant(stress)
    stuff.tII = second_invariant(stress)
    strain = des.read_field(frame, 'strain')
    stuff.sI = first_invariant(strain)
    stuff.sII = second_invariant(strain)
    strain_rate = des.read_field(frame, 'strain-rate')
    stuff.srI = first_invariant(strain_rate)
    stuff.srII = second_invariant(strain_rate)

    marker_data = des.read_markers(frame, markersetname)
    field = marker_data[markersetname + '.coord']
    stuff.m_x = field[:,0]
    stuff.m_z = field[:,1]
    stuff.m_id = marker_data[markersetname + '.id']
    stuff.m_mat = marker_data[markersetname + '.mattype']
    stuff.m_time = marker_data[markersetname + '.time']
    return stuff


def reldiff(oldf, newf):
    m = np.abs(oldf).max()
    diff = np.abs(newf - oldf)
    if m == 0.:
        return diff.max(), diff.std()
    else:
        return diff.max()/m, diff.std()/m


def compare(old, new):
    max, sigma = reldiff(old.T, new.T)
    print('  Temperature:\t\t', max, sigma)

    max, sigma = reldiff(old.x, new.x)
    print('  X coordinate:\t\t', max, sigma)

    max, sigma = reldiff(old.z, new.z)
    print('  Z coordinate:\t\t', max, sigma)

    max, sigma = reldiff(old.vx, new.vx)
    print('  X velocity:\t\t', max, sigma)

    max, sigma = reldiff(old.vz, new.vz)
    print('  Z velocity:\t\t', max, sigma)

    max, sigma = reldiff(old.pls, new.pls)
    print('  Pl. strain:\t\t', max, sigma)

    max, sigma = reldiff(old.tI, new.tI)
    print('  Stress I:\t\t', max, sigma)

    max, sigma = reldiff(old.tII, new.tII)
    print('  Stress II:\t\t', max, sigma)

    max, sigma = reldiff(old.sI, new.sI)
    print('  Strain I:\t\t', max, sigma)

    max, sigma = reldiff(old.sII, new.sII)
    print('  Strain II:\t\t', max, sigma)

    max, sigma = reldiff(old.srI, new.srI)
    print('  Strain rate I:\t', max, sigma)

    max, sigma = reldiff(old.srII, new.srII)
    print('  Strain rate II:\t', max, sigma)

    max, sigma = reldiff(old.m_x, new.m_x)
    print('  Marker X:\t\t', max, sigma)

    max, sigma = reldiff(old.m_z, new.m_z)
    print('  Marker Z:\t\t', max, sigma)

    max, sigma = reldiff(old.m_mat, new.m_mat)
    print('  Marker Material:\t', float(max), sigma)

    max, sigma = reldiff(old.m_time, new.m_time)
    print('  Marker Time:\t\t', max, sigma)

    return


olddir = sys.argv[1]
frame = int(sys.argv[2])

curdir = os.getcwd()
newdir = curdir

# name holder
old = 0
new = 0

modelname = 'benchmark'
markersetname = 'markerset'

try:
    # read old and new results

    os.chdir(olddir)
    des = Dynearthsol(modelname)
    old = read_data(des, frame)

    os.chdir(newdir)
    des = Dynearthsol(modelname)
    new = read_data(des, frame)

    # compare results
    print()
    print('Relative difference (max, stddev) of frame =', frame,
          ' step =', int(des.steps[frame]))
    compare(old, new)

finally:
    # restort to original directory
    os.chdir(curdir)
