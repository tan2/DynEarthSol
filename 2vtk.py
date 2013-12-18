#!/usr/bin/env python
'''Convert the binary output of DynEarthSol3D to VTK files.

usage: 2vtk.py [-2 -3 -a -t -h] modelname [start [end]]

options:
    -2          save 2D data (default)
    -3          save 3D data
    -a          save data in ASCII format (default: binary)
    -t          save all tensor components (default: no component)
    -h,--help   show this help
'''

from __future__ import print_function
import sys, os
import base64, zlib
import numpy as np

# 2D or 3D data?
ndims = 2

# Save in ASCII or encoded binary.
# Some old VTK programs cannot read binary VTK files.
output_in_binary = True

# Save indivisual components?
output_tensor_components = False


def main(modelname, start, end):
    path = os.path.dirname(modelname)
    prefix = os.path.basename(modelname)
    if path:
        os.chdir(path)
    tmp = np.fromfile(prefix+'.info', dtype=float, sep=' ')
    tmp.shape = (-1, 8)

    frames = tmp[:,0].astype(int)
    nnode_list = tmp[:,5].astype(int)
    nelem_list = tmp[:,6].astype(int)

    if end == -1:
        end = len(frames)

    if ndims == 2:
        nstr = 3
        component = ('XX', 'ZZ', 'XZ')
    else:
        nstr = 6
        component = ('XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ')

    for i, rec in enumerate(frames[start:end]):
        # convert from numpy.int to python int
        rec = int(rec)
        suffix = '{0:0=6}'.format(rec)
        print('Converting frame #{0}'.format(suffix))
        nnode = nnode_list[i+start]
        nelem = nelem_list[i+start]

        filename = '{0}.{1}.vtu'.format(prefix, suffix)
        fvtu = open(filename, 'w')

        try:
            vtu_header(fvtu, nnode, nelem)

            #
            # node-based field
            #
            fvtu.write('  <PointData>\n')

            vel = np.fromfile(prefix+'.vel.'+suffix, dtype=np.float64, count=ndims*nnode)
            vel.shape = (nnode, ndims)
            if ndims == 2:
                # VTK requires vector field (velocity, coordinate) has 3 components.
                # Allocating a 3-vector tmp array for VTK data output.
                tmp = np.zeros((nnode, 3), dtype=vel.dtype)
                tmp[:,:ndims] = vel
            else:
                tmp = vel
            vtk_dataarray(fvtu, tmp, 'velocity', 3)

            force = np.fromfile(prefix+'.force.'+suffix, dtype=np.float64, count=ndims*nnode)
            force.shape = (nnode, ndims)
            if ndims == 2:
                tmp = np.zeros((nnode, 3), dtype=force.dtype)
                tmp[:,:ndims] = force
            else:
                tmp = force
            vtk_dataarray(fvtu, tmp, 'force', 3)

            temperature = np.fromfile(prefix+'.temperature.'+suffix, dtype=np.float64, count=nnode)
            vtk_dataarray(fvtu, temperature, 'temperature')

            #bcflag = np.fromfile(prefix+'.bcflag.'+suffix, dtype=np.int32, count=nnode)
            #vtk_dataarray(fvtu, bcflag, 'BC flag')

            # node number for debugging
            vtk_dataarray(fvtu, np.arange(nnode, dtype=np.int32), 'node#')

            fvtu.write('  </PointData>\n')

            #
            # element-based field
            #
            fvtu.write('  <CellData>\n')

            quality = np.fromfile(prefix+'.meshquality.'+suffix, dtype=np.float64, count=nelem)
            vtk_dataarray(fvtu, quality, 'mesh quality')
            plstrain = np.fromfile(prefix+'.plstrain.'+suffix, dtype=np.float64, count=nelem)
            vtk_dataarray(fvtu, plstrain, 'plastic strain')

            strain_rate = np.fromfile(prefix+'.strain-rate.'+suffix, dtype=np.float64, count=nstr*nelem)
            strain_rate.shape = (nelem, nstr)
            srII = second_invariant(strain_rate)
            vtk_dataarray(fvtu, np.log10(srII+1e-45), 'strain-rate II log10')
            if output_tensor_components:
                for d in range(nstr):
                    vtk_dataarray(fvtu, strain_rate[:,d], 'strain-rate ' + component[d])

            strain = np.fromfile(prefix+'.strain.'+suffix, dtype=np.float64, count=nstr*nelem)
            strain.shape = (nelem, nstr)
            sI = first_invariant(strain)
            sII = second_invariant(strain)
            vtk_dataarray(fvtu, sI, 'strain I')
            vtk_dataarray(fvtu, sII, 'strain II')
            if output_tensor_components:
                for d in range(nstr):
                    vtk_dataarray(fvtu, strain[:,d], 'strain ' + component[d])

            stress = np.fromfile(prefix+'.stress.'+suffix, dtype=np.float64, count=nstr*nelem)
            stress.shape = (nelem, nstr)
            tI = first_invariant(stress)
            tII = second_invariant(stress)
            vtk_dataarray(fvtu, tI, 'stress I')
            vtk_dataarray(fvtu, tII, 'stress II')
            if output_tensor_components:
                for d in range(nstr):
                    vtk_dataarray(fvtu, stress[:,d], 'stress ' + component[d])

            effvisc = tII / (srII + 1e-45)
            vtk_dataarray(fvtu, effvisc, 'effective viscosity')

            density = np.fromfile(prefix+'.density.'+suffix, dtype=np.float64, count=nelem)
            vtk_dataarray(fvtu, density, 'density')

            viscosity = np.fromfile(prefix+'.viscosity.'+suffix, dtype=np.float64, count=nelem)
            vtk_dataarray(fvtu, viscosity, 'viscosity')

            volume = np.fromfile(prefix+'.volume.'+suffix, dtype=np.float64, count=nelem)
            vtk_dataarray(fvtu, volume, 'volume')

            volume_old = np.fromfile(prefix+'.volume_old.'+suffix, dtype=np.float64, count=nelem)
            vtk_dataarray(fvtu, 1 - volume_old/volume, 'dvol')

            # element number for debugging
            vtk_dataarray(fvtu, np.arange(nelem, dtype=np.int32), 'elem#')

            fvtu.write('  </CellData>\n')

            #
            # node coordinate
            #
            fvtu.write('  <Points>\n')
            coord = np.fromfile(prefix+'.coord.'+suffix, dtype=np.float64, count=ndims*nnode)
            coord.shape = (nnode, ndims)
            if ndims == 2:
                # VTK requires vector field (velocity, coordinate) has 3 components.
                # Allocating a 3-vector tmp array for VTK data output.
                tmp = np.zeros((nnode, 3), dtype=coord.dtype)
                tmp[:,:ndims] = coord
            else:
                tmp = coord
            vtk_dataarray(fvtu, tmp, '', 3)
            fvtu.write('  </Points>\n')

            #
            # element connectivity & types
            #
            fvtu.write('  <Cells>\n')
            conn = np.fromfile(prefix+'.connectivity.'+suffix, dtype=np.int32, count=(ndims+1)*nelem)
            conn.shape = (nelem, ndims+1)
            vtk_dataarray(fvtu, conn, 'connectivity')
            vtk_dataarray(fvtu, (ndims+1)*np.array(range(1, nelem+1), dtype=np.int32), 'offsets')
            if ndims == 2:
                # VTK_ TRIANGLE == 5
                celltype = 5
            else:
                # VTK_ TETRA == 10
                celltype = 10
            vtk_dataarray(fvtu, celltype*np.ones((nelem,), dtype=np.int32), 'types')
            fvtu.write('  </Cells>\n')

            vtu_footer(fvtu)
            fvtu.close()

        except:
            # delete partial vtu file
            fvtu.close()
            os.remove(filename)
            raise
    return


def vtk_dataarray(f, data, data_name=None, data_comps=None):
    if data.dtype in (np.int32,):
        dtype = 'Int32'
    elif data.dtype in (np.single, np.float32):
        dtype = 'Float32'
    elif data.dtype in (np.double, np.float64):
        dtype = 'Float64'
    else:
        raise Error('Unknown data type: ' + name)

    name = ''
    if data_name:
        name = 'Name="{0}"'.format(data_name)

    ncomp = ''
    if data_comps:
        ncomp = 'NumberOfComponents="{0}"'.format(data_comps)

    if output_in_binary:
        fmt = 'binary'
    else:
        fmt = 'ascii'
    header = '<DataArray type="{0}" {1} {2} format="{3}">\n'.format(
        dtype, name, ncomp, fmt)
    f.write(header)
    if output_in_binary:
        header = np.zeros(4, dtype=np.int32)
        header[0] = 1
        a = data.tostring()
        header[1] = len(a)
        header[2] = len(a)
        b = zlib.compress(a)
        header[3] = len(b)
        f.write(base64.standard_b64encode(header))
        f.write(base64.standard_b64encode(b))
    else:
        data.tofile(f, sep=' ')
    f.write('\n</DataArray>\n')
    return


def vtu_header(f, nnode, nelem):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<UnstructuredGrid>
<Piece NumberOfPoints="{0}" NumberOfCells="{1}">
'''.format(nnode, nelem))
    return


def vtu_footer(f):
    f.write(
'''</Piece>
</UnstructuredGrid>
</VTKFile>
''')
    return


def first_invariant(t):
    return np.sum(t[:,:ndims], axis=1) / ndims


def second_invariant(t):
    '''The second invariant of the deviatoric part of a symmetric tensor t,
    where t[:,0:ndims] are the diagonal components;
      and t[:,ndims:] are the off-diagonal components.'''

    # second invariant: sqrt(0.5 * t_ij**2)
    if ndims == 2:
        return np.sqrt(0.25 * (t[:,0] - t[:,1])**2 + t[:,2]**2)
    else:
        a = (t[:,0] + t[:,1] + t[:,2]) / 3
        return np.sqrt( 0.5 * ((t[:,0] - a)**2 + (t[:,1] - a)**2 + (t[:,2] - a)**2) +
                        t[:,3]**2 + t[:,4]**2 + t[:,5]**2)


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    else:
        for arg in sys.argv[1:]:
            if arg.lower() in ('-h', '--help'):
                print(__doc__)
                sys.exit(0)

    if '-2' in sys.argv:
        ndims = 2
    if '-3' in sys.argv:
        ndims = 3
    if '-a' in sys.argv:
        output_in_binary = False
    if '-t' in sys.argv:
        output_tensor_components = True

    # delete options
    for i in range(len(sys.argv)):
        if sys.argv[1][0] == '-':
            del sys.argv[1]
        else:
            # the rest of argv cannot be options
            break

    modelname = sys.argv[1]

    if len(sys.argv) < 3:
        start = 0
    else:
        start = int(sys.argv[2])

    if len(sys.argv) < 4:
        end = -1
    else:
        end = int(sys.argv[3]) + 1

    main(modelname, start, end)
