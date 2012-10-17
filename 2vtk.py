#!/usr/bin/env python
'''Convert the binary output of DynEarthSol3D to VTK files.

usage: 2vtk.py modelname [start [end]]
'''

import sys, os
import base64, zlib
import numpy as np

ndims = 2

# Save in ASCII or encoded binary.
# Some old VTK programs cannot read binary VTK files.
output_in_binary = 0


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

    for i, rec in enumerate(frames[start:end]):
        # convert from numpy.int to python int
        rec = int(rec)
        suffix = '{0:0=6}'.format(rec)
        print 'Converting frame #{}'.format(suffix)
        nnode = nnode_list[i+start]
        nelem = nelem_list[i+start]

        fvtu = open('{}.{}.vtu'.format(prefix, suffix, rec), 'w')
        vtu_header(fvtu, nnode, nelem)


        # node-based field
        fvtu.write('  <PointData>\n')
        # node number for debugging
        vtk_dataarray(fvtu, np.arange(nnode, dtype=np.int32), 'node#')

        fvtu.write('  </PointData>\n')

        # element-based field
        fvtu.write('  <CellData>\n')
        # element number for debugging
        vtk_dataarray(fvtu, np.arange(nelem, dtype=np.int32), 'elem#')

        fvtu.write('  </CellData>\n')

        # coordinate
        coord = np.fromfile(prefix+'.coord.'+suffix, dtype=np.float64, count=ndims*nnode)
        coord.shape = (nnode, ndims)
        fvtu.write('  <Points>\n')
        # VTK requires vector field (velocity, coordinate) has 3 components.
        # Allocating a 3-vector tmp array for VTK data output.
        tmp = np.zeros((nnode, 3), dtype=coord.dtype)
        tmp[:,:ndims] = coord
        vtk_dataarray(fvtu, tmp, '', 3)
        fvtu.write('  </Points>\n')

        # cell connectivity
        conn = np.fromfile(prefix+'.connectivity.'+suffix, dtype=np.int32, count=(ndims+1)*nelem)
        conn.shape = (nelem, ndims+1)
        fvtu.write('  <Cells>\n')
        vtk_dataarray(fvtu, conn, 'connectivity')
        vtk_dataarray(fvtu, (ndims+1)*np.array(range(1, nelem+1), dtype=np.int32), 'offsets')
        vtk_dataarray(fvtu, 5*np.ones((nelem,), dtype=np.int32), 'types')
        fvtu.write('  </Cells>\n')

        vtu_footer(fvtu)
        fvtu.close()
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


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print __doc__
        sys.exit(1)
    else:
        for arg in sys.argv[1:]:
            if arg.lower() in ('-h', '--help'):
                print __doc__
                sys.exit(0)
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
