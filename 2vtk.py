#!/usr/bin/env python
# encoding: utf-8
'''Convert the binary output of DynEarthSol3D to VTK files.

usage: 2vtk.py [-a -c -t -h] modelname [start [end]]

options:
    -a          save data in ASCII format (default: binary)
    -c          save files in current directory
    -t          save all tensor components (default: no component)
    -h,--help   show this help
'''

from __future__ import print_function, unicode_literals
import sys, os
import base64, zlib
import numpy as np

# 2D or 3D data?
ndims = 2

# Save in ASCII or encoded binary.
# Some old VTK programs cannot read binary VTK files.
output_in_binary = True

# Save the resultant vtu files in current directory?
output_in_cwd = False

# Save indivisual components?
output_tensor_components = False


class Dynearthsol:
    '''Read output file of 2D/3D DynEarthSol'''

    def __init__(self, modelname):
        self.modelname = modelname
        self.read_info()
        self.read_header(self.frames[0])


    def read_info(self):
        tmp = np.fromfile(self.modelname + '.info', dtype=float, sep=' ')
        tmp.shape = (-1, 8)
        self.frames = list(tmp[:,0].astype(int))
        self.nnode_list = tmp[:,5].astype(int)
        self.nelem_list = tmp[:,6].astype(int)
        return


    def get_fn(self, frame):
        return '{0}.save.{1:0=6}'.format(self.modelname, frame)


    def read_header(self, frame):
        headerlen = 4096
        fname = self.get_fn(frame)
        with open(fname) as f:
            header = f.read(headerlen).splitlines()
            #print(header)

        # parsing 1st line
        first = header[0].split(' ')
        if (first[0] != '#' or
            first[1] != 'DynEarthSol' or
            first[2].split('=')[0] != 'ndims' or
            first[3].split('=')[0] != 'revision'):
            print('Error:', fname, 'is not a valid DynEarthSol output file!')
            sys.exit(1)

        self.ndims = int(first[2].split('=')[1])
        self.revision = int(first[3].split('=')[1])
        if self.ndims == 2:
            self.nstr = 3
            self.component_names = ('XX', 'ZZ', 'XZ')
        else:
            self.nstr = 6
            self.component_names = ('XX', 'YY', 'ZZ', 'XY', 'XZ', 'YZ')

        # parsing other lines
        self.field_pos = {}
        for line in header[1:]:
            if line[0] == '\x00': break  # end of record
            name, pos = line.split('\t')
            self.field_pos[name] = int(pos)

        #print(self.field_pos)
        return


    def read_field(self, frame, name):
        pos = self.field_pos[name]
        i = self.frames.index(frame)
        nnode = self.nnode_list[i]
        nelem = self.nelem_list[i]

        dtype = np.float64 if name != 'connectivity' else np.int32
        count = 0
        shape = (-1,)
        if name in set(['strain', 'strain-rate', 'stress', 'stress averaged']):
            count = self.nstr * nelem
            shape = (nelem, self.nstr)
        elif name in set(['density', 'material', 'mesh quality',
                          'plastic strain', 'plastic strain-rate',
                          'viscosity']):
            count = nelem
        elif name in set(['connectivity']):
            count = (self.ndims + 1) * nelem
            shape = (nelem, self.ndims+1)
        elif name in set(['coordinate', 'velocity', 'velocity averaged', 'force']):
            count = self.ndims * nnode
            shape = (nnode, self.ndims)
        elif name in set(['temperature']):
            count = nnode

        fname = self.get_fn(frame)
        with open(fname) as f:
            f.seek(pos)
            field = np.fromfile(f, dtype=dtype, count=count).reshape(shape)
        return field


def main(modelname, start, end):
    prefix = modelname
    if output_in_cwd:
        output_prefix = os.path.basename(modelname)
    else:
        output_prefix = modelname

    des = Dynearthsol(modelname)
    ndims = des.ndims

    if end == -1:
        end = len(des.frames)

    for i, frame in enumerate(des.frames[start:end]):
        des.read_header(frame)
        suffix = '{0:0=6}'.format(frame)
        print('Converting frame #{0}'.format(suffix))

        filename = '{0}.{1}.vtu'.format(output_prefix, suffix)
        fvtu = open(filename, 'wb')

        nnode = des.nnode_list[i+start]
        nelem = des.nelem_list[i+start]
        try:
            vtu_header(fvtu, nnode, nelem)

            #
            # node-based field
            #
            fvtu.write(b'  <PointData>\n')

            # averaged velocity is more stable and is preferred
            try:
                convert_field(des, frame, 'velocity averaged', fvtu)
            except KeyError:
                convert_field(des, frame, 'velocity', fvtu)

            convert_field(des, frame, 'force', fvtu)

            convert_field(des, frame, 'temperature', fvtu)
            #convert_field(des, frame, 'bcflag', fvtu)

            # node number for debugging
            vtk_dataarray(fvtu, np.arange(nnode, dtype=np.int32), 'node number')

            fvtu.write(b'  </PointData>\n')
            #
            # element-based field
            #
            fvtu.write(b'  <CellData>\n')

            convert_field(des, frame, 'mesh quality', fvtu)
            convert_field(des, frame, 'plastic strain', fvtu)
            convert_field(des, frame, 'plastic strain-rate', fvtu)

            strain_rate = des.read_field(frame, 'strain-rate')
            srII = second_invariant(strain_rate)
            vtk_dataarray(fvtu, np.log10(srII+1e-45), 'strain-rate II log10')
            if output_tensor_components:
                for d in range(des.nstr):
                    vtk_dataarray(fvtu, strain_rate[:,d], 'strain-rate ' + des.component_names[d])

            strain = des.read_field(frame, 'strain')
            sI = first_invariant(strain)
            sII = second_invariant(strain)
            vtk_dataarray(fvtu, sI, 'strain I')
            vtk_dataarray(fvtu, sII, 'strain II')
            if output_tensor_components:
                for d in range(des.nstr):
                    vtk_dataarray(fvtu, strain[:,d], 'strain ' + des.component_names[d])

            # averaged stress is more stable and is preferred
            try:
                stress = des.read_field(frame, 'stress averaged')
            except KeyError:
                stress = des.read_field(frame, 'stress')
            tI = first_invariant(stress)
            tII = second_invariant(stress)
            vtk_dataarray(fvtu, tI, 'stress I')
            vtk_dataarray(fvtu, tII, 'stress II')
            if output_tensor_components:
                for d in range(des.nstr):
                    vtk_dataarray(fvtu, stress[:,d], 'stress ' + des.component_names[d])

            convert_field(des, frame, 'density', fvtu)
            convert_field(des, frame, 'material', fvtu)
            convert_field(des, frame, 'viscosity', fvtu)
            effvisc = tII / (srII + 1e-45)
            vtk_dataarray(fvtu, effvisc, 'effective viscosity')

            # element number for debugging
            vtk_dataarray(fvtu, np.arange(nelem, dtype=np.int32), 'elem number')

            fvtu.write(b'  </CellData>\n')

            #
            # node coordinate
            #
            fvtu.write(b'  <Points>\n')
            convert_field(des, frame, 'coordinate', fvtu)
            fvtu.write(b'  </Points>\n')

            #
            # element connectivity & types
            #
            fvtu.write(b'  <Cells>\n')
            convert_field(des, frame, 'connectivity', fvtu)
            vtk_dataarray(fvtu, (des.ndims+1)*np.array(range(1, nelem+1), dtype=np.int32), 'offsets')
            if des.ndims == 2:
                # VTK_ TRIANGLE == 5
                celltype = 5
            else:
                # VTK_ TETRA == 10
                celltype = 10
            vtk_dataarray(fvtu, celltype*np.ones((nelem,), dtype=np.int32), 'types')
            fvtu.write(b'  </Cells>\n')

            vtu_footer(fvtu)
            fvtu.close()

        except:
            # delete partial vtu file
            fvtu.close()
            os.remove(filename)
            raise
    return


def convert_field(des, frame, name, fvtu):
    field = des.read_field(frame, name)
    if name in ('coordinate', 'velocity', 'velocity averaged', 'force'):
        if des.ndims == 2:
            # VTK requires vector field (velocity, coordinate) has 3 components.
            # Allocating a 3-vector tmp array for VTK data output.
            i = des.frames.index(frame)
            tmp = np.zeros((des.nnode_list[i], 3), dtype=field.dtype)
            tmp[:,:des.ndims] = field
        else:
            tmp = field

        # Rename 'velocity averaged' to 'velocity'
        if name == 'velocity averaged': name = 'velocity'

        vtk_dataarray(fvtu, tmp, name, 3)
    else:
        vtk_dataarray(fvtu, field, name)
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
    f.write(header.encode('ascii'))
    if output_in_binary:
        header = np.zeros(4, dtype=np.int32)
        header[0] = 1
        a = data.tostring()
        header[1] = len(a)
        header[2] = len(a)
        b = zlib.compress(a)
        header[3] = len(b)
        f.write(base64.standard_b64encode(header.tostring()))
        f.write(base64.standard_b64encode(b))
    else:
        data.tofile(f, sep=b' ')
    f.write(b'\n</DataArray>\n')
    return


def vtu_header(f, nnode, nelem):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<UnstructuredGrid>
<Piece NumberOfPoints="{0}" NumberOfCells="{1}">
'''.format(nnode, nelem).encode('ascii'))
    return


def vtu_footer(f):
    f.write(
b'''</Piece>
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

    if '-a' in sys.argv:
        output_in_binary = False
    if '-c' in sys.argv:
        output_in_cwd = True
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
