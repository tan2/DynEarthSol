#!/usr/bin/env python

from __future__ import print_function, unicode_literals
import sys
import numpy as np

# 2D or 3D data?
ndims = 2

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
        self.steps = list(tmp[:,1].astype(int))
        self.time = list(tmp[:,2].astype(float))
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
                          'viscosity', 'edvoldt', 'volume']):
            count = nelem
        elif name in set(['connectivity']):
            count = (self.ndims + 1) * nelem
            shape = (nelem, self.ndims+1)
        elif name in set(['coordinate', 'velocity', 'velocity averaged', 'force']):
            count = self.ndims * nnode
            shape = (nnode, self.ndims)
        elif name in set(['temperature', 'mass', 'tmass', 'volume_n']):
            count = nnode
        else:
            raise NameError('uknown field name: ' + name)

        fname = self.get_fn(frame)
        with open(fname) as f:
            f.seek(pos)
            field = np.fromfile(f, dtype=dtype, count=count).reshape(shape)
        return field


    def read_markers(self, frame, markername):
        'Read and return marker data'
        fname = self.get_fn(frame)
        with open(fname) as f:

            pos = self.field_pos[markername+' size']
            f.seek(pos)
            nmarkers = np.fromfile(f, dtype=np.int32, count=1)[0]

            marker_data = {'size': nmarkers}

            # floating point
            for name in (markername+'.coord',):
                pos = self.field_pos[name]
                f.seek(pos)
                tmp = np.fromfile(f, dtype=np.float64, count=nmarkers*self.ndims)
                marker_data[name] = tmp.reshape(-1, self.ndims)
                #print(marker_data[name].shape, marker_data[name])

            # int
            for name in (markername+'.elem', markername+'.mattype', markername+'.id'):
                pos = self.field_pos[name]
                f.seek(pos)
                marker_data[name] = np.fromfile(f, dtype=np.int32, count=nmarkers)
                #print(marker_data[name].shape, marker_data[name])

        return marker_data


