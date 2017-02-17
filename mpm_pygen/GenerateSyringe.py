#!/usr/bin/env python
import sys
import math

import Primitives3d
import CSGTree3d
import Grid3d

if (len(sys.argv) > 1):
    particle_filename = sys.argv[-1] + '.particles'
    grid_filename = sys.argv[-1] + '.grid'
else:
    print sys.argv[0], "outprefix"
    exit(0)
print "files named"

#grid properties
#Ly = Lx = Lz = 0.4
Lx = 0.005
Ly = 0.014
Lz = 0.0001
#Ne = 40
Nx = 50
Ny = 140
Nz = 1
lmpp = 2
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81

# bulk properties
bulk_properties = { 'rho': 940.0 }
grain_width = Lx*10/Nx
grain_depth = 0.45*Ly
grain_height = 0.5*Lz

bulk_primitive = Primitives3d.Box(0,grain_width,
                                 (0.95*Ly)-grain_depth,(0.95*Ly),
                                 0,grain_height,
                                 )
bulk_body = CSGTree3d.Node(bulk_primitive)
print "bulk created"

# add points to arrays
bulk_point_array = []
grid.point_array = []
grid.generate_point_array(bulk_primitive)
bulk_point_array = grid.point_array

# free block properties
block_properties = { 'rho': 2700.0 }
block_width = grain_width#*Lx
block_height = 0.5*Lz
block_depth = 0.05*Ly

block_primitive = Primitives3d.Box(#(Lx-block_width)/2, (Lx+block_width)/2,
                                 0,block_width,
                                 Ly-block_depth, Ly,
                                 #0,block_height,
                                 0, grain_height,
                                 )
block_body = CSGTree3d.Node(block_primitive)
print "block created"

# add points to arrays
block_point_array = []
grid.generate_point_array(block_primitive)
block_point_array = grid.point_array


# add arrays to mpm dictionary
mpm_points = []
nb1 = len(bulk_point_array) #number of material points in body 1
nb2 = len(block_point_array)

grid.write_grid_file(grid_filename)
with open(particle_filename, 'w') as f:
    f.write("%d\n" % 2)
    f.write("%d\n" % nb1)
    f.write("%d\n" % nb2)
    pID = 0
    for p in bulk_point_array:
    #    mpm_point = {
    #                    'body': 1,
    #                    'm': bulk_properties['rho']*grid.material_point_volume,
    #                    'v': grid.material_point_volume,
    #                    'x': p.x, 'y': p.y, 'z': p.z,
    #                    'x_t':0, 'y_t':0, 'z_t':0 #,
    #                    #'sxx': 0, 'sxy': 0, 'syy': 0
    #                }
    #    if (mpm_point['body'] == 1):
    #        nb1 += 1
    #    mpm_points.append(mpm_point)
        f.write("%g %g %g %g %g %g %g %g %g %g\n" % (0, pID, bulk_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, p.z, 0, 0, 0))
        pID += 1
    pID = 0
    for p in block_point_array:
        f.write("%g %g %g %g %g %g %g %g %g %g\n" % (1, pID, block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, p.z, 0, 0, 0))
        pID += 1
print "file written"
