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
Lx = 1.0
Ly = 0.2
Lz = 1.0
#Ne = 40
Nx = 50
Ny = 10
Nz = 50
lmpp = 2
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
grid2 = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81

# free block properties
block_properties = { 'rho': 1500.0 }
block_width = 0.1#*Lx
block_height = 0.1#*Lz
block_depth = 0.1#*Ly

block_primitive = Primitives3d.Box(#(Lx-block_width)/2, (Lx+block_width)/2,
                                 block_width,block_width+block_width,
                                 (Ly-block_depth)/2, (Ly+block_depth)/2,
                                 #0,block_height,
                                 (Lz-block_height)/2, (Lz+block_height)/2,
                                 )
block_body = CSGTree3d.Node(block_primitive)
print "body 1 created"

# add points to arrays
block_point_array = []
grid.generate_point_array(block_primitive)
block_point_array = grid.point_array

print "point array 1 created"

# free block properties
block2_properties = { 'rho': 1500.0 }
block2_width = 0.1#*Lx
block2_height = 0.1#*Lz
block2_depth = 0.1#*Ly

block2_primitive = Primitives3d.Box(Lx-2*block2_width,Lx-block2_width,
                                 (Ly-block2_depth)/2, (Ly+block2_depth)/2,
                                 (Lz-block2_height)/2, (Lz+block2_height)/2,
                                 )
block2_body = CSGTree3d.Node(block2_primitive)
print "body 2 created"

# add points to arrays
block2_point_array = []
grid2.generate_point_array(block2_primitive)
block2_point_array = grid2.point_array

print "point array 2 created"

# add arrays to mpm dictionary
mpm_points = []
nb1 = len(block_point_array) #number of material points in body 1
nb2 = len(block2_point_array)

grid.write_grid_file(grid_filename)
with open(particle_filename, 'w') as f:
    f.write("%d %d %d\n" % (nb1+nb2,nb1,nb2))
    for p in block_point_array:
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
        f.write("%g %g %g %g %g %g %g %g %g\n" % (1, block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, p.z, 1, 0, 1))
    for p in block2_point_array:
        f.write("%g %g %g %g %g %g %g %g %g\n" % (2, block2_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, p.z, -1, 0, 1))
    

print "file written"
