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
Ly = Lx = Lz = 1.0
Ne = 10
lmpp = 2
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Ne, lmpp)
print "grid created"

# global properties
g = -9.81

# free block properties
block_properties = { 'rho': 700.0 }
block_width = 0.1*Lx
block_height = 0.1*Ly
block_depth = 0.1*Lz

block_primitive = Primitives3d.Box(Lx/2, Lx/2+block_width,
                                 Ly/2, Ly/2+block_width,
                                 Lz/2, Lz/2+block_width,
                                 )
block_body = CSGTree3d.Node(block_primitive)
print "body created"

# add points to arrays
block_point_array = []
grid.generate_point_array(block_primitive)
block_point_array = grid.point_array

#for p in grid.point_array:
#    if block_body.evaluate(p):
#        block_point_array.append(p)
print "point array created"

# add arrays to mpm dictionary
mpm_points = []
nb1 = len(block_point_array) #number of material points in body 1
nb2 = 0

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
        f.write("%g %g %g %g %g %g %g %g %g\n" % (1, block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, p.z, 0, 0, 0))

print "file written"
