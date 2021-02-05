#!/usr/bin/env python
import sys
import math

import Primitives3d
import CSGTree3d
import Grid3d

if (len(sys.argv) > 1):
    particle_filename = sys.argv[-1] + '.points'
    grid_filename = sys.argv[-1] + '.grid'
else:
    print sys.argv[0], "outprefix"
    exit(0)
print "files named"

#grid properties
#Ly = Lx = Lz = 0.4
Lx = 1.0
Ly = 1.0
Lz = 1.0
#Ne = 40
Nx = 50
Ny = 50
Nz = 1
lmpp = 8
Lz = 1.0*lmpp
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81

# free block properties
block_properties = { 'rho': 2500}
block_width = 0.2
block_height = 0.2
block_depth = 1.0
hx = Lx/Nx
block_primitive = Primitives3d.Box(0.5-block_width/2.0, 0.5+block_width/2.0,
                                 0.5-block_height/2.0, 0.5+block_height/2.0,
                                 0, 1.0,
                                 )
block_body = CSGTree3d.Node(block_primitive)
print "body created"

# add points to arrays
block_point_array = []
grid.generate_point_array(block_primitive)
block_point_array = grid.point_array

print "point array created"

# add arrays to mpm dictionary
mpm_points = []
nb = len(block_point_array) #number of material points in body

grid.write_grid_file(grid_filename)
with open(particle_filename, 'w') as f:
    f.write("%d\n" % nb)
    for p in block_point_array:
        vx = -(p.y-0.5)
        vy = (p.x-0.5)
        f.write("%g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, vx, vy, 1))

print "file written"
