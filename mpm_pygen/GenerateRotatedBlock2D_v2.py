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
R = 30.0
Lx = 0.2
Ly = 0.2
Lz = 1.0
#Ne = 40
Nx = 100
Ny = 100
Nz = 1
lmpp = 2
Lz = 1.0*lmpp
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81

# free block properties
block_properties = { 'rho': 1000 }
block_width = 0.1
block_height = 0.08
block_depth = 1.0
hx = Lx/Nx
block_primitive = Primitives3d.Box(0, block_width,
                                   0, block_height,
                                 0, 1.0,
                                 )
block_body = CSGTree3d.Node(block_primitive)
print "body created"

# empty space properties
empty_primitive = Primitives3d.Box(0, Lx,
                                   0, Ly,
                                 0, 1.0,
                                 )
empty_body = CSGTree3d.Node(empty_primitive)

# add points to arrays
block_point_array = []
#grid.generate_point_array(block_primitive)
grid.generate_point_array(empty_primitive)
block_point_array = grid.point_array

print "point array created"

# add arrays to mpm dictionary
mpm_points = []
nb = len(block_point_array) #number of material points in body

grid.write_grid_file(grid_filename)
with open(particle_filename, 'w') as f:
    f.write("%d\n" % nb)
    for p in block_point_array:
        #rotate point by R
        x = p.x
        y = p.y
        #first point
        p.x = x*math.cos(R*math.pi/180.0) + y*math.sin(R*math.pi/180.0)
        p.y = -x*math.sin(R*math.pi/180.0) + y*math.cos(R*math.pi/180.0)
        if block_primitive.encompasses(p):
            f.write("%g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, 0, 0, 1))
        #second point
        p.x = -x*math.cos(R*math.pi/180.0) + y*math.sin(R*math.pi/180.0)
        p.y = x*math.sin(R*math.pi/180.0) + y*math.cos(R*math.pi/180.0)
        if block_primitive.encompasses(p):
            f.write("%g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, 0, 0, 1))

print "file written"
