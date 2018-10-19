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
Nx = 40
Ny = 40
Nz = 1
lmpp = 3
Lz = 1.0*lmpp
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81
ri = 0.25;
ro = 1.0;

# free block properties
block_properties = { 'rho': 1000 }
block_width = 0.4
block_height = 0.4
block_depth = 1.0
hx = Lx/Nx
block_include = Primitives3d.Cylinder(Primitives3d.Point(0, 0, 0.0),
                                      Primitives3d.Point(0, 0, 1.0),
                                      ro,
                                      )
block_exclude = Primitives3d.Cylinder(Primitives3d.Point(0,0,0.0),
                                      Primitives3d.Point(0,0,1.0),
                                      ri,
                                      )
block_body_i = CSGTree3d.Node(block_include)
block_body_e = CSGTree3d.Node(block_exclude)
print "body created"

# add points to arrays
block_point_array = []
grid.generate_point_array(block_include)
block_point_array = grid.point_array

print "point array created"

# add arrays to mpm dictionary
mpm_points = []
nb = len(block_point_array) #number of material points in body

grid.write_grid_file(grid_filename)
with open(particle_filename, 'w') as f:
    f.write("%d\n" % (4*nb))
    for p in block_point_array:
        if not block_exclude.encompasses(p):
            f.write("%g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, 0, 0, 1))
            f.write("%g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, -p.x, p.y, 0, 0, 1))
            f.write("%g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, -p.y, 0, 0, 1))
            f.write("%g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, -p.x, -p.y, 0, 0, 1))

print "file written"
