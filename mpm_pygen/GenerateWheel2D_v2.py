#!/usr/bin/env python
import sys
import math

import Primitives3d
import CSGTree3d
import Grid3d

if (len(sys.argv) > 1):
    wheel_filename = sys.argv[-2] + '.points'
    bulk_filename = sys.argv[-1] + '.points'
else:
    print sys.argv[0], "outprefix"
    exit(0)
print "files named"

#grid properties
#Ly = Lx = Lz = 0.4
Lx = 2.0
Ly = 1.0
Lz = 1.0
#Ne = 40
Nx = 100
Ny = 50
Nz = 1
lmpp = 3
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81

# free block properties
block_properties = { 'rho': 2000.0 }
block_width = 0.2#*Lx
block_height = 0.2
block_depth = 0.2#*Ly
spike_height = 0.15
wheel_frac = 0.5
nspike = 4

block_primitive = Primitives3d.WheelOne(
                        Primitives3d.Point(Lx*(1.0-2.0/Nx) - 3*block_width, 0.3 + block_height/2 + spike_height, 0.0), 
                        Primitives3d.Point(Lx*(1.0-2.0/Nx) - 3*block_width, 0.3 + block_height/2 + spike_height, 1.0),
                        block_width/2,
                        spike_height,
                        wheel_frac,
                        nspike)
                  #Primitives3d.Cylinder(Primitives3d.Point(Lx-block_width,0.3+block_depth/2,0),
                  #                      Primitives3d.Point(Lx-block_width,0.3+block_depth/2,block_height),
                  #                      block_width/2)
                  #Primitives3d.Box(#(Lx-block_width)/2, (Lx+block_width)/2,
                  #               block_width,block_width+block_width,
                  #               (Ly-block_depth)/2, (Ly+block_depth)/2,
                  #               #0,block_height,
                  #               (Lz-block_height)/2, (Lz+block_height)/2,
                  #               )
block_body = CSGTree3d.Node(block_primitive)
print "body 1 created"

# add points to arrays
block_point_array = []
grid.generate_point_array(block_primitive)
block_point_array = grid.point_array

print "point array 1 created"

# free block properties
block2_properties = { 'rho': 1500.0 }
block2_width = Lx
block2_height = 0.3
block2_depth = Ly

block2_primitive = Primitives3d.Box(0,block2_width,
                                 0,block2_height,
                                 0,1.0,
                                 )
block2_body = CSGTree3d.Node(block2_primitive)
print "body 2 created"

# add points to arrays
block2_point_array = []
grid.point_array = []
grid.generate_point_array(block2_primitive)
block2_point_array = grid.point_array

print "point array 2 created"

# add arrays to mpm dictionary
mpm_points = []
nb1 = len(block_point_array) #number of material points in body 1
nb2 = len(block2_point_array)

with open(wheel_filename, 'w') as w, open(bulk_filename, 'w') as b:
    w.write("%d\n" % nb1)
    b.write("%d\n" % nb2)
    pID = 0
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
        w.write("%g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, 0, 0, 1))
    for p in block2_point_array:
        b.write("%g %g %g %g %g %g %i\n" % (block2_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, 0, 0, 1))

print "file written"
