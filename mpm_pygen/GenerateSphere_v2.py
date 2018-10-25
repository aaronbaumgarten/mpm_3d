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
Lx = 2 / 39.37
Ly = 7 / 39.37
Lz = 4 / 39.37
#Ne = 40
dx = 3.0/8 / 10 / 39.37
Nx = int(math.ceil(Lx / dx))
Ny = int(math.ceil(Ly / dx))
Nz = int(math.ceil(Lz / dx))
lmpp = 2
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81

# free block properties
block_properties = { 'rho': 7800.0 }
block_radius = 3.0 / 8 / 39.37
block_center = Primitives3d.Point(0, 2 * dx + block_radius, Lz / 2 + 2 * dx + block_radius)
block_primitive = Primitives3d.Sphere(block_center, block_radius)
block_body = CSGTree3d.Node(block_primitive)
block_angle = 20
block_speed = 500
block_velocity = (block_speed * math.cos(block_angle), -block_speed * math.sin(block_angle))
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
    #    mpm_point = {
    #                    'body': 1,
    #                    'm': bulk_properties['rho']*grid.material_point_volume,
    #                    'v': grid.material_point_volume,
    #                    'x': p.x, 'y': p.y, 'z': p.z,
    #                    'x_t':0, 'y_t':0, 'z_t':0 #,
    #                    #'sxx': 0, 'sxy': 0, 'syy': 0
    #                    'active':1
    #                }
    #    if (mpm_point['body'] == 1):
    #        nb1 += 1
    #    mpm_points.append(mpm_point)
        f.write("%g %g %g %g %g %g %g %g %i\n" % (block_properties['rho']*grid.material_point_volume, grid.material_point_volume, p.x, p.y, p.z, 0, block_velocity[0], block_velocity[1], 1))

print "file written"
