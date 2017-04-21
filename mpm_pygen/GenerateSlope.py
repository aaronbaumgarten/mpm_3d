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
Lx = 10.0
Ly = 3.0
Lz = 0.2
#Ne = 40
Nx = 50
Ny = 15
Nz = 1
lmpp = 3
hx = Lx/Nx
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81
phi = 0.54

# grain properties
bulk_properties = { 'rho': 2700*phi }
grain_width = Lx
grain_depth = 0.5 * Lx
grain_height = Lz/lmpp
fluid_depth = 0.2 * Lx

bulk_primitive = Primitives3d.Slope2D(0.3*Lx, 0.7*Lx,
                                      0.1*Lx, 0.3*Lx,
                                      0, grain_height
                                      )
bulk_body = CSGTree3d.Node(bulk_primitive)
print "bulk created"

# add points to arrays
bulk_point_array = []
grid.point_array = []
grid.generate_point_array(bulk_primitive)
bulk_point_array = grid.point_array


# fluid properties
fluid_properties = { 'rho': 1000.0 }
fluid_width = Lx
fluid_height = Lz/(lmpp)
fluid_depth = fluid_depth#Ly

fluid_primitive = Primitives3d.Box(0,Lx,
                                 0, fluid_depth,
                                 0, fluid_height
                                 )
fluid_body = CSGTree3d.Node(fluid_primitive)
print "fluid created"

# add points to arrays
fluid_point_array = []
grid.point_array = []
grid.generate_point_array(fluid_primitive)
fluid_point_array = grid.point_array


# add arrays to mpm dictionary
mpm_points = []
nb1 = len(bulk_point_array) #number of material points in body 1
nb2 = len(fluid_point_array)

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
    for p in fluid_point_array:
        pointMass = fluid_properties['rho']*grid.material_point_volume
        pointVolume = grid.material_point_volume
        if bulk_primitive.encompasses(p):
            pointMass *= 1-phi
            pointVolume *= 1-phi
        f.write("%g %g %g %g %g %g %g %g %g %g\n" % (1, pID, pointMass, pointVolume, p.x, p.y, p.z, 0, 0, 0))
        pID += 1

print "file written"
