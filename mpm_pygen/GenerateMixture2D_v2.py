#!/usr/bin/env python
import sys
import math

import Primitives3d
import CSGTree3d
import Grid3d

if (len(sys.argv) > 2):
    bulk_filename = sys.argv[-2] + '.points'
    fluid_filename = sys.argv[-1] + '.points'
else:
    print sys.argv[0], "bulkprefix fluidprefix"
    exit(0)
print "files named"

#grid properties
#Ly = Lx = Lz = 0.4
Lx = 2.0
Ly = 1.0
#Ne = 40
Nx = 150
Ny = 75
Nz = 1
lmpp = 1
Lz = Nz*lmpp
hx = Lx/Nx
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81
phi = 0.50

# grain properties
bulk_properties = { 'rho': 1620*phi }
grain_width = 0.06
grain_depth = 0.048
grain_height = 1.0
fluid_depth = 0.10#0.1 * Lx

bulk_primitive = Primitives3d.Box(0.0, Lx,
                                  0.0, 0.5, #0.1*Lx, 0.2*Lx,
                                  0, 1.0
                                  )
bulk_body = CSGTree3d.Node(bulk_primitive)
print "bulk created"

# add points to arrays
bulk_point_array = []
grid.point_array = []
grid.generate_point_array(bulk_primitive)
bulk_point_array = grid.point_array


# fluid properties
fluid_properties = { 'rho': 1620.0 }
fluid_width = Lx
fluid_height = 0.08
#fluid_depth = fluid_depth#Ly

fluid_primitive = Primitives3d.Box(0.0, Lx,
                                 0.0, 0.5,
                                 0, 1.0
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

with open(bulk_filename, 'w') as b, open(fluid_filename, 'w') as f:
    b.write("%d\n" % nb1)
    f.write("%d\n" % nb2)
    for p in bulk_point_array:
        b.write("%g %g %g %g %g %g %g\n" % (bulk_properties['rho']*grid.material_point_volume, grid.material_point_volume, (p.x), p.y, 0, 0, 1))
    for p in fluid_point_array:
        pointMass = fluid_properties['rho']*grid.material_point_volume
        pointVolume = grid.material_point_volume
        if bulk_primitive.encompasses(p):
            pointMass *= 1.0-phi
            #pointVolume *= 1.0-phi
        f.write("%g %g %g %g %g %g %g\n" % (pointMass, pointVolume, (p.x), p.y, 0, 0, 1))

print "file written"
