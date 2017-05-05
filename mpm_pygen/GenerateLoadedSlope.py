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
Lx = 15.0
Ly = 4.0
Lz = 0.1
#Ne = 40
Nx = 60
Ny = 16
Nz = 1
Lz = Lx/Nx * Nz
lmpp = 3
hx = Lx/Nx
grid = Grid3d.CartesianPointGrid(Lx, Ly, Lz, Nx, Ny, Nz, lmpp)
print "grid created"

# global properties
g = -9.81
phi = 1800.0/2700.0

# grain properties
bulk_properties = { 'rho': 1800 } # 2700*phi }
grain_width = Lx
grain_depth = 0.1 * Lx
grain_height = Lz/lmpp
fluid_depth = 2.0#0.1 * Lx

bulk_primitive = Primitives3d.Slope2D(0.2*Lx, 0.4*Lx,
                                      1.0, 3.0, #0.1*Lx, 0.2*Lx,
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

# load properties
load_properties = { 'rho': 1100.0 }
load_width = 1.0
load_height = Lz/(lmpp)
load_depth = 0.5

load_primitive = Primitives3d.Box(0.2*Lx - load_width, 0.2*Lx,
                                 3.0, 3.0+load_depth,
                                 0, load_height
                                 )
load_body = CSGTree3d.Node(load_primitive)
print "load created"

# add points to arrays
load_point_array = []
grid.point_array = []
grid.generate_point_array(load_primitive)
load_point_array = grid.point_array

# add arrays to mpm dictionary
mpm_points = []
nb1 = len(bulk_point_array) #number of material points in body 1
nb2 = len(fluid_point_array)
nb3 = len(load_point_array)

grid.write_grid_file(grid_filename)
with open(particle_filename, 'w') as f:
    f.write("%d\n" % 3)
    f.write("%d\n" % nb1)
    f.write("%d\n" % nb2)
    f.write("%d\n" % nb3)
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
            #pointVolume *= 1-phi
        f.write("%g %g %g %g %g %g %g %g %g %g\n" % (1, pID, pointMass, pointVolume, p.x, p.y, p.z, 0, 0, 0))
        pID += 1
    pID = 0
    for p in load_point_array:
        pointMass = load_properties['rho']*grid.material_point_volume
        pointVolume = grid.material_point_volume
        f.write("%g %g %g %g %g %g %g %g %g %g\n" % (2, pID, pointMass, pointVolume, p.x, p.y, p.z, 0, 0, 0))
        pID += 1

print "file written"
