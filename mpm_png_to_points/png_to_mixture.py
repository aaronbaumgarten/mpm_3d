#!/usr/bin/env python
import sys
import os
try:
    import Image
except:
    from PIL import Image
import numpy as np

# linear material points per pixel (linear, so 2 -> 4 points in 1 pixel)
lmpp = 2
rho = 1500
rho_f = 1000
phi = 0.6
# cell spacing for material point
cs = 1.0 / lmpp
s = map(lambda k: (0.5 + k) * cs, range(0, lmpp))
xy_s = [(x,y) for x in s for y in s]

if (len(sys.argv) > 3):
    imgfile = Image.open(sys.argv[-3])
    outfile = open(sys.argv[-2], "w")
    fluidfile = open(sys.argv[-1], "w")
else:
    print sys.argv[0], "imagefile outfile fluidfile"
    exit(0)

# make image black/white
im = imgfile.convert("1")

print "File", sys.argv[-2], "has size", im.size, "pixels."

#assign length and lower left position.
Length = 1.0;
X_LL = 0.0;
Y_LL = 0.0;

dx = float(Length) / im.size[0]
dy = dx

# pixel spacing
#dx = float(1) / im.size[0]
#dy = float(1) / im.size[1]

ij_array = [(i, j) for i in xrange(0, im.size[0]) for j in xrange(0, im.size[1])]
pixel_array = map(im.getpixel, ij_array)

# if pixel is black, we have material points in that location
filled_array = [ij_array[idx] for idx, pixel in enumerate(pixel_array) if pixel == 0]
fluid_array = [ij_array[idx] for idx, pixel in enumerate(pixel_array) if pixel == 255]

#print pixel_array
#exit(0)

#print filled_array

# i and j are measured from top left of image (i column, j row).
# mpm_2d measures from bottom left, so convert coordinates.
# xy nodes now has a tuple containing the coordinates of the bottom left node
# of the filled element
xy_nodes = map(lambda ij:
                ((float(Length)*float(ij[0])/im.size[0]),(float(Length)/im.size[0])*im.size[1] * (1.0-(float(ij[1])+1)/im.size[1])),
                filled_array)
xy_nodes_fluid = map(lambda ij:
                ((float(Length)*float(ij[0])/im.size[0]),(float(Length)/im.size[0])*im.size[1] * (1.0-(float(ij[1])+1)/im.size[1])),
                fluid_array)

xy_material_points = [(n[0]+dx*sp[0], n[1]+dy*sp[1]) for n in xy_nodes for sp in xy_s]
xy_fluid_points = [(n[0]+dx*sp[0], n[1]+dy*sp[1]) for n in xy_nodes_fluid for sp in xy_s]

#print xy_nodes
#print xy_material_points
#sys.exit(0)

material_points = map(lambda xy: {
                        'm':cs*cs*rho*dx*dy, 'v':cs*cs*dx*dy,
                        'x':xy[0]+X_LL, 'y':xy[1]+Y_LL,
                        'x_t':0, 'y_t':0,
                        'sxx':0, 'sxy':0, 'syy':0}, xy_material_points)

fluid_points =      map(lambda xy: {
                        'm':cs*cs*rho_f*dx*dy, 'v':cs*cs*dx*dy,
                        'x':xy[0]+X_LL, 'y':xy[1]+Y_LL,
                        'x_t':0, 'y_t':0,
                        'sxx':0, 'sxy':0, 'syy':0}, xy_fluid_points)

print "Have", len(material_points), "solid points and", (len(material_points)+len(fluid_points)), "fluid points."
print "Writing to file", sys.argv[-2], " and", sys.argv[-1]

outfile.write("%d\n" % len(material_points))
for point in material_points:
    outfile.write("%g %g %g %g %g %g 1\n" % (point['m'], point['v'], point['x'], point['y'], point['x_t'], point['y_t']))

fluidfile.write("%d\n" % (len(material_points)+len(fluid_points)))
for point in material_points:
    fluidfile.write("%g %g %g %g %g %g 1\n" % (point['m']*(float(rho_f)/float(rho))*(1.0-phi), point['v'], point['x'], point['y'], point['x_t'], point['y_t']))

for point in fluid_points:
    fluidfile.write("%g %g %g %g %g %g 1\n" % (point['m'], point['v'], point['x'], point['y'], point['x_t'], point['y_t']))
