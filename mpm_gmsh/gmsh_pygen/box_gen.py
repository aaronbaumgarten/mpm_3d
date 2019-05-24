#!/usr/bin/env python

from point import Point
from line import Line
from surface import Surface
from volume import Volume
from tetrahedron import Tetrahedron
from box import Box
from gmsh import Gmsh
import math

g = Gmsh()

# Taken from Soliman's paper
radius = 0.0047625 
L_impression = 0.381
W_impression = 0.1905
H_impression = 0.0508

# Refinements
ref_air = 2 * radius
ref_sand = radius
ref_med = radius / 3
ref_high = radius / 8

# Domain size definitions
L = L_impression
W = radius * 10 # Note: this is < W_impression / 2, so we may have problems
H = H_impression * 4

# Bed size definitions
# L_bed = L
# W_bed = W
H_bed = H_impression * 2

# Refinement area definition
L_ref = L / 2
W_ref = W / 2
H_ref = H_impression

# Define less refined sand boxes
o = Point(0, 0, 0, ref_sand)
v2 = Box(o, L, W, H_bed - H_ref)

p6 = Point(0, W_ref, H_bed - H_ref, ref_sand)
v6 = Box(p6, L, W - W_ref, H_ref)

p9 = Point(0, 0, H_bed - H_ref, ref_sand)
v9 = Box(p9, (L - L_ref) / 2, W_ref, H_ref)

p10 = Point( (L + L_ref) / 2, 0, H_bed - H_ref, ref_sand)
v10 = Box(p10, (L - L_ref) / 2, W_ref, H_ref)

# Define less refined air boxes
p1 = Point(0, 0, H_bed + H_ref, ref_air)
v1 = Box(p1, L, W, H - H_bed - H_ref)

p5 = Point(0, W_ref, H_bed, ref_air)
v5 = Box(p5, L, W - W_ref, H_ref)

p7 = Point(0, 0, H_bed, ref_air)
v7 = Box(p7, (L - L_ref) / 2, W_ref, H_ref)

p8 = Point( (L + L_ref) / 2, 0, H_bed, ref_air)
v8 = Box(p8, (L - L_ref) / 2, W_ref, H_ref)

# Define refined sand and air boxes
p4 = Point( (L - L_ref) / 2, 0, H_bed - H_ref, ref_high)
v4 = Box(p4, L_ref, W_ref, H_ref)

p3 = Point( (L - L_ref) / 2, 0, H_bed, ref_high)
v3 = Box(p3, L_ref, W_ref, H_ref)

# Add volumes

g.addvolume(v1)
g.addvolume(v2)
g.addvolume(v3)
g.addvolume(v4)
g.addvolume(v5)
g.addvolume(v6)
g.addvolume(v7)
g.addvolume(v8)
g.addvolume(v9)
g.addvolume(v10)

g.generatefile("sandcollision boxy.geo")
