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
L_impression = 0.25
W_impression = 0.15
H_impression = 0.05

# Refinements
ref_1 = 5 * radius
ref_2 = 2 * radius 
ref_3 = radius / 5

# Domain size definitions
L = 2*L_impression
W = W_impression
H = H_impression * 4

# Bed size definitions
# L_bed = L
# W_bed = W
H_bed = H_impression * 2

# Refinement area definition
L_ref = L_impression
W_ref = radius
W_med = W_impression/2.0
H_ref = H_impression
L_0 = L_ref / 2.0

# Define most refined box along impacting direction
p212 = Point(L_0, 0, H_bed - H_ref, ref_3)
b212 = Box(p212, L_ref, W_ref, H_ref)

# Define second box
p222 = Point(L_0, W_ref, H_bed - H_ref, ref_2)
b222 = Box(p222, L_ref, W_med - W_ref, H_ref)

# Define remaining space
p111 = Point(0, 0, 0, ref_1)
b111 = Box(p111, L_0, W_ref, H_bed - H_ref)

p112 = Point(0, 0, H_bed - H_ref, ref_1)
b112 = Box(p112, L_0, W_ref, H_ref)

p113 = Point(0, 0, H_bed, ref_1)
b113 = Box(p113, L_0, W_ref, H - H_bed)

p121 = Point(0, W_ref, 0, ref_1)
b121 = Box(p121, L_0, W_med - W_ref, H_bed - H_ref)

p122 = Point(0, W_ref, H_bed - H_ref, ref_1)
b122 = Box(p122, L_0, W_med - W_ref, H_ref)

p123 = Point(0, W_ref, H_bed, ref_1)
b123 = Box(p123, L_0, W_med - W_ref, H - H_bed)

p131 = Point(0, W_med, 0, ref_1)
b131 = Box(p131, L_0, W - W_med, H_bed - H_ref)

p132 = Point(0, W_med, H_bed - H_ref, ref_1)
b132 = Box(p132, L_0, W - W_med, H_ref)

p133 = Point(0, W_med, H_bed, ref_1)
b133 = Box(p133, L_0, W - W_med, H - H_bed)

p211 = Point(L_0, 0, 0, ref_1)
b211 = Box(p211, L_ref, W_ref, H_bed - H_ref)

p213 = Point(L_0, 0, H_bed, ref_1)
b213 = Box(p213, L_ref, W_ref, H - H_bed)

p221 = Point(L_0, W_ref, 0, ref_1)
b221 = Box(p221, L_ref, W_med - W_ref, H_bed - H_ref)

p223 = Point(L_0, W_ref, H_bed, ref_1)
b223 = Box(p223, L_ref, W_med - W_ref, H - H_bed)

p231 = Point(L_0, W_med, 0, ref_1)
b231 = Box(p231, L_ref, W - W_med, H_bed - H_ref)

p232 = Point(L_0, W_med, H_bed - H_ref, ref_1)
b232 = Box(p232, L_ref, W - W_med, H_ref)

p233 = Point(L_0, W_med, H_bed, ref_1)
b233 = Box(p233, L_ref, W - W_med, H - H_bed)

p311 = Point(L_0 + L_ref, 0, 0, ref_1)
b311 = Box(p311, L - L_0 - L_ref, W_ref, H_bed - H_ref)

p312 = Point(L_0 + L_ref, 0, H_bed - H_ref, ref_1)
b312 = Box(p312, L - L_0 - L_ref, W_ref, H_ref)

p313 = Point(L_0 + L_ref, 0, H_bed, ref_1)
b313 = Box(p313, L - L_0 - L_ref, W_ref, H - H_bed)

p321 = Point(L_0 + L_ref, W_ref, 0, ref_1)
b321 = Box(p321, L - L_0 - L_ref, W_med - W_ref, H_bed - H_ref)

p322 = Point(L_0 + L_ref, W_ref, H_bed - H_ref, ref_1)
b322 = Box(p322, L - L_0 - L_ref, W_med - W_ref, H_ref)

p323 = Point(L_0 + L_ref, W_ref, H_bed, ref_1)
b323 = Box(p323, L - L_0 - L_ref, W_med - W_ref, H - H_bed)

p331 = Point(L_0 + L_ref, W_med, 0, ref_1)
b331 = Box(p331, L - L_0 - L_ref, W - W_med, H_bed - H_ref)

p332 = Point(L_0 + L_ref, W_med, H_bed - H_ref, ref_1)
b332 = Box(p332, L - L_0 - L_ref, W - W_med, H_ref)

p333 = Point(L_0 + L_ref, W_med, H_bed, ref_1)
b333 = Box(p333, L - L_0 - L_ref, W - W_med, H - H_bed)

# Add volumes

g.addvolume(b111)
g.addvolume(b112)
g.addvolume(b113)
g.addvolume(b121)
g.addvolume(b122)
g.addvolume(b123)
g.addvolume(b131)
g.addvolume(b132)
g.addvolume(b133)

g.addvolume(b211)
g.addvolume(b212)
g.addvolume(b213)
g.addvolume(b221)
g.addvolume(b222)
g.addvolume(b223)
g.addvolume(b231)
g.addvolume(b232)
g.addvolume(b233)

g.addvolume(b311)
g.addvolume(b312)
g.addvolume(b313)
g.addvolume(b321)
g.addvolume(b322)
g.addvolume(b323)
g.addvolume(b331)
g.addvolume(b332)
g.addvolume(b333)

g.generatefile("subdomain_3on8in.geo")
