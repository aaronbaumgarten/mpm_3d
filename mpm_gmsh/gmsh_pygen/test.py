#!/usr/bin/env python

from point import Point
from line import Line
from surface import Surface
from volume import Volume
from tetrahedron import Tetrahedron
from box import Box
from gmsh import Gmsh

g = Gmsh()
o = Point(0, 0, 0, 2)#.25)
sand = Box(o, 1, 2, 1)
p0 = Point(0, 0, 1, 2)#.75)
air = Box(p0, 1, 2, 1)
p1 = Point(.25, .5, .75, .5)#.1)
impact = Box(p1, .5, 1, .25)

g.addvolume(air)
g.addvolume(sand)
g.addvolume(impact)

g.generatefile("out.gmsh")
