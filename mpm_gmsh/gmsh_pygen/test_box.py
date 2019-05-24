#!/usr/bin/env python

from point import Point
from line import Line
from surface import Surface
from volume import Volume
from tetrahedron import Tetrahedron
from box import Box
from gmsh import Gmsh

g = Gmsh()
p = Point(0, 0, 0, .1)
v = Box(p, 1, 1, 1)
g.addvolume(v)

g.generatefile("box.gmsh")
