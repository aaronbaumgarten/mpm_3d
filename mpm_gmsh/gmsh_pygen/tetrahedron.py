#!/usr/bin/env python

from point import Point
from line import Line
from surface import Surface
from volume import Volume
from gmsh import Gmsh

class Tetrahedron(Volume):
    def __init__(self, corners, vid=None):
        if len(corners) != 4:
            raise ValueError("A tetrahedron has four corners.")
        p0 = corners[0]
        p1 = corners[1]
        p2 = corners[2]
        p3 = corners[3]
        s0 = Surface.frompoints([p0, p1, p2])
        s1 = Surface.frompoints([p0, p1, p3])
        s2 = Surface.frompoints([p0, p2, p3])
        s3 = Surface.frompoints([p1, p2, p3])
        Volume.__init__(self, [s0, s1, s2, s3])
        self.corners = corners
        
