#!/usr/bin/env python

from point import Point
from line import Line
from surface import Surface
from volume import Volume
from gmsh import Gmsh

class Box(Volume):
    def __init__(self, corner, L, W, H, vid=None):
        p = corner
        corners = []
        corners.append(p)
        corners.append(Point(p.x + L, p.y, p.z, p.ref))
        corners.append(Point(p.x + L, p.y + W, p.z, p.ref))
        corners.append(Point(p.x, p.y + W, p.z, p.ref))
        corners.append(Point(p.x, p.y + W, p.z + H, p.ref))
        corners.append(Point(p.x + L, p.y + W, p.z + H, p.ref))
        corners.append(Point(p.x + L, p.y, p.z + H, p.ref))
        corners.append(Point(p.x, p.y, p.z + H, p.ref))
        surfaces = []
        surfaces.append(Surface.frompoints(corners[:4]))
        surfaces.append(Surface.frompoints(corners[4:]))
        surfaces.append(Surface.frompoints(corners[:2] + corners[6:]))
        surfaces.append(Surface.frompoints(corners[2:6]))
        surfaces.append(Surface.frompoints(corners[1:3] + corners[5:7]))
        surfaces.append(Surface.frompoints(corners[:1] + corners[3:5] + \
                corners[-1:]))
        Volume.__init__(self, surfaces)
        self.corners = corners
        self.L = L
        self.W = W
        self.H = H

