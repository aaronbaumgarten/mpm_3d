#!/usr/bin/env python
import math

class Point:
    """A class for points in gmsh"""
    def __init__(self, x, y, z, ref):
        self.x = x
        self.y = y
        self.z = z
        if ref <= 0:
            raise ValueError("Point refinement must be positive")
        self.ref = ref
        self.pid = None
    
    def __eq__(self, p):
        return (self.x == p.x) and (self.y == p.y) and (self.z == p.z)

    def __ne__(self, p):
        return not (self == p)
    
    def __str__(self):
        return "Point({}) = {{{}, {}, {}, {}}};\n".format(self.pid, \
                self.x, self.y, self.z, self.ref)
    
    def dist(p1, p2):
        dx = p1.x - p2.x
        dy = p1.y - p2.y
        dz = p1.z - p2.z
        return sqrt(dx*dx + dy*dy + dz*dz)

    @staticmethod
    def coplanar(p1, p2, p3, p4, eps=1e-10):
        d = [ [p2.x - p1.x, p2.y - p1.y, p2.z - p1.z], \
            [p3.x - p1.x, p3.y - p1.y, p3.z - p1.z], \
            [p4.x - p1.x, p4.y - p1.y, p4.z - p1.z] ]
        
        return abs(d[0][0] * ( d[1][1]*d[2][2] - d[2][1]*d[1][2] ) \
                - d[0][1] * ( d[1][0]*d[2][2] - d[2][0] * d[1][2] ) + \
                d[0][2] * ( d[1][0]*d[2][1] - d[2][0]*d[1][1] ) ) < eps
