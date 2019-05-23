#!/usr/bin/env python

from point import Point
import math

class Line:
    """A class for lines in gmsh"""
    def __init__(self, p1, p2, lid=None):
        self.p1 = p1
        self.p2 = p2
        if lid == 0:
            raise ValueError("Line id cannot be 0")
        self.lid = lid
    
    def __eq__(self, l):
        return (self.p1 == l.p1) and (self.p2 == l.p2)

    def __ne__(self,l):
        return not (self == l)
    
    def __neg__(self):
        if self.lid == None:
            return Line(self.p2, self.p1, None)
        return Line(self.p2, self.p1, -self.lid)
    
    def __str__(self):
        if self.lid > 0:
            return "Line({}) = {{{}, {}}};\n".format(self.lid, \
                    self.p1.pid, self.p2.pid)
        else:
            return "Line({}) = {{{}, {}}};\n".format(-self.lid, \
                    self.p2.pid, self.p1.pid) 

    def __contains__(self, p):
        return self.p1 == p or self.p2 == p
    
    def length(self):
        return (self.p1).dist(self.p2)
    
    def connected(l1, l2):
        return l1.p2 == l2.p1
    
    def shares_point(l1, l2):
        return l1.p1 in l2 or l1.p2 in l2

    @staticmethod
    def order(l1, l2):
        if not Line.shares_point(l1,l2):
            raise ValueError("Cannot order lines which don't share a point")
        if l1.p2 == l2.p1:
            return (l1, l2)
        if l1.p1 == l2.p1:
            return (-l1, l2)
        if l1.p2 == l2.p2:
            return (l1, -l2)
        return (-l1, -l2)
    
    @staticmethod
    def coplanar(l1, l2):
        return Point.coplanar(l1.p1, l1.p2, l2.p1, l2.p2)

    def crossprod(v1, v2):
        ox = v1.y * v2.z - v2.y * v1.z
        oy = v1.z * v2.x - v2.z * v1.x
        oz = v1.x * v2.y - v2.x * v1.y
        return (ox, oy, oz)
