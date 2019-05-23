#!/usr/bin/env python

from point import Point
from line import Line

class Surface:
    """A class for surfaces in gmsh"""
    def __init__(self, lines, sid=None):
        """Precondition: the lines are in order"""
        self.lines = lines[:]
        self.lines[0], self.lines[1] = Line.order(lines[0], lines[1])
        for i in range(1, len(lines)):
            prev = self.lines[i-1]
            curr = self.lines[i]
            if not prev.shares_point(curr):
                raise ValueError("Lines do not connect")
            if not prev.connected(curr):
                self.lines[i] = -curr
        for i in range(2, len(self.lines)-1):
            if not Line.coplanar(self.lines[i], self.lines[0]):
                raise ValueError("Lines are not coplanar")
        
        self.points = []
        for l in self.lines:
            self.points.append(l.p1)
        if sid == 0:
            raise ValueError("Surface id cannot be 0")
        self.sid = sid
    
    def __neg__(self):
        lines = []
        n = len(self.lines)
        for i in range(n):
            lines.append(-self.lines[n - i - 1])
        if self.sid == None:
            return Surface(lines, None)
        return Surface(lines, -self.sid)
    
    def __eq__(self, s):
        if not (self.points[0] in s.points) or len(self.points)!=len(s.points):
            return False
        else:
            shift = s.points.index(self.points[0])
            return all(self.points[i] == s.points[(i + shift) % len(s.points)] \
                    for i in range (1, len(s.points)))
    
    def __iter__(self):
        return iter(self.lines)
    
    def __str__(self):
        if self.sid < 0:
            return str(-self)
        return "Curve Loop({0}) = {{{1}}};\nPlane Surface({0}) = {{{0}}};\n".\
                format(self.sid,", ".join([str(l.lid) for l in self.lines]))
    
    @classmethod
    def frompoints(cls, pts, sid=None):
        lines = []
        for i in range(1, len(pts)):
            lines.append(Line(pts[i-1], pts[i]))
        lines.append(Line(pts[-1], pts[0]))
        return cls(lines, sid)
    
    def __contains__(self, line):
        return any(l == line for l in self.lines)
    
    def sharesedge(self, other):
        return any (l in other or -l in other for l in self.lines)
        
    def normal(self):
        v1 = self.lines[0]
        v2 = self.lines[1]
        return v1.crossprod(v2)
    
    def centroid(self):
        sx, sy, sz = 0, 0, 0
        for p in self.points:
            sx += p.x
            sy += p.y
            sz += p.z
        n = len(self.points)
        return (sx/n, sy/n, sz/n)
