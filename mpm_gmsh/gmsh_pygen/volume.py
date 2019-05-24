#!/usr/bin/env python

from point import Point
from line import Line
from surface import Surface

class Volume:
    """A class for volumes in gmsh"""
    def __init__(self, surfaces, vid=None):
        def validsurfaces(surf):
            for s1 in surf:
                for e in s1:
                    a = 0
                    for s2 in surf:
                        if s1 is s2:
                            continue
                        if e in s2 or -e in s2:
                            a += 1
                    if a != 1:
                        print(a)
                        return False
            return True
        
        if not validsurfaces(surfaces):
            raise ValueError("Surfaces do not form a valid volume.")
        surfaces = surfaces[:]
        
        checked = [False] * len(surfaces)
        oriented = [0]
        
        while len(oriented) > 0:
            i = oriented.pop()
            if checked[i]:
                continue
            checked[i] = True
            s1 = surfaces[i]
            for j in range(len(surfaces)):
                if checked[j] or j in oriented:
                    continue
                s2 = surfaces[j]
                for e1 in s1.lines:
                    for e2 in s2.lines:
                        if e1 == e2:
                            surfaces[j] = -s2
                            oriented.append(j)
                        elif e1 == -e2:
                            oriented.append(j)
        self.surfaces = surfaces
        if vid != None and vid <= 0:
            raise ValueError("Volume id must be positive.")
        self.vid = vid
    
    def __eq__(self, v):
        return len(self.surfaces) == len(v.surfaces) and all(s in v.surfaces \
                for s in self.surfaces)
    
    def __str__(self):
        return "Surface Loop({0}) = {{{1}}};\nVolume({0}) = {{{0}}};\n".\
                format(self.vid, ", ".join([str(s.sid) for s in self.surfaces]))
