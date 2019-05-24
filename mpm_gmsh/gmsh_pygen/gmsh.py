#!/usr/bin/env python

from point import Point
from line import Line
from surface import Surface
from volume import Volume

class Gmsh:
    """A class for making gmsh files"""
    def __init__(self):
        self.points = []
        self.lines = []
        self.surfaces = []
        self.volumes = []
    
    def addpoint(self, point):
        """Adds point to this gmsh file.
        If point already exists, updates the refinement of the gmsh point
        to the minimum of the two refinements, and the point id (pid) of point.
        If point does not exist and has no pid, updates the pid and adds
        point to this gmsh file.
        If point does not exist and has a used pid, overwrites current 
        point with the same pid
        If point does not exist and has an unused pid, raises an error.
        """
        if point in self.points:
            i = self.points.index(point)
            ref = min(self.points[i].ref, point.ref)
            self.points[i].ref = ref
            point.ref = ref
            point.pid = i + 1
        elif point.pid == None:
            point.pid = len(self.points) + 1
            self.points.append(point)
        elif point.pid <= 0:
            raise ValueError("Point id must be positive")
        elif point.pid <= len(self.points):
            self.points[point.pid - 1] = point
        else:
            raise ValueError("If a point id is specified it must already exist")

    def addline(self, line):
        """Adds line to this gmsh file.
        If the line has a negative line id (lid), adds -line.
        If the line has a positive, used lid, overwrites current line with 
        the same lid and adds points in line to this gmsh file (does not 
        remove points in overwritten line).
        If line has a nonnegative, unused lid, raises an error.
        If line has no lid and already exists, updates refinement of points
        in line and lid of line
        If line has no lid and does not exist, adds line (and its points) to
        this gmsh file and adds lid to line.
        """
        if line.lid == None:
            self.addpoint(line.p1)
            self.addpoint(line.p2)
            if line in self.lines:
                i = self.lines.index(line)
                line.lid = i + 1
            elif -line in self.lines:
                i = self.lines.index(-line)
                line.lid = -(i + 1)
            else:
                line.lid = len(self.lines) + 1
                self.lines.append(line)
        elif line.lid == 0:
            raise ValueError("Line id cannot be 0")
        elif line.lid < 0:
            self.addline(-line)
        elif line.lid <= len(self.lines):
            self.addpoint(line.p1)
            self.addpoint(line.p2)
            self.lines[line.lid - 1] = line
        else:
            raise ValueError("If a line id is specified it must already exist")

    def addsurface(self, surface):
        """Adds surface to this gmsh file.
        If the surface has a negative surface id (sid), adds -surface.
        If the surface has a positive, used sid, overwrites current surface
        with the same sid and adds lines (and therefore points) in surface 
        to this gmsh file (but does not remove lines and points in 
        overwritten surface).
        If surface has a nonnegative, unused sid, raises an error.
        If surface has no sid and does not exist, adds surface (with its 
        lines and points) to this gmsh file and adds sid to surface.
        If surface has no sid and already exists, updates refinement of 
        points in surface.
        """
        if surface.sid == None:
            for l in surface.lines:
                self.addline(l)
            if surface in self.surfaces:
                i = self.surfaces.index(surface)
                surface.sid = i + 1
            elif -surface in self.surfaces:
                i = self.surfaces.index(-surface)
                surface.sid = -(i + 1)
            else:
                surface.sid = len(self.surfaces) + 1
                self.surfaces.append(surface)
        elif surface.sid < 0:
            self.addsurface(-surface)
        elif surface.sid == 0:
            raise ValueError("Surface id cannot be 0")
        elif surface.sid > 0 and surface.sid <= len(self.surfaces):
            for l in surface.lines:
                self.addline(l)
            self.surfaces[surface.sid - 1] = surface
        else:
            raise ValueError("If a surface id is specified it must already exist")
    
    def addvolume(self, volume):
        """Adds volume to this gmsh file.
        If the volume has a positive, used vid, overwrites current volume
        with the same vid and adds surfaces (and therefore points and lines)
        in volume to this gmsh file (but does not remove surfaces, lines, 
        and points in overwritten volume).
        If volume has a nonpositive or unused vid, raises an error.
        If volume has no vid and does not exist, adds volume (with its 
        surfaces, lines, and points) to this gmsh file and adds vid to 
        volume.
        If volume has no vid and already exists, updates refinement of 
        points in volume.
        """
        if volume.vid == None:
            if volume in self.volumes in self.volumes:
                i = self.volumes.index(volume)
                volume.vid = i + 1
                for s in volume.surfaces:
                    for p in s.points:
                        self.addpoint(p)
            else:
                for s in volume.surfaces:
                    self.addsurface(s)
                volume.vid = len(self.volumes) + 1
                self.volumes.append(volume)
        elif volume.vid <= 0:
            raise ValueError("Volume id must be positive")
        elif volume.vid <= len(self.volumes):
            for s in volume.surfaces:
                self.addsurface(s)
            self.volumes[volume.vid - 1] = volume
        else:
            raise ValueError("If a volume id is specified it must already exist")
    
    def generatefile(self, filename):
        with open(filename, 'w') as f:
            f.write("// Generated gmsh file\n\n")
             
            if len(self.points) > 0:
                f.write("// Points\n\n")
                for p in self.points:
                    f.write(str(p))
                f.write("\n")

            if len(self.lines) > 0:
                f.write("// Lines\n\n")
                for l in self.lines:
                    f.write(str(l))
                f.write("\n")
            
            if len(self.surfaces) > 0:
                f.write("// Surfaces\n\n")
                for s in self.surfaces:
                    f.write(str(s)) 
                f.write("\n")
           
            if len(self.volumes) > 0:
                f.write("// Volumes\n\n")
                for v in self.volumes:
                    f.write(str(v)) 
                f.write("\n") 
            
    """
    def deprecatedgeneratefile(self, filename):
        with open(filename, 'w') as f:
            f.write("// Generated gmsh file\n\n")
            
            p_i = 1
            l_i = 1
            s_i = 1
            v_i = 1
            
            if len(self.points) > 0:
                f.write("// Stand-alone points\n\n")
                for p in self.points:
                    f.write("Point({}) = {};\n".format(p_i, str(p)))
                    p_i += 1
                f.write("\n")

            if len(self.lines) > 0:
                f.write("// Stand-alone lines\n\n")
                for l in self.lines:
                    f.write("// Line {}\n".format(l_i))
                    
                    p1_i = p_i
                    f.write("Point({}) = {};\n".format(p_i, str(l.p1)))
                    p_i += 1
                    
                    p2_i = p_i
                    f.write("Point({}) = {};\n".format(p_i, str(l.p2)))
                    p_i += 1
                    
                    f.write("Line({}) = {{{}, {}}};\n".format(l_i, p1_i, \
                            p2_i))
                    l_i += 1
                f.write("\n")
            
            if len(self.surfaces) > 0:
                f.write("// Stand-alone surfaces\n\n")
                for s in self.surfaces:
                    f.write("// Surface {}\n".format(s_i))
                    f.write("Point({}) = {};\n".format(p_i, \
                            str(s.points[0])))
                    n = len(s.points)
                    for i in range(1, n):
                        f.write("Point({}) = {};\n".format(p_i + i, \
                                str(s.points[i])))
                        f.write("Line({}) = {{{}, {}}};\n".format(l_i + \
                                i - 1, p_i + i - 1, p_i + i))
                    f.write("Line({}) = {{{}, {}}};\n".format(l_i + n - 1, \
                            p_i + n - 1, p_i))
                    p_i += n
                    f.write("Curve Loop({}) = {{{}}};\n".format(s_i, \
                            ", ".join([ str(l_i + x) for x in range(n) ])))
                    l_i += n
                    f.write("Plane Surface({0}) = {{{0}}};\n".format(s_i))
                    s_i += 1
                f.write("\n")
            
            if len(self.volumes) > 0:
                f.write("// Stand-alone volumes\n\n")
                for v in self.volumes:
                    f.write("// Volume {}\n".format(v_i))
                    m = len(v.surfaces)
                    for j in range(m):
                        s = v.surfaces[j]
                        f.write("Point({}) = {};\n".format(p_i, \
                                str(s.points[0])))
                        n = len(s.points)
                        for i in range(1, n):
                            f.write("Point({}) = {};\n".format(p_i + i, \
                                    str(s.points[i])))
                            f.write("Line({}) = {{{}, {}}};\n".format(l_i \
                                    + i - 1, \ p_i + i - 1, p_i + i))
                        f.write("Line({}) = {{{}, {}}};\n".format(l_i + n \
                                - 1, p_i + n - 1, p_i))
                        p_i += n
                        f.write("Curve Loop({}) = {{{}}};\n".format(s_i + \
                                j, ", ".join([ str(l_i + x) for x in \
                                range(n) ])))
                        l_i += n
                        f.write("Plane Surface({0}) = {{{0}}};\n".format(s_i+j))
                    f.write("Surface Loop({}) = {{{}}};\n".format(v_i, \
                            ", ".join([ str(s_i + x) for x in range(m) ])))
                    s_i += m
                    f.write("Volume({0}) = {{{0}}};\n".format(v_i))
                    v_i += 1
     """
