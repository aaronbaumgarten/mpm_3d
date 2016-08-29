#!/usr/bin/env python
import math

class Point:
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __repr__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'

    def __eq__(self, rhs):
        return (self.x == rhs.x and self.y == rhs.y and self.z == rhs.z)
    def __req__(self, lhs):
        return (self == lhs)

    def __neg__(self):
        return Point(-self.x, -self.y, -self.z)
    def __pos__(self):
        return Point(self.x, self.y, self.z)

    def __add__(self, rhs):
        return Vector(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    def __radd__(self, lhs):
        return (self + lhs)

    def __sub__(self, rhs):
        return (self + (-rhs))
    def __rsub__(self, lhs):
        return (lhs + (-self))
    

class Vector:
    def __init__(self, x, y, z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __repr__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'

    def magnitude(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)
    def __pos__(self):
        return Vector(self.x, self.y, self.z)

    def __add__(self, rhs):
        return Vector(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    def __radd__(self, lhs):
        return (self + lhs)

    def __sub__(self, rhs):
        return (self + (-rhs))
    def __rsub__(self, lhs):
        return (lhs + (-self))

class Circle:
    def __init__(self, center, radius):
        self.center = center
        self.radius = float(radius)
    def __repr__(self):
        return '[center: ' + str(self.center) + ', radius: ' + str(self.radius) + ']'
    def encompasses(self, point):
        return (((point - self.center).magnitude() - self.radius) <= 0)

class Box:
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.ymin = float(ymin)
        self.ymax = float(ymax)
        self.zmin = float(zmin)
        self.zmax = float(zmax)
    
    def encompasses(self, point, tol=1e-8):
        inside = True
        return (point.x>=self.xmin and point.x<=self.xmax and point.y>=self.ymin and point.y<=self.ymax and point.z>=self.zmin and point.z<=self.zmax)

#class ConvexPolygon:
#    def __init__(self, ccw_vertices):
#        self.ccw_vertices = ccw_vertices
#
#    def orientation(self, p1,p2,p3):
#        theta = (p2.y - p1.y) * (p3.x - p2.x) - (p3.y - p2.y) * (p2.x - p1.x)
#        return theta
#
#    def signed_angle(self, ap, bp):
#        theta_ap = math.atan2(ap.y, ap.x)
#        theta_bp = math.atan2(bp.y, bp.x)
#        theta = theta_ap - theta_bp
#        while (theta >= (math.pi)):
#            theta -= (2*math.pi)
#        while (theta < (-math.pi)):
#            theta += (2*math.pi)
#        return theta
#
#    def encompasses(self, point, tol=1e-8):
#        inside = True
#        anglesum = 0
#        for i in range(1, len(self.ccw_vertices)):
#            a = self.ccw_vertices[i-1]
#            b = self.ccw_vertices[i]
#            bp = point - b
#            ap = point - a
#            anglesum += self.signed_angle(ap, bp)
#
        # return (abs(anglesum % (2 * math.pi)) <= tol)
#        return (abs(anglesum) >= math.pi)