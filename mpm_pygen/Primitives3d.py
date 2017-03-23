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

    def dot(self,vec):
        return (self.x*vec.x + self.y*vec.y + self.z*vec.z)

    def cross(self,vec):
        return Vector(self.y*vec.z-self.z*vec.y,-self.x*vec.z+self.z*vec.x,self.x*vec.y-self.y*vec.x)

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

    def __mul__(self, rhs):
        return (Vector(self.x*rhs,self.y*rhs,self.z*rhs))
    def __rmul__(self,lhs):
        return (self*lhs)

class Sphere:
    def __init__(self, center, radius):
        self.center = center
        self.radius = float(radius)
    def __repr__(self):
        return '[center: ' + str(self.center) + ', radius: ' + str(self.radius) + ']'
    def encompasses(self, point):
        return (((point - self.center).magnitude() - self.radius) <= 0)

class Peanut:
    def __init__(self, center_start, center_end, radius):
        self.center_start = center_start
        self.center_end = center_end
        self.radius = float(radius)
    def __repr__(self):
        return '[start center-line: ' + str(self.center_start) + ', end center-line: ' + str(self.center_end) + ', radius: ' + str(self.radius) + ']'
    def encompasses(self, point):
        avec = self.center_end - self.center_start
        bvec = point - self.center_start
        cval = avec.dot(bvec) * 1.0/(avec.magnitude())
        rvec = bvec - (avec*(cval/avec.magnitude()))
        dx = cval
        dy = rvec.magnitude()
        if (dy>=0 and dx>=0):
            theta = math.atan(dy/dx)
        elif (dy>=0 and dx<0):
            theta= (math.atan(dy/dx) + math.pi)
        elif (dy<0 and dx<=0):
            theta= math.atan(dy/dx) + math.pi
        else:
            theta= math.atan(dy/dx) + math.pi + math.pi
        return (bvec.magnitude() <= self.radius*(2+math.cos(2*theta)))

class Cylinder:
    def __init__(self, center_start, center_end, radius):
        self.center_start = center_start
        self.center_end = center_end
        self.radius = float(radius)
    def __repr__(self):
        return '[start center-line: ' + str(self.center_start) + ', end center-line: ' + str(self.center_end) + ', radius: ' + str(self.radius) + ']'
    def encompasses(self, point):
        avec = self.center_end - self.center_start
        bvec = point - self.center_start
        cval = avec.dot(bvec) * 1.0/(avec.magnitude())
        rvec = bvec - (avec*(cval/avec.magnitude()))
        return ((rvec.magnitude() - self.radius) <= 0 and cval >= 0 and cval <= avec.magnitude())

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

class WheelOne:
    def __init__(self, center_start, center_end, radius, height, shaft_frac, N_tread):
        self.center_start = center_start
        self.center_end = center_end
        self.radius = float(radius)
        self.height = float(height)
        self.shaft_frac = float(shaft_frac)
        self.N_tread = float(N_tread)
        
    def __repr__(self):
        return '[center_start: ' + str(self.center_start) + ', center_end: ' + str(self.center_end) + ', radius: ' + str(self.radius) + ', height: ' + str(self.height) + ', shaft_frac: ' + str(self.shaft_frac) + ', N_tread: ' + str(self.N_tread) + ']'
    def encompasses(self, point):
        avec = self.center_end - self.center_start
        bvec = point - self.center_start
        cval = avec.dot(bvec) * 1.0/(avec.magnitude())
        rvec = bvec - (avec*(cval/avec.magnitude()))
        #return ((rvec.magnitude() - self.radius) <= 0 and cval >= 0 and cval <= avec.magnitude())

        radius=self.radius
        height= self.height
        shaft_frac= self.shaft_frac
        N_tread= self.N_tread
        #define normal lying on z-plane
        if (avec.x == 0 and avec.y == 0):
            nvec = Vector(1,0,0)
        else:
            nvec = Vector(-avec.y,avec.x,0) * (1.0/math.sqrt(avec.x*avec.x + avec.y*avec.y))
        nvec2 = nvec.cross(avec)
        nvec2 *= 1.0/nvec2.magnitude()
        dx = rvec.dot(nvec)
        dy = rvec.dot(nvec2);
        is_inside=0
        p1x=p1y=p2x=p2y=p3x=p3y=p4x=p4y=0
        a1=a2=a3=a4=b1=b2=b3=b4=c1=c2=c3=c4=d1=d2=dist1=dist2=0


        if (rvec.magnitude() <= (radius+height) and cval >= 0 and cval <= avec.magnitude()):
            if (dy>=0 and dx>=0):
       	        theta = math.atan(dy/dx)
            elif (dy>=0 and dx<0):
                theta= (math.atan(dy/dx) + math.pi)
       	    elif (dy<0 and dx<=0):
                theta= math.atan(dy/dx) + math.pi
            else:
       	        theta= math.atan(dy/dx) + math.pi + math.pi

            if (math.sqrt(dy*dy + dx*dx) <=radius):
                is_inside = 1

            elif ( (theta % (2*math.pi/N_tread)) >= ((2*math.pi/N_tread)*shaft_frac) ):
	           
                TreadStart = theta -( theta % (2*math.pi/N_tread) ) + (2*math.pi/N_tread)*shaft_frac
                TreadEnd   = theta -( theta % (2*math.pi/N_tread) ) + (2*math.pi/N_tread)
                TreadMid   = theta -( theta % (2*math.pi/N_tread) ) + (2*math.pi/N_tread)*(shaft_frac+1)*0.5
                p1x= radius * math.cos(TreadStart)
                p1y= radius * math.sin(TreadStart)

                p4x= radius * math.cos(TreadEnd)
                p4y= radius * math.sin(TreadEnd)

                p2x= p1x + height * math.cos(TreadMid)
                p2y= p1y + height * math.sin(TreadMid)

                p3x= p4x + height*math.cos(TreadMid)
                p3y= p4y + height*math.sin(TreadMid)

        	## LINE 12 AND 34
                a1= p2y-p1y
                b1= p1x-p2x
                c1= p1y*p2x- p1x*p2y
                a2= p4y-p3y
                b2= p3x-p4x
                c2= p3y*p4x- p3x*p4y   
                d1= abs( (a1*dx + b1*dy + c1) / ( math.sqrt(a1*a1+b1*b1) ) )
                d2= abs((a2*dx + b2*dy + c2) / ( math.sqrt(a2*a2+b2*b2) ) )
                dist1= math.sqrt( (p2y-p3y)*(p2y-p3y) + (p2x-p3x)*(p2x-p3x) )
        	## LINE 23 AND 14
                a3= p2y-p3y
                b3= p3x-p2x
                c3= p3y*p2x- p3x*p2y
                a4= p4y-p1y
                b4= p1x-p4x
                c4= p1y*p4x- p1x*p4y   
                d3= abs(( a3*dx+b3*dy+c3 )/( math.sqrt(a3*a3+b3*b3) ) ) 
                d4= abs(( a4*dx+b4*dy+c4  )/( math.sqrt(a4*a4+b4*b4) ) )
                dist2= math.sqrt( (p1y-p2y)*(p1y-p2y) + (p1x-p2x)*(p1x-p2x) )
	          
                if ( ((d1+d2) <= dist1*1.01) and ((d3+d4) <= dist2*1.01)) : 
                    is_inside = 1
                else:
            	    is_inside = 0
            else:
                is_inside = 0
        else:
            is_inside = 0
	        
    ##################################

        return (is_inside) 
            # checking if point is within the wheel

    ##################################
    ##################################

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
