""" Module for various geometry entities
to be used in Civil Engineering
makes use of python's complex numbers as
x,y pairs.
"""

import cmath
from dataclasses import dataclass
import itertools
import math
from math import e, pi
import matplotlib.pyplot as plt
from matplotlib import patches

P = {} #Points
S = {} #Segments
C = {} #Curves
Ch ={} #Chains

def sign( in_val )->int:
    if in_val == 0: return 0
    if in_val >  0: return 1
    if in_val <  0: return -1


def normalize( angle:float) ->float:
    """ Return radians between 0 and 2pi
    for any input.
    """
    remain = math.fmod(angle,(2*pi))
    if remain < 0:
        return remain + 2*pi
    return remain

def norm_as_delta ( angle:float) ->float:
    """ normalize to between -pi and pi
    this is better for a delta in bearing
    not the bearing itself.  Note all bearings
    are still reachable.
    """
    return normalize( angle +pi ) - pi

def change_in_bearing( start:float, end:float)->float:
    return norm_as_delta(normalize(end) - normalize(start))

#class Geo:
#    pass

class Point(complex):
    ''' At new Point called with  Point(x,y) or as Point.from_complex( complex_number )
        internally the number is stored as a complex base class'''
    is_Point = True
    def __init__(self, x:float, y:float )->float:
        self.val = complex(x, y)
        self.X = x
        self.Y = y

    @classmethod
    def from_complex( cls, complex_numb:complex ):
        """Use this to define a Point you already have a complex number"""
        X = complex_numb.real
        Y = complex_numb.imag
        return cls(X,Y)

    def phase(self):
        """ returns the CCW  angle from Y=0 (East). This is not the
        same as the Bearing! """
        return cmath.phase(self)

    def __repr__(self):
        return f"Point x={self.val.real}, y={self.val.imag}"

    def __str__(self):
        return f"Point x={self.val.real}, y={self.val.imag}"

class Bearing(float):
    """ Angle from Y=0 (East) in the range of 0 to 2*pi
    we will want to add and subtract with changes in angles that will be in float"""

    is_Bearing = True

    def __init__(self, brg:float):
        self.val = brg

    def Opposite(self):
        """ Pi radians from the Bearing good for backsights in Civil
        Engineering"""
        return normalize( self.val + pi )

    def __repr__(self):
        return f"Bearing={self.val}"

@dataclass
class Ray:
    """Class for when we have a know Point and only one particular Bearing
    """
    is_Ray = True
    def __init__(self, pt:Point, brg:Bearing):
        self.Point = pt
        self.Bearing = Bearing(brg)

    def move_to(self, distance:float)->Point:
        ''' return a point on the ray, distance away from the ray origen
        '''
        _new_point = self.Point + cmath.rect(distance, self.Bearing)
        return Point.from_complex(_new_point)

    def equal(self, ob):
        """ Are two Rays equal"""
        if isinstance(ob, Ray) and \
           self.Point == ob.Point and \
           self.Bearing == ob.Bearing:
           return True
        return False

    def angle_to(self, inBearing:float )->float:
        start = inBearing
        end = self.Bearing
        return norm_as_delta(start - end )

    def distance_to(self, inPoint:Point)->float:
        ''' What is the distance along the Ray till the 
        point in quesition is at a right angle
        to the Ray '''
        #make a Segment to the new point
        #segments stored as 2 points Bearing is a method()
        mySeg = Segment(self.Point, inPoint)
        theta = self.angle_to(mySeg.Bearing() )
        return math.cos(theta) * mySeg.Len()

    def offset_to(self, inPoint:Point)->float:
        ''' Return the perpendicular distance from ray 
        to inPoint.'''
        mySeg = Segment(self.Point, inPoint)
        theta = self.angle_to(mySeg.Bearing() )
        return math.sin(theta) * mySeg.Len()
         
        


    def __repr__(self):
        _s = ""
        _s +=  "Ray:\n"
        _s += f"       {self.Point}\n"
        _s += f"       {self.Bearing}\n"
        _s += f"-- Ray End --"
        return _s


class Delta(float):
    """ Delta is a difference in two Bearings
    we are going to limit it to +/- pi/2
    """
    is_Delta = True
    def __init__(self, ang):
        if ang< -pi/2 or ang >pi/2:
            raise ValueError(
                "%s Delta must be between +/- pi/2" % ang)
        self.val = float(ang)

    #def __repr__(self):
    #    return f"Delta= {self.val} rad"

class Line:
    is_Line = True

    def __init__(self, pt:Point, bearing:Bearing):
        self.Point = pt
        self.Bearing = Bearing(bearing)
    #need to be able to accept an offset

    def __repr__(self):
        _s = ""
        _s +=  "Line:\n"
        _s += f"       Point= {self.Point}"
        _s += f"       Bearing= {self.Bearing}\n"
        _s += f"--"
        return _s

class Segment(Ray):
    is_Segment = True
    def __init__(self, Pt1:Point, Pt2:Point):
        self.Pt1 = Pt1
        self.Pt2 = Pt2
    def Movement(self)->complex:
        """Subtract destination from start to get movement"""
        return self.Pt2 - self.Pt1

    def Len(self)->float:
        """Length of Segment """
        return abs( self.Movement() )

    def Bearing(self)->Bearing:
        """Bearing of Segment (phase is the Bearing but it is stored
        as a Bearing"""
        return  Bearing( cmath.phase( self.Movement() ))

    def inRay(self)->Ray:
        """Point and Bearing for the first Point """
        return Ray(self.Pt1, self.Bearing() )

    def outRay(self)->Ray:
        """Point and Bearing for the second Point """
        return Ray(self.Pt2, self.Bearing() )

    def delta_to(self, inBearing)->float:
        return norm_as_delta(self.Bearing() - inBearing)



    #need to be able to accept an offset

    def __repr__(self):
        s = ""
        s +=  "Segment:\n"
        s += f"        Pt1= {self.Pt1}\n"
        s += f"        Pt2= {self.Pt2}\n"
        s += f"        L()= {self.Len()}\n"
        s += f"--"
        return s

    def as_Line(self)->Line:
        """ Make a Line from a Point and a Bearing """
        return Line(self.Pt1, self.Bearing() )


    def patch(self, **kwargs):
        """Use matplotlib.patches for consistency, method patches.Arc
        works best for class Curve so we will also use method patches.Polygon
        for the Segemnts. """
        return patches.Polygon( ( (self.Pt1.X, self.Pt1.Y), \
                                  (self.Pt2.X, self.Pt2.Y), \
                                ), fill=False, closed=False, **kwargs )

class Curve:
    is_Curve = True

    """Curve is defined by a (PC, CC, Delta),
    PC != CC and -pi< Delta < pi and not zero
    """

    def __init__(self, point1:Point, point2:Point,  angle , *args, **kwargs):
        if abs( point1 - point2) < 0.1:
            raise ValueError(
                "%s.__new__ radius is too small" % self)
        if abs(angle == 0):
            raise ValueError(
                "%s.__new__ Delta can not be zero" % self)

        self.PC = point1
        self.CC = point2
        self.PC_to_CC = self.CC - self.PC
        self.PC_to_CC_asSeg = Segment(self.CC, self.PC)

        self.R = abs( point1 - point2 )
        self.Delta = angle
        #self.Len = self.R * abs(self.Delta)
        # E distance from curve to PI
        self.E = self.R * (1/math.cos(self.Delta/2.0) - 1)
        bearing_CC_to_PI = cmath.phase(self.PC - self.CC) + self.Delta/2.0
        CC_to_PI = cmath.rect(self.R+self.E, bearing_CC_to_PI)
        self.PI = Point.from_complex( self.CC + CC_to_PI )

        self.CC_to_PC = -1 * self.PC_to_CC
        #rotation of the CC_to_PC
        #same as  mult by   cos(angle) + i*(sin(angle))
        self.CC_to_PT =  self.CC_to_PC * e **complex(0, angle )

        # two PC is global other two are movements
        self.PT = Point.from_complex( self.PC + self.PC_to_CC + self.CC_to_PT )
        self.Chord = Segment( self.PC, self.PT )
        self.T =  self.R * math.tan( self.Delta / 2.0)
        self.inBearing =  self.Chord.Bearing() - self.Delta /2.0
        self.outBearing = self.inBearing + self.Delta
    #need to be able to accept an offset

    def Len(self)->float:
        return self.R * abs(self.Delta)


    @classmethod
    def from_PC_bearing_R_Delta( cls, pc:Point, brg:float, R:float, delta:float ):
        #pc given
        #find cc
        cc = Point.from_complex( pc + cmath.rect( R, brg +sign(delta)*pi/2 ) )
        return cls(pc, cc, delta)

    @classmethod
    def from_Ray_T_Delta( cls, _inRay:Ray, _T:float, delta:float ):
        _pc = _inRay.Point
        _start_bearing= _inRay.Bearing
        _R = math.cos( delta/2.0 ) * _T
        return cls.from_PC_bearing_R_Delta(_pc, _start_bearing, _R, delta)





    def inRay(self)->Ray:
        return Ray(self.PC, Bearing(self.inBearing))
    def outRay(self)->Ray:
        return Ray(self.PT, Bearing(self.outBearing))


    def patch(self,**kwargs):
        x = self.CC.X
        y = self.CC.Y
        brg_CC_to_PC = math.degrees( cmath.phase(self.CC_to_PC) )
        brg_CC_to_PT = math.degrees( cmath.phase(self.CC_to_PT) )
        r = 2 * self.R #actually the axis
        if self.Delta > 0:
            brg_start = brg_CC_to_PC
            brg_end   = brg_CC_to_PT
        else:
            brg_start = brg_CC_to_PT
            brg_end   = brg_CC_to_PC
        return patches.Arc( (x, y),r, r, 0, brg_start, brg_end, **kwargs )


    def __repr__(self):
        rep =  "Curve:\n"
        rep += f"PC= {self.PC}\n"
        rep += f"CC= {self.CC}\n"
        rep += f"PT= {self.PT}\n"
        rep += f" R= {self.R}\n"
        rep += f"Delta= {self.Delta}\n"
        rep += f"Len()= {self.Len()}\n"
        rep += f"Chord= {self.Chord}\n"
        #rep += f"-- Curve end --\n"
        # the segment finish up with a -- so last line not needed
        return rep

class Chain:
    '''
    '''
    is_Curve = True
    def __init__(self, name,*args, StartSta=None, **kwargs):
        self.Name = name
        self.Routes = []
        self.RoutesSta = []

        if StartSta == None:
            self.StartSta=0.0
        else:
            self.StartSta=StartSta

        for ob in args:
            self.addRoute(ob)

    def addRoute( self, ob):
        if self.validate_add(ob):
            self.Routes.append(ob)
            #start station of each segment

            x = []
            for i in self.Routes:
                x.append( i.Len() )
            start = self.StartSta
            self.RoutesSta = list(itertools.accumulate(x, initial=start))
        else:
            raise ValueError(
                "%s This route can't be added here" % self )

    def forward( self, distance:float)->Segment:
        '''Add a new segment to the Chain using the outRay of the
        previously last route that was defined
        '''
        if len(self.Routes) == 0:
            raise ValueError(
                "forward() method not allowd for Chain with no routes"  )

        lastRay = self.Routes[-1].outRay()
        pt2 = lastRay.move_to( distance )
        self.addRoute( Segment( lastRay.Point, pt2 ))
        return self.Routes[-1].outRay


        

    def validate_add( self, ob)->bool:
        if not isinstance(ob, (Segment,Curve)):
            raise ValueError(
                "%s not a Segment or Curve" % self )

        if len(self.Routes)==0: #first add needs no checking
            return True
        else:
            lastRay = self.Routes[-1].outRay()
            if ob.inRay() == lastRay:
                return True
            else:
                raise ValueError(
                    "%s inRay needs to equal outRay" % self)
        return False





def main():
    points = [ Point(0,0), Point(100,0)]

    curves = [ Curve( Point(11,10), Point(11,15), Delta(pi/4) ) ]

    PIs = [ Point(0,0), Point(10,0), Point(20,3), Point(30,4), \
            Point(40,4), Point(50,5), Point(60,5), ]

    #import X30_alignment_dict

    p1= Point( 100, 0)
    print(p1)
    p2 = Point( 90, 0)
    print(p2)
    p3 = Point.from_complex( complex(80,25) )
    print(p3)
    b = Bearing( pi/4 )
    print(b)
    d = Delta( -pi/4 )
    print(d)
    c = Curve(p1, p2, d)
    print(f"\n\n{c}")
    s1 = Segment( c.PT, p3)
    print("\n------\n")
    L1 =  s1.as_Line()
    print (L1)

    def as_XY( list_of_points ):
        x = list()
        y = list()
        for i in list_of_points:
            x.append(i.X)
            y.append(i.Y)
        return (x,y)

    pts = [ p1, p2, p3, c.PC, c.CC, c.PT]
    fig, ax = plt.subplots()
    ax.scatter( as_XY(pts)[0],  as_XY(pts)[1] )

    segPatch = Segment( c.CC, c.PT).patch()
    print(f"{segPatch=}")
    circlePatch = c.patch()
    print(f"{circlePatch=}")

    ax.add_patch( Segment( c.PC, c.CC).patch(linestyle=':') )
    ax.add_patch( Segment( c.CC, c.PT).patch(linestyle=':') )
    ax.add_patch( c.patch(linestyle='-') )
    ax.add_patch( s1.patch(linestyle='-') )
    plt.axis('scaled')
    plt.show()

    mychain = Chain("Ch_1", s1)
    mychain.StartSta = 12044
    #print(f"{mychain=}")
    print(f"{mychain.Name}")
    print(f"{mychain.Routes}")

    print(  s1 )
    s2 = Segment( s1.Pt2, s1.Pt2 + s1.Movement() )
    print(  s2 )
    mychain.addRoute( s2 )
    print( s1.outRay() ,  s2.inRay() )
    s3 = Segment( s2.Pt2, s2.Pt2 + s2.Movement() )
    print(  s3 )
    mychain.addRoute( s3 )

    print(f"{mychain.Routes=}")
    print(f"{mychain.RoutesSta=}")




if __name__ == "__main__":
    main()



