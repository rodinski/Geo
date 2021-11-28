import  math
from math import e, pi
from dataclasses import dataclass
import itertools 

import cmath
import matplotlib.pyplot as plt
from   matplotlib.path import Path
import matplotlib.patches as patches 
#use Arc     for Curves and 
#    Polygon for Segments 

pi2 = pi / 2.0
pi4 = pi / 4.0

P = dict() #Points
S = dict() #Segments
C = dict() #Curves
Ch = dict() #Chains

def normalize( angle:float) ->float:
    """ Return radians between 0 and 2pi
    """
    remain = math.fmod(angle,(2*pi))
    if remain < 0:
        return remain + 2*pi
    else:
        return remain


class Geo:
    pass

class Point(complex):
    ''' At new Point called with  Point(x,y) or as Point.from_complex( complex_number )
        internally the number is stored as a complex base class'''
    is_Point = True
    def __init__(self, x:float, y:float,  *args, **kwargs)->float:
        self.val = complex(x, y)
        self.X = x
        self.Y = y
 
    @classmethod
    def from_complex( cls, complex_numb:complex ):
        X = complex_numb.real
        Y = complex_numb.imag
        return cls(X,Y)

    def phase(self):
        return cmath.phase(self)

    def __repr__(self):
        return(f"Point x={self.val.real}, y={self.val.imag}")
    
    def __str__(self):
        return(f"Point x={self.val.real}, y={self.val.imag}")

class Bearing(float):
    is_Bearing = True

    def normalRange(self, phase ):
        return normalize( phase )

    def __init__(self, brg, *args, **kwargs):
        self.val = self.normalRange(brg)  

    def Opposite(self):
        return normalize( self.val + pi )

    def __repr__(self):
        return f"Bearing={self.val}"

@dataclass
class Ray:
    is_Ray = True
    def __init__(self, pt:Point, brg:Bearing):
        self.Point = pt
        self.Bearing = Bearing(brg)
    def equal(self, ob):
        if isinstance(ob, Ray) and \
           self.Point == ob.Point and \
           self.Bearing == ob.Bearing:
           return True
        else:
           return False
    def __repr__(self):
        s = ""
        s +=  "Ray:\n"
        s += f"       {self.Point}\n"
        s += f"       {self.Bearing}\n"
        s += f"--\n"
        return s 


class Delta(Bearing):
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
    
    def __init__(self, pt:Point, bearing:Bearing, *args, **kwargs):
        self.Point = pt
        self.Bearing = Bearing(bearing)
    #need to be able to accept an offset 

    def __repr__(self):
        s = ""
        s +=  "Line:\n"
        s += f"       Point= {self.Point}"
        s += f"       Bearing= {self.Bearing}\n"
        s += f"--\n"
        return s 

class Segment:
    is_Segment = True
    def __init__(self, Pt1:Point, Pt2:Point, *args, **kwargs):
        self.Pt1 = Pt1 
        self.Pt2 = Pt2 
    def Movement(self)->complex:
        return self.Pt2 - self.Pt1
    def Len(self)->float:
        return abs( self.Movement() )
    def Bearing(self)->Bearing:
        return  Bearing( cmath.phase( self.Movement() ))

    def inRay(self)->Ray:
        return Ray(self.Pt1, self.Bearing() )

    def outRay(self)->Ray:
        return Ray(self.Pt2, self.Bearing() )

    #need to be able to accept an offset 

    def __repr__(self):
        s = ""
        s +=  "Segment:\n"
        s += f"        Pt1= {self.Pt1}\n"
        s += f"        Pt2= {self.Pt2}\n"
        s += f"          L= {self.Len()}\n"
        s += f"--\n"
        return s

    def as_Line(self)->Line:
        return Line(self.Pt1, self.Bearing() )


    def patch(self,**kwargs):
        '''Use matplotlib.patches for consistency 
        patches.Arc works best for clase Curve so you Polygon for the segemnts
        '''
        return patches.Polygon( ( (self.Pt1.X, self.Pt1.Y), \
                                  (self.Pt2.X, self.Pt2.Y), \
                                ), fill=False, closed=False, **kwargs ) 

class Curve:
    is_Curve = True
    ''' Curve is defined by a (PC, CC, Delta), since the arc will always tend to the CC, 
    the notion of a negative Delta makes no sense so we will only allow 0 < Delta < 2pi


    '''
    def __init__(self, point1:Point, point2:Point,  angle:Delta , *args, **kwargs):
        if abs( point1.val - point2.val) < 0.1:
            raise ValueError(
                "%s.__new__ radius is too small" % self.__name__)
        if abs(angle == 0):
            raise ValueError(
                "%s.__new__ Delta can not be zero" % self.__name__)
 
        self.PC = point1
        self.CC = point2
        self.PC_to_CC = self.CC - self.PC
        self.PC_to_CC_asSeg = Segment(self.CC, self.PC)

        self.R = abs( point1 - point2 )
        self.Delta = angle.val
        self.L = self.R * abs(self.Delta)
        # E distance from curve to PI
        self.E = self.R * (1/math.cos(self.Delta/2.0) - 1)
        bearing_CC_to_PI = cmath.phase(self.PC - self.CC) + self.Delta/2.0
        CC_to_PI = cmath.rect(self.R+self.E, bearing_CC_to_PI)
        self.PI = Point.from_complex( self.CC + CC_to_PI )

        self.CC_to_PC = -1 * self.PC_to_CC 
        #rotation of the CC_to_PC 
        #same as  mult by   cos(angle) + i*(sin(angle))
        self.CC_to_PT =  self.CC_to_PC * e **complex(0, angle.val )

        # two PC is global other two are movements
        self.PT = Point.from_complex( self.PC + self.PC_to_CC + self.CC_to_PT )
        self.Chord = Segment( self.PC, self.PT )
        self.T =  self.R * math.tan( self.Delta / 2.0) 
        self.inBearing =  self.Chord.Bearing() - self.Delta /2.0
        self.outBearing = self.inBearing + self.Delta
    #need to be able to accept an offset 

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
        rep += f"   L= {self.L}\n"
        rep += f"Chord= {self.Chord}\n"
        rep += f"--\n"
        return rep

class Chain:
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
            for i in self.Routes[1:]:
                x.append( i.Len() )
            start = self.StartSta
            self.RoutesSta = list(itertools.accumulate(x, initial=start))
        else:
            raise ValueError(
                "%s.__new__ This route can't by added here" % self.__name__ )

    def validate_add( self, ob)->bool:
        if not isinstance(ob, (Segment,Curve)):
            raise ValueError(
                "%s.__new__ not a Segment or Curve" % self.__name__ )
            return False

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

if __name__ == "__main__":

    points = [ Point(0,10), Point(10,10)]
    curves = [ Curve( Point(11,10), Point(11,15), Delta(pi/4) ) ] 
    PIs = [ Point(0,0), Point(10,0), Point(20,3), Point(30,4), \
            Point(40,4), Point(50,5), Point(60,5), ] 

    p1= Point( 100, 0)
    print(p1)
    p2 = Point( 90, 0)
    print(p2)
    p3 = Point.from_complex( complex(80,25) )
    print(p3)
    b = Bearing( pi4 ) 
    print(b)
    d = Delta( -pi4 ) 
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
    #plt.show()

    mychain = Chain("Ch_1", s1)
    mychain.StartSta = 12044
    print(f"{mychain=}")
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
    

    
