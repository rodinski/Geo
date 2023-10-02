""" Module for various geometry entities
to be used in Civil Engineering
makes use of python's complex numbers as
x, y pairs.
"""

import cmath
from dataclasses import dataclass
import math
from math import e, pi, cos, sin, tan
import matplotlib.pyplot as plt
from matplotlib import patches
from collections import namedtuple
import itertools
import IPython

P = {}  # Points
S = {}  # Segments
C = {}  # Curves
Ch = {}  # Chains
Dist_os_tup = namedtuple("Dist_os_tup", "distance offset")


def sign(in_val) -> int:
    ''' Return a value based only on the numerical sign or zero of
    the input value'''
    if in_val == 0:
        return int(0)
    if in_val > 0:
        return int(1)
    if in_val < 0:
        return int(-1)

def normalize(angle: float) -> float:
    """ Return radians between 0 and 2pi
    for any input.
    """
    remain = math.fmod(angle, (2*pi))
    if remain < 0:
        return remain + 2*pi
    return remain

def norm_as_delta(angle: float) -> float:
    """ normalize to between -pi and pi
    this is better for a delta in bearing
    not the bearing itself.  Note, all bearings
    are still reachable."""
    return normalize(angle + pi) - pi

def conjugate_angle(inAngle):
    if abs(inAngle) > 2*pi:
        ValueError: print("Angle out of range in congugate_angle function")
    if inAngle == 0:
        return 2*pi
    # change the sign b/c we would going other way around the circle
    return -1.0*sign(inAngle)*(2*pi-abs(inAngle))

def change_in_bearing(start: float, end: float) -> float:
    ''' Returns the change in bearing needed to get from start to end'''
    return norm_as_delta(normalize(end) - normalize(start))

def xy(points_list: list) -> tuple:
    """Return a list of .x and a list of .y properites
    for the given points_list.
    Useful for use in matplotlib plts and scatter.
    Such as  ax.scatter(*xy(pts), marker='+', color='g')
    """
    x = [arg.X for arg in points_list]
    y = [arg.Y for arg in points_list]
    return(x, y)

def as_XY(list_of_points):
    """From a list of Points return a tuple containing
    two list   ([x's], [y's])"""
    x = list()
    y = list()
    for i in list_of_points:
        x.append(i.X)
        y.append(i.Y)
    return (x, y)

class Point(complex):
    """ At new Point called with  Point(x, y) or as 
    Point.from_complex(complex_number)
    internally the number is stored as a complex base class
    """
    is_Point = True
 

    def __init__(self, x: float, y: float) -> float:
        Decimals = 4
        x = round(x, Decimals)
        y = round(y, Decimals)
        self.val = complex(x, y)
        self.X = x
        self.Y = y

    @classmethod
    def from_complex(cls, complex_numb:complex):
        """Use this to define a Point you already have a complex number
        """
        X = complex_numb.real
        Y = complex_numb.imag
        return cls(X, Y)

    def phase(self):
        ''' returns the CCW  angle from Y=0 (East). This is not
        same as the Bearing! '''
        return cmath.phase(self)

    # def move_to(self, inPoint):
    #    return cmath.polar(inPoint - self.val)

    def __repr__(self):
        return f"Point x={self.val.real:.4f}, y={self.val.imag:0.4f}"

    def __str__(self):
        return f"Point x={self.val.real:.4f}, y={self.val.imag:0.4f}"

class Angle(float):
    '''Angle will be stored as a float with a nice __str__().
    All math funtions will return a float without the __str__()
    if unit='deg' is used the input should be in degrees and a
    conversion will happen so the internal storage will be radians
    '''
    is_float = True


    def __new__(cls, real, unit='rad'):
        if unit.lower() == 'rad':
            real = float(real)
        if unit.lower() == 'deg':
            real = float(math.radians(real))
        return float.__new__(cls, real)


    def __init__(self, real: float, unit='rad'):
        # .real attibute should have been set in __new__
        #
        self.unit = unit

    def __str__(self):
        return '%g rad  %g deg' % (self, math.degrees(self))

    def __repr(self):
        return "Angle({self:,f})"

class Distance(float):
    is_float = True

    def __init__(self, flt: float):
        Decimals = int(6)
        self.val = round(flt, Decimals)

    def __repr__(self):
        return f"Distance({self.val:,f})"

class Bearing(float):
    """ Angle from Y=0 (East) in the range of 0 to 2*pi
    we will want to add and subtract with changes in angles that will be in float"""

    is_Bearing = True

    def __init__(self, brg: float):
        self.val = brg

    def Opposite(self):
        """Pi radians from the Bearing good for backsights in Civil
        Engineering.
        """
        return normalize(self.val + pi)

    def lt_90(self):
        """Return the bearing that is 90 CCW"""
        return self.val + pi/2

    def rt_90(self):
        """Return the bearing that is 90 CW """
        return self.val - pi/2

    def __repr__(self):
        return f"Bearing={self.val}"

@dataclass
class Ray:
    """Class for when we have a know Point and only one particular Bearing
    """
    is_Ray = True

    def __init__(self, pt: Point, brg: Bearing):
        self.Point = pt
        self.Bearing = Bearing(brg)

    def move_to(self, distance: float, offset=0.0) -> Point:
        ''' return a point on the ray, distance away from the ray origen
        '''
        _new_point = self.Point + cmath.rect(distance, self.Bearing) \
                     + cmath.rect(offset, normalize(self.Bearing - pi/2.0))
        return Point.from_complex(_new_point)

    def multipoint(self, start: float,
                    spacing: float = 0, 
                    n_steps: int = 0,
                   point_at_start: bool=False) -> tuple:
        """Return a tuple of points input
        start distance from Pt1 along the bearing, return of a point here
        is controlled by point_at_start=T/F (default=False)
        spacing is a distance to the next point
        n_steps integer number of points to produce
        """
        ret = ()
        if point_at_start:
            ret.append(self.move_to(start))
        for i in range(n_steps):
            distance = start + (i+1)*spacing
            ret.append(self.move_to(distance))

    def equal(self, ob):
        ''' Are two Rays equal?'''
        if isinstance(ob, Ray) and \
           self.Point == ob.Point and \
           self.Bearing == ob.Bearing:
           return True
        return False

    def angle_to(self, destBearing: float) -> float:
        ''' Return the turn angle needed to get at some destination bearing.
        Output limited to +/- pi '''
        start = destBearing
        end = self.Bearing
        return norm_as_delta(start - end)

    def distance_to(self, inPoint: Point) -> float:
        ''' What is the distance along the Ray till the
        point in quesition is at a right angle
        to the Ray '''
        # make a Segment to the new point
        # segments stored as 2 points Bearing is a method()
        mySeg = Segment(self.Point, inPoint)
        theta = self.angle_to(mySeg.Bearing())
        ret = math.cos(theta) * mySeg.Len()
        return ret
    
    def offset_to(self, inPoint: Point) -> float:
        #Return the perpendicular distance from ray 
        #to inPoint.
        mySeg = Segment(self.Point, inPoint)
        theta = self.angle_to(mySeg.Bearing()) * -1 #(+) offsets have negative theta
        #-1 above needed to get sign of offset correct.
        offset = math.sin(theta) * mySeg.Len()
        return offset

    def copy_parallel(self, offset):
        return Ray( self.move_to(0, offset), self.Bearing)

    def patch(self, scale=1, width=1,  **kwargs):
        """Use matplotlib.patches for consistency, method patches.Arc
        works best for class Curve so we will also use method patches.Polygon
        for the Segemnts. """

        return patches.Arrow (self.Point.X, self.Point.Y, \
                             scale *cos(self.Bearing), scale *sin(self.Bearing), \
                             width,  **kwargs)

        '''to inPoint. +(pos) = left, -(neg) = right'''
        #mySeg = Segment(self.Point, inPoint)
        #theta = self.angle_to(mySeg.Bearing())
        #return math.sin(theta) * mySeg.Len()

    def offset(self, dist:float) -> Point:
        ''' Return a new point perpendicular to the bearing a distance
        of D from defining Point of the Ray. Do this my turning +90 deg
        and moving the distance noted. -(neg) will move backwards.'''
        newPoint = self.Point + cmath.rect(dist, self.Bearing.lt_90())
        return Point.from_complex(newPoint)


    def __repr__(self):
        _s = ""
        _s +=  "Ray:\n"
        _s += f"       {self.Point}\n"
        _s += f"       {self.Bearing}\n"
        _s += f"-- Ray End --"
        return _s

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

    def Movement(self) -> complex:
        '''Subtract destination from start to get movement'''
        return self.Pt2 - self.Pt1

    def Len(self) -> float:
        '''Length of Segment '''
        return abs(self.Movement())

    def Bearing(self) -> Bearing:
        """Bearing of the segment """
        return  Bearing(cmath.phase(self.Movement()))

    def inRay(self) -> Ray:
        """Point and Bearing for the first Point """
        return Ray(self.Pt1, self.Bearing())

    def outRay(self) -> Ray:
        """Point and Bearing for the second Point """
        return Ray(self.Pt2, self.Bearing())

    def delta_angle_to_some_bearing(self, inBearing) -> float:
        """ Reurn the turning angle needed from the inRay to some 
        given bearing"""
        return norm_as_delta(self.Bearing() - inBearing)

    def distance_offset(self, obPoint: Point) -> (float, float):
        if not (0.0 <= self.inRay().distance_to(obPoint) <= self.Len()):
           return Dist_os_tup(None, None)
        else:
           distance = self.inRay().distance_to(obPoint)
           offset = self.inRay().offset_to(obPoint)
           return Dist_os_tup(distance, offset)


    def move_to(self, distance: float, offset=0.0) -> Point:
        if distance > self.Len():
            raise ValueError(
                    f"Distance of {distance} is greater than length of Segment={self.Len()}")
        if distance < 0.0:
            raise ValueError(
                    f"Distance of {distance} is negative")
        _new_point = self.Pt1 + cmath.rect(distance, self.Bearing()) \
                     + cmath.rect(offset, normalize(self.Bearing() - pi/2.0))
        return Point.from_complex(_new_point)

    def copy_parallel(self, offset):
        """Return copy of the Segment, offset has no limits
        """
        _pt1 = self.move_to(0, offset)
        _pt2 = self.move_to(self.Len(), offset)
        return Segment(_pt1, _pt2)

    def __repr__(self):
        s = ""
        s +=  "Segment:\n"
        s += f"        Pt1= {self.Pt1}\n"
        s += f"        Pt2= {self.Pt2}\n"
        s += f"        L()= {self.Len()}\n"
        s += f"--"
        return s

    def as_Line(self) -> Line:
        """ Make a Line from a Point and a Bearing """
        return Line(self.Pt1, self.Bearing())


    def patch(self, **kwargs):
        """Use matplotlib.patches for consistency, method patches.Arc
        works best for class Curve so we will also use method patches.Polygon
        for the Segemnts. """
        return patches.Polygon(((self.Pt1.X, self.Pt1.Y), \
                                  (self.Pt2.X, self.Pt2.Y), \
                                ), fill=False, closed=False, **kwargs)

class Curve:
    is_Curve = True

    """
    Curve is defined by a (PC, CC, Delta), PC != CC and
    -2pi< Delta < 2pi and not zero
    """

    def __init__(self, point1:Point, point2:Point,  Delta:Angle , *args, **kwargs):
        if abs(point1 - point2) < 0.1:
            raise ValueError(
                "%s.__new__ radius is too small" % self)
        if abs(Delta == 0):
            raise ValueError(
                "%s.__new__ Delta can not be zero" % self)


        self.PC = point1   #Point
        self.CC = point2   #Point
        self.PC_to_CC = self.CC - self.PC   #Complex
        self.CC_to_PC = self.PC - self.CC   #Complex
        self.PC_to_CC_asSeg = Segment(self.CC, self.PC) #Segment

        self.R = Distance(abs(point1 - point2))
        self.Delta = Delta
        # E distance from curve to PI, value error on cos(pi/2.0))
        self.E = self.R * (1.0 / math.cos(Delta/2.0) - 1.0) #
        bearing_CC_to_PI = cmath.phase(self.PC - self.CC) + self.Delta/2.0
        CC_to_PI = cmath.rect(self.R+self.E, bearing_CC_to_PI)
        self.PI = Point.from_complex(self.CC + CC_to_PI)

        # rotation of the CC_to_PC using imag part of a complex number
        # same as  mult by   cos(Delta) + i*(sin(Delta))
        self.CC_to_PT =  self.CC_to_PC * e **complex(0, self.Delta)

        # Point.from_complex is still a point, PC->to_CC->to_PT is two steps
        # could have been CC->to_PT
        self.PT = Point.from_complex(self.PC + self.PC_to_CC + self.CC_to_PT) #Point

        self.Chord = Segment(self.PC, self.PT)
        self.T =  Distance(self.R * math.tan(self.Delta / 2.0))
        self.inBearing =  self.Chord.Bearing() - self.Delta /2.0
        self.outBearing = self.inBearing + self.Delta
    #need to be able to accept an offset

    def Len(self) -> float:
        return Distance(self.R * abs(self.Delta))


    @classmethod
    def from_PC_bearing_R_Delta(cls, pc:Point, brg:float, R:float, delta:float):
        #pc given
        #find cc
        cc = Point.from_complex(pc + cmath.rect(R, brg + sign(delta)*pi/2))
        return cls(pc, cc, delta)

    @classmethod
    def from_Ray_T_Delta(cls, _inRay:Ray, _T:float, delta:float):
        _pc = _inRay.Point
        _start_bearing= _inRay.Bearing
        _R = math.cos(delta/2.0) * _T
        return cls.from_PC_bearing_R_Delta(_pc, _start_bearing, _R, delta)

    def inRay(self) -> Ray:
        return Ray(self.PC, Bearing(self.inBearing))

    def outRay(self) -> Ray:
        return Ray(self.PT, Bearing(self.outBearing))

    def is_offset_valid(self, offset: float) -> bool:
        """
        When copying a paralell offset the offset has limits.
        CCW (positive) curves must have an offset larger than -R.
        CW (negative) curves must have an offset smaller than R.
        -offsets are to the left as traveling in direction of the curve.
        +offsets are to the right.
        """
        if sign(self.Delta) > 0:
            return offset > -self.R
        else:
            return offset < self.R


    #  fix this and inversings return need more consistency
    def has_dist_os(self, obPoint:Point) -> str:
        """
        For the input point determin if it is within the sweep of the curve
        A curve sweep may cross over the -/+Pi line and the logic becomes difficutly.
        Howerver all curves have a non-curve part of the whold cirle.  Alwasy one side
        sweep crosses the -/+Pi line and the other side does not.  Try to work with the
        side that does not.
        """

        # logic is easy to see if point is within the wedge
        pc_brg = cmath.phase(self.CC_to_PC)
        defined_curve = -pi <= (pc_brg + self.Delta)  <= +pi
        ob_phase = cmath.phase(obPoint-self.CC)

        if defined_curve:  #curve does not cross West axis
            domain = sorted([ pc_brg , pc_brg+self.Delta ]) #put in numerical order
            #check this wedge, if present return True
            #import pdb; pdb.set_trace()
            if domain[0] <= ob_phase <= domain[1]:
                return True
            else:
                return False
        else:  #ok so it will be easier to deal with the other (negative) wedge
            #switch the returns from above
            #need value of opposite curve AND need to switch the sign
            opposite_delta = (2.0*pi - abs(self.Delta)) * (-1*sign(self.Delta))
            domain2 = sorted([ pc_brg , pc_brg +opposite_delta ]) #put in numerical order
            if domain2[0] <= ob_phase <= domain2[1]:
                return False
            else:
                return True

    def distance_offset(self, obPoint:Point) -> (float, float):
        '''
        Return a tuple of  (distance, os) for a given point.
        Return (None, None) when the point is not within the domain of the curve
        The curves domain may be numerically discontinuous over the
        +pi / -pi line
        
        Use Curve.has_dist_os()  to check if we are in the curves domain.

        Note all curves are one part of a total circle. Of the two curves
        that complete a circle one will (and one will not) cross the 
        +pi/-pi line. This logic is used to simplify finding the angle
        from the PC_CC_point of interest angle.

        The sign of Delta is important

        Start of the curve domain is the phase of the CC to PC line.
        the domain increases for left (ccw) turns
        the domain decreases for right (cw) turns
        Delta has a domain of -2pi to 2pi 

        helper functions: has_dist_os() and conjugate_angle
        '''

        mySegment = Segment(self.CC, obPoint)
        # print(f"{Curve.has_dist_os(self, obPoint)=}")
        if not Curve.has_dist_os(self, obPoint):
            return Dist_os_tup(None, None)
        # find the angle  PC_CC_obpoint  and the compliment of this angle.
        # not which one has the same sign as the self.Delta
        brg_pc = cmath.phase(self.CC_to_PC)
        brg_obPoint = cmath.phase(obPoint-self.CC)
        # find the angle that matches the sign of self.Delta
        if sign(self.Delta) == sign(brg_obPoint - brg_pc):
            angle = brg_obPoint - brg_pc
        else:
            # need the explementary (or conjugate) angle
            angle = conjugate_angle(brg_obPoint - brg_pc)
        dist = Distance(abs(angle) * self.R)
        offset = Distance((self.R - mySegment.Len())* sign(self.Delta) * -1) # -1 needed to get offset sign correct
            # R and Len() do not have sign, the only sign has to do with Delta
            # we need to use correctly use sign(Delta) to get Lt(-) vs Rt(+) sedt correctly
        return Dist_os_tup(dist, offset)

    def normal_at_point(self, obPoint: Point) -> Bearing:
        if not self.has_dist_os(obPoint):
            return None
        return Bearing(cmath.phase(obPoint - self.CC))

    def tangent_at_point(self, obPoint:Point) -> Bearing:
        if not self.has_dist_os(obPoint):
            return None
        norm = self.normal_at_point(obPoint)
        if sign(self.Delta > 0):
            turn = pi/2
        else:
            turn = -pi/2
        return Bearing(norm + turn)

    def copy_parallel(self, offset):
        """Return a curve with a changed radius"""
        # check if offset results in a valid R
        # do we add or subtact the offset

        if not self.is_offset_valid(offset):
            raise ValueError(
                    f"This {offset=} will result in a negative R")
        _D = self.Delta
        _D_sign = sign(self.Delta)
        if _D_sign == 1:
            _new_R = self.R + offset
        else:
            _new_R = self.R - offset


        # make that new curve
        # print(f"{self.R=}  {_new_R=}")
        _cc = self.CC
        _pc = self.inRay().copy_parallel(offset).Point


        # from_PC_bearing_R_Delta(cls, pc:Point, brg:float, R:float, delta:float)
        return Curve(_pc, _cc, _D)


    '''def distance_and_offset(self, obPoint:Point)->(float, float):

        pc_brg = cmath.phase(self.CC_to_PC)
        start1 = pc_brg

        #Wwedges that cross the -/+pi boundary have two seperate domains
        #that need checking. Maybe we can avoid this check . . .

        #Two pie-shaped wedge pieces, only one of them can cross the -/+pi line
        #only work with the pieces that does not cross this line

        #did we stay in the boundary    -Pi < PC+Delta < +Pi
        
       
        if  -pi <=  (pc_brg + self.Delta) <= +pi: #logic is easy to see if point is within the wedge
            domain = sorted([ start1 , pc_brg+self.Delta ])
            if   not domain[0] <= obPoint.phase() <= domain[1]:
                return Dist_os_tup(None, None)
        
        else:
            #intead of setting two valid ranges can we test the not condition
            #the outer wedge? 

            #set two valid ranges domain
            end1 = pi * sign(self.Delta)  #-pi for cw and +pi for ccw
            start2 = -end1
            remaining_turn = (start1 + self.Delta) - end1  
            end2 = start2 +  remaining_turn
            domain = sorted([start1, end1]) + sorted([start2, end2])
            if not (domain[0] <= obPoint.phase() <= domain[1] 
               or 
               domain[2] <= obPoint.phase() <= domain[3]):
               return Dist_os_tup(None, None)



        # Work form the curve CC, create a Segment to obPoint
        mySegment = Segment(self.CC, obPoint)

        #start from two extremes and work back to the test point, if they 
        # travel in seperate ways the test must have been beween them
        mid = mySegment.Bearing()
        a = Segment(self.CC, self.PC).Bearing()  
        #b = Segment(self.CC, self.PT).Bearing()
        #if  norm_as_delta(mid-a) * norm_as_delta(mid-b) <= 0: 
            # one must have been (pos) one (neg) 
        distance = Distance(abs(mid-a)*self.R)
            # difference in lenght then l/r correction based on
            # the sign of the Delta angle
            #print(f"{self.R=}")
            #print(f"{mySegment.Len()=}")
            #print(f"{sign(self.Delta)=}")
            
        offset = Distance((self.R - mySegment.Len())* sign(self.Delta) * -1) # -1 needed to get offset sign correct
            #R and Len() do not have sign, the only sign has to do with Delta
            #we need to use correctly use sign(Delta) to get Lt(-) vs Rt(+) sedt correctly
        return Dist_os_tup(distance, offset)

        #else:
        #    return (None, None)
    '''

    def move_to(self, distance:float, offset=0.0) -> Point:
        '''Return a point on the curve at distance from the PC'''
        if distance > self.Len():
            raise ValueError(
                    f"Distance of {distance} is greater than length of Arc={self.Len()}")
        if distance < 0:
            raise ValueError(
                    f"Distance of {distance} is negative")
        #if Delta>0 -> CCW movement with cc on the left  so min offset is -R
        #if Delta<0 ->  CW movement with cc on the right so max offset is +R
        #check both
        if self.Delta > 0 and offset < -self.R:
            raise ValueError(
                    f"{offset=} is too negative, goes left more than {self.R=}")
        if self.Delta < 0 and offset > +self.R:
            raise ValueError(
                    f"{offset=} is too possitive, goes right more than {self.R=}")

        #just move to cc turn and move back
        if abs(distance) <= self.Len():
            #start a the cc and look a that pc
            #change your bearing by the arc of the movement
            brg_from_cc = cmath.phase(self.CC_to_PC) + sign(self.Delta) * distance / self.R
            #move from the CC exactly one R in this new direction
            _new_point = self.CC + cmath.rect(self.R, brg_from_cc)
            #find new bearing by staring at the the inBearin and adjusting by the arc of movement 
            _new_bearing = self.inBearing + sign(self.Delta) * distance / self.R
            #adjust this new point on the curve by turning right 90deg and going forward or backward
            _new_point += cmath.rect(offset, normalize(_new_bearing - pi/2.0))
            return Point.from_complex(_new_point)
        return "Error"

    def patch(self, **kwargs):
        x = self.CC.X
        y = self.CC.Y
        brg_CC_to_PC = math.degrees(cmath.phase(self.CC_to_PC))
        brg_CC_to_PT = math.degrees(cmath.phase(self.CC_to_PT))
        r = 2 * self.R #actually the axis
        if self.Delta > 0:
            brg_start = brg_CC_to_PC
            brg_end   = brg_CC_to_PT
        else:
            brg_start = brg_CC_to_PT
            brg_end   = brg_CC_to_PC
        return patches.Arc((x, y), r, r, angle=0, theta1=brg_start, theta2=brg_end, **kwargs)


    def patch_all(self, ** kwargs):
       ret = []
       ret.append(self.patch(**kwargs))
       ret.append(Segment(self.CC, self.PC).patch(linestyle=':'))
       ret.append(Segment(self.CC, self.PT).patch(linestyle=':'))
       ret.append(Segment(self.PC, self.PT).patch(linestyle=':'))
       #ret.append(Segment(self.PC, self.PI).patch(linestyle=':'))
       #ret.append(Segment(self.PI, self.PT).patch(linestyle=':'))
       return ret

    def __repr__(self):
        rep =  "Curve:\n"
        rep += f"PC= {self.PC}\n"
        rep += f"CC= {self.CC}\n"
        rep += f"PT= {self.PT}\n"
        rep += f" R= {self.R}\n"
        rep += f" E= {self.E}\n"
        rep += f"Delta= {self.Delta}\n"
        rep += f"Len()= {self.Len()}\n"
        rep += f"Chord= {self.Chord}\n"
        rep += f"{cmath.phase(self.CC_to_PC)=}\n"
        rep += f"{cmath.phase(self.CC_to_PT)=}\n"
        rep += f"{self.CC_to_PC - self.CC_to_PT=}\n"
        #rep += f"-- Curve end --\n"
        # the segment finish up with a -- so last line not needed
        return rep

class Chain:
    '''A list of connected segment and or curves
    '''
    is_Curve = True
    def __init__(self, name, *args, StartSta=None, **kwargs):
        self.Name = name
        self.Routes = []
        self.RoutesSta = []

        if StartSta == None:
            self.StartSta=0.0
        else:
            self.StartSta=StartSta

        for ob in args:
            self.addRoute(ob)

    def addRoute(self, ob):
        if self.validate_add(ob):
            self.Routes.append(ob)
            #start station of each segment

            x = []
            for i in self.Routes:
                x.append(i.Len())
            start = self.StartSta
            self.RoutesSta = list(itertools.accumulate(x, initial=start))
        else:
            raise ValueError(
                "%s This route can't be added here" % self)

    def forward(self, distance:float, R=float('inf'), Delta=float(0)):
        '''Add a new segment to the Chain using the outRay of the
        previously last route that was defined
        '''
        # Need to get inputs for curves to work
        if len(self.Routes) == 0 and (math.isinf(R) or Delta==0):
            raise ValueError(
                "forward() method not allowed for Chain with no routes")

        lastRay = self.Routes[-1].outRay()
        #segment add
        if math.isinf(R) and Delta ==0:
            lastRay = self.Routes[-1].outRay()
            pt2 = lastRay.move_to(distance)
            self.addRoute(Segment(lastRay.Point, pt2))
            return self.Routes[-1].outRay
        #curve add
        if math.isinf(R) and Delta != 0:
            R = distance / Delta
            self.addRoute(Curve.from_PC_bearing_R_Delta(lastRay.Point, lastRay.Bearing, R, Delta))
            return self.Routes[-1].outRay


    def validate_add(self, ob) -> bool:
        if not isinstance(ob, (Segment, Curve)):
            raise ValueError(
                "%s not a Segment or Curve" % self)

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

    def move_to(self, sta:float, offset=0.0) -> Point:
        '''Return a new point along a chain defined by station and offset'''
        #print(f"{sta}")
        #print(f"{self.RoutesSta}")
        #Lens =[]
        #for r in self.Routes:
        #    Lens.append(r.Len())
        #print(f"{Lens=}")

        if sta < self.RoutesSta[0]:
            raise ValueError("Sta should be greater than the chain start station")
        if sta > self.RoutesSta[-1]:  #RoutesSta has end of chain station
            raise ValueError("Sta should be less than the chain ending station")

        for i in range(len(self.RoutesSta[:-1])):
            if self.RoutesSta[i] <= sta and sta <= self.RoutesSta[i+1]:
                distance = sta - self.RoutesSta[i]
                return self.Routes[i].move_to(distance, offset)

    def copy_parallel(self, offset:float, StartSta=None):
        """
        Return a new chain at the given offset
        """
        # check that all curves have a valid offset
        _new_routes = ['new_name']
        for r in self.Routes:
            if isinstance(r, Curve):                # test the curves
                r.is_offset_valid(offset) # will stop with a valueError
            _new_routes.append(r.copy_parallel(offset))

        _new_chain = Chain(*_new_routes, StartSta=StartSta)
        return _new_chain


    def outRay(self) -> Ray:
        return self.Routes[-1].outRay()

    def point_list(self) -> list:
        ret = []
        for r in self.Routes:
            if isinstance(r, Curve):
                ret.extend([r.PC, r.CC, r.PT])
            else:
                ret.extend([r.Pt1, r.Pt2])
        return ret


    def patch_list(self) -> list:
        ret = []
        for r in self.Routes:
            if isinstance(r, Curve):
                for curve_patch in r.patch_all():
                    ret.append(curve_patch)
            else:
                ret.append(r.patch())
        return ret

    def inverse(self, point:Point, decimals=2) -> tuple:
        '''For this Chain(self) return the station and offset of the input point
        there might be multiple valid station offset pairs, return 
        the one with the shortest absolue offset distance from the chain.
        decimals=2 is the default decimal precision. 
        '''
        # might be better to return a named tuple or of list of named tubles.
        # ('station':???, 'offset':???)

        #Fix this and the inverse for curves!
        subList = []
        for r, s_sta in zip(self.Routes, self.RoutesSta[:-1]):
            dist, offset = r.distance_offset(point)
            if dist is not None and offset is not None: # skip Nones
                sta    = round(dist + s_sta, decimals)
                offset = round(offset,       decimals)
                subList.append((sta, offset))
        #print(subList)
            #import pdb; pdb.set_trace()

        subList = sorted(subList, key = lambda x:abs(x[1]))  #sort by the ascending abs of the offset
        #print(subList)
        if len(subList) == 0:
            return Dist_os_tup(None, None)
        return Dist_os_tup(subList[0][0], subList[0][1])

    def __repr__(self) -> str:
        _s = ""
        _s += "Chain: " + self.Name +"\n"
        for sta, rte in zip(self.RoutesSta[:-1], self.Routes):
            _s += f"<< Station: {sta} >>" + "\n"
            _s += rte.__str__() + "\n"

        _s += f"<< Station: {self.RoutesSta[-1]} >>" + "\n"
        _s += "End Chain: " + self.Name
        return _s


def main():
    import random
    def randomPoints(pt, noise=5, number=10) -> list:
      point_list = []
      for i in range(number):
          x=pt.X + random.uniform(-noise, +noise)
          y=pt.Y + random.uniform(-noise, +noise)
          point_list.append(Point(x, y))
      return point_list
    ru = lambda a, b: random.uniform(a, b) # b/w a and b
    ru_10 = lambda : random.uniform(-10, 10)
    ru_100 = lambda : random.uniform(90, 110)

    R= 15
    points = []
    x=5 # ru_10()
    y= 2 #ru_10()
    points.append(Point(x, y))
    #points.append(Point(ru_100(), points[-1].Y + ru_100()))
    points.append(Point(102, points[-1].Y + 105))

    segments = []
    segments.append(Segment(points[-2], points[-1]))
    seg_length = segments[-1].Len()
    myChain = Chain("default", segments[-1], StartSta=1000)
    #points.append(myChain.outRay().move_to(20, 0))
    #rand_R = ru(seg_length/10, seg_length/2)
    rand_R = seg_length/5
    #rand_delta_n = ru(-90, -90-180)
    rand_delta_n = -70

    myChain.addRoute(Curve.from_PC_bearing_R_Delta( 
                            pc=myChain.Routes[-1].outRay().Point, #pc=points[-1], 
                            brg=myChain.Routes[-1].outRay().Bearing,
                            R=rand_R,
                            delta=Angle(rand_delta_n, unit='deg')
                            )
                     )
    #points.append(myChain.outRay().move_to(20, 0))
    myChain.addRoute(Curve.from_PC_bearing_R_Delta( 
                            pc=myChain.Routes[-1].outRay().Point, #pc=points[-1], 
                            brg=myChain.Routes[-1].outRay().Bearing,
                            R=rand_R,
                            delta=Angle(-rand_delta_n, unit='deg')
                            )
                     )
    myChain.forward(distance=100)

    domain = {
    'low_x': min([item.X for item in myChain.point_list() ]), 
    'high_x': max([item.X for item in myChain.point_list() ]),
    'low_y': min([item.Y for item in myChain.point_list() ]),
    'high_y': max([item.Y for item in myChain.point_list() ]) }
    print(domain)

    rand_points = []
    bad_points = []
    rand_sta_off = []

    print(myChain.RoutesSta)

    # for i in range(3):
    d_os = []
    y = 10 
    # for i in  [Point(0, 0), Point(95, y), Point(100, y), Point(110, y), Point(120, y), Point(130, y), Point(140, y), Point(150, y), ]:
    for loop in range(30):
        i = Point(random.uniform(domain['low_x']-20, 2*domain['high_x']), random.uniform(domain['low_y']-20, 2*domain['high_y']))
        inver = myChain.inverse(i, )
        if inver[0] is None:
            bad_points.append(i)
        else:
            rand_points.append(i)
            segments.append(Segment(myChain.move_to(inver[0], 0), i))  #note there is no offset

    fig, ax = plt.subplots()
    #not plot vs scatter
    ax.scatter(as_XY(points)[0],  as_XY(points)[1],           marker='+', color='r')
    ax.scatter(as_XY(rand_points)[0],  as_XY(rand_points)[1], marker='.', color='g')
    ax.plot(as_XY(bad_points)[0],  as_XY(bad_points)[1],      marker='o', color='r', markersize=2, linestyle="None")
    ax.scatter(as_XY(d_os)[0],  as_XY(d_os)[1],               marker='+', color='g')


    myChain2 = myChain.copy_parallel(28)
    for rp in rand_points:
        pass
        #print(myChain.inverse(rp))
        for c in [1, 2]:
            if myChain.Routes[c].normal_at_point(rp):
                n1 = myChain.Routes[c].normal_at_point(rp)
                t1 = myChain.Routes[c].tangent_at_point(rp)
                print(c, n1, t1, n1-t1)

    # input()
    for s in segments:
        ax.add_patch(s.patch(linestyle=':', linewidth=0.2))
    for p in myChain.patch_list():
        #print(p)
        ax.add_patch(p)
    for p in myChain2.patch_list():
        #print(p)
        ax.add_patch(p)
    plt.axis('scaled')
    plt.show()
    IPython.embed()


if __name__ == "__main__":
    main()

