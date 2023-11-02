""" Module for various geometry entities to be used in Civil 
Engineering makes use of python's complex numbers as x, y pairs.

Every complex number has a phase = cmath.phase( complex) it is
the angle from the real axis it only goes from -pi to pi.  This
is good for finding one and only one direction.

It has proplems when adds or subracts will cross over the +/- pi
break.  

Let "phase" method determin a global _direction_

Let  Angle(floats) determin changes in a direction.  Angles have to
have  a direction +/- (ccw/cw) 

Additions of Bearing and Angle are possible but they result in
a float that must be converted back to a Bearing.
(Future work would be to overload the Bearing's +, - methods)

"""

import math
from math import e, pi, cos, sin, tan, isclose
import cmath
from cmath import phase, rect, polar
from dataclasses import dataclass
from collections import namedtuple
import itertools
import IPython
import matplotlib.pyplot as plt
from matplotlib import patches

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

# normalize_as_phase = -pi pi
#def float_to_phase(inTheta: float):
#    return phase(rect(1.0,inTheta)) #(-π, π]

# turns should go from (-2π, 2π] AND KEEP there sign
def normalize_angle(theta:float):
    # keep chiping away as needed
    while theta < 2*-pi or theta > 2*pi:
        theta = theta - sign(theta)*2*pi
    return theta

def normalize_bearing(theta:float):
    # keep chiping away as needed
    while theta < -pi or theta > pi:
        theta = theta - sign(theta)*2*pi
    if isclose( abs(theta), pi):
        theta = pi
    return theta


def conjugate_angle(inAngle):
    if abs(inAngle) > 2*pi:
        ValueError: print("Angle out of range in congugate_angle function")
    if inAngle == 0:
        return 2*pi
    # change the sign b/c we would going other way around the circle
    return -1.0*sign(inAngle)*(2*pi-abs(inAngle))


def xy(points_list: list) -> tuple:
    """Return a list of .x and a list of .y properites
    for the given points_list.
    Useful for use in matplotlib plts and scatter.
    Such as  ax.scatter(*xy(pts), marker='+', color='g')
    """
    x = [arg.X for arg in points_list]
    y = [arg.Y for arg in points_list]
    return(x, y)


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
        """ returns the angle of a complex number from X=0(East). This is not
        same as the Bearing! 
        -pi < phase(complex) < pi
        """
        return Bearing(cmath.phase(self))  #technically this is a bearing

    def equal(self, ob):
        ''' Are two Points equal?'''
        if isinstance(ob, Point) and \
                      math.isclose(self.X, ob.X) and \
                      math.isclose(self.Y, ob.Y):
           return True
        return False

    def __repr__(self):
        return f"Point({self.val.real:.4f}, {self.val.imag:0.4f})"

    def __str__(self):
        return f"Point x={self.val.real:.4f}, y={self.val.imag:0.4f}"
 
 
class Angle(float):
    '''Angle will be stored as a float with a nice __str__().
    All math funtions will return a float without the __str__()
    if unit='deg' is used the input should be in degrees and a
    conversion will happen so the internal storage will be radians

    Angles are relative and have no global directoin buth they do
    have a significant sign (ccw vs cw).

    normalize_angle return a float in range (-2π, 2π)
    '''
    is_float = True

    def __new__(cls, _numb, unit='rad'):
        if unit.lower() == 'rad':
            _numb = float(_numb)
        if unit.lower() == 'deg':
            _numb = float(math.radians(_numb))
        _numb = normalize_angle(_numb)
        return float.__new__(cls, _numb)

    def __init__(self, _numb: float, unit='rad'):
        self.unit = unit
        self.sign = sign(self) 

    def __str__(self):
        return '%g rad  %g deg' % (self, math.degrees(self))

    def __repr__(self):
        return f"Angle({self:,f})"

class Distance(float):
    is_float = True

    def __init__(self, flt: float):
        Decimals = int(4)
        self.val = round(flt, Decimals)

    def __repr__(self):
        return f"Distance({self.val:,f})"

class Bearing(float):
    """ Angle from real axis,(Y=0; East)) this is the same as 
    cmath.phase(complex_number) we will want to add and subtract with 
    changes in angles (delta_phase) that will be in float"""

    def __new__(cls, _ang, unit='rad'):
        if unit.lower() == 'rad':
            _ang = float(_ang)
        if unit.lower() == 'deg':
            _ang = float(math.radians(_ang))
        _ang = normalize_bearing(_ang)
        return float.__new__(cls, _ang)


    def __init__(self, real: float, unit='rad'):
        self.unit = unit


    def opposite(self):
        """Pi radians from the Bearing good for backsights in Civil
        Engineering.
        """
        return Bearing(cmath.phase(cmath.rect(1.00, self + pi)))

    def lt_90(self):
        """Return the bearing that is 90 CCW"""
        # cmath.rect( r, theta)  returns a complex
        return Bearing(cmath.phase(cmath.rect(1.00,  self + pi/2.0)))

    def rt_90(self):
        """Return the bearing that is 90 CW """
        # cmath.rect( r, theta)  returns a complex
        return Bearing(cmath.phase(cmath.rect(1.00,  self - pi/2.0)))

    def __str__(self):
        return '%g rad  %g deg' % (self, math.degrees(self))

    def __repr__(self):
        return f"Bearing({self:f})"

@dataclass
class Ray:
    """Class for when we have a know Point and a Bearing
    """
    is_Ray = True

    def __init__(self, pt: Point, brg: Bearing):

        if not isinstance(pt, Point):
            raise TypeError(
                    f"Provide a Point object")
        # print(f"{isinstance(brg, Bearing)=}")
        if not isinstance(brg, Bearing):
            raise TypeError(
                    f"Provide a Bearing (not an Angle or float)")

        self.Point = pt
        self.bearing = brg

    def set_point(self, distance: float, offset=0.0) -> Point:
        ''' Return a point at the given distance and offset from the Ray origin 
        '''

        _bearing = self.bearing
        _bearing_rt_90 = self.bearing.rt_90()
        _new_point = self.Point \
                     + cmath.rect(distance, self.bearing) \
                     + cmath.rect(offset, _bearing_rt_90)
                     # start + move along arrow + ,move to offset 
        return Point.from_complex(_new_point)
    

    def set_multipoint(self, start: float = 0.,
                    spacing: float = 0., 
                    n_steps: int = 0,
                   point_at_start: bool=False) -> tuple:
        """Return a tuple of points input start distance from Pt1 along 
        the bearing, return of a point here is controlled by 
        point_at_start=T/F (default=False) spacing is a distance to the 
        next point n_steps integer number of points to produce
        """
        ret = ()
        if point_at_start:
            ret.append(self.set_point(start))
        for i in range(n_steps):
            distance = start + (i+1)*spacing
            ret.append(self.set_point(distance))

    def equal(self, ob):
        ''' Are two Rays equal?'''
        if isinstance(ob, Ray) and \
           self.Point == ob.Point and \
           self.bearing == ob.bearing:
           return True
        return False

    def get_angle_to_bearing(self, destBearing: float) -> float:
        """ Return the turn angles needed to get at some destination bearing.
        Output limited to +/- pi. 
        Rotate both rays by the negative of the self.bearing, then the phase 
        of the destination is the angle_to_bearing.
        rotation of destination = start.bearing - end.bearing and
        rotation of self bearing = self.bearing - self.bearing = 0
        """
        start = destBearing
        end = self.bearing
        return start - end  # -pi < x < pi

    def get_distance_to(self, inPoint: Point) -> float:
        ''' What is the distance along the Ray till the
        point in quesition is at a right angle
        to the Ray '''
        # make a Segment to the new point
        # segments stored as 2 points Bearing is a method()

        # complex that takes you from ray.Point to inPoint
        _movement_to_point = inPoint - self.Point 
        theta = cmath.phase(_movement_to_point) - self.bearing #angle = diff in phase
        ret = math.cos(theta) * abs(_movement_to_point)
        return ret
    
    def get_offset_to(self, inPoint: Point) -> float:
        """Return the perpendicular distance from ray to inPoint."""

        # complex that takes you from ray.Point to inPoint
        _movement_to_point = inPoint - self.Point 
        theta = cmath.phase(_movement_to_point) - self.bearing #angle = diff in phase
        ret = math.sin(theta) * abs(_movement_to_point)
        return ret
    
        #mySeg = Segment(self.Point, inPoint)
        #theta = self.get_angle_to(mySeg.bearing) * -1 #(+) offsets have negative theta
        ##-1 above needed to get sign of offset correct.
        #offset = math.sin(theta) * mySeg.length
        #return offset

    def get_distance_offset(self, inPoint: Point) -> tuple:
        _distance = self.get_distance_to(inPoint)
        _offset = self.get_offset_to(inPoint)
        return Dist_os_tup(_distance, _offset)

    def copy_parallel(self, offset):
        #new offset the point and use same Bearing
        return Ray( self.set_point(0, offset), self.bearing)

    def patch(self, scale=1, width=1,  **kwargs):
        """Use matplotlib.patches for consistency, method patches.Arc
        works best for class Curve so we will also use method patches.Polygon
        for the Segemnts. """

        return patches.Arrow (self.Point.X, self.Point.Y, \
                             scale *cos(self.bearing), scale *sin(self.bearing), \
                             width,  **kwargs)

    def __repr__(self):
        return f"Ray(self.Point.__repr__(), self.bearing.__repr__())"

    def __str__(self):
        _s = "\n"
        _s +=  "Ray:\n"
        _s += f"       .{self.Point}\n"
        _s += f"       .{self.bearing}\n"
        _s += f"-- Ray End --"
        return _s

#class Line:
#    is_Line = True
#
#    def __init__(self, pt:Point, bearing:Bearing):
#        self.Point = pt
#        self.bearing = Bearing(bearing)
#    #need to be able to accept an offset
#
#    def __repr__(self):
#        _s = ""
#        _s +=  "Line:\n"
#        _s += f"       Point= {self.Point}"
#        _s += f"       Bearing= {self.bearing}\n"
#        _s += f"--"
#        return _s


@dataclass
#class Segment(Ray):
class Segment:

    is_Segment = True
    def __init__(self, Pt1:Point, Pt2:Point):

        if not isinstance(Pt1, Point):
            raise TypeError(
                    f"Provide a Point object")
        if not isinstance(Pt2, Point):
            raise TypeError(
                    f"Provide a Point object")
        self.Pt1 = Pt1
        self.Pt2 = Pt2
        self.length = Distance(abs(self.Pt2 - self.Pt1))
        self.bearing = Bearing(cmath.phase(self.Pt2 - self.Pt1))

    #def Len(self) -> float:
    #    '''Length of Segment '''
    #    return self.length

    # remove eventually and just use attrabute .bearing
    #def Bearing(self) -> Bearing:
    #    """Bearing of the segment """
    #    return  Bearing(cmath.phase(self.Pt2 - self.Pt1))

    def inRay(self) -> Ray:
        """Point and Bearing for the first Point """
        return Ray(self.Pt1, self.bearing)

    def outRay(self) -> Ray:
        """Point and Bearing for the second Point """
        return Ray(self.Pt2, self.bearing)
    
    # needs redone b/c either CW or CCW will get you there!
    # finding the angle to a point is ok but delta_angle implize a direction
    #def delta_angle_to_some_bearing(self, inBearing) -> float:
    #    """ Reurn the turning angle needed from the inRay to some 
    #    given bearing"""
    #    return norm_as_delta(self.bearing - inBearing)

    def get_distance_offset(self, obPoint: Point) -> (float, float):
        if not (0.0 <= self.inRay().get_distance_to(obPoint) <= self.length):
           return Dist_os_tup(None, None)
        else:
           distance = self.inRay().get_distance_to(obPoint)
           offset = self.inRay().get_offset_to(obPoint)
           return Dist_os_tup(distance, offset)

    # remove eventually
    def distance_offset(self, obPoint: Point) -> (float, float):
        return self.get_distance_offset(obPoint)


    def set_point(self, distance: float, offset=0.0) -> Point:
        if distance > self.length:
            raise ValueError(
                    f"Distance of {distance} is greater than length of Segment={self.length}")
        if distance < 0.0:
            raise ValueError(
                    f"Distance of {distance} is negative")
        _new_point = self.Pt1 + cmath.rect(distance, self.bearing) \
                     + cmath.rect(offset, self.bearing - pi/2.0)
        return Point.from_complex(_new_point)

    # remove this
    #def move_to(self, distance: float, offset=0.0) -> Point:
    #    return self.set_point(distance, offset)



    def copy_parallel(self, offset):
        """Return copy of the Segment, offset has no limits
        """
        _pt1 = self.set_point(0, offset)
        _pt2 = self.set_point(self.length, offset)
        return Segment(_pt1, _pt2)


    def __repr__(self):
        s = ""
        s += f"Segment({self.Pt1.__repr__()}, {self.Pt2.__repr__()})"
        return s

    def patch(self, **kwargs):
            """Use matplotlib.patches for consistency, method patches.Arc
            works best for class Curve so we will also use method patches.Polygon
            for the Segemnts. """
            return patches.Polygon(((self.Pt1.X, self.Pt1.Y), \
                                      (self.Pt2.X, self.Pt2.Y), \
                                    ), fill=False, closed=False, **kwargs)
@dataclass
class Curve:
    is_Curve = True

    """
    Curve is defined by a (PC, CC, Delta), PC != CC and
    -2pi< Delta < 2pi and not zero
    """

    def __init__(self, Pt1:Point, Pt2:Point,  Delta:Angle , *args, **kwargs):

        if not isinstance(Pt1, Point):
            raise TypeError(
                    f"Provide a Point object")

        if abs(Pt1 - Pt2) < 0.1:
            raise ValueError(
                "%s.__new__ radius is too small" % self)
        if abs(Delta == 0):
            raise ValueError(
                "%s.__new__ Delta can not be zero" % self)


        self.PC = Pt1   #Point
        self.CC = Pt2   #Point
        self.PC_to_CC = self.CC - self.PC   #Complex
        self.CC_to_PC = self.PC - self.CC   #Complex
        self.PC_to_CC_asSeg = Segment(self.CC, self.PC) #Segment

        self.R = Distance(abs(Pt1 - Pt2))
        self.Delta = Delta
        self.sign_delta = sign(Delta)  #ccw=+   cw=-
        self.length = abs(self.Delta * self.R)
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
        self.inBearing =  Bearing(self.Chord.bearing - self.Delta /2.0)
        self.outBearing = Bearing(self.inBearing + self.Delta)
    #need to be able to accept an offset

    #def Len(self) -> float:
    #    return Distance(self.R * abs(self.Delta))


    @classmethod
    def from_PC_bearing_R_Delta(cls, pc:Point, brg:float, R:float, delta:float):
        #pc given
        #find cc
        cc = Point.from_complex(pc + cmath.rect(R, brg + sign(delta)*pi/2))
        return cls(pc, cc, delta)

    @classmethod
    def from_Ray_T_Delta(cls, _inRay:Ray, _T:float, delta:float):
        _pc = _inRay.Point
        _start_bearing= _inRay.bearing
        _R = math.cos(delta/2.0) * _T
        return cls.from_PC_bearing_R_Delta(_pc, _start_bearing, _R, delta)

    def inRay(self) -> Ray:
        return Ray(self.PC, self.inBearing)  #should already be Point and Bearing

    def outRay(self) -> Ray:
        return Ray(self.PT, self.outBearing)  #should already be Point and Bearing

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
        Howerver all curves are one part of the whold cirle.  Alwasy one wedge
        sweep crosses the -/+Pi line and the other side does not.  Try to work with the
        side that does not.
        """

        # logic it is easy to see if point is within the wedge
        pc_brg = cmath.phase(self.CC_to_PC)
        defined_curve = -pi <= (pc_brg + self.Delta)  <= +pi #set a bool
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

    def get_distance_offset(self, obPoint:Point) -> (float, float):
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
        offset = Distance((self.R - mySegment.length)* sign(self.Delta) * -1) # -1 needed to get offset sign correct
            # R and length do not have sign, the only sign has to do with Delta
            # we need to use correctly use sign(Delta) to get Lt(-) vs Rt(+) sedt correctly
        return Dist_os_tup(dist, offset)

    # eventually remove this
    def distance_offset(self, obPoint:Point) -> (float, float):
        return self.get_distance_offset(obPoint)

    def normal_at_point(self, obPoint: Point) -> Bearing:
        if not self.has_dist_os(obPoint):
            return None
        #phase of a movement
        return Bearing(cmath.phase(obPoint - self.CC)) 

    def tangent_at_point(self, obPoint:Point) -> Bearing:
        if not self.has_dist_os(obPoint):
            return None
        norm = self.normal_at_point(obPoint)
        if sign(self.Delta > 0):
            turn = pi/2
        else:
            turn = -pi/2
        return Bearing(norm + turn)  # good bearing+angle=float

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


    def set_point(self, distance:float, offset=0.0) -> Point:
        '''Return a point on the curve at distance from the PC'''
        if distance > self.length:
            raise ValueError(
                    f"Distance of {distance} is greater than length of Arc={self.length}")
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
        if abs(distance) <= self.length:
            #start a the cc and look a that pc
            #change your bearing by the arc of the movement
            #distance has no negative so sign(self.Delta) is needed 
            brg_from_cc = cmath.phase(self.CC_to_PC) + sign(self.Delta) * distance / self.R
            #move from the CC exactly one R in this new direction
            _new_point = self.CC + cmath.rect(self.R, brg_from_cc)
            #find new bearing by staring at the the inBearin and adjusting by the arc of movement 
            _new_bearing = self.inBearing + sign(self.Delta) * distance / self.R
            #adjust this new point on the curve by turning right 90deg and going forward or backward
            _new_point += cmath.rect(offset, _new_bearing - pi/2.0) 
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

    def __str__(self):
        rep =  "Curve:\n"
        rep += f"  PC= {self.PC}\n"
        rep += f"  CC= {self.CC}\n"
        rep += f"  PT= {self.PT}\n"
        rep += f"   R= {self.R}\n"
        rep += f"   E= {self.E}\n"
        rep += f"Delta= {self.Delta}\n"
        rep += f"length= {self.length}\n"
        rep += f"Chord= {self.Chord}\n"
        rep += f"{cmath.phase(self.CC_to_PC)=}\n"
        rep += f"{cmath.phase(self.CC_to_PT)=}\n"
        rep += f"{self.CC_to_PC - self.CC_to_PT=}"
        #rep += f"-- Curve end --\n"
        # the segment finish up with a -- so last line not needed
        return rep

    def __repr__(self):
        rep =  "Curve:(\n"
        rep += f"  {self.PC.__repr__()}, \n"
        rep += f"  {self.CC.__repr__()}, \n"
        rep += f"  {self.Delta.__repr__()})"
        return rep
class Chain:
    """A list of connected segments and/or curves
    """
    is_Curve = True
    def __init__(self,  *args, name="no_name", start_station=None, **kwargs):
        self.Name = name
        self.Routes = []
        self.RoutesSta = []

        if start_station == None:
            self.start_station=0.0
        else:
            self.start_station=float(start_station)

        for ob in args:
            self.addRoute(ob)

    def addRoute(self, ob):
        if self.validate_add(ob):
            self.Routes.append(ob)
            #start station of each segment

            x = []
            for i in self.Routes:
                x.append(i.length)
            start = self.start_station
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
            pt2 = lastRay.set_point(distance)
            self.addRoute(Segment(lastRay.Point, pt2))
            return self.Routes[-1].outRay
        #curve add
        if math.isinf(R) and Delta != 0:
            R = distance / Delta
            self.addRoute(Curve.from_PC_bearing_R_Delta(lastRay.Point, lastRay.bearing, R, Delta))
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

    def set_point(self, sta:float, offset=0.0) -> Point:
        '''Return a new point along a chain defined by station and offset'''
        #print(f"{sta}")
        #print(f"{self.RoutesSta}")
        #Lens =[]
        #for r in self.Routes:
        #    Lens.append(r.length)
        #print(f"{Lens=}")

        if sta < self.RoutesSta[0]:
            raise ValueError("Sta should be greater than the chain start station")
        if sta > self.RoutesSta[-1]:  #RoutesSta has end of chain station
            raise ValueError("Sta should be less than the chain ending station")

        for i in range(len(self.RoutesSta[:-1])):
            if self.RoutesSta[i] <= sta and sta <= self.RoutesSta[i+1]:
                distance = sta - self.RoutesSta[i]
                return self.Routes[i].set_point(distance, offset)

    def copy_parallel(self, offset:float, start_station=None, name='parallel_copy'):
        """
        Return a new chain at the given offset
        """
        # check that all curves have a valid offset
        _new_routes = [] #pre-populate with some name change this with a new __int__ function putting 'name' as only keyward arg
        for r in self.Routes:
            if isinstance(r, Curve):                # test the curves
                r.is_offset_valid(offset) # will stop with a valueError
            _new_routes.append(r.copy_parallel(offset))

        #import pdb; pdb.set_trace()
        _new_chain = Chain( *_new_routes, start_station=start_station, name=name)
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
            dist, offset = r.get_distance_offset(point)
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


    def __str__(self) -> str:
        _s = ""
        _s += f"Chain:\n   {self.Name=}\n"
        for sta, rte in zip(self.RoutesSta[:-1], self.Routes):
            _s += f"<< Station: {sta} >>" + "\n"
            _s += rte.__str__() + "\n"

        _s += f"<< Station: {self.RoutesSta[-1]} >>" + "\n"
        _s += f"End Chain:  name={self.Name}"
        return _s

    def __repr__(self) -> str:  #needs to have list(Segment( ) and Curve( ))
        _s = ""
        _s += "Chain(\n" 
        for rte in  self.Routes:
            _s += f"{rte.__repr__()}, \n"
        _s += f"  name='{self.Name}', \n"
        _s += f"  start_station={self.start_station}\n"
        _s += f")"
        return _s


def main():
    import random

    R= 15
    points = []
    x=5 
    y= 2 
    points.append(Point(x, y))
    points.append(Point(102, points[-1].Y + 105))

    segments = []
    segments.append(Segment(points[-2], points[-1]))
    seg_length = segments[-1].length
    myChain = Chain(segments[-1], name="default",  StartSta=1000)
    rand_R = seg_length/5
    rand_delta_n = -70

    myChain.addRoute(Curve.from_PC_bearing_R_Delta( 
                            pc=myChain.Routes[-1].outRay().Point, #pc=points[-1], 
                            brg=myChain.Routes[-1].outRay().bearing,
                            R=rand_R,
                            delta=Angle(rand_delta_n, unit='deg')
                            )
                     )
    myChain.addRoute(Curve.from_PC_bearing_R_Delta( 
                            pc=myChain.Routes[-1].outRay().Point, #pc=points[-1], 
                            brg=myChain.Routes[-1].outRay().bearing,
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

    d_os = []
    y = 10 
    for loop in range(300):
        i = Point(random.uniform(domain['low_x']-20, 2*domain['high_x']), random.uniform(domain['low_y']-20, 2*domain['high_y']))
        inver = myChain.inverse(i, )
        if inver[0] is None:
            bad_points.append(i)
        else:
            rand_points.append(i)
            segments.append(Segment(myChain.set_point(inver[0], 0), i))  #note there is no offset

    fig, ax = plt.subplots()
    #not plot vs scatter
    ax.scatter(*xy(points),      marker='+', color='r')
    ax.scatter(*xy(rand_points), marker='.', color='g')
    ax.plot(*xy(bad_points),     marker='o', color='r', markersize=2, linestyle="None")
    ax.scatter(*xy(d_os),        marker='+', color='g')


    myChain2 = myChain.copy_parallel(28)
    for rp in rand_points:
        pass
        #print(myChain.inverse(rp))
        for c in [1, 2]:
            if myChain.Routes[c].normal_at_point(rp):
                n1 = myChain.Routes[c].normal_at_point(rp)
                t1 = myChain.Routes[c].tangent_at_point(rp)
                # print(c, n1, t1, n1-t1)

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
    print("\n", myChain2, "\n")
    print( myChain2.__repr__() )
    #IPython.embed()

if __name__ == "__main__":
    main()

