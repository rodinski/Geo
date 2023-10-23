import geo_
from geo_ import  Point, Ray, Bearing, Curve, Angle, Segment, xy
from math import pi, tan , atan2, atan, cos,sqrt
import matplotlib.pyplot as plt
from matplotlib import patches
import random



def r_c_intersection_count( r1:Ray, c1:Curve )->int:

    if not isinstance(r1, Ray):
        raise TypeError( f"Provide a Ray object")
    if not isinstance(c1, Curve):
        raise TypeError( f"Provide a Curve object")

    dist_to_CC = abs(r1.get_offset_to(c1.CC))
    if dist_to_CC < c1.R:
        return  2
    if dist_to_CC == c1.R:
        return  1
    return 0
    
def c_c_intersection_count( c1:Curve, c2:Curve )->int:

    if not isinstance(c1, Curve):
        raise TypeError( f"Provide a Curve object")
    if not isinstance(c2, Curve):
        raise TypeError( f"Provide a Curve object")

    cc_to_cc = abs( c1.CC - c2.CC ) 
    if  cc_to_cc > c1.R + c2.R :
        return 0
    if cc_to_cc + c1.R  < c2.R:
        return  0
    if cc_to_cc + c2.R  < c1.R:
        return  0

    if cc_to_cc + c1.R  == c2.R:
        return  1
    if cc_to_cc + c2.R  == c1.R:
        return  1
    return 2

def as_XY( list_of_points ):
    x = list()
    y = list()
    for i in list_of_points:
        x.append(i.X)
        y.append(i.Y)
    return (x,y)

def f(myStr):
    retList = []
    for i in myStr:
        retList.append( 'f"{' + i +'=}"' )
        #print(i)
    return retList


def curve_curve_intersect( c1:Curve, c2:Curve)->((Point, Point)):
    # assume c1.R > c2.R  then check
    b_C= c1
    l_c= c2
    if b_C.R < l_c.R:
        b_C= c2
        l_c= c1
    d= abs( l_c.CC - b_C.CC)
    r= l_c.R
    R= b_C.R
    x= (d**2 - r**2 + R**2 ) / (2* d)

    a=sqrt( 4*d**2*R**2-(d**2-r**2+R**2)**2) / d
    y= a/2
    for i in  ['d', 'r', 'R', 'x', 'a', 'y'  ]:
        print( f'{i}=', eval(i) ) 
    return [ Point(x, y), Point(x, -y) ]

def ray_curve_intersect( r1:Ray, c1:Curve)->tuple:
    if not r_c_intersection_count(r1,c1):
        return (None)
    distance_to_cc = r1.get_distance_to(c1.CC)
    offset_to_cc = r1.get_offset_to(c1.CC)
    a = r1.get_offset_to(c1.CC)
    c = c1.R
    try:
        b = sqrt(c**2 - a**2)
    except:
        raise ValueError(f"Can not calucutlate sqrt(c**2-a**2)")
    greater_point= r1.set_point(distance_to_cc - b)
    lesser_point= r1.set_point(distance_to_cc + b)
    if not b:
        return(greater_point)
     
    return tuple(sorted([lesser_point, greater_point], 
        key=lambda x: r1.get_distance_to(x)))
     

def main():
    pts = [] 
    myC1 = Curve( Point( 8,0), Point( 0,0), Angle(2.0*pi))
    myC2 = Curve( Point(19,0), Point(13,0), Angle(2.0*pi))

    rayPt=Point(-3,-9)
    pts.append(rayPt)
    
    myRay1 = Ray( rayPt,  Bearing(pi/12 ))
    if c_c_intersection_count(myC1, myC2) > 0:
        pts.extend(curve_curve_intersect( myC1, myC2 ))

    near, far= ray_curve_intersect(myRay1, myC1)
    pts.extend((near,far)) 
    near, far= ray_curve_intersect(myRay1, myC2)
    pts.extend((near,far)) 


    fig, ax = plt.subplots()
    ax.scatter( *xy(pts) )

    ax.add_patch( myC1.patch(linestyle='-'))
    ax.add_patch( myC2.patch(linestyle='-'))
    ax.add_patch( Segment(pts[-4],pts[-1]).patch())
    plt.axis('scaled')
    plt.show()

if __name__ == "__main__":
    main()
