import Geo
from Geo import  Point, Ray, Curve
from math import pi, tan , atan2, atan, cos,sqrt
import matplotlib.pyplot as plt
from matplotlib import patches
import random


def r_c_intersections( r1:Ray, c1:Curve )->int:
    dist_to_CC = r1.offset_to(c1.CC)
    if dist_to.CC < c1.R:
        return  2
    if dist_to.CC == c1.R:
        return  1
    return 0
    
def c_c_intersections( c1:Curve, c2:Curve )->int:
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


def c_c_2 ( c1:Curve, c2:Curve)->( ( Point, Point) ):
    # assume c1.R > c2.R  then check
    b_C = c1
    l_c = c2
    if b_C.R < l_c.R:
        b_C = c2
        l_c = c1
    d = abs( l_c.CC - b_C.CC)
    r = l_c.R
    R = b_C.R
    x = (d**2 - r**2 + R**2 ) / (2* d)

    a =sqrt( 4*d**2*R**2-(d**2-r**2+R**2)**2) / d
    y = a/2
    for i in  ['d', 'r', 'R', 'x', 'a', 'y'  ]:
        print( f'{i}=', eval(i) ) 
    return [ Point(x, y), Point(x, -y) ]


def main():
    
    myC1 = Curve( Point( 8,0), Point( 0,0), 2.0*pi) 
    myC2 = Curve( Point(19,0), Point(13,0), 2.0*pi) 
    if c_c_intersections(myC1, myC2) > 0:
        pts =  c_c_2( myC1, myC2 )

    fig, ax = plt.subplots()
    ax.scatter( as_XY(pts)[0],  as_XY(pts)[1] )

    ax.add_patch( myC1.patch(linestyle='-') )
    ax.add_patch( myC2.patch(linestyle='-') )
    plt.axis('scaled')
    plt.show()

if __name__ == "__main__":
    main()
