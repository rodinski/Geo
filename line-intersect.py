import Geo
from Geo import  Point, Ray
from math import pi, tan , atan2, atan, cos
import matplotlib.pyplot as plt
from matplotlib import patches
import random

def bearing_as_slope( brg1 )->float:
    ''' good for Rays that are to be treated as lines
    always return change in y as x INCREASES'''
    #deal with a  +/- pi only 
    norm_brg = Geo.norm_as_delta(brg1)
    if  abs(norm_brg)  >  pi/2 :
        # as pi and renorm to a delta from 0
        norm_brg  = Geo.norm_as_delta( norm_brg + pi) 
    if norm_brg in { -pi/ 2 , pi/2}:
        return float( 'inf' )   
    return tan( norm_brg ) 

def intersect_lines( r1:Ray, r2:Ray)-> Point:
    # same point?
    if r1.Point == r2.Point:
        return r1.Point
    # same bearing?
    m1 = bearing_as_slope(r1.Bearing)
    m2 = bearing_as_slope(r2.Bearing)

    if m1 == m2:
        return None
    # set x
    if m1 == float('inf'):
        x = r1.Point.X
    elif m2 == float('inf'):
        x = r2.Point.X
    else:
        x = (-r2.Point.X*m2 -r1.Point.Y  +r1.Point.X*m1 +r2.Point.Y ) / (m1-m2) 
    #set y
    if m1 == 0.0:
        y = r1.Point.Y
    elif m2 == 0.0:
        y = r2.Point.Y
    else:
        if m1 != float('inf'):
            y = m1*(x - r1.Point.X) +r1.Point.Y
        else: 
            y = m2*(x - r2.Point.X) +r2.Point.Y
    return Point(x, y)


def as_XY( list_of_points ):
    x = list()
    y = list()
    for i in list_of_points:
        x.append(i.X)
        y.append(i.Y)
    return (x,y)

def rand_bearing( ):
    return random.random() * 4*pi - 2*pi 
    

def rb(a=-10,b=10):
    return random.randint(a,b)


def main():

    fig, ax = plt.subplots()
    A = Point(rb(), rb() )
    B = Point(rb(), rb() )
    pts = []
    pts.append( A)
    pts.append( B)
    ax.scatter( as_XY(pts)[0],  as_XY(pts)[1] ,color='red',  )
    Aray = Ray ( A, rand_bearing() )
    ax.add_patch( Aray.patch( scale = 50) )

    pts = [ ]

    for i in range(5):
        B_bearing = rand_bearing()
        Bray = Ray ( B, B_bearing)
        pts.append( intersect_lines( Aray, Bray))
        ax.scatter( as_XY(pts)[0],  as_XY(pts)[1] , color='blue',  )
        ax.add_patch( Bray.patch(scale = 10) )

    plt.axis('scaled')
    plt.show()

if __name__ == "__main__":
    main()
