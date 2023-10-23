import geo_
from geo_ import  Point, Ray, Segment, Bearing, xy
from math import pi, tan , atan2, atan, cos, degrees, radians, isclose
import matplotlib.pyplot as plt
from matplotlib import patches
import random

def bearing_as_slope( theta:float )->float:
    """ Good for Rays that are to be treated as lines
    always return change in y as x INCREASES
    tan(theta) does this with a check of vertical
    """
    if isclose( theta, -pi/2) or isclose( theta, pi/2):
    #if theta in { Bearing(-pi/ 2) , Bearing(pi/2) }:
        return float( 'inf' )   
    return tan(theta)

def intersect_lines( r1:Ray, r2:Ray)-> Point:
    # same point?
    if r1.Point == r2.Point:
        return r1.Point
    # same bearing?
    m1 = bearing_as_slope(r1.bearing)
    m2 = bearing_as_slope(r2.bearing)

    if isclose(m1, m2):
        return None
    # set x
    if m1 == float('inf'):
        x = r1.Point.X
    elif m2 == float('inf'):
        x = r2.Point.X
    else:
        x = (-r2.Point.X*m2 -r1.Point.Y  +r1.Point.X*m1 +r2.Point.Y ) / (m1-m2) 
    #set y
    if isclose(m1, 0.0):
        y = r1.Point.Y
    elif isclose(m2, 0.0):
        y = r2.Point.Y
    else:
        if m1 != float('inf'):
            y = m1*(x - r1.Point.X) +r1.Point.Y
        else: 
            y = m2*(x - r2.Point.X) +r2.Point.Y
    return Point(x, y)



def rand_bearing( ):
    return random.random() * 4*pi - 2*pi 

def rb(a=-10,b=10):
    return random.randint(a,b)

def main():

    fig, ax = plt.subplots()
    #pts = []
    #rays = []
    #xtheta = pi/5
    #ytheta = pi/10
    #lines = 5
    #for x in range(0, lines):
    #    xRay = Ray(Point(x,0), xtheta)
    #    ax.add_patch( xRay.patch( scale = 5) )
    #    for y in range(0, lines):
    #        yRay = Ray(Point(0,y), ytheta)
    #        ax.add_patch( yRay.patch( scale = 5) )
    #        intersect = intersect_lines(xRay,yRay)
    #        if intersect:
    #            pts.append(intersect)
    #ax.scatter( *xy(pts) ,color='red',marker='+'  )

    #See if rays that should cross at right angle inscribe
    #a half-circle
    D=15
    c_pts = []
    segments = []
    for x in range(100):
        theta_btm = random.random()*pi/2
        ray_btm = Ray(Point(0,0), Bearing(theta_btm) )

        theta_top = theta_btm-pi/2
        ray_top = Ray(Point(0,D), Bearing(theta_top))
        intersection_pt = intersect_lines(ray_top, ray_btm) 

        c_pts.append(intersection_pt)
        segments.append(Segment(Point(0,0), intersection_pt))
        segments.append(Segment(Point(0,D), intersection_pt))

    for s in segments:
        ax.add_patch(s.patch(lw=0.1))

    ax.scatter( *xy(c_pts) ,color='blue',marker='.'  )

    plt.axis('scaled')
    plt.show()

if __name__ == "__main__":
    main()
