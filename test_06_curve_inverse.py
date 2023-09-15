import Geo
from Geo import Point, Angle
import matplotlib.pyplot as plt
import math
import cmath
import random


myCurve = Geo.Curve(Point(-.0001,.5), Point(0,0), Angle(-90, unit='deg'))

points =[]
bad_pts=[]
for theta in range(0,370,10):
    t_rad = math.radians(theta) 
    #print( cmath.rect(1, t_rad))
    pt = Point.from_complex(cmath.rect(1,t_rad))
    ret = myCurve.distance_and_offset(pt)
    d = ret[0]
    os = ret[1]
    if d and os:
        points.append(pt)
    else:
        bad_pts.append(pt)
    print( theta, ret )

x_min = myCurve.CC.X - myCurve.R
x_max = myCurve.CC.X + myCurve.R
y_min = myCurve.CC.Y - myCurve.R
y_max = myCurve.CC.Y + myCurve.R
for i in range(1000):
    pt = Point(random.uniform(x_min,x_max), random.uniform(y_min,y_max))
    ret = myCurve.distance_and_offset(pt)
    d = ret[0]
    os = ret[1]
    print(i, ret, myCurve.has_dist_os(pt), myCurve.distance_offset(pt) ) 
    if myCurve.has_dist_os(pt):
        points.append(pt)
    else:
        bad_pts.append(pt)


pt = Point(0.1,10)
print( myCurve.distance_and_offset(pt), myCurve.distance_offset(pt), myCurve.Len() ) 
pt = Point(10,0.1)
print( myCurve.distance_and_offset(pt), myCurve.distance_offset(pt), myCurve.Len() ) 

fig, ax = plt.subplots()

print (myCurve.PC.X)
print (myCurve.PC.Y)
ax.scatter( [myCurve.PC.X], [myCurve.PC.Y], marker='o', color='g')
ax.scatter( [myCurve.PT.X], [myCurve.PT.Y], marker='o', color='r')

ax.scatter( Geo.as_XY(points)[0], Geo.as_XY(points)[1], marker='+', color='g' )
ax.scatter( Geo.as_XY(bad_pts)[0], Geo.as_XY(bad_pts)[1], marker='.', color='r' )
for p in myCurve.patch_all():
    ax.add_patch(p)
plt.axis('scaled')
plt.show()
