import random
from Geo import Point, Segment, Ray, Curve, Chain
import Geo
import random
import cmath
from math import pi
import matplotlib.pyplot as plt
from matplotlib import patches

pts = []
n_pts = []
rays = []
segs = []
curves= []
chains = []

os_x = 555
os_y = 666

for x in range(10):
    for y in range(10):
        pts.append( Point( x+os_x, y+os_y) )
        n_pts.append( Point( -x+os_x, -y+os_y))
        #n_pts.append( Point.from_complex( pts[-1] ))

rays.append( Ray(Point( os_x, os_y ), random.uniform( -pi, pi)))

all_pts = pts + n_pts
x_min = min( all_pts, key=lambda p: p.X).X
x_max = max( all_pts, key=lambda p: p.X).X

print( x_max - x_min) 

y_min = min( all_pts, key=lambda p: p.Y).Y
y_max = max( all_pts, key=lambda p: p.Y).Y

def xy(inlist:list)->tuple:
    x = [ arg.X for arg in inlist ]
    y = [ arg.Y for arg in inlist ]
    return(x,y)

fig, ax = plt.subplots()
curves.append( Curve(Point(os_x, y_max), Point(os_x,os_y),-0.75*pi))
curves.append( Curve(Point(x_min, os_y), Point(os_x,os_y), 0.75*pi))
for c in curves:
    for p in c.patch_all():
        ax.add_patch( p )

ax.scatter( *xy(pts), marker='+', color='r' )
ax.scatter( *xy(n_pts), marker='o', color='g' )
plt.axis('scaled')
plt.show()

