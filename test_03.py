import cmath
from math import pi
from Geo import *

import matplotlib.pyplot as plt
from   matplotlib.path import Path
import matplotlib.patches as patches 

def E(R, delta)->float:
    return  R*(1/math.cos(delta/2.0) - 1)

def as_XY( list_of_points ):
    x = list() 
    y = list() 
    for i in list_of_points:
        x.append(i.X)
        y.append(i.Y)
    return (x,y)


segments = []
curves = []
segments.append( Segment( Point(0,0), Point(100,0)))

pc_bearing = segments[-1].Bearing()
pc =  segments[-1].Pt2
R = 150

#cc = pc + cmath.rect(R, pc_bearing + pi/2)
#print(f"{cc=}")
print(f"{pc_bearing=}")
print(f"{R=}")

c1 = Curve.from_PC_bearing_R_Delta( pc, pc_bearing, R, pi/2*1.85)

myChain = Chain("given name is oval", segments[0] )
myChain.addRoute(c1)
myChain.forward(100)
print(myChain.RoutesSta) 

lastRay = myChain.Routes[-1].outRay()
c2 = Curve.from_Ray_T_Delta( lastRay, 375, 0.2*pi )
myChain.addRoute(c2)




fig, ax = plt.subplots()
#ax.scatter( as_XY(pts)[0],  as_XY(pts)[1] )

for r in myChain.Routes:
    ax.add_patch( r.patch(linestyle=':'))
    print(r)
    pass
plt.axis('scaled')
plt.show()


print(myChain.Routes[2].inRay().distance_to( Point(135,295)))

