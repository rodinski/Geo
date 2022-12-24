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



PIs = [ Point(0,0), Point(5,0), Point(10,3), Point(15,4), \
        Point(20,0), Point(25,5), Point(30,20) ] 
PIs = []


bigR =5
n = 60
for i in range(0, 360+n, n):
    PIs.append( Point.from_complex(cmath.rect(bigR, math.radians(i)) ))
R=4.5
Segments = []
pts = []
fig, ax = plt.subplots()

pts.append(PIs[0])
pts.append(PIs[-1])

for i in range(1, len(PIs)-1):

    print( PIs[i-1], PIs[i], PIs[i+1],)
    PI = PIs[i]
    pts.append(PI)
    bk = Segment(PIs[i],PIs[i-1], )
    ah = Segment(PIs[i],PIs[i+1], )
 
    print(f"{ah.Bearing()=}")
    print(f"{bk.Bearing().Opposite()=}")
    print(f"{norm_as_delta(ah.Bearing()-bk.Bearing().Opposite())=}")
    #input("wait")
    delta = norm_as_delta(ah.Bearing() - bk.Bearing().Opposite() )

    print(f"{delta=}")
    print(f"{E(R,delta)=}")
    toCC = E(R,delta) + R
    a = normalize(ah.Bearing() )
    b = normalize(bk.Bearing() )

    delta_a_to_b = norm_as_delta(b-a)

    #need to make this delta part of the Bearing class
    #note one can not just take average of two bearings. 
    toCCbrg = normalize( a + delta_a_to_b/2 )


    CC = Point.from_complex(PI + cmath.rect(toCC, toCCbrg) )
    PC = Point.from_complex(CC + cmath.rect(R, toCCbrg -pi -delta/2.0) )
    PT = Point.from_complex(CC + cmath.rect(R, toCCbrg -pi +delta/2.0) )
    pts.append( CC )
    pts.append( PC )
    pts.append( PT )
    c = Curve( PC,CC, delta) 
    ax.add_patch( Segment(PC,PT).patch(linestyle=':') )
    ax.add_patch( c.patch() )
    
    
    print(f"{PI=}")
    print(f"{CC=}")
    print()
print(pts)
ax.scatter( as_XY(pts)[0],  as_XY(pts)[1] )
plt.axis('scaled')
plt.show()
