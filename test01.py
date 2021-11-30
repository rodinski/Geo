from Geo import *


P[1] = Point( 0, 3)
P[2] = Point( 5, 3)
P[3] = Point(-5, 3)

C[1] = Curve(P[2], P[1], Delta(pi/6) )
P[4] = C[1].PI
C[2] = Curve(P[2], P[1], Delta(-pi/6) )
C[3] = Curve(P[3], P[1], Delta(pi/3) )
P[5] = C[3].PI
#C[4] = Curve(P[3], P[1], Delta(-pi/3) )

fig, ax = plt.subplots()
for k,v in P.items():

    print( k, v)
    ax.scatter(v.X,v.Y)
    
for k,c in C.items():
    print( k, c)
    ax.scatter(c.PT.X,c.PT.Y)
    ax.add_patch( c.patch(linestyle='-' ) ) 
    ax.add_patch( Segment( c.PC, c.PT).patch(linestyle=':') )

plt.axis('scaled')
plt.show()
