import itertools
from Geo import Point, Segment, Curve, Chain
import Geo 
import IPython
import matplotlib.pyplot as plt

Rbase= 1234
R= Rbase
GDSPC= 8.5


pts = []
d = {}
d['mdx_1'] = {'pts': [], 'ch':None, 'seg':[] }
pts = d['mdx_1']['pts'] #just a refference for easy typing

pts.append( Point(R, -100) )
pts.append( Point(R,    0) )
pts.append( Point(0,    0) )
delta_curve = Geo.Angle( 60, unit='deg')

d['mdx_1']['ch'] =Chain ( "mdx_G_1" , Segment(pts[0], pts[1]), 
        Curve(pts[1], pts[2], delta_curve))



for name in range(2,6):
    i = name -1
    key = 'mdx_' +  str(name) 
    d[key] = {'pts': [], 'ch':None, 'seg':[]}
    R = Rbase - i * GDSPC
    pts = d[key]['pts']  # just a refference for easy typing
    pts.append( Point(R, -100) )
    pts.append( Point(R,    0) )
    pts.append( Point(0,    0) )
    d_1curve = Geo.Angle( 90, unit='deg')
    d[key]['ch'] = Chain ( "mdx_G_1" , Segment(pts[0], pts[1]), 
        Curve(pts[1], pts[2], delta_curve))

fig, ax = plt.subplots()

for R in d.keys():
    print( R )
    myChain = d[R]['ch']
    #ax.scatter( [myCurve.PC.X], [myCurve.PC.Y], marker='o', color='g')
    #ax.scatter( [myCurve.PT.X], [myCurve.PT.Y], marker='o', color='r')
    #ax.scatter( Geo.as_XY(points)[0], Geo.as_XY(points)[1], marker='+', color='g' )
    #ax.scatter( Geo.as_XY(bad_pts)[0], Geo.as_XY(bad_pts)[1], marker='.', color='r' )
    for rt in myChain.Routes:
        ax.add_patch( rt.patch() )

bts = itertools.accumulate( [ 180,240,180], initial=200)
for bt in bts:
    ref = d['mdx_1']
    ref['pts'].append( ref['ch'].move_to(bt,100) )
    ref['pts'].append( ref['ch'].move_to(bt,-100) )
    ref['seg'].append( Segment (ref['pts'][-2], ref['pts'][-1] ))
    ax.add_patch( ref['seg'][-1].patch() )

for key in d.keys():
    pt = d[key]['pts']
    ax.scatter( Geo.as_XY(pt)[0], Geo.as_XY(pt)[1], marker='o', color='b' )

print(bts)
plt.axis('scaled')
plt.show()

IPython.embed()
