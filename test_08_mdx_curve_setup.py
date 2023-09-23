import itertools
from Geo import Point, Segment, Curve, Chain, Ray, Angle, xy
import Geo 
import IPython
import matplotlib.pyplot as plt
import cmath
import collections
import pprint as pp

#setup for defaultdict to always produce a dict
tree = lambda: collections.defaultdict(tree)
t = tree()

R= 900
GDSPC= 8.0
n_girders = 4
delta_curve = Geo.Angle( -45, unit='deg')
spans = itertools.accumulate([ 200, 200, 150 ], initial=10)
sweeps = [90, 60, 80, 30 ]

cc = Point(1000,0)
t['cc'] = cc
pc = Point.from_complex( cc + complex(-R,0))
curve_1 = Curve( pc, cc, delta_curve)

pts = [pc,  curve_1.PT ]
brg = 0
for span, sweep in zip( spans, sweeps):
    brg += 1
    print(f"{span=}    {sweep=}")
    print(span)
    pt =  curve_1.move_to( span, 0)
    t['pts']['bent'][brg] = pt
    pts.append(pt)

pp.pprint(t)
pp.pprint(pts)
fig, ax = plt.subplots()
ax.scatter( *xy(pts) )
plt.axis('scaled')
plt.show()
    




"""
fig, ax = plt.subplots()
for i,R in enumerate([900, 908, 916, 924]): #radui
    mdx = "mdx_" + str(i + 1)
    pc = Point(1000-R,0)
    curve=Curve(pc, cc, delta_curve)
    d[mdx]['curves'].append(curve)
    d[mdx]['ch'] = Chain( mdx, curve )
    ax.add_patch( curve.patch() )


#for R in d.keys():
#    print( R )
#    myChain = d[R]['ch']
#    #ax.scatter( [myCurve.PC.X], [myCurve.PC.Y], marker='o', color='g')
#    #ax.scatter( [myCurve.PT.X], [myCurve.PT.Y], marker='o', color='r')
#    #ax.scatter( Geo.as_XY(points)[0], Geo.as_XY(points)[1], marker='+', color='g' )
#    #ax.scatter( Geo.as_XY(bad_pts)[0], Geo.as_XY(bad_pts)[1], marker='.', color='r' )
#    for rt in myChain.Routes:
#        ax.add_patch( rt.patch() )

bts = itertools.accumulate( [ 200,200,], initial=10)
sweeps = [90, 60, 80]
for bt,sweep in zip(bts, sweeps):
    ref = d['mdx_1']
    ref['pts'].append( ref['ch'].move_to(bt,0) )
    tangent = Segment(cc, ref['pts'][-1] ).Bearing().rt_90()
    brg_line = tangent + Angle(sweep, unit='deg')
    ref['pts'].append(Point.from_complex( ref['pts'][-1] + cmath.rect(100, brg_line)))
    print(ref['pts'])
    
    ref['seg'].append( Segment (ref['pts'][-2], ref['pts'][-1] ))
    print(ref['seg'])
    ax.add_patch( ref['seg'][-1].patch() )

for key in d.keys():
    pt = d[key]['pts']
    ax.scatter( Geo.as_XY(pt)[0], Geo.as_XY(pt)[1], marker='o', color='b' )

print(bts)
plt.axis('scaled')
plt.show()
"""
#IPython.embed()
