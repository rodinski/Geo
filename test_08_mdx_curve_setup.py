import itertools
from Geo import Point, Segment, Curve, Chain, Ray, Angle, xy
import Geo
import IPython
import matplotlib.pyplot as plt
import cmath
from math import pi
from collections import defaultdict
import collections
import pprint as pp
from curve_intersection import ray_curve_intersect


# setup for defaultdict to always produce a dict

nested_dict = lambda: defaultdict(nested_dict)
t = nested_dict()

def walk_tree(inTree):
    for k, v in inTree.items():
        if isinstance(v, dict):
            print(f"got on {v=}")
            walk_tree(v)


R = 900
GDSPC = 8.0
n_girders = 4
delta_curve = Geo.Angle(-45, unit='deg')
spans = itertools.accumulate([200, 200, 150], initial=10)
sweeps = [90, 60, 80, 30]

cc = Point(1000, 0)
t['cc'] = cc
pc = Point.from_complex(cc + complex(-R, 0))
curve_1 = Curve(pc, cc, delta_curve)

pts = [pc, curve_1.PT, cc]
brg = 0
for span, sweep in zip(spans, sweeps):
    brg += 1
    print(f"{span=}    {sweep=}")
    print(span)
    pt = curve_1.move_to(span, 0)
    ref = t['pts']['bent'][brg]
    ref['pt'][1] = pt
    brg_brg = curve_1.tangent_at_point(pt) + \
              Angle(sweep, unit='deg')
    ref['ray']  = Ray(pt, brg_brg)
    ref['pt'][2] = ref['ray'].move_to(40)
    pts.append(pt)
    for i in [1, 2]:
        pts.append(ref['pt'][i])
brace=itertools.accumulate([20, 20, 20, 20 ,20,
                           20, 20, 20, 20 , 40,
                           20, 20, 20, 20 , 20,
                           20, 20, 20   ], initial=10)
for b in brace:
    print(b)
    pts.append(curve_1.move_to(b, 0))


GDSPC = 8.0
for i in range(1, 5):
    mdx_g = i + 1
    delta_R = i*GDSPC
    R_i = R + delta_R

    t[mdx_g][i]['circle'] = Curve(cc+R_i, cc, -2*pi)
    c_ref = t[mdx_g][i]['circle']
    r_ref = t['pts']['bent'][2]['ray']
    print(ray_curve_intersect(r_ref, c_ref) )
    pts.append( ray_curve_intersect(r_ref, c_ref)[1])


#print(f"\n\n{dict(t)=}")

fig, ax = plt.subplots()
ax.scatter(*xy(pts))
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
# IPython.embed()
