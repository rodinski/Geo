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
class NestedDict(dict):    
    def __missing__(self, key):  
        self[key] = NestedDict()
        return self[key]

t = NestedDict()


R = 900
girder_spacing = 8.0
n_girders = 4
delta_curve = Geo.Angle(-45, unit='deg')
spans = itertools.accumulate([200, 200, 150], initial=10)
sweeps = [90, 60, 60, 30]

cc = Point(1000, 0)
pc = Point.from_complex(cc + complex(-R, 0))
t['cc'] = cc
t['pc'] = pc
curve_1 = Curve(pc, cc, delta_curve)
t['curve_1'] = curve_1

pts = [cc, pc, curve_1.PT]

bent = 0   # will increment
for span, sweep in zip(spans, sweeps):
    bent += 1
    print(f"{span=}    {sweep=}")
    print(span)
    pt = curve_1.move_to(span, 0)  # locate bents on curve_1 
    ref = t['pts']['bent'][bent] 
    ref['pt']['mdx'][1] = pt              # far point on Ray
    pts.append(pt)
    
    # bearing of brg line = tang + sweep 
    brg_brg = curve_1.tangent_at_point(pt) + Angle(sweep, unit='deg')
    ref['ray']  = Ray(pt, brg_brg)
    ref['pt']['left'] = ref['ray'].move_to(100) 
    ref['pt']['right'] = ref['ray'].move_to(-100) 
    pts.append(pt)
    for item in [ 'left', 'right']:            
        pts.append(ref['pt'][item])
        
    girder_spacing = 8.0
    for i in range(1, 5):  # loop over girders 1,2,3,4
        mdx_g = i + 1
        delta_R = i*girder_spacing
        R_i = R + delta_R

        c_ref = t['mdx_g'][mdx_g]['circle'] = Curve(Point.from_complex(cc+R_i), cc, -2*pi)
        r_ref = t['pts']['bent'][bent]['ray']
        pts.append( ray_curve_intersect(r_ref, c_ref)[1])    # second of two intersections






brace=itertools.accumulate([20, 20, 20, 20 ,20,
                            20, 20, 20, 40, 20,
                            20, 20, 20, 20, 20,
                            20, 20              ], initial=10 + 20)

t['int_braces']['mdx'][1] = []
brace_pts = []
for b in brace:
    print(b)
    b_pt = curve_1.move_to(b, 0)
    brace_pts.append(b_pt)
    t['int_braces']['mdx'][1].append(b_pt)

for mdx_g in [2, 3, 4, 5]:
    t['int_braces']['mdx'][mdx_g] = list()
    c_ref = t['mdx_g'][mdx_g]['circle']
    for pt in t['int_braces']['mdx'][1]:   #each brace on mdx_1
        r_ref = Ray(cc, Segment(cc, pt).Bearing())
        pt = ray_curve_intersect(r_ref, c_ref)[1]
        brace_pts.append(pt)
        t['int_braces']['mdx'][mdx_g].append(pt)


for mdx_g in [2, 3, 4, 5]:
    for i, brace in enumerate(t['int_braces']['mdx'][mdx_g]):
        outer = brace
        inner = t['int_braces']['mdx'][mdx_g-1][i]
        print( f"{mdx_g=}\t{i=}{abs(outer-inner)=}")


#print(f"\n\n{dict(t)=}")

fig, ax = plt.subplots()
ax.scatter(*xy(pts))
ax.scatter(*xy(brace_pts), marker='+', color='r')
plt.axis('scaled')
plt.show()


IPython.embed()
