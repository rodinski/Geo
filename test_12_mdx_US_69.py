import itertools
from geo_ import Point, Bearing, Angle, Segment, Curve, Chain, Ray, Angle, xy
import Geo
import IPython
import matplotlib.pyplot as plt
import cmath
from math import pi
from collections import defaultdict
import collections
import pprint as pp
from curve_intersection_ import ray_curve_intersect


# setup for defaultdict to always produce a dict
class NestedDict(dict):    
    def __missing__(self, key):  
        self[key] = NestedDict()
        return self[key]

t = NestedDict()
pts = []
brg_pts = []
brace_pts = []
FS_pts = []

n_girders = 5
girder_spacing = 8.167
t['CL_R'] = 1200
cc = Point( t['CL_R'], 0)
t['CL_Bent_1'] = Point.from_complex((cc + complex(-t['CL_R'], 0)))
delta_curve = Geo.Angle(-30, unit='deg')
spans = itertools.accumulate([200, 270], initial=0)
sweeps = [90, 90, 90]

pc = Point.from_complex(cc + complex(-t['CL_R'], 0))
t['cc'] = cc
t['CL_pc'] = pc
t['CL_curve'] = Curve(t['CL_Bent_1'], cc, delta_curve)
t['curve_1'] = Curve(pc, cc, delta_curve)
for s in spans:
    brg_pts.append( t['CL_curve'].set_point(s,0))
    



s1 = 200 / 10
s2 = 270 /13
s3 = s1 + s2


brace=itertools.accumulate([s1, s1, s1, s1 ,s1,
                            s1, s1, s1, s3, s2,
                            s2, s2, s2, s2, s2,
                            s2, s2, s2, s2, s2   ], initial= s1 )


fig, ax = plt.subplots()
ax.add_patch(t['CL_curve'].patch())
ax.scatter(*xy(pts))
ax.scatter(*xy(brg_pts), marker='+', color='r')
ax.scatter(*xy(brace_pts), marker='+', color='r')
ax.scatter(*xy(FS_pts), marker='^', color='g')
plt.axis('scaled')
plt.show()


#IPython.embed()
