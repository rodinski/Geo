import civil.utility  as civil  #move this
import civil.vertical as vertical   #move this
from geo_ import Point, Bearing, Angle, Segment, Curve, Chain, Ray, Angle, xy, Distance, translate_rotate
import itertools
import IPython
import matplotlib.pyplot as plt
import cmath
from math import pi, degrees, sqrt, sin, cos, tan
from collections import defaultdict
import collections
import pprint as pp
from curve_intersection_ import ray_curve_intersect
from line_intersect_ import intersect_lines
from itertools import accumulate
from collections import namedtuple

import curve_intersection
# setup for defaultdict to always produce a dict
class NestedDict(dict):    
    def __missing__(self, key):  
        self[key] = NestedDict()
        return self[key]
    

def walkDict( inDict, depth=0):
     
    pre = "\t"*depth
    for k, v in inDict.items():
        if not isinstance(v, (NestedDict, dict)):
            print(f"{pre}{k} -> {type(v)} ")
        else:
            print(f"{pre}{k}")
            walkDict(v, depth=depth +1 )
    depth -= 1


def find_key(d:dict, value):
    """walks a nested dict looking for a value. Returns an order list of
    of keys that references the value"""
    for k,v in d.items():
        if isinstance(v, dict):
            p = find_key(v, value)
            if p:
                return [k] + p
        elif v == value:
            return [k]


def walkDict_values( inDict, depth=0, retDict={}):
    """walk a dict and returns a dict of lists. The members of a list share the same type. 
    The keys are all the unique types in the inDict"""
    retStr = ""
    """walk all the NestedDict and dict"""
    pre = "\t"*depth
    for k, v in inDict.items():
        if not isinstance(v, (NestedDict, dict)):
            retStr += (f"{pre}{k} -> {type(v)}")
            if type(v) not in retDict:
                retDict[type(v)] = list()
            retDict[type(v)].append(v)
        else:
            #print(f"{pre}{k}")
            walkDict_values(v, depth=depth +1, retDict=retDict )
    depth -= 1
    return retDict


t = NestedDict()
c='curves'
s='segments'
Ch='chains'
Ra = 'rays'
p='pts'
t[s]= []  #list of segments

pts = []


#set CL
t["CL-K96"] = ( Segment(Point(0, 0), Point(1000,0)) )

#make Chain of road
t[Ch]["Road"] = Chain( t["CL-K96"] , name="Road", start_station=0.00)

ref = t[Ch]["Road"]

#set ProGr
t[Ch]["ProGr"] =  ref.copy_parallel( 32, start_station=0, name="ProGr")
#t[Ch]["GL_A"]

t["ProGr"][1] = ref.set_point(136.64, 32) 
t["ProGr"][2] = ref.set_point(185.64, 32) 
t["ProGr"][3] = ref.set_point(261.64, 32) 
t["ProGr"][4] = ref.set_point(337.64, 32) 
t["ProGr"][5] = ref.set_point(409.64, 32) 
for p in t["ProGr"].values():
    pts.append(p)

myPierBrg = Bearing( -(49 + 21/60), unit='deg') 
#set Ray at each pier
print( t[Ra] )

for i in range(1,6):   #5 rays
    t[Ra][i] = Ray( t["ProGr"][i], myPierBrg) 
print(t[Ra])

#set GL A B C

t["GL"]["A"] = t["CL-K96"].copy_parallel( 14 - 1.5 + 4 )

t["GL"]["P_1"]["G_A"] = intersect_lines(t["GL"]["A"].inRay(), t[Ra][1] )
t["GL"]["P_5"]["G_A"] = intersect_lines(t["GL"]["A"].inRay(), t[Ra][5] )

#Graph GL A 
t[s].append(Segment(t["GL"]["P_1"]["G_A"], t["GL"]["P_5"]["G_A"]))

t["GL"]["B"] = t["GL"]["A"].copy_parallel( 7 )
t["GL"]["C"] = t["GL"]["A"].copy_parallel( 14 )
t["GL"]["D"] = t["GL"]["A"].copy_parallel( 21 )

t["Bearing"]["P_1"]["G_D"] = intersect_lines(t["GL"]["D"].inRay(), t[Ra][1] )
t["Bearing"]["P_2"]["G_D"] = intersect_lines(t["GL"]["D"].inRay(), t[Ra][2] )
t["Bearing"]["P_3"]["G_D"] = intersect_lines(t["GL"]["D"].inRay(), t[Ra][3] )
t[s].append(Segment(t["Bearing"]["P_1"]["G_D"], t["Bearing"]["P_3"]["G_D"]))
print(t[s])


#t['segments'].append(t["GL"]["A"])
#t['segments'].append(t["GL"]["B"])
#t['segments'].append(t["GL"]["C"])

#set several points at intersection of Pier and GLs
#actually intersection of Ray with Ray
for GL in ["A", "B", "C"]:
    for i in range(1,6):
        #set many of the GL Pier Points
        pts.append( intersect_lines(t["GL"][GL].inRay(), t[Ra][i] ))
        t["Bearing"]["P_" + str(i)]["G_" + GL] = pts[-1]



#t["Bearing"]["P_5"]["G_C"] = pts[-1]
#t["Bearing"]["P_4"]["G_C"] = pts[-2]
#t["Bearing"]["P_3"]["G_C"] = pts[-3]

pts.append( intersect_lines(t["GL"]["C"].copy_parallel(7).inRay(), t[Ra][2] ))
pts.append( intersect_lines(t["GL"]["C"].copy_parallel(7).inRay(), t[Ra][1] ))
pts.append( intersect_lines(t["GL"]["C"].copy_parallel(7).inRay(), t[Ra][2] ))
pts.append( intersect_lines(t["GL"]["C"].copy_parallel(7).inRay(), t[Ra][3] ))

# next 3 all start form _D  
pts.append( Point.from_complex( t["Bearing"]["P_1"]["G_D"]+ cmath.rect( 30+3.75/12, myPierBrg )))
t["Bearing"]["P_1"]["G_G"] = pts[-1]

pts.append( Point.from_complex( t["Bearing"]["P_2"]["G_D"]+ cmath.rect( 28+2.625/12, myPierBrg )))
t["Bearing"]["P_2"]["G_G"] = pts[-1]

pts.append( Point.from_complex( t["Bearing"]["P_3"]["G_D"]+ cmath.rect( 24+11.25/12, myPierBrg )))
t["Bearing"]["P_3"]["G_G"] = pts[-1]

for i in range(1,6):
    print(t["Bearing"]["P_" + str(i)].keys())

# next 2 all start form _C  
pts.append( Point.from_complex( t["Bearing"]["P_4"]["G_C"]+ cmath.rect( 30+11/12, myPierBrg )))
t["Bearing"]["P_4"]["G_G"] = pts[-1]

pts.append( Point.from_complex( t["Bearing"]["P_5"]["G_C"]+ cmath.rect( 27+10/12, myPierBrg )))
t["Bearing"]["P_5"]["G_G"] = pts[-1]

t[s].append(Segment(t["Bearing"]["P_1"]["G_G"], t["Bearing"]["P_5"]["G_G"]))


#set new beams
pts.append( Point.from_complex( t["Bearing"]["P_1"]["G_G"]+ cmath.rect( 5 + 7.75/12, myPierBrg )))
t["Bearing"]["P_1"]["G_J"] = pts[-1]
pts.append( Point.from_complex( t["Bearing"]["P_2"]["G_G"]+ cmath.rect( 7 + 8.25/12, myPierBrg )))
t["Bearing"]["P_2"]["G_J"] = pts[-1]

pts.append( Point.from_complex( t["Bearing"]["P_3"]["G_G"]+ cmath.rect( 5 + 5.75/12, myPierBrg )))
t["Bearing"]["P_3"]["G_H"] = pts[-1]
pts.append( Point.from_complex( t["Bearing"]["P_3"]["G_G"]+ cmath.rect( 10 + 11.5/12, myPierBrg )))
t["Bearing"]["P_3"]["G_J"] = pts[-1]

pts.append( Point.from_complex( t["Bearing"]["P_4"]["G_G"]+ cmath.rect( 7+( 5+ 5/16) /12, myPierBrg )))
t["Bearing"]["P_4"]["G_H"] = pts[-1]
pts.append( Point.from_complex( t["Bearing"]["P_4"]["G_G"]+ cmath.rect( 14+(2+ 5/8 )/12, myPierBrg )))
t["Bearing"]["P_4"]["G_J"] = pts[-1]

pts.append( Point.from_complex( t["Bearing"]["P_5"]["G_G"]+ cmath.rect( 8+( 7+ 7/8) /12, myPierBrg )))
t["Bearing"]["P_5"]["G_H"] = pts[-1]
pts.append( Point.from_complex( t["Bearing"]["P_5"]["G_G"]+ cmath.rect( 17+(3+ 3/4) /12, myPierBrg )))
t["Bearing"]["P_5"]["G_J"] = pts[-1]


print( t[Ch]["Road"].inverse( t["Bearing"]["P_1"]["G_J"]) )
print( t[Ch]["Road"].inverse( t["Bearing"]["P_5"]["G_J"]) )


fig, ax = plt.subplots()
fig.set_size_inches((16,12))

#plot everythinkg in the list of pts
#print(f"{pts=}")

ax.scatter( *xy(pts), marker='+', color='b' )

for k,v  in t[Ch].items():
#    print(v)
    for p in v.patch_list():
        ax.add_patch(p)

for segment in t['segments']:
    pass
    ax.add_patch(segment.patch())

#listed stations at supports
sta_list = sorted( [ ])    # list of stations to look at
for sta in sta_list: 
    print( f"{civil.format_STA(sta)},   El {pgl.profile_grade(sta)}")

plt.axis('scaled')
ax.set_xlim([ 100, 600])
ax.set_ylim([ -100, 0])

#plt.show()
#import pdb; pdb.set_trace()
t[Ch]["Road"].inverse( Point(99, -3))

print( walkDict_values(t) )

for k,mylist in walkDict_values(t).items():
    for obj in mylist:
        print( obj.__repr__(), find_key(t, obj))

import IPython
IPython.embed()
