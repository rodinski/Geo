import itertools
from geo_ import Point, Bearing, Angle, Segment, Curve, Chain, Ray, Angle, xy
#import Geo
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
s = NestedDict()
t['curves'] = []
t['segments'] =[]
pts = []

myChain = Chain( Segment( Point(2764400.2656, 1087103.3073), Point(2764401.7092, 1087261.0125)), start_station=16500.00 )
myChain.addRoute( Curve( Point(2764401.7092, 1087261.0125), Point(2762551.2067, 1087277.9524), Angle(+0.231232738031)) )
myChain.addRoute( Curve( Point(2764356.3396, 1087685.5572), Point(2763185.8096, 1087421.2477), Angle(+0.370495973988)) )
myChain.addRoute( Segment( Point(2764181.2163, 1088091.446), Point(2764059.354, 1088272.4409)) )
myChain.addRoute( Curve( Point(2764059.354, 1088272.4409), Point(2765054.7605, 1088942.6393), Angle(-0.405282539789)) )
myChain.addRoute( Segment( Point(2763875.746, 1088719.2003), Point(2763733.4805, 1089469.889)) )

t["BentSta"][1] = 17348.18
t["BentSta"][2] = 17473.18
t["BentSta"][3] = 17598.18
t["BentSta"][4] = 17708.18
t["BentSta"][5] = 17818.18
t["BentSta"][6] = 17923.56
t["BentSta"][7] = 18016.52
t["BentSta"][8] = 18169.83
t["BentSta"][9] = 18329.83
t["BentSta"][10] = 18362.33


t["BentRot"][1] =  Angle(0, unit='deg')
t["BentRot"][2] =  Angle(0, unit='deg')
t["BentRot"][3] =  Angle(0, unit='deg')
t["BentRot"][4] =  Angle(0, unit='deg')
t["BentRot"][5] =  Angle(0, unit='deg')
t["BentRot"][6] =  Angle(-30, unit='deg')
t["BentRot"][7] =  Angle(-54, unit='deg')
t["BentRot"][8] =  Angle(-55, unit='deg')
t["BentRot"][9] =  Angle(-50, unit='deg')
t["BentRot"][10] = Angle(-15, unit='deg')


#for i in range(1,11):
#    point = myChain.set_point(t["BentSta"][i], 0.0)
#    print( point )
#    print( myChain.inverse(point) )
#    print( myChain.normal_at_point( point ) )
#    input("pause")

bearing_length = 50
# Point and CL of each bent

for i in range(1,11):
    sta =t["BentSta"][i] 
    angle = t["BentRot"][i] 
    
    t["Bent_CL_Pt"][i]= point_CL = myChain.set_point(t["BentSta"][i], 0.0)
    pts.append(point_CL)
    bearing = Bearing( myChain.normal_at_point(point_CL)+ angle )
    t["Bearing_Ray"][i] = ray_rt = Ray(t["Bent_CL_Pt"][i],  bearing )

    t["Bent_Right_Point"][i] = point_right = ray_rt.set_point( bearing_length, 0)
    t["Bent_Left_Point"][i] = point_left = ray_rt.set_point(-bearing_length, 0)

    t["Bent_CL"][i] = Segment(point_left, point_right)
    t["segments"].append( Segment(point_left, point_right))

    # set bearing lines at int_bents
    if i > 1 and i < 11:
        bent_CL_offset = 1.0
        if i in [ 6,7]:  # high skew 
            bent_CL_offset = 1.25
        t["Bent_Bearing"][i]["ahead"] = ahead = t["Bent_CL"][i].copy_parallel(-bent_CL_offset)
        t["Bent_Bearing"][i]["back"]  = back =  t["Bent_CL"][i].copy_parallel( bent_CL_offset)
        t["segments"].append( ahead )
        t["segments"].append( back )

# Chords b/w Piers
for i in range(1,10):
    t["Span_Chord"][i] = Segment(t["Bent_CL_Pt"][i], t["Bent_CL_Pt"][i+1] )
    t['segments'].append( t["Span_Chord"][i] ) 
    #t['segments'].append( t["Span_Chord"][i].copy_parallel(10) ) 

t["G3_os"][1]["ahead"] = 0.75
t["G3_os"][2]["back"] =  0.75

t["G3_os"][2]["ahead"] = 0.333
t["G3_os"][3]["back"] =  0.333

t["G3_os"][3]["ahead"] = 0.0
t["G3_os"][4]["back"] =  2.0

t["G3_os"][4]["ahead"] = 1.666
t["G3_os"][5]["back"] =  1.666

t["G3_os"][5]["ahead"] = 1.5
t["G3_os"][6]["back"] =  1.5

t["G3_os"][6]["ahead"] = 2.0
t["G3_os"][7]["back"] =  2.0

t["G3_os"][7]["ahead"] = 0.0
t["G3_os"][8]["back"] =  0.0

t["G3_os"][8]["ahead"] = 1.75
t["G3_os"][9]["back"] =  1.75


# Girder 3
bearing_pts = []
t["segments_G3"] = []
for i in range(1,9):
    ref = t["Beam_span"][i]["G"][3]
    #print(t["Bearing_Ray"][i  ],   t["G3_os"][i]["ahead"])
    #use the in/outRay of each Span_Chord
    ref["start"] = t["Span_Chord"][i].inRay().set_point( 0.0, t["G3_os"][i]["ahead"] )
    ref["end"]  =  t["Span_Chord"][i].outRay().set_point( 0.0, t["G3_os"][i+1]["back"] )
    t["segments_G3"].append( Segment( ref["start"], ref["end"] ) ) 

    bearing_pts.append( ref["start"] )
    bearing_pts.append( ref["end"  ] )


#right edge of deck
#t["BentSta"][1] = 17348.18
#t["BentSta"][10] = 18362.33
my_Right_EOD = myChain.copy_parallel( 21.333 ).split(start_sta=800, end_sta=800+900)
my_Left_EOD =  myChain.copy_parallel( -21.333).split(start_sta=800, end_sta=800+700)

   

#plt.figure(figsize=(8, 8))
fig, ax = plt.subplots( figsize=(9,9) )

ax.scatter(*xy(pts))
ax.scatter(*xy(bearing_pts), marker='+', color='r')

for curve in t['curves']:
    for p in curve.patch_all():
        ax.add_patch( p)

for segment in t['segments']:
    ax.add_patch(segment.patch() )
    pts.append(segment.Pt1)
    pts.append(segment.Pt2)

for segment in t['segments_G3']:
    ax.add_patch(segment.patch( linewidth=3, linestyle='--'))
ax.scatter( *xy(pts) )



for chain in [ myChain, my_Right_EOD, my_Left_EOD ]:
    for patch in chain.patch_list():
        ax.add_patch(patch)

plt.axis('scaled')
ax.set_xlim( 2763600, 2763600+800 )
ax.set_ylim(  1087600, 1087600+1500  )
plt.show()


def walkDict_values( inDict, depth=0, retDict={}):
    retStr = ""
    """walk all the NestedDict and dict"""
    pre = "\t"*depth
    for k, v in inDict.items():
        if not isinstance(v, (NestedDict, dict)):
            retStr += (f"{pre}{k} -> {type(v)}\n")
            if type(v) not in retDict:
                retDict[type(v)] = list()
            retDict[type(v)].append(v)
        else:
            print(f"{pre}{k}\n")
            walkDict_values(v, depth=depth +1, retDict=retDict )
    depth -= 1
    return retDict


def walkDict( inDict, depth=0):
    """walk all the NestedDict and dict"""
    pre = "\t"*depth
    keyDict = {}
    for k, v in inDict.items():
        if not isinstance(v, (NestedDict, dict)):
            print(f"{pre}{k} -> {type(v)} ")
            keyDict[k] ={}
        else:
            print(f"{pre}{k}")
            walkDict(v, depth=depth +1)
    depth -= 1
    return keyDict

#walkDict(t) 
IPython.embed()
