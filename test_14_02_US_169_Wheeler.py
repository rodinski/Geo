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
from line_intersect_ import intersect_lines
import yaml


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
    t["Bent_CL_Ray"][i] = ray_rt = Ray(t["Bent_CL_Pt"][i],  bearing )

    t["Bent_Right_Point"][i] = point_right = ray_rt.set_point( bearing_length, 0)
    t["Bent_Left_Point"][i] = point_left = ray_rt.set_point(-bearing_length, 0)

    t["Bent_CL"][i] = Segment(point_left, point_right)
    t["segments"].append( Segment(point_left, point_right))

    # set bearing lines at end_bents
    t["Bent_Bearing"][1]["ahead"] = t["Bent_CL"][1]     
    t["Bent_Bearing"][11]["back"] = t["Bent_CL"][11]

    if i > 1 and i < 11:
        bent_CL_offset = 1.0
        if i in [ 5]:  # special support 5
            bent_CL_offset = 1.4166
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

t["G3_os"][1] = 0.75
t["G3_os"][2] = 0.333
t["G3_os"][3] = None
t["G3_os"][4] = 1.666
t["G3_os"][5] = 1.5
t["G3_os"][6] = 2.0
t["G3_os"][7] = None
t["G3_os"][8] = 1.75


# Girder 3
bearing_pts = []
t["segments_G3"] = []


#not able to set girders in spans 3 or 7 'continue' on these cases
for span in range(1,9):
    if span == 3 or span == 7:
        continue 
    ref = t["Beam_span"][i]["G"][3]
    #make off set the chord then intersect it with the two bearing lines
    myBmWeb = t["Span_Chord"][span].copy_parallel( t["G3_os"][span] )

    ref["start"] = intersect_lines(myBmWeb.inRay(), t["Bent_Bearing"][span]["ahead"].inRay() )
    ref["end"] = intersect_lines(myBmWeb.inRay(), t["Bent_Bearing"][span+1]["back"].inRay() )

    bearing_pts.append( ref["start"] )
    bearing_pts.append( ref["end"] )

parallel = [8.9166, 9.1666, None, 8.0833, 7.9166, 7.9166, None, 8.000]
for span, spac in zip(range(1,9), parallel):
    if spac == None:
        continue
    print(span, spac)

    #already have zero loacation = G3
    for n in [-2, -1, 0,  1, 2]:
        ref = t["Beam_span"][span]["G"][n+3]
        #make off set the chord then intersect it with the two bearing lines
        myBmWeb = t["Span_Chord"][span].copy_parallel( t["G3_os"][span]  + spac * n )

        ref["start"] = intersect_lines(myBmWeb.inRay(), t["Bent_Bearing"][span]["ahead"].inRay() )
        ref["end"] = intersect_lines(myBmWeb.inRay(), t["Bent_Bearing"][span+1]["back"].inRay() )

        bearing_pts.append( ref["start"] )
        bearing_pts.append( ref["end"] )


# how to input span3 and span 7 girders??
# for span 7 the CL of all the gireders is defined by a movement along the CL of 
# piers 7 and 8.  Starting with GL 3, from Bent_CL_Pt_7 move 1.75' SE along the CL
# of the bent 7  place a point.   From Bent_CL_Pt_8 move 1.50' SE along the CL of bent 8
# and place a point.  Girder 3 lies on this line. 


#set span 3 bearings   Bearin_Ray != Bent_Bearing
CL_cap3_os = 0.000 
spacing_CL_3 = 9.1666

CL_cap4_os = 2.000
spacing_CL_4 = 8.1666
for i in range(1,6):
    t['span_3']['web_CL']['start'][i] = ref = t["Bent_CL_Ray"][3].set_point( CL_cap3_os + (i - 3) * spacing_CL_3 )
    pts.append(ref)

    t['span_3']['web_CL']['end'][i] = ref = t["Bent_CL_Ray"][4].set_point( CL_cap4_os + (i - 3) * spacing_CL_4 )
    pts.append(ref)

    web_cl = Segment(t['span_3']['web_CL']['start'][i] , t['span_3']['web_CL']['end'][i] )

    t["Beam_span"][3]['G'][i]['start'] = ref =intersect_lines( web_cl.inRay(), t["Bent_Bearing"][3]['ahead'].inRay() ) 
    bearing_pts.append(ref)
    t["Beam_span"][3]['G'][i]['end'] = ref =intersect_lines( web_cl.inRay(), t["Bent_Bearing"][4]['back'].inRay() ) 
    bearing_pts.append(ref)

#span 3 bearings set
del(ref) 




#set span 7 bearings   Bearin_Ray != Bent_Bearing
CL_cap7_os = 1.75 
spacing_CL_7 = 13.250

CL_cap8_os = 1.50
spacing_CL_8 = 13.375
for i in range(1,6):
    t['span_7']['web_CL']['start'][i] = ref = t["Bent_CL_Ray"][7].set_point( CL_cap7_os + (i - 3) * spacing_CL_7 )
    pts.append(ref)

    t['span_7']['web_CL']['end'][i] = ref = t["Bent_CL_Ray"][8].set_point( CL_cap8_os + (i - 3) * spacing_CL_8 )
    pts.append(ref)

    web_cl = Segment(t['span_7']['web_CL']['start'][i] , t['span_7']['web_CL']['end'][i] )

    t["Beam_span"][7]['G'][i]['start'] = ref =intersect_lines( web_cl.inRay(), t["Bent_Bearing"][7]['ahead'].inRay() ) 
    bearing_pts.append(ref)
    t["Beam_span"][7]['G'][i]['end'] = ref =intersect_lines( web_cl.inRay(), t["Bent_Bearing"][8]['back'].inRay() ) 
    bearing_pts.append(ref)
#span 7 bearings set
del(ref) 


    

#right edge of deck
#t["BentSta"][1] = 17348.18
#t["BentSta"][10] = 18362.33
my_Right_EOD = myChain.copy_parallel( 21.333 ).split(start_sta=800, end_sta=800+900)
my_Left_EOD =  myChain.copy_parallel( -21.333).split(start_sta=800, end_sta=800+700)

   

#plt.figure(figsize=(8, 8))
fig, ax = plt.subplots( figsize=(9,9) )

ax.scatter(*xy(pts))
ax.scatter(*xy(bearing_pts), marker='x', color='r', s=80)

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
#IPython.embed()


'''
In [6]: walkDict(t)
'''
# input()


for b in range(1,10):
    print()
    for g in range(1, 6):
        ref = t["Beam_span"][b]["G"]
        val = ref.get(g, None)
        if val is None:
            continue
        if g in ref.keys():
          # do both start and end
          if ref[g].get('start', None) and ref[g].get('end', None):
              s_val = myChain.inverse( ref[g]['start'] )
              e_val = myChain.inverse( ref[g]['end'] )
              print( f't["Beam_span"][{b}]["G"][{g}]\t{s_val.distance:.2f}\t{s_val.offset:6.2f}\t{e_val.distance:.2f}\t{e_val.offset:6.2f} ' )
          else:
              continue
              
import civil.vertical
'''
class ProG:
    def __init__(self, name, pointList ):
'''
''' ProG are defined with a list.  
The list must start and end with a "named point tuple"
midpoints must be "named pi tuples" 
Once the input is validated two lists are created
self.pog = list of "named pc and pt tuples"
self.set = list of "named segment tuples"
'''

import civil.vertical as vc
myPG = vc.ProG( name = 'PGL', 
    pointList=[

    vc.point( 16905.00,764.91), 
    vc.pi( 17175.00, 765.06,540),
    vc.pi( 17740.00, 791.00, 450), 
    vc.pi( 17965.00, 792.46, 450), 
    vc.pi( 18190.00, 793.92, 450),  
    vc.point( 18415.00,787.01) 
    ])

#print(myPG)
#print( dir(myPG) )
import cross_slopes_US_169_Wheeler as xs
#IPython.embed()

BI = dict()
BI["Span"] = dict()
for b in range(1,10):
    BI["Span"][b] = dict()
    BI["Span"][b]["G"] = dict()

    print()
    for g in range(1, 6):
        BI["Span"][b]["G"][g] = dict()
        BI["Span"][b]["G"][g]["Nth"] = dict()

        ref = t["Beam_span"][b]["G"]
        val = ref.get(g, None)
        if val is None:
            continue
        if g in ref.keys():
          # do both start and end
          if ref[g].get('start', None) and ref[g].get('end', None):
              s_val = myChain.inverse( ref[g]['start'] )
              e_val = myChain.inverse( ref[g]['end'] )

              s_ProG = myPG.profile_grade(s_val.distance)
              e_ProG = myPG.profile_grade(e_val.distance)
              # print(s_ProG)
              # print(s_val.distance, s_val.offset)
              # input()

              s_xs_correction = xs.xs_correction(s_val.distance, s_val.offset)
              s_Deck = round( s_ProG  + s_xs_correction, 3 )

              e_xs_correction = xs.xs_correction(e_val.distance, e_val.offset)
              e_Deck = round( e_ProG  + e_xs_correction, 3 )

              res = [ f'"Beam_span"][{b}]["G"][{g}]', 
                      f'{s_val.distance:.2f}', 
                      f'{s_ProG:.2f}', 
                      f'{s_xs_correction:.2f}', 
                      f'{s_Deck:.2f}', 

                      f'{e_val.distance:.2f}', 
                      f'{e_ProG:.2f}', 
                      f'{e_xs_correction:.2f}', 
                      f'{e_Deck:.2f}' ]

              BI["Span"][b]["G"][g]["Nth"][0]  = { "Sta": s_val.distance, "offset": s_val.offset,  "ProG": s_ProG, "xs_corr": s_xs_correction, "el_Deck": s_Deck }
              BI["Span"][b]["G"][g]["Nth"][10] = { "Sta": e_val.distance, "offset": e_val.offset,  "ProG": e_ProG, "xs_corr": e_xs_correction, "el_Deck": e_Deck }


              #print( f't["Beam_span"][{b}]["G"][{g}]\t{s_val.distance:.2f}\t{s_ProG:.2f}\t{s_xs_correction:.2f}\t{s_Deck:.2f}  |\t{e_val.distance:.2f}\t{e_ProG:.2f}\t{e_xs_correction:.2f}\t{e_Deck:.2f}')

              print( "\t".join( res ) )
          else:
              continue
with open( "test_14_beam_info.yaml", 'w') as fh:
    yaml.dump( BI, fh,  sort_keys=False )

#print(dir(xs)) 
print( BI )
IPython.embed()
