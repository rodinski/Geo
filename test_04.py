from Geo import *
import matplotlib.pyplot as plt
from   matplotlib.path import Path
import matplotlib.patches as patches 


fig, ax = plt.subplots()

PIs = [ Point(0,0), Point(5,0), Point(10,3), Point(15,4), \
        Point(20,0), Point(25,5), Point(30,20) ] 


myseg = Segment( Point(10,10), Point(20,12))
myCurve = Curve( myseg.Pt2, Point(18,40), -0.2)

mychain = Chain( "ChaneName", myseg ) 
mychain.addRoute(myCurve) 

print(mychain.__str__() )
