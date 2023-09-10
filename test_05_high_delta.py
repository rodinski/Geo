import Geo
import matplotlib.pyplot as plt
from   matplotlib.path import Path
import matplotlib.patches as patches 

fig, ax = plt.subplots()

pc =  Geo.Point(0,0)
cc =  Geo.Point(-0.3,0.8)
delta = Geo.Angle( 170, unit='deg')
print(delta)
curves = []

#curves.append( Geo.Curve(pc,cc,delta) )
curves.append( Geo.Curve(pc,cc, Geo.Angle(160, unit='deg')))
curves.append( Geo.Curve(pc,cc, Geo.Angle(220, unit='deg')))
curves.append( Geo.Curve(pc,cc, Geo.Angle(300, unit='deg')))

for curve in curves:
    print("\n",curve)
    for r in curve.patch_all():
        ax.add_patch( r ) 
plt.axis('scaled')
plt.show()

