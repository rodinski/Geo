from Geo import *
import matplotlib.pyplot as plt
from   matplotlib.path import Path
import matplotlib.patches as patches 
import math
import line_intersect

fig, ax = plt.subplots()

R = 980
Delta = math.radians(55 + ( 35 + 59.2/60 ) / 60 )
print(Delta)

pc_sta = 2196.18
pt_sta = 3147.17
arc_len = pt_sta - pc_sta

pt_pc =  Point(0,0)
pt_cc =  Point(0,R)
print(pt_pc,pt_cc)

curve = Curve( pt_pc, pt_cc, Delta )

my_chain = Chain( "BL", curve, StartSta=2196.18 )
my_chain.forward(100)
    
print( pc_sta, pt_sta, arc_len ) 
print(my_chain) 
b1 =  my_chain.move_to(2951.25)
b2 =  my_chain.move_to(3026.00)
b3 =  my_chain.move_to(3106.00)
b4 =  my_chain.move_to(3179.00)

b1_brg = Segment(b1, curve.CC).Bearing()
s1_C_brg = Bearing( b1_brg - math.radians(87 +(45 + 50/60)/60))

b1_C_point = Ray( b1       , b1_brg).move_to( 4 + (0 + 3/16)/16  )
b1_A_point = Ray(b1_C_point, b1_brg).move_to( 9.5 +11  )
b1_B_point = Ray(b1_C_point, b1_brg).move_to( 9.5  )
b1_D_point = Ray(b1_C_point, b1_brg).move_to( -9.5  )
b1_F_point = Ray(b1_C_point, b1_brg).move_to( -9.5 -9.5 )

print(f"    {b1_brg=}  ", type(b1_brg) )
print(f"  {b1_C_point=}  ", type(b1_C_point) )

b1_left = Ray(b1, b1_brg).move_to(30, 21/12)
b1_right = Ray(b1, b1_brg).move_to(-30, 21/12)
b1_brg_line = Segment(b1_left, b1_right)
ax.add_patch(b1_brg_line.patch())

b1_C = Ray(b1, b1_brg).move_to( 4+ (0 + 3/16)/12 )

b2_brg =  Segment(b2, curve.CC).Bearing()
b2_bk_left  = Ray(b2, b2_brg).move_to( 25, -10/12)
b2_bk_right = Ray(b2, b2_brg).move_to(-20, -10/12)
b2_ah_left  = Ray(b2, b2_brg).move_to( 25, +10/12)
b2_ah_right = Ray(b2, b2_brg).move_to(-20, +10/12)
b2_bk_line = Segment(b2_bk_left, b2_bk_right)
b2_ah_line = Segment(b2_ah_left, b2_ah_right)
ax.add_patch(b2_bk_line.patch())
ax.add_patch(b2_ah_line.patch())

b2_C_point = Ray(b2, b2_brg).move_to(4 +1/12 )
b2_C_point_2 = line_intersect.intersect_lines( Ray(b1_C, s1_C_brg), Ray(b2, b2_brg)  )
print( abs(b2_C_point_2 - b2_C_point) )
b2_A_point = Point.from_complex(b2_C_point + cmath.rect(  9.5+ 11, b2_brg ) )
b2_B_point = Point.from_complex(b2_C_point + cmath.rect(  9.5, b2_brg ) )
b2_D_point = Point.from_complex(b2_C_point + cmath.rect( -9.5, b2_brg ) )
b2_F_point = Point.from_complex(b2_C_point + cmath.rect( -19 , b2_brg ) )

s1_c_seg = Segment(b1_C, b2_C_point)
ax.add_patch( s1_c_seg.patch() )

span_1 = Segment(b1, b2)
span_2 = Segment(b2, b3)
span_3 = Segment(b3, b4)
ax.add_patch( span_1.patch() )

b3_brg =  Segment(b3, curve.CC).Bearing()
b3_A_point = Ray(b3, b3_brg).move_to(4 +1/12 +9.5 +11)
b3_B_point = Ray(b3, b3_brg).move_to(4 +1/12 +9.5    )
b3_C_point = Ray(b3, b3_brg).move_to(4 +1/12         )
b3_D_point = Point.from_complex(b3_C_point + cmath.rect( -9.5, b3_brg ) )
b3_F_point = Point.from_complex(b3_D_point + cmath.rect( -9.5, b3_brg ) )

b4_brg =  Segment(b4, curve.CC).Bearing()
b4_A_point  = Ray(b4, b4_brg).move_to(4 +1/12 +9.5 +11, -2.5)
b4_B_point  = Ray(b4, b4_brg).move_to(4 +1/12 +9.5    , -2.5)
b4_C_point  = Ray(b4, b4_brg).move_to(4 +1/12         , -2.5)
b4_D_point = Point.from_complex(b4_C_point + cmath.rect( -9.5, b4_brg ) )
b4_F_point = Point.from_complex(b4_D_point + cmath.rect( -9.5, b4_brg ) )


print(b2_brg, type(b2_brg) )
point_list =  [ b1, b2, b3, b4, b1_left, b1_right, ] 
point_list += [ b1_C, ] 
point_list += [ b1_C_point, b1_B_point, b1_A_point,  b1_D_point, b1_F_point, ]
point_list += [ b2_C_point_2, b2_C_point, b2_B_point, b2_A_point,  b2_D_point, b2_F_point, ]
point_list += [ b3_C_point, b3_B_point, b3_A_point,  b3_D_point, b3_F_point, ] 
point_list += [ b4_C_point, b4_B_point, b4_A_point,  b4_D_point, b4_F_point, ] 

xs = [ i.X for i in point_list]
ys = [ i.Y for i in point_list]
ax.scatter(xs, ys )


#b3_bk_left = Ray(b3, b3_brg).move_to(100,-10/12)
#b3_bk_right = Ray(b3, b3_brg).move_to(-100,-10/12)
#b3_ah_left = Ray(b3, b3_brg).move_to(100,+10/12)
#b3_ah_right = Ray(b3, b3_brg).move_to(-100,+10/12)
#b3_bk_line = Segment(b3_bk_left, b3_bk_right)
#b3_ah_line = Segment(b3_ah_left, b3_ah_right)
#ax.add_patch(b3_bk_line.patch())
#ax.add_patch(b3_ah_line.patch())

#s3_C_brg = Bearing( b4_brg + math.radians(87 +(45 + 50/60)/60))
#b3_C_point = line_intersect.intersect_lines( Ray(b1_C,s3_C_brg), Ray(b2, b2_brg)  )
#b3_B_point = Point.from_complex(b2_C_point + cmath.rect( 9.5, b2_brg ) )
#b3_A_point = Point.from_complex(b2_B_point + cmath.rect( 11, b2_brg ) )
#b3_D_point = Point.from_complex(b2_C_point + cmath.rect( -9.5, b2_brg ) )
#b3_F_point = Point.from_complex(b2_C_point + cmath.rect( -19 , b2_brg ) )






for p in my_chain.patch_list():
    ax.add_patch( p )
plt.axis('scaled')
plt.show()
