from angles import Angle, normalize 
import angles 
import math
from math import pi

def AsBearing( obj:float)-> float:
    return normalize( obj, 0.0, 2*pi)

def AsDelta( obj:float)-> float:
    return normalize( obj, -pi, pi)

if __name__ == "__main__" :
   from random import uniform

   low  =   0
   high =   2* pi

   for i in range(10):
       a_rad = round( uniform(low, high),2 )
       b_rad = round( uniform(low, high),2 )
       c =  a_rad - b_rad
       c_Delta =  round(AsDelta(a_rad - b_rad), 2) 

       fstr = f"{i=} {AsBearing(a_rad):.2f}  {AsBearing(b_rad):.2f}  {c_Delta:6.2f} "
       print( fstr )
