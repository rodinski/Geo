import Geo
from pprint import pprint 

def pi_from_d( in_tuple: tuple):
    outList = []
    outList.append( (0,0) )
    for t in in_tuple:
        x =outList[-1][0] + t[0]
        y =outList[-1][1] + t[1]
        outList.append( (x,y) )
    return outList


deltas = [ (-1000,   0),
           (   50, 100),
           ( -400, -100),
           ( -600,    0),
           (-200, 800),
           (-100,  50),
           (-300, 500),
           (500, 150),
           (600, -500),
           (600, -400),
           (200,   25),
           (200, -150),
           (1000,   0), 
           (   0, -300) 
           ]

a = pi_from_d(deltas)
print(f"{a=}")
pprint(a)

'''
points = [ Point(0    ,     0  ), 
           Point(-5000,     0  )
           Point(-4950,   300  )
           Point(-6000,     0  )
           Point(-9000,     0  )
           Point(-10000, 3000  )
           Point(-10500, 3050  )
           Point(-10500, 3250  )
           Point(-5000,     0  )
           Point(-5000,     0  )
           Point(-5000,     0  )
           Point(-5000,     0  )]
'''


#mychain = Chain("Monza", s1)
