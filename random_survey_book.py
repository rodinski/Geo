""" Constuct random inputs that basically conform 
to geopak coco commands will use this output to 
test a parser.  Certain deviations from strict 
MX cogo sintax have been allowed the hope is that 
any legacy command will be understood yet new 
naming with "_" and multiple alpha digit 
combinations are allowed
"""
from random import  randrange, sample, uniform
from math import modf

def r_from( obj, k=1 ) ->str:
    return   str(sample(list( obj ), k )[0] )
def digit()->str: 
    return samp(list("0123456789"))

def r_float(min=-85,max=85)->str: 
    return "%.2f" % uniform(min,max) 

def dms(min=-85,max=85)->str:
    ''' 
    '''
    angle = uniform(min,max)
    mf = modf(angle)
    d = mf[1]
    m = int( abs(mf[0])  * 60)
    s = ( abs(mf[0]) * m/60  ) 

    #print(angle, '\t', end='')
    d=  " %d" % d
    m = " %d" % m 
    s = " %.2f" % s
    
    #randomly drop some of the numbers
    one_in = 3
    if randrange(one_in)<1:
        s = s[:randrange(len(s))] 
    return d+m+s

def float_or_dms(min=-85,max=85)->str:
    myList = [ r_float(min,max), dms(min, max) ]
    return sample(myList,1)[0]

def samp( inList )->str:
    return sample( inList, 1)[0]

#recrate this each time otherwis r_from doesn't get rerun
def pt_id()->str:
      myList = [ r_from('123'),  "P"+r_from('456'), "p_"+r_from('789') ]
      return sample(myList, 1)[0] 

def bear()->str:
      myList = [ r_float(), dms(), r_from('NSns') + dms()+ ' ' + r_from('EWew') ]
      return sample(myList,1)[0]

def ne(min=0,max=10000, dec=2)->str:     #for now don't use negatives
    
    return f"N {uniform(min,max):.{dec}f} E {uniform(min,max):.{dec}f}"

def location()->str:
    valid_pt = ["Pt1", "23", "Point_4", "PT_5", ne()]
    return sample(valid_pt, 1)[0]

def btr()->str:
    btrList = [  
             " ".join(["PC", location(), "DB", bear() ]),
             " ".join(["PI", location(), "DB", bear() ]),
             " ".join(["PB", location(), "DB", bear(), "TL", r_float(0,999) ] )
              ]
    return sample(btrList,1)[0]


def element()->str:
    elementList= [ 
        " ".join([samp(["RAD", "radius", "T", "tangent", "L"]), r_float(0,1000)] ),
        " ".join([samp(["DEG", "DEGREE"]), dms()] ),
        " ".join(["POC", pt_id()] ),
        ]
    return sample(elementList,1)[0]

def atr()->str:
     atrList = [ 
          " ".join(["DA " + bear()  ] ),
          "PA " + pt_id(),           # more logic need here!
          " ".join([r_from("MP ") 
              + samp([" DEF "," DEFLECTION "]) 
              + float_or_dms() ] )
          ]
     return sample(atrList,1)[0]

class Rand_curve:
    """Each of these objects should be readabel this is for testing"""

    def __init__(self):
        s = ["store", "sto"]   #"s" is valid for storing points not curves
        cur = ["curve", "cur" ]
        n_id = ["abc",  "abc123",  ]
        self.item = [  
                       samp(s),
                       samp(cur),
                       samp(n_id),
                       btr(),
                       element(),
                       atr(),
                    ]
                
    def __str__(self)->str:
        return " | ".join(self.item)


def main():
    """only ran for testing"""
    print("Example stor curve:")
    for i in range(20):
        myRand_curve = Rand_curve()
        print( myRand_curve)
    for i in range(400):
        print( dms(-2, 2),  )
        #print( "%8s %20s %50s  %26s %26s" % (pt_id(),location(), btr(),element(),atr() ))
        pass

if __name__ == "__main__":
    main()
