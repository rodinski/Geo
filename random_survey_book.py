from random import  randrange, sample

def r_from( obj, k=1 ) ->str:
    return   str(sample(list( obj ), k )[0] )


digit = lambda : samp(list("0123456789"))
r_float = lambda: str(randrange(0,10000)/ 100) 

def dms()->str:
    ''' 
    '''
    d = " " + r_from("- 012345") + digit()
    m = " " + r_from(" 012346")  + digit()
    s = " " +str( randrange(0,10000)/100 )
    return d+m+s

def float_or_dms()->str:
    myList = [ r_float(), dms() ]
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

def ne()->str:     #for now don't use negatives
    return f"N {randrange(0, 10000)/100} E {randrange(0, 10000)/100}"

def location()->str:
    valid_pt = [  "Pt1", "23", "Point_4", "PT_5", ne()  ]
    return sample(valid_pt, 1)[0]

def btr()->str:
    btrList = [  " ".join(["PC", location(), "DB", bear() ])
                ," ".join(["PI", location(), "DB", bear() ])
                ," ".join(["PB", location(), "DB", bear(), "TL", r_float() ] )]
    return sample(btrList,1)[0]

def element()->str:
    elementList= [ 
        " ".join([ samp( ["RAD", "radius", "T", "tangent", "L"]), r_float()] )
        ," ".join([ samp( ["DEG", "DEGREE"]), dms()] ) 
        ," ".join([ "POC",  pt_id()  ] ) 
        ]
    return sample(elementList,1)[0]


def atr()->str:
     atrList = [ " ".join(["DA " + bear()  ] )
          ,"PA " +  pt_id()            # more logic need here!
          ," ".join([r_from("MP ") +samp([" DEF "," DEFLECTION "]) + float_or_dms() ] )
         ]
     return sample(atrList,1)[0]

class Rand_curve:


    def __init__(self):
        s = ["store", "sto"]   #"s" is valid for storing points not curves
        cur = ["curve", "cur" ]
        n_id = ["abc",  "abc123",  ]

        self.item = [   samp(s) 
                      , samp(cur) 
                      , samp(n_id)
                      , btr()
                      , element()
                      , atr()
                    ]
                
    def __str__(self)->str:
        return " | ".join(self.item)

if __name__ == "__main__":
    print("Example stor curve:")
    for i in range(20):
        myRand_curve = Rand_curve()
        print( myRand_curve)
    for i in range(400):
        pass
        print( bear(), atr() )
        #print( "%8s %20s %50s  %26s " % (pt_id(),location(), btr(),element()))
