from lark import  Lark, Transformer, Tree
import random_survey_book as rsb 
import random
mygrammer = open( "store_curve.lark", 'r' )


l = Lark(mygrammer, debug=True, parser='earley', )

def rand_case(instr:str) -> str:
    ret_str = ''
    for i in instr:
        i = i.upper() 
        if random.choice([True,False]):
            i = i.lower()
        ret_str += i
    return ret_str

def main( ):
    instr = [ 
     "STORE CURVE CURVE_01 PC 1002 DB 1001 TO 1002 DEGREE 3 45 DA 20 38 23.7 STATION PC 20291.93"
    ,"STORE CURVE CURVE_02 PI 1002 DB 12.3  DEGREE 3 45 DA 20 38 23.7 STATION PC 20291.93"
    ,"STO CUR STER4 PC N 107.43 E 26.01 DB 4.1432 DEGREE 3.83 DEL -62.87 STA PC 1866.0"
    ,"STO CUR STER5 PC N 107.57 E 85.46 DB 01 02 34.56 DEGREE 7.84 DEL 36.661 STA PC 3769.8"
    ]

    mystring="store | cur | abc123 | PB n 1 e 2 DB 123.45 TLl 80.8 | DEGREE  15 44 34.64 | M DEF 123.4"
    print( rand_case(mystring) )
    myTree = l.parse(mystring)
    print( myTree.pretty(), type(myTree) )

    #print(type(myTree), type(myTree.iter_subtrees_topdown() ) )
    for item in myTree.iter_subtrees_topdown():
        print(type(item)) #, item) 

    if False:
        for n,mystring in enumerate(range(10)):
            mySC = rsb.Rand_curve()
            #print( mySC )
            print( l.parse( str(mySC) ).pretty() )

        


if __name__ == "__main__":
    main()


