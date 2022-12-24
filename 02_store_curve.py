from lark import Lark
import random_survey_book as rsb
myGrammer = open( "store_curve.lark", 'r' )

l = Lark(myGrammer, debug=True, parser='earley')

def main( ):
    instr = [ 
     "STORE CURVE CURVE_01 PC 1002 DB 1001 TO 1002 DEGREE 3 45 "#DA 20 38 23.7 STATION PC 20291.93"
    ,"STORE CURVE CURVE_02 PI 1002 DB 12.3  DEGREE 3 45 "#DA 20 38 23.7 STATION PC 20291.93"
    ,"STO CUR STER4 PC N 107.43 E 26.01 DB 4.1432 DEGREE 3.83 "#DEL -62.87 STA PC 1866.0"
    ,"STO CUR STER5 PC N 107.57 E 85.46 DB 01 02 34.56 DEGREE 7.84 "#DEL 36.661 STA PC 3769.8"
    ]

    test ="store | cur | abc123 | PB 4 DB 123.45 TLl 80.8 | DEGREE  15 44 34.64 | M DEF 123.4"
    print(test)
    print( l.parse(test).pretty())
    if True:
        for item in range(100):
            mySC = rsb.Rand_curve()
            print( mySC )
            print( l.parse( str(mySC) ).pretty() )
            


if __name__ == "__main__":
    main()

"""
ABOUT        ADD           ADJUSTMENT    AH          AHEAD        ALIGNMENT
ALIGNMENTS   ALL           AND           ANG         ANGLE        APPEND
ARC          AREA          AS            ASEC        AT           AZIMUTH
BACK         BEARING       BEG           BEGIN       BK           BLD
BY           CATALOG       CC            CCL         CEAL         CELL
CENTER       CHA           CHAIN         CHORD       CIRCLE       CL
CLOSURE      COMMANDS      COMPASS       CON         CONCENTRIC   COORDINATE
COPY         CS            CURVE         CURVES      A            DATA
DB           DEF           DEFINE        DEFLECTION  DEGREE       DELETE
DELTA        DESCRIBE      DESCRIPTION   DIRECTION   DISTANCE     DIV
DRIVE        DTM           E             EAS         EASEMENT     EDIT
ELEMENT      ELEV          ELEVATION     END         EQN          EQUATE
EVEN         EXIT          EXTERIOR      EXTERNAL    FEATURE      FILE
FIX          FOR           FROM          G           GEOMETRY     GEOPAK
GRADE        H             HP            IMPROVEMEN  IN           INC
INCLUDE      INCOMPLETE    INCREMENT     INPUT       INSERT       INTERIOR
INTERSECT    INTERSECTION  INVERSE       IS          JOB          K
L            LAYOUT        LC            LEAST       LEGAL        LENGTH
LESS         LINE          LIST          LOAD        LOCATE       LP
LS           LT            M             MACROS      MAKE         MEASUREMENT
MODIFY       MUL           N             NE          NEAR         NEXT
NUMBER       OFF           OFFSET        ON          OPEN         OPERATOR
OUTPUT       OWNER         P             PA          PARALLEL     PARCEL
PATH         PB            PC            PI          PLUS         POC
POINT        POINTS        POT           PRINT       PRINTER      PROFILE
PROFILES     PROJECT       PT            R           RADIUS       READ
RECOVERY     REFERENCE     REGION        RESTORE     ROTATE       ROTATION
RT           SAVE          SC            SCALE       SCS          SECTION
SEGMENT      SEL           SET           SL          SPI          SPIRAL
ST           STA           STATION       STEP        STORE        STUFF
SUBDIVIDE    SUBJECT       SURVEY        T           TAKE         TAN
TANGENT      TD            THRU          TL          TO           TRANSFORM
TRANSLATION  TRANSPOSE     TRAVERSE      TS          TTL          TYPE
USING        VC            VPI           X           XY           Y
"""

"""
$ CHAIN ELM_BRANCH
SET FEATURE ADD 2002
STORE POINT 1001 N 0000.00   E 0.000 STA 20000.00
STORE POINT 1002 N  291.93   E 0.000 STA 20291.93
STORE LINE LINE01 1001 TO 1002
STORE CURVE CURVE_01 PC 1002 DB 1001 TO 1002 DEGREE   3 45 DA 20 38 23.7 STATION PC 20291.93
STORE CHAIN ELM_BRANCH 1001 1002 CUR CURVE_01  AS IS
$
$ FIND THE CENTER "CC" OF CURVE 01 AND CALL IT CC_01
EQUATE CC01 TO CC CURV01
SET FEATURE ADD 3002
LOCATE BENT20 ON CHAIN ELM_BRANCH STA 20203.00
$STORE SOME VARIABLE LENGTHS AND RECALL THEM WITH %NAME
STORE DISTANCE SPAN2 60
LOCATE BENT30 ON CURVE CURVE_01 %SPAN2 FROM BENT20
STORE ANGLE BSKEW 315 00 00
LOCATE ANGLE 1001 BENT20 2001 %BSKEW 100
"""
