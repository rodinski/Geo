import yaml
import IPython

filename= "test_14_beam_info.yaml"
with open( filename, 'r') as fh:
     BI = yaml.safe_load(fh) 



#print( BI )



def get_skips( span:int, girder:int ) -> int:
    '''return the coorect number of blank rows for a given span girder combination
    default to "\n"'''
    skips = { "S4_G1": 2, "S4_G2": 1, "S4_G3": 2, "S4_G4": 2, "S4_G5": 0, 
              "S6_G1": 2, "S6_G2": 4, "S6_G3": 4, "S6_G4": 4, "S6_G5": 0, 
              "S8_G1": 2, "S8_G2": 4, "S8_G3": 4, "S8_G4": 4, "S8_G5": 0,  }

    if span in { 4, 6, 8}:
        mykey = "S" + str(span) + "_G" + str(girder)
        return  skips[mykey]
    else:
        return 0

def loop_match_excel_shape( unit ):
    '''minic the rows that amk has used in her three tabs
    unit, s, and g are all int '''
    
    for g in [1, 2, 3, 4, 5]:
     
        if unit not in { 1, 2, 3 }:
            raise ValueError( f"unit value of {unit}, is out of range" )
        if unit == 1:
            span_array = [ 1, 2, 3, 4 ]
        if unit == 2:
            span_array = [ 5, 6]
        if unit == 3:
            span_array = [ 7, 8]

        for s in span_array:

              for Nth in [0, 10]:
                  ref = BI["Span"][s]["G"][g]["Nth"][Nth]

                  tag = f"Span_{s}.G_{g}.Nth_{Nth}"

                  print(f"{tag}\t{ref['Sta']}\t{ref['offset']}")   # \t{ref['ProG']}
                  if Nth == 0:
                        print( "\n"*9, end = None )  
              print( "\n" * get_skips( s, g), end = None)


def loop_for_equations( unit ):
    base_tab_name = "Haunch Calcs Unit "   
    tab_ref = f"='{base_tab_name}{unit}'!"   #to look like excel equation for a general tab
    myrow = 25

    for g in [1, 2, 3, 4, 5]:
     
        if unit not in { 1, 2, 3 }:
           raise ValueError( f"unit value of {unit}, is out of range" )
        if unit == 1:
           span_array = [ 1, 2, 3, 4 ]
        if unit == 2:
           span_array = [ 5, 6]
        if unit == 3:
           span_array = [ 7, 8]

        for s in span_array:
            for Nth in [0, 10]:
               ref = BI["Span"][s]["G"][g]["Nth"][Nth]

               tag = f"Span_{s}.G_{g}.Nth_{Nth}"

               print(f"{tag}\t{tab_ref}A{myrow}\t{tab_ref}F{myrow}")   # \t{ref['ProG']}
               myrow += 1                 # always increament 1
               if Nth == 0:
                        myrow += 9        # increment to get to the next bearin
               myrow += get_skips( s, g)   # amk might have had inconsistent gaps b/w girders
    pass


if __name__ == "__main__":
    loop_match_excel_shape( 1 )
    loop_match_excel_shape( 2 )
    loop_match_excel_shape( 3 )
    loop_for_equations( 3 )



    
#IPython.embed()

