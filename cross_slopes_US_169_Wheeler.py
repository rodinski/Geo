# setup a data structure of ordered stations
# with os_el as an ordered list of  ( offset , delta_el ) pairs
# offsets are listed from left to right


lft =  -9.333       # left edge of deck FROM BREAK
rht =  33.333       # right edge of deck FROM BREAK

edge_lt_pos_6= round ( lft * 0.06 , 5)  # lef edge at +6% slope   left_turn_curve
edge_rt_pos_6= round ( rht * 0.06 , 5) # right edge at +6% slope

edge_lt_0= round ( lft * -0.00 , 5)  # left edge at 0.0% slope    zero
edge_rt_0= round ( rht * -0.00 , 5) # right edge at 0.0% slope

edge_lt_neg_6= round ( lft * -0.06 , 5)  # lef edge at -6% slope   right_turn_curve
edge_rt_neg_6= round ( rht * -0.06 , 5) # right edge at -6% slope

edge_lt_neg_2= round ( lft * -0.02 , 5)  # lef edge at -6% slope   right_turn_curve
edge_rt_neg_2= round ( rht * -0.02 , 5) # right edge at -6% slope


#os_el is a list of (off_set_from_CL, delta_el_from_PGL) 

s = [ 
      {"sta":17103.63, "os_el": [ (-9.3333 - 12.0, edge_lt_pos_6 ), (0, 12 * 0.06), ( 21.33, edge_rt_pos_6 ) ] },
      {"sta":17423.32, "os_el": [ (-9.3333 - 12.0, edge_lt_pos_6 ), (0, 12 * 0.06), ( 21.33, edge_rt_pos_6 ) ] },
      {"sta":17639.32, "os_el": [ (-9.3333 - 12.0, 0.0           ), (0,         0), ( 21.33, 0.0) ] },  
      {"sta":17855.32, "os_el": [ (-9.3333 - 12.0, edge_lt_neg_6 ), (0, 12 * -0.06), ( 21.33, edge_rt_neg_6 ) ] },
      {"sta":18186.77, "os_el": [ (-9.3333 - 12.0, edge_lt_neg_6 ), (0, 12 * -0.06), ( 21.33, edge_rt_neg_6 ) ] },

      {"sta":18330.77, "os_el": [ (-9.3333 - 12.0, edge_lt_neg_2 ), (0, 12 * -0.02), ( 21.33, edge_rt_neg_2 ) ] },
      {"sta":18533.65, "os_el": [ (-9.3333 - 12.0, edge_lt_neg_2 ), (0, 12 * -0.02), ( 21.33, edge_rt_neg_2 ) ] },

      ]
def Average(lst):
    return round( sum(lst) / len(lst),5)

def  slopes_from_sta( station ):
    if station < s[0]["sta"]:
        return "sta out of range too small"
        
    if station > s[-1]["sta"]:
        return "sta out of range too large"

    for (bk,ah) in zip( s[:-1],s[1:] ):
        #print (bk["sta"],  ah["sta"])
        if bk["sta"] <= station and station  <= ah["sta"]:
            L_sta = ah["sta"] - bk["sta"] 
            x_sta = station - bk["sta"]
            pct = x_sta / L_sta
            
            #print(f"\nStation: {station}    pct: {pct} ")
            #print( bk ) 
            #print( "" )
            #print( ah ) 
            os_el_returnList = list()
            for (os_el_bk, os_el_ah ) in zip(bk["os_el"], ah["os_el"]):
                #for all of the listed breaks do both a linear interperlation of
                # both offset and el   (offsets should have no change)
                offset_start = os_el_bk[0]
                offset_end = os_el_ah[0]
                offset_delta = offset_end - offset_start
                offset_return = offset_start + offset_delta * pct
                #
                el_start = os_el_bk[1]
                el_end = os_el_ah[1]
                el_delta = el_end - el_start
                el_return = round( el_start + el_delta * pct, 4)

                os_el_returnList.append( ( offset_return, el_return) )
            # return a list of offset, el pairs    
            return os_el_returnList 
    return False    #something went wrong

def cross_slope_from_breaks( offset, breakList ):
    if offset < breakList[0][0]:
        return "offset out of range, too small"
    if offset > breakList[-1][0]:
        return "offset out of range, too large"

    for (lft_break, rht_break) in zip( breakList[:-1], breakList[1:]):
        os_start = lft_break[0]
        el_start = lft_break[1]
        os_end = rht_break[0]
        el_end = rht_break[1]
        os_L = os_end - os_start
        el_delta = el_end - el_start
        if os_start <= offset and offset <= os_end:
            os_percent = (offset - os_start) / os_L
            el_return = el_start + el_delta * os_percent
            # return a single elevation
            return round(el_return,4)
    return False  #something went wrong

def xs_correction( my_sta, my_offset):
    tuple_of_breaks = slopes_from_sta( my_sta)
    return cross_slope_from_breaks( my_offset , tuple_of_breaks)
    


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    x = list()
    y1 = list()
    y2 = list()
    y3 = list()
    y4 = list()
    y5 = list()


    # step thru all the available stations 10 at a time plot the following offsets.
    for my_sta in range( int(s[0]["sta"]+1), int(s[-1]["sta"]-1), 10):
        x.append(  my_sta)
        y1.append( xs_correction(my_sta, -20.00 ) )
        y2.append( xs_correction(my_sta, -12.00 ) )
        y3.append( xs_correction(my_sta,  -4.00 ) )
        y4.append( xs_correction(my_sta,   0.00 ) )
        y5.append( xs_correction(my_sta,  12.00  ) )

    plt.plot(x,y1)
    plt.plot(x,y2)
    plt.plot(x,y3)
    plt.plot(x,y4)
    plt.plot(x,y5)
    plt.show()
