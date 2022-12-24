"""
This is a Docstring
"""
from math import pi, radians, fmod
from  angles import normalize


class Bearing(float):
    """
    Bearing as a float normally from 0 to pi
    with the ability to be a Delta(Angle)  bewteen -pi and pi
    """
    
    def __init__(self, brg):
        self.val = float( normalize(brg) )
      
    def Opposite(self):
        return normalize( self.val + pi )
    
    def as_Delta(self):
        return normalize( self.val, -pi, pi )

    def Delta_to(self, other=None ):
        """ Bearing.Delta_to( x )  = bearing.self + x
        """
        if other == None:
            other = self.val
        return normalize( other - self.val, -pi, pi )

    def __str__(self):
        return f"Bearing={self.val}"
    
if __name__ == "__main__":

    start = Bearing(radians(45))
    end   = Bearing(radians(315))
    print( start)
    print( f"{    end=}")
    print( Bearing(end-start).as_Delta() )
    print( start.Delta_to() )
