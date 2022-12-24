import cmath

class Point(complex):
    ''' At new Point called with  Point(x,y) or as Point.from_complex( complex_number )
        internally the number is stored as a complex base class'''
    is_Point = True

    def __init__(self, x:float, y:float,  *args, **kwargs)->float:
        self.val = complex(x, y)
        self.X = self.val.real
        self.Y = self.val.imag


    @classmethod
    def from_complex( cls, complex_numb:complex ):
        X = complex_numb.real
        Y = complex_numb.imag
        return cls(X,Y)

    def phase(self):
        return cmath.phase(self)

    def __str__(self):
        return(f"Point x={self.val.real}, y={self.val.imag}")




if __name__ == "__main__" :
   from random import uniform
   pts = []
   for i in range(10):
       x = uniform(0,5)
       y = uniform(0,5)
       pts.append( Point(x,y) ) 
       print(pts[-1],  pts[-1].__repr__  )


