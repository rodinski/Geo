import unittest
from unittest import TestCase
from geo_ import Point, Bearing, Angle
import math
import cmath
from math import pi, isclose

import random

class test_Point(TestCase):

    def test_Point_from_complex_x(self):
        """
        Test that Points are made from complex numbers
        """
        print("")
        p_list = [ (3, 4), (4, 3), (-3, -4), (-4,-3) ]
        n = 0
        for x, y in p_list:
            n = n + 1
            with self.subTest(f"test_Point_from_complex_x{n}"):
                self.assertEqual(Point(x, y).X, float(x) )
            with self.subTest(f"test_Point_from_complex_y{n}"):
                self.assertEqual(Point(x, y).Y, float(y) )

    def test_Point_phase(self):
        """
        Test Point.phase() 
        """
        print("")
        p_list = [ (9, 0), (0, 9), (-9, 0),
                   (0,-9), (9, 9), (3, 3),
                    (3, 4),(-3, 4),
                  ]
        n = 0
        for x, y in p_list:
            n = n + 1
            with self.subTest(f"Test Point.phase() {n}"):
                self.assertEqual(Point(x,y).phase(), Bearing(math.atan2(y,x)) )
        for n in range(n+1, n+1000):
            x = random.random() *2 -1
            y = random.random() *2 -1
            with self.subTest(f"Test Point.phase() {n}"):
                self.assertEqual(Point(x,y).phase(), Bearing(math.atan2(y,x)) )
                self.assertEqual(Point(x,y), Point.from_complex( complex(x,y)) )

class test_Direction(TestCase):
    def test_Angle(self):
        """
        Test Angle(x, unit='deg' ) 
        """
        print("")
        for n in range(5):
            rad = random.random()*2*pi - pi
            with self.subTest(f"Test Angle()_{n}"):
                self.assertEqual(Angle(rad), rad)
                # see if wokring with deg gives same answer
                deg = math.degrees(rad)
                self.assertTrue(math.isclose(Angle(deg, unit='deg'), Angle(rad)))


    def test_Bearing(self):
        """
        Test Bearing(x, unit='deg' ) 
        """
        print("")
        deg_list = []
        for n in range(18):
            deg = n*45 - 360.
            deg_list.append(deg)
            rad = math.radians(deg)

            _a = Bearing(rad)
            _b = Bearing( deg, unit='deg')
            with self.subTest(f"Test Bearing()_{n}"):
                #all bearing are normalized (-π, +π) so we do same to the rad
                self.assertTrue(math.isclose( _a,_b,))
               
                a = _b.rt_90()
                b = Bearing(cmath.phase(cmath.rect(1.000, rad - math.pi/2.0)))
                #print(n, a, b,  a-b)
                self.assertTrue(math.isclose( a,b, rel_tol=0.000001, abs_tol=0.000001))
                    

                self.assertTrue(math.isclose(
                    Bearing(deg, unit='deg').opposite(), 
                    Bearing(cmath.phase(cmath.rect(1.00, rad + math.pi))),
                    rel_tol=0.000001, abs_tol=0.000001))

                #print(f"{deg.opposite()=:.4f}\t{rad=:.4f}\t\t{_a=}\t{_b=}\t{_a-_b}")
                _a = Bearing(deg, unit='deg').opposite()
                _b = Bearing(cmath.phase(cmath.rect(1.00, rad + math.pi)))
                _c = Bearing(cmath.phase(cmath.rect(1.00,  rad - math.pi)))
                #print(f"{n}\t{_a=:.6f}\t{_b:.6f} ")
                self.assertTrue(math.isclose(_a, _b,
                    rel_tol=0.000001, abs_tol=0.000001))
                self.assertTrue(math.isclose(_a, _c,
                    rel_tol=0.000001, abs_tol=0.000001))


    def test_Angle(self):
        """
        Test Angle(x, unit='deg' ) 
        """
        print("")
        for n in range(18):
            deg = n*45 - 360.
            rad = math.radians(deg)
            _bearing = random.random() * 2*pi - pi

            a_rad = Angle(rad)
            a_deg = Angle( deg, unit='deg')

            self.assertEqual( a_rad, a_deg)
            self.assertEqual( a_rad + pi, a_deg + pi)
            self.assertEqual( a_rad + 2*pi, a_deg + 2 * pi)
            self.assertEqual( a_rad + 2*pi, a_deg + 2 * pi)
            self.assertEqual( _bearing + a_rad , _bearing +  a_deg )


if __name__ == '__main__':
    unittest.main(verbosity=2)
