import unittest
from unittest import TestCase
from geo_ import Point, Bearing, Angle, Ray, xy
import math
import cmath
from math import pi, isclose
import matplotlib.pyplot as plt
from matplotlib import patches

import random

class test_01_Point(TestCase):

    def test_Ray(self):
        """
        Test Rays
        """
        print("")
        pts = []
        for n in range(18):
            deg = n*45 - 360.
            rad = math.radians(deg)

            _a = Bearing(rad)
            _b = Bearing( deg, unit='deg')
            pt_0 = Point(0, 0)
            with self.subTest(f"Test Ray()_{n}"):
                R = random.random() * 10 + 100
                _x = math.cos(rad)
                _y = math.sin(rad)
                pt_1 = Point(R*_x, R*_y)
                pts.append(pt_1)
                ray = Ray( pt_0, _a )
                pt_2 = ray.set_point( R, 0)
                self.assertTrue(math.isclose(pt_1.X, pt_2.X) and math.isclose(pt_1.Y, pt_2.Y))
                self.assertTrue(pt_1.equal(pt_2) and pt_2.equal(pt_1))
        
        fig, ax = plt.subplots()
        #not plot vs scatter
        ax.scatter(*xy(pts),      marker='+', color='r')
        plt.show()

if __name__ == '__main__':
    unittest.main(verbosity=2)
