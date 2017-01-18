#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
  $Id: test_pnumeric.py,v 1.1 2007/01/21 01:38:35 jp Exp $

  test_pnumeric.py
     Tests for pnumeric Python module.

  Copyright (c) 2007 Jiří Popek <jiri.popek@gmail.com>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""

import unittest
from pnumeric import *
from math import pi, sin


class TestSequenceFunctions(unittest.TestCase):

    def test_creating(self):
        '''creating matrixes'''
        a2 = Matrix([[.1, .2], [-.2, .1]])
        a3 = Matrix([[.2, .4, .2], [-.2, .2, .0], [.2, .2, -.2]])
        a43 = Matrix([[.2, .4, .2], [-.2, .2, .0], [.2, .2, -.2], [3., 4., 5.]])

    def test_compare(self):
        '''matrix comparison'''
        a3  = Matrix([[.134345674, .4, .2], [-.2, .2, .0], [.2, .2, -.2]])
        tmp = Matrix([[.134345675, .4, .2], [-.2, .2, .0], [.2, .2, -.2]])
        self.assertEqual(a3, tmp)

        # note EPS is set to 1e-8
        a3  = Matrix([[213, .4, .2], [-.2, .2, .0], [.2, .13434564, -.2]])
        tmp = Matrix([[213, .4, .2], [-.2, .2, .0], [.2, .13434565, -.2]])
        self.assert_(a3 != tmp)

    def test_inverse(self):
        '''matrix inversion'''
        a3 = Matrix([[.2, .4, .2], [-.2, .2, .0], [.2, .2, -.2]])
        i_a3 = a3.inv()
        # should return
        # 1       -3        1
        # 1        2        1
        # 2       -1       -3
        tmp = Matrix([[1, -3, 1], [1, 2, 1], [2, -1, -3]])
        self.assertEqual(i_a3, tmp)

    def test_shape(self):
        '''shape'''
        a3 = Matrix([[.2, .4, .2], [-.2, .2, .0], [.2, .2, -.2], [3., 4., 5.]])
        self.assertEqual((a3.rows, a3.cols), (4, 3))
        self.assertEqual(a3.shape, (4, 3))

    def test_rmultiply(self):
        '''test rmultiply'''
        a = Matrix([[1, 2], [3, 4]])
        am = a*.5 # Float
        res = Matrix([[0.5, 1], [1.5, 2]])
        self.assertEqual(am, res)
        am = a*5 # Int
        res = Matrix([[5, 10], [15, 20]])
        self.assertEqual(am, res)

    def test_lmultiply(self):
        '''test lmultiply'''
        a = Matrix([[1, 2], [3, 4]])
        am = .5*a # Float
        res = Matrix([[0.5, 1], [1.5, 2]])
        self.assertEqual(am, res)
        am = 5*a # Int
        res = Matrix([[5, 10], [15, 20]])
        self.assertEqual(am, res)

    def test_multiply(self):
        '''test multiply'''
        a = Matrix([[1, 2], [3, 4]])
        b = Matrix([[5, 6], [7, 8]])
        ab= Matrix([[19, 22], [43, 50]])
        res = a*b
        self.assertEqual(ab, res)

    def test_add(self):
        '''test add'''
        a = Matrix([[1, 2], [3, 4]])
        b = Matrix([[5, 6], [7, 8]])
        ab= Matrix([[6, 8], [10, 12]])
        res = a+b
        self.assertEqual(ab, res)
        res = b+a
        self.assertEqual(ab, res)
        sa = Matrix([[1.5, 2.5], [3.5, 4.5]])
        res = a+0.5
        self.assertEqual(sa, res)

    def test_sub(self):
        '''test subtract'''
        a = Matrix([[1, 2], [3, 4]])
        b = Matrix([[5, 6], [7, 8]])
        ab= Matrix([[6, 8], [10, 12]])
        res = ab-b
        self.assertEqual(a, res)
        sa = Matrix([[1.5, 2.5], [3.5, 4.5]])
        res = ab-3
        self.assertEqual(Matrix([[3, 5], [7, 9]]), res)
        res = 3-ab
        self.assertEqual(Matrix([[-3, -5], [-7, -9]]), res)

    def test_matrix_row(self):
        a3 = Matrix([[.2, .4, .2], [-.2, .2, .0], [.2, .2, -.2]])
        r = a3[1]
        tmp = Vector([-.2, .2, .0])
        self.assertEqual(r, tmp)

    def test_setitem(self):
        '''matrix set item'''
        a3 = Matrix([[.2, .4, .2], [-.2, .2, .0], [.2, .2, -.2]])
        a3[1][2] = 33
        tmp = Matrix([[.2, .4, .2], [-.2, .2, 33], [.2, .2, -.2]])
        self.assertEqual(a3, tmp)

    def test_vector(self):
        v = Vector([.2, .4, .2, -3333, 41, 89.444444])

    def test_vector_mul(self):
        v = Vector([1, 2, 3])

        self.assertEqual(v*5.0, Vector([5, 10, 15]))
        self.assertEqual(5*v, Vector([5, 10, 15]))
        self.assertEqual(5.0*v, Vector([5.0, 10.0, 15.0]))
        self.assertEqual(Vector([2])*v, Vector([2.0, 4.0, 6.0]))
        self.assertEqual(Vector([1, 2, 3])*v, Vector([1.0, 4.0, 9.0]))
        self.assertRaises(ValueError, lambda: Vector([1, 2])*v)

    def test_vector_abs(self):
        v = abs(Vector([.2, -3333, 41, -89.444444]))
        self.assertEqual(v, Vector([.2, 3333, 41, 89.444444]))

    def test_slice(self):  #TODO
        v = Vector([0, 1, 2, -4, 6, 7, 8, 9])
        s = v[2:5]
        self.assertEqual(s, Vector([2.0, -4.0, 6.0]))

    def test_vrange(self):
        v = vrange(10)
        self.assertEqual(v, Vector([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))
        v = vrange(5, 8)
        self.assertEqual(v, Vector([5, 6, 7]))
        v = vrange(0, 10, 3)
        self.assertEqual(v, Vector([0, 3, 6, 9]))
        v = vrange(3, 0, -2)
        self.assertEqual(v, Vector([3, 1]))
        v = vrange(3, 0)
        self.assertEqual(v, None)
        v = vrange(0.5, 1.1, 0.1)
        v2 = Vector([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
        self.assertEqual(v, v2)

    def x_test_fft(self):
        v = Vector([0, 1, 0, -1, 0, 1, 0, -1])
        out = fft(v)

    def x_test_fft2(self):
        from matplotlib import pylab
        from math import pi, sin
        v = [2*pi*sin(.1*y) for y in range(128)]
        vv = Vector(v)
        f = fft(vv)
        pylab.plot(f[:len(f)/2])
        #pylab.plot(pylab.fft(v))
        pylab.show()

    def test_rms(self):
        v = [sin(2*pi*y/1023.) for y in range(1024)]
        v = Vector(v)
        r = rms(v)
        self.assert_(r>0.7067 and r<0.7068)

    def test_rect(self):
        v = Vector(100 * [1.0])
        self.assertEquals(v, Vector(rect(100)))

    def test_hann(self):
        v = Vector([0.00000, 0.11698, 0.41318, 0.75000, 0.96985, 0.96985, 0.75000, 0.41318, 0.11698, 0.00000])
        h = Vector([round(x, 5) for x in hann(10)])
        self.assertEquals(v, h)

    def test_hamm(self):
        v = Vector([0.080000, 0.187620, 0.460122, 0.770000, 0.972259, 0.972259, 0.770000, 0.460122, 0.187620, 0.080000])
        h = Vector([round(x, 6) for x in hamming(10)])
        self.assertEquals(v, h)


if __name__ == '__main__':
    unittest.main()
