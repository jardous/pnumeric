#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension

module1 = Extension('pnumeric', sources = ['pnumeric.c', 'vector.c', 'matrix.c', 'cgensupport.c', 'm2/m2.c', 'kf.c', 'fft.c', 'window.c', 'hpspectrum.c'])

setup (name = 'pNumeric',
       version = '1.0',
       description = 'Light matrix manipulation package with Kalman filter',
       ext_modules = [module1])
