#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

An example showing how to use pyefd for fitting points along a closed curve.

The demo curves are taken from:
"Optimized Fourier representations for three-dimensional magnetic surfaces"
by S. P. Hirshman and H. K. Meier,
Phys. Fluids 28 (5) 1985 (https://doi.org/10.1063/1.864972)


Created by Jonathan Schilling <jonathan.schilling@mail.de> on 2021-12-09.

"""

import numpy as np

# hack to get the import from within module directory to work
# see also: https://stackoverflow.com/a/16985066
import sys
import os
SCRIPT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
if not SCRIPT_DIR in sys.path:
    sys.path.append(SCRIPT_DIR)

import pyefd

# bean
rCos = np.array([3.0, 1.042, 0.502, -0.0389])
zSin = np.array([0.0, 1.339, 0.296, -0.0176])

# diamond
# rCos = np.array([0.5, 0.427, 0.0,  0.0732])
# zSin = np.array([0.0, 0.427, 0.0, -0.0732])

# D
# rCos = np.array([3.0, 0.991,  0.136])
# zSin = np.array([0.0, 1.409, -0.118])

# belt
# rCos = np.array([3.0, 0.453, 0.0, 0.0  ])
# zSin = np.array([0.0, 0.6  , 0.0, 0.196])

# ellipse
# rCos = np.array([3.0, 1.0])
# zSin = np.array([0.0, 3.0])

# evaluate curve geometry given as Fourier coefficients above
n = 300
contour = np.zeros([n,2])
theta = np.linspace(0.0, 2.0*np.pi, n)
mMax = len(rCos)
for m in range(mMax):
    contour[:,0] += rCos[m] * np.cos(m*theta)
    contour[:,1] += zSin[m] * np.sin(m*theta)

# apply pyefd to get elliptic Fourier descriptors
coeffs = pyefd.elliptic_fourier_descriptors(contour)
a0, c0 = pyefd.calculate_dc_coefficients(contour)
pyefd.plot_efd(coeffs, locus=(a0,c0), contour=contour)
