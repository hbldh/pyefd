PyEFD
=====

[![image](https://travis-ci.org/hbldh/pyefd.svg?branch=master)](https://travis-ci.org/hbldh/pyefd)
[![Documentation Status](https://readthedocs.org/projects/pyefd/badge/?version=latest)](http://pyefd.readthedocs.org/en/latest/?badge=latest)
[![image](http://img.shields.io/pypi/v/pyefd.svg)](https://pypi.python.org/pypi/pyefd/)
[![image](http://img.shields.io/pypi/l/pyefd.svg)](https://pypi.python.org/pypi/pyefd/)
[![image](https://coveralls.io/repos/github/hbldh/pyefd/badge.svg?branch=master)](https://coveralls.io/github/hbldh/pyefd?branch=master)

An Python/NumPy implementation of a method for approximating a contour with a Fourier series, as described in [1].

Installation
------------

``` {.sourceCode .bash}
$ pip install pyefd
```

Usage
-----

Given a closed contour of a shape, generated by e.g. [scikit-image](http://scikit-image.org/) or
 [OpenCV](http://opencv.org/), this package can fit a 
 [Fourier series](https://en.wikipedia.org/wiki/Fourier_series) approximating the shape of the contour.

### General usage examples

This section describes the general usage patterns of `pyefd`.

```python
from pyefd import elliptic_fourier_descriptors
coeffs = elliptic_fourier_descriptors(contour, order=10)
```

The coefficients returned are the `a_n`, `b_n`, `c_n` and `d_n` of the following Fourier series 
representation of the shape.

The coefficients returned are by default normalized so that they are rotation and size-invariant. 
This can be overridden by calling:

```python
from pyefd import elliptic_fourier_descriptors
coeffs = elliptic_fourier_descriptors(contour, order=10, normalize=False)
```

Normalization can also be done afterwards:

```python
from pyefd import normalize_efd
coeffs = normalize_efd(coeffs)
```

### OpenCV example

If you are using [OpenCV](http://opencv.org/) to generate contours, this example shows how to 
connect it to `pyefd`.

```python
import cv2 
import numpy
from pyefd import elliptic_fourier_descriptors

# Find the contours of a binary image using OpenCV.
contours, hierarchy = cv2.findContours(
    im, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

# Iterate through all contours found and store each contour's 
# elliptical Fourier descriptor's coefficients.
coeffs = []
for cnt in contours:
    # Find the coefficients of all contours
    coeffs.append(elliptic_fourier_descriptors(
        numpy.squeeze(cnt), order=10))
```

### Using EFD as features

To use these as features, one can write a small wrapper function:

```python
from pyefd import elliptic_fourier_descriptors

def efd_feature(contour):
    coeffs = elliptic_fourier_descriptors(contour, order=10, normalize=True)
    return coeffs.flatten()[3:]
```

If the coefficients are normalized, then `coeffs[0, 0] = 1.0`, `coeffs[0, 1] = 0.0` and 
`coeffs[0, 2] = 0.0`, so they can be disregarded when using the elliptic Fourier descriptors as features.

See [1] for more technical details.

Testing
-------

Run tests with with [Pytest](http://pytest.org/latest/):

```bash
$ py.test tests.py
```

The tests include a single image from the MNIST dataset of handwritten digits ([2]) as a contour to use for testing.

Documentation
-------------

See [ReadTheDocs](http://pyefd.readthedocs.org/).

References
----------

[1]: [Frank P Kuhl, Charles R Giardina, Elliptic Fourier features of a closed contour, Computer Graphics and Image Processing, Volume 18, Issue 3, 1982, Pages 236-258, ISSN 0146-664X, <http://dx.doi.org/10.1016/0146-664X(82)90034-X>.](http://www.sci.utah.edu/~gerig/CS7960-S2010/handouts/Kuhl-Giardina-CGIP1982.pdf)

[2]: [LeCun et al. (1999): The MNIST Dataset Of Handwritten Digits](http://yann.lecun.com/exdb/mnist/)