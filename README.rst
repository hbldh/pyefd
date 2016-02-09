PyEFD
=====

.. image:: https://travis-ci.org/hbldh/pyefd.svg?branch=master
    :target: https://travis-ci.org/hbldh/pyefd
.. image:: http://img.shields.io/pypi/v/pyefd.svg
    :target: https://pypi.python.org/pypi/pyefd/
.. image:: http://img.shields.io/pypi/dm/pyefd.svg
    :target: https://pypi.python.org/pypi/pyefd/
.. image:: http://img.shields.io/pypi/l/pyefd.svg
    :target: https://pypi.python.org/pypi/pyefd/

An Python/NumPy implementation of the method described in \[1\].

Usage
-----
::

    $ pip install pyefd

Given a closed contour of a shape, this package can fit a 
[Fourier series](https://en.wikipedia.org/wiki/Fourier_series) 
approximating the shape of the contour::

    from pyefd import elliptic_fourier_descriptors
    coeffs = elliptic_fourier_descriptors(contour, order=10)
   
The coefficients returned are the :math:`a_n`, :math:`b_n`, :math:`c_n` and :math:`d_n` of
the following Fourier series representation of the shape:

.. math::
    \begin{align*}
        \hat{x}(t) & = A_0 + \sum_{n=1}^N\left( a_n \cos \frac{2n\pi t}{T} + b_n \sin \frac{2n\pi t}{T} \right)
        \hat{y}(t) & = C_0 + \sum_{n=1}^N\left( c_n \cos \frac{2n\pi t}{T} + d_n \sin \frac{2n\pi t}{T} \right) 
    \end{align*}

See \[1\] for more technical details.

Testing
-------

Run tests with::

    $ python setup.py test

or with [Pytest](http://pytest.org/latest/)::

    $ py.test tests.py

The tests includes a single image from the MNIST dataset of handwritten digits (\[2\]) as a contour to use
for testing.

Documentation
-------------

See the [Github pages](http://hbldh.github.io/pyefd).

References
----------

\[1\] [Frank P Kuhl, Charles R Giardina, Elliptic Fourier features of a closed contour, 
Computer Graphics and Image Processing, Volume 18, Issue 3, 1982, Pages 236-258, 
ISSN 0146-664X, http://dx.doi.org/10.1016/0146-664X(82)90034-X.](http://www.sciencedirect.com/science/article/pii/0146664X8290034X)

\[2\] [LeCun et al. (1999): The MNIST Dataset Of Handwritten Digits](http://yann.lecun.com/exdb/mnist/)
