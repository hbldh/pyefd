#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`efd`
==========

Created by hbldh <henrik.blidh@nedomkull.com>
Created on 2016-01-30

A Python implementation of the method described in [1]_ for
calculating Fourier coefficients for characterizing
closed contours.

References:
-----------

.. [1] F. P. Kuhl and C. R. Giardina, “Elliptic Fourier Features of a
   Closed Contour," Computer Vision, Graphics and Image Processing,
       Vol. 18, pp. 236-258, 1982.

.. [2] Oivind Due Trier, Anil K. Jain and Torfinn Taxt, “Feature Extraction
   Methods for Character Recognition - A Survey”, Pattern Recognition
   Vol. 29, No.4, pp. 641-662 (1996)

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

try:
    _range = xrange
except NameError:
    _range = range


def elliptic_fourier_descriptors(contour, order=10, normalize=False):
    """Calculate elliptical Fourier descriptors for a contour.

    :param contour: A contour array of size [M x 2].
    :type contour: :py:class:`numpy.ndarray`
    :param order: The order of Fourier coefficients to calculate.
    :type order: int
    :param normalize: If the coefficients should be normalized
                      as is described in [1]_ and [2]_.
    :type normalize: bool
    :return: A [n x 4] array of Fourier coefficients.
    :rtype: :py:class:`numpy.ndarray`

    """
    dxy = np.diff(contour, axis=0)
    dt = np.sqrt((dxy ** 2).sum(axis=1))
    t = np.concatenate([([0., ]), np.cumsum(dt)])
    T = t[-1]

    phi = (2 * np.pi * t) / T

    coeffs = np.zeros((order, 4))
    for n in _range(1, order + 1):
        const = T / (2 * n * n * np.pi * np.pi)
        phi_n = phi * n
        d_cos_phi_n = np.cos(phi_n[1:]) - np.cos(phi_n[:-1])
        d_sin_phi_n = np.sin(phi_n[1:]) - np.sin(phi_n[:-1])
        a_n = const * np.sum((dxy[:, 0] / dt) * d_cos_phi_n)
        b_n = const * np.sum((dxy[:, 0] / dt) * d_sin_phi_n)
        c_n = const * np.sum((dxy[:, 1] / dt) * d_cos_phi_n)
        d_n = const * np.sum((dxy[:, 1] / dt) * d_sin_phi_n)
        coeffs[n - 1, :] = a_n, b_n, c_n, d_n

    if normalize:
        coeffs = normalize_efd(coeffs)

    return coeffs


def normalize_efd(coeffs, size_invariant=True):
    """Normalizes an array of Fourier coefficients.

    See details in [1]_ or [2]_.

    :param coeffs: A [n x 4] Fourier coefficient array.
    :type coeffs: :py:class:`numpy.ndarray`
    :param size_invariant: If size invariance normalizing should be done as well.
                           Default is `True`
    :type size_invariant: bool
    :return: The normalized [n x 4] Fourier coefficient array.
    :rtype: :py:class:`numpy.ndarray`

    """
    # Make the coefficients have a zero phase shift from
    # the first major axis. Theta_1 is that shift angle.
    theta_1 = 0.5 * np.arctan2(
        2 * ((coeffs[0, 0] * coeffs[0, 1]) + (coeffs[0, 2] * coeffs[0, 3])),
        ((coeffs[0, 0] ** 2) - (coeffs[0, 1] ** 2) + (coeffs[0, 2] ** 2) - (coeffs[0, 3] ** 2)))
    # Rotate all coefficients by theta_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs[n - 1, :] = np.dot(
            np.array([[coeffs[n - 1, 0], coeffs[n - 1, 1]],
                      [coeffs[n - 1, 2], coeffs[n - 1, 3]]]),
            np.array([[np.cos(n * theta_1), -np.sin(n * theta_1)],
                      [np.sin(n * theta_1), np.cos(n * theta_1)]])).flatten()

    # Make the coefficients rotation invariant by rotating so that
    # the semi-major axis is parallel to the x-axis.
    psi_1 = np.arctan2(coeffs[0, 2], coeffs[0, 0])
    psi_rotation_matrix = np.array([[np.cos(psi_1), np.sin(psi_1)],
                                    [-np.sin(psi_1), np.cos(psi_1)]])
    # Rotate all coefficients by -psi_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs[n - 1, :] = psi_rotation_matrix.dot(
            np.array([[coeffs[n - 1, 0], coeffs[n - 1, 1]],
                      [coeffs[n - 1, 2], coeffs[n - 1, 3]]])).flatten()

    if size_invariant:
        # Obtain size-invariance by normalizing.
        coeffs /= np.abs(coeffs[0, 0])

    return coeffs


def calculate_dc_coefficients(contour):
    """Calculate the A0 and C0 coefficients of the elliptic Fourier series

    :param contour: A contour array of size [M x 2].
    :type contour: :py:class:`numpy.ndarray`
    :return: The A0 and C0 coefficients.
    :rtype: tuple

    """
    dxy = np.diff(contour, axis=0)
    dt = np.sqrt((dxy ** 2).sum(axis=1))
    t = np.concatenate([([0., ]), np.cumsum(dt)])
    T = t[-1]

    xi = np.cumsum(dxy[:, 0]) - (dxy[:, 0] / dt) * t[1:]
    A0 = (1 / T) * np.sum(((dxy[:, 0] / (2 * dt)) * np.diff(t ** 2)) + xi * dt)
    delta = np.cumsum(dxy[:, 1]) - (dxy[:, 1] / dt) * t[1:]
    C0 = (1 / T) * np.sum(((dxy[:, 1] / (2 * dt)) * np.diff(t ** 2)) + delta * dt)

    # A0 and CO relate to the first point of the contour array as origin.
    # Adding those values to the coefficients to make them relate to true origin.
    return contour[0, 0] + A0, contour[0, 1] + C0


def plot_efd(coeffs, locus=(0., 0.), image=None, contour=None, n=300):
    """Plot a [2 x (n/2)] grid of successive truncations of the series.

    :param coeffs:  [n x 4] Fourier coefficient array.
    :type coeffs: :py:class:`numpy.ndarray`
    :param locus: The A_0 and C_0 elliptic locus in [1]_ and [2]_.
    :type locus: list, tuple or :py:class:`numpy.ndarray`
    :param n: Number of points to use for plotting of Fourier series.
    :type n: int

    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Cannot plot: matplotlib was not installed.")
        return

    N = coeffs.shape[0]
    N_half = int(np.ceil(N / 2))
    n_rows = 2

    t = np.linspace(0, 1.0, n)
    xt = np.ones((n,)) * locus[0]
    yt = np.ones((n,)) * locus[1]

    for n in _range(coeffs.shape[0]):
        xt += (coeffs[n, 0] * np.cos(2 * (n + 1) * np.pi * t)) + \
              (coeffs[n, 1] * np.sin(2 * (n + 1) * np.pi * t))
        yt += (coeffs[n, 2] * np.cos(2 * (n + 1) * np.pi * t)) + \
              (coeffs[n, 3] * np.sin(2 * (n + 1) * np.pi * t))
        ax = plt.subplot2grid((n_rows, N_half), (n // N_half, n % N_half))
        ax.set_title(str(n + 1))
        if contour is not None:
            ax.plot(contour[:, 1], contour[:, 0], 'c--', linewidth=2)
        ax.plot(yt, xt, 'r', linewidth=2)
        if image is not None:
            ax.imshow(image, plt.cm.gray)

    plt.show()


