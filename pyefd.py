#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

A Python implementation of the method described in [#a]_ and [#b]_ for
calculating Fourier coefficients for characterizing
closed contours.

References
----------

.. [#a] F. P. Kuhl and C. R. Giardina, “Elliptic Fourier Features of a
   Closed Contour," Computer Vision, Graphics and Image Processing,
   Vol. 18, pp. 236-258, 1982.

.. [#b] Oivind Due Trier, Anil K. Jain and Torfinn Taxt, “Feature Extraction
   Methods for Character Recognition - A Survey”, Pattern Recognition
   Vol. 29, No.4, pp. 641-662, 1996

Created by hbldh <henrik.blidh@nedomkull.com> on 2016-01-30.

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


def elliptic_fourier_descriptors(
    contour, order=10, normalize=False, return_transformation=False
):
    """Calculate elliptical Fourier descriptors for a contour.

    :param numpy.ndarray contour: A contour array of size ``[M x 2]``.
    :param int order: The order of Fourier coefficients to calculate.
    :param bool normalize: If the coefficients should be normalized;
        see references for details.
    :param bool return_transformation: If the normalization parametres should be returned.
        Default is ``False``.
    :return: A ``[order x 4]`` array of Fourier coefficients and optionally the
        transformation parametres ``scale``, ``psi_1`` (rotation) and ``theta_1`` (phase)
    :rtype: ::py:class:`numpy.ndarray` or (:py:class:`numpy.ndarray`, (float, float, float))

    """
    dxy = np.diff(contour, axis=0)
    dt = np.sqrt((dxy ** 2).sum(axis=1))
    t = np.concatenate([([0.0]), np.cumsum(dt)])
    T = t[-1]

    phi = (2 * np.pi * t) / T

    orders = np.arange(1, order + 1)
    consts = T / (2 * orders * orders * np.pi * np.pi)
    phi = phi * orders.reshape((order, -1))

    d_cos_phi = np.cos(phi[:, 1:]) - np.cos(phi[:, :-1])
    d_sin_phi = np.sin(phi[:, 1:]) - np.sin(phi[:, :-1])

    a = consts * np.sum((dxy[:, 0] / dt) * d_cos_phi, axis=1)
    b = consts * np.sum((dxy[:, 0] / dt) * d_sin_phi, axis=1)
    c = consts * np.sum((dxy[:, 1] / dt) * d_cos_phi, axis=1)
    d = consts * np.sum((dxy[:, 1] / dt) * d_sin_phi, axis=1)

    coeffs = np.concatenate(
        [
            a.reshape((order, 1)),
            b.reshape((order, 1)),
            c.reshape((order, 1)),
            d.reshape((order, 1)),
        ],
        axis=1,
    )

    if normalize:
        coeffs = normalize_efd(coeffs, return_transformation=return_transformation)

    return coeffs


def normalize_efd(coeffs, size_invariant=True, return_transformation=False):
    """Normalizes an array of Fourier coefficients.

    See [#a]_ and [#b]_ for details.

    :param numpy.ndarray coeffs: A ``[n x 4]`` Fourier coefficient array.
    :param bool size_invariant: If size invariance normalizing should be done as well.
        Default is ``True``.
    :param bool return_transformation: If the normalization parametres should be returned.
        Default is ``False``.
    :return: The normalized ``[n x 4]`` Fourier coefficient array and optionally the
        transformation parametres ``scale``, :math:`psi_1` (rotation) and :math:`theta_1` (phase)
    :rtype: :py:class:`numpy.ndarray` or (:py:class:`numpy.ndarray`, (float, float, float))

    """
    # Make the coefficients have a zero phase shift from
    # the first major axis. Theta_1 is that shift angle.
    theta_1 = 0.5 * np.arctan2(
        2 * ((coeffs[0, 0] * coeffs[0, 1]) + (coeffs[0, 2] * coeffs[0, 3])),
        (
            (coeffs[0, 0] ** 2)
            - (coeffs[0, 1] ** 2)
            + (coeffs[0, 2] ** 2)
            - (coeffs[0, 3] ** 2)
        ),
    )
    # Rotate all coefficients by theta_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs[n - 1, :] = np.dot(
            np.array(
                [
                    [coeffs[n - 1, 0], coeffs[n - 1, 1]],
                    [coeffs[n - 1, 2], coeffs[n - 1, 3]],
                ]
            ),
            np.array(
                [
                    [np.cos(n * theta_1), -np.sin(n * theta_1)],
                    [np.sin(n * theta_1), np.cos(n * theta_1)],
                ]
            ),
        ).flatten()

    # Make the coefficients rotation invariant by rotating so that
    # the semi-major axis is parallel to the x-axis.
    psi_1 = np.arctan2(coeffs[0, 2], coeffs[0, 0])
    psi_rotation_matrix = np.array(
        [[np.cos(psi_1), np.sin(psi_1)], [-np.sin(psi_1), np.cos(psi_1)]]
    )
    # Rotate all coefficients by -psi_1.
    for n in _range(1, coeffs.shape[0] + 1):
        coeffs[n - 1, :] = psi_rotation_matrix.dot(
            np.array(
                [
                    [coeffs[n - 1, 0], coeffs[n - 1, 1]],
                    [coeffs[n - 1, 2], coeffs[n - 1, 3]],
                ]
            )
        ).flatten()

    size = coeffs[0, 0]
    if size_invariant:
        # Obtain size-invariance by normalizing.
        coeffs /= np.abs(size)

    if return_transformation:
        return coeffs, (size, psi_1, theta_1)
    else:
        return coeffs


def calculate_dc_coefficients(contour):
    """Calculate the :math:`A_0` and :math:`C_0` coefficients of the elliptic Fourier series.

    :param numpy.ndarray contour: A contour array of size ``[M x 2]``.
    :return: The :math:`A_0` and :math:`C_0` coefficients.
    :rtype: tuple

    """
    dxy = np.diff(contour, axis=0)
    dt = np.sqrt((dxy ** 2).sum(axis=1))
    t = np.concatenate([([0.0]), np.cumsum(dt)])
    T = t[-1]

    xi = np.cumsum(dxy[:, 0]) - (dxy[:, 0] / dt) * t[1:]
    A0 = (1 / T) * np.sum(((dxy[:, 0] / (2 * dt)) * np.diff(t ** 2)) + xi * dt)
    delta = np.cumsum(dxy[:, 1]) - (dxy[:, 1] / dt) * t[1:]
    C0 = (1 / T) * np.sum(((dxy[:, 1] / (2 * dt)) * np.diff(t ** 2)) + delta * dt)

    # A0 and CO relate to the first point of the contour array as origin.
    # Adding those values to the coefficients to make them relate to true origin.
    return contour[0, 0] + A0, contour[0, 1] + C0


def reconstruct_contour(coeffs, locus=(0, 0), num_points=300):
    """Returns the contour specified by the coefficients.

    :param coeffs: A ``[n x 4]`` Fourier coefficient array.
    :type coeffs: numpy.ndarray
    :param locus: The :math:`A_0` and :math:`C_0` elliptic locus in [#a]_ and [#b]_.
    :type locus: list, tuple or numpy.ndarray
    :param num_points: The number of sample points used for reconstructing the contour from the EFD.
    :type num_points: int
    :return: A list of x,y coordinates for the reconstructed contour.
    :rtype: numpy.ndarray

    """
    t = np.linspace(0, 1.0, num_points)
    # Append extra dimension to enable element-wise broadcasted multiplication
    coeffs = coeffs.reshape(coeffs.shape[0], coeffs.shape[1], 1)

    orders = coeffs.shape[0]
    orders = np.arange(1, orders + 1).reshape(-1, 1)
    order_phases = 2 * orders * np.pi * t.reshape(1, -1)

    xt_all = coeffs[:, 0] * np.cos(order_phases) + coeffs[:, 1] * np.sin(order_phases)
    yt_all = coeffs[:, 2] * np.cos(order_phases) + coeffs[:, 3] * np.sin(order_phases)

    xt_all = xt_all.sum(axis=0)
    yt_all = yt_all.sum(axis=0)
    xt_all = xt_all + np.ones((num_points,)) * locus[0]
    yt_all = yt_all + np.ones((num_points,)) * locus[1]

    reconstruction = np.stack([xt_all, yt_all], axis=1)
    return reconstruction


def plot_efd(coeffs, locus=(0.0, 0.0), image=None, contour=None, n=300):
    """Plot a ``[2 x (N / 2)]`` grid of successive truncations of the series.

    .. note::

        Requires `matplotlib <http://matplotlib.org/>`_!

    :param numpy.ndarray coeffs: ``[N x 4]`` Fourier coefficient array.
    :param list, tuple or numpy.ndarray locus:
        The :math:`A_0` and :math:`C_0` elliptic locus in [#a]_ and [#b]_.
    :param int n: Number of points to use for plotting of Fourier series.

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
        xt += (coeffs[n, 0] * np.cos(2 * (n + 1) * np.pi * t)) + (
            coeffs[n, 1] * np.sin(2 * (n + 1) * np.pi * t)
        )
        yt += (coeffs[n, 2] * np.cos(2 * (n + 1) * np.pi * t)) + (
            coeffs[n, 3] * np.sin(2 * (n + 1) * np.pi * t)
        )
        ax = plt.subplot2grid((n_rows, N_half), (n // N_half, n % N_half))
        ax.set_title(str(n + 1))

        if image is not None:
            # A background image of shape [rows, cols] gets transposed
            # by imshow so that the first dimension is vertical
            # and the second dimension is horizontal.
            # This implies swapping the x and y axes when plotting a curve.
            if contour is not None:
                ax.plot(contour[:, 1], contour[:, 0], "c--", linewidth=2)
            ax.plot(yt, xt, "r", linewidth=2)
            ax.imshow(image, plt.cm.gray)
        else:
            # Without a background image, no transpose is implied.
            # This case is useful when (x,y) point clouds
            # without relation to an image are to be handled.
            if contour is not None:
                ax.plot(contour[:, 0], contour[:, 1], "c--", linewidth=2)
            ax.plot(xt, yt, "r", linewidth=2)
            ax.axis("equal")

    plt.show()
