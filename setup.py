# -*- coding: utf-8 -*-
"""
:mod:`setup`
============

.. module:: setup
    :platform: Unix, Windows
    :synopsis: The Python Packaging setup file.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-09-17, 12:39

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import os
from setuptools import setup, find_packages

# Get the long description from the README file
try:
    with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.rst')) as f:
        long_description = f.read()
except:
    with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'README.md')) as f:
        long_description = f.read()

setup(
    name='pyefd',
    version='0.1.0.dev1',
    author='Henrik Blidh',
    author_email='henrik.blidh@nedomkull.com',
    description='Python implementation of "Elliptic Fourier Features of a Closed Contour"',
    long_description=long_description,
    license='MIT',
    url='https://github.com/hbldh/pyefd',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: MIT',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords=["elliptical fourier descriptors", "shape descriptors", "image analysis"],
    packages=find_packages(exclude=('tests', )),
    install_requires=[],
    package_data={},
    dependency_links=[],
    ext_modules=[],
    entry_points={},
)
