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
import sys
from codecs import open
from setuptools import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py register')
    os.system('python setup.py sdist upload')
    os.system('python setup.py bdist_wheel upload')
    sys.exit()


version = '0.1.2'
requires = ["numpy>=1.7.0"]


def read(f):
    return open(f, encoding='utf-8').read()


setup(
    name='pyefd',
    version=version,
    author='Henrik Blidh',
    author_email='henrik.blidh@nedomkull.com',
    description='Python implementation of "Elliptic Fourier Features of a Closed Contour"',
    long_description=read('README.rst') + '\n\n' + read('HISTORY.rst'),
    license='MIT',
    url='https://github.com/hbldh/pyefd',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
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
    keywords=["elliptic fourier descriptors", "shape descriptors", "image analysis"],
    py_modules=['pyefd'],
    test_suite="tests",
    zip_safe=False,
    include_package_data=True,
    platforms='any',
    install_requires=requires,
    setup_requires=['pytest-runner', ],
    tests_require=['pytest', ],
    package_data={},
    dependency_links=[],
    ext_modules=[],
    entry_points={},
)
