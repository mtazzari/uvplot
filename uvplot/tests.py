#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_allclose

from uvplot import UVTable
from .uvtable import COLUMNS_V0


def create_sample_uvtable(uvtable_filename):
    u = np.random.randn(10000)*1e4  # m
    v = np.random.randn(10000)*1e4  # m
    re = np.random.randn(10000)  # Jy
    im = np.random.randn(10000)  # Jy
    w = np.random.rand(10000)*1e4

    wle = 0.88e-3  # m

    np.savetxt(uvtable_filename, np.stack((u, v, re, im, w), axis=-1))

    return (u, v, re, im, w), wle


def test_init_uvtable():

    uvtable_filename = "/tmp/uvtable.txt"
    uvtable, wle = create_sample_uvtable(uvtable_filename)
    u, v, re, im, w = uvtable

    # test reading from file
    # format='uvtable'
    uvt_file = UVTable(filename=uvtable_filename, wle=wle, columns=COLUMNS_V0)

    # test providing uvtable tuple
    # takes u, v in units of observing wavelength
    uvt_uvtable1 = UVTable(uvtable=(u, v, re, im, w), wle=wle, columns=COLUMNS_V0)
    uvt_uvtable2 = UVTable(uvtable=(u/wle, v/wle, re, im, w), columns=COLUMNS_V0)

    reference = np.hypot(u/wle, v/wle)

    assert_allclose(uvt_uvtable1.uvdist, reference, atol=0, rtol=1e-16)
    assert_allclose(uvt_uvtable2.uvdist, reference, atol=0, rtol=1e-16)
    assert_allclose(uvt_file.uvdist, reference, atol=0, rtol=1e-16)


def test_deproject():

    uvtable_filename = "/tmp/uvtable.txt"
    uvtable, wle = create_sample_uvtable(uvtable_filename)

    uv = UVTable(filename=uvtable_filename, wle=wle, columns=COLUMNS_V0)

    inc = np.radians(30)
    uv_30 = uv.deproject(inc, inplace=False)

    assert_allclose(uv_30.u, uv.u*np.cos(inc))

    uv.deproject(inc)  # inplace=True by default
    assert_allclose(uv_30.u, uv.u)


def test_uvcut():

    uvtable_filename = "/tmp/uvtable.txt"
    uvtable, wle = create_sample_uvtable(uvtable_filename)

    uv = UVTable(filename=uvtable_filename, wle=wle, columns=COLUMNS_V0)
    uv.header = dict(test='test_header')

    maxuv = 5e3
    uvt = uv.uvcut(maxuv)

    # manual uvcut
    uvcut = uv.uvdist <= maxuv

    assert_allclose(uvt.u, uv.u[uvcut])
    assert_allclose(uvt.v, uv.v[uvcut])
    assert_allclose(uvt.re, uv.re[uvcut])
    assert_allclose(uvt.im, uv.im[uvcut])
    assert_allclose(uvt.weights, uv.weights[uvcut])

    assert uv.header == uvt.header
