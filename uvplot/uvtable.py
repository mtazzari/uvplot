#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np

__all__ = ["UVTable"]


class UVTable(object):
    """
    UV table.

    """
    def __init__(self, filename, wle, format='uvtable'):
        """
        Init the UVTable object by importing the visibilities from file.

        Parameters
        ----------
        filename : str
            Filename of the uv table.
        wle : float
            Observing wavelength.
            **units**: typically meters, same as the u, v coordinates in the uv table.

        Notes
        -----
        The (u, v) coordinates are stored in units of the observing wavelength `wle`.

        """
        self.filename = filename
        self.wle = wle

        if format == 'uvtable':
            uvdata = np.loadtxt(filename)
        else:
            raise NotImplementedError

        self.u = uvdata[:, 0].copy()/wle
        self.v = uvdata[:, 1].copy()/wle
        self.re = uvdata[:, 2].copy()
        self.im = uvdata[:, 3].copy()
        self.w = uvdata[:, 4].copy()

        self.ndat = len(self.u)

    @property
    def u(self):
        return self.u

    @u.setter
    def u(self, value):
        self.u = value

    @property
    def v(self):
        return self.v

    @v.setter
    def v(self, value):
        self.v = value
        
    @property
    def Re(self):
        return self.Re

    @Re.setter
    def Re(self, value):
        self.Re = value
        
    @property
    def Im(self):
        return self.Im

    @Im.setter
    def Im(self, value):
        self.Im = value
        
    @property
    def w(self):
        return self.w

    @w.setter
    def w(self, value):
        self.w = value

    def uvbin(self, uvbin_size, inc=0, PA=0):
        raise NotImplementedError

    def apply_phase(self, dRA=0, dDec=0):
        raise NotImplementedError

    def uvplot(self):
        raise NotImplementedError

