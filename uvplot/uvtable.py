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

        self._u = uvdata[:, 0].copy()/wle
        self._v = uvdata[:, 1].copy()/wle
        self._re = uvdata[:, 2].copy()
        self._im = uvdata[:, 3].copy()
        self._w = uvdata[:, 4].copy()

        self.ndat = len(self.u)
        self._uvdist = None

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, value):
        self._u = value

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        self._v = value
        
    @property
    def re(self):
        return self._re

    @re.setter
    def re(self, value):
        self._re = value
        
    @property
    def im(self):
        return self._im

    @im.setter
    def im(self, value):
        self._im = value
        
    @property
    def w(self):
        return self.w

    @w.setter
    def w(self, value):
        self.w = value

    @property
    def uvdist(self):
        if self._uvdist is None:
            self._uvdist = np.hypot(self._u, self._v)

        return self._uvdist

    def uvbin(self, uvbin_size, inc=0, PA=0):
        raise NotImplementedError

    def apply_phase(self, dRA=0, dDec=0):
        """
        Apply a phase change in the uv space, corresponding to a translation by
        angular offsets (dRA, dDec) in the plane of the sky.

        Parameters
        ----------
        dRA : float, optional
            Offset along the Right Ascension direction. dRA>0 translates image towards East.
            **units**: rad
        dDec : float, optional
            Offset along the Declination direction. dDec>0 translates image towards North.
            **units**: rad

        """
        dRA *= 2. * np.pi
        dDec *= 2. * np.pi

        phi = self.u * dRA + self.v * dDec
        vis = (self._re + 1j*self._im) * (np.cos(phi) + 1j*np.sin(phi))

        self._re = vis.real
        self._im = vis.imag

    def uvplot(self):
        raise NotImplementedError

