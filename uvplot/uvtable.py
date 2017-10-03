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
    def __init__(self, uvtable=None, filename="", wle=1, **kwargs):
        """
        Init the UVTable object by importing the visibilities from file.

        Parameters
        ----------
        filename : str
            Name of the file containing the uv table.
        wle : float, optional
            Observing wavelength. Default is 1, i.e. the `(u,v)` coordinates are
            assumed to be expressed in units of wavelength.
            **units**: same as the u, v coordinates in the uv table.

        Notes
        -----
        The (u, v) coordinates are stored in units of the observing wavelength `wle`.

        """
        if filename != "":
            self.filename = filename

            uvdata = self.read_uvtable(self.filename, kwargs.get('format', 'uvtable'))

            u = uvdata[:, 0]
            v = uvdata[:, 1]
            re = uvdata[:, 2]
            im = uvdata[:, 3]
            weights = uvdata[:, 4]

        if uvtable is not None:
            u, v, re, im, weights = uvtable

        self._u = u/wle
        self._v = v/wle
        self._re = re
        self._im = im
        self._weights = weights

        self.ndat = len(self.u)
        self._uvdist = None

    @staticmethod
    def read_uvtable(filename, format):
        """ Read uvtable from file, given a specific format. """
        if format == 'uvtable':
            uvdata = np.loadtxt(filename)

        else:
            raise NotImplementedError

        return uvdata

    @property
    def u(self):
        return self._u

    @u.setter
    def u(self, value):
        self._u = np.ascontiguousarray(value)

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, value):
        self._v = np.ascontiguousarray(value)
        
    @property
    def re(self):
        return self._re

    @re.setter
    def re(self, value):
        self._re = np.ascontiguousarray(value)
        
    @property
    def im(self):
        return self._im

    @im.setter
    def im(self, value):
        self._im = np.ascontiguousarray(value)
        
    @property
    def weights(self):
        return self._weights

    @weights.setter
    def weights(self, value):
        self._weights = np.ascontiguousarray(value)

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

    @staticmethod
    def rotate(x, y, theta):
        """
        Rotate `(x,y)` coordinates counter-clockwise by an angle `theta`.

        Parameters
        ----------
        x : array-like, float
            X coordinates.
        y : array-like, float
            Y coordinates.
        theta : float
            Angle of rotation.
            **units**: rad

        Returns
        -------
        x_r, y_r : array-like
            X and Y coordinates rotated by an angle `theta`.

        """
        if theta == 0:
            return x, y

        cos_t = np.cos(theta)
        sin_t = np.sin(theta)

        x_r = x*cos_t - y*sin_t
        y_r = x*sin_t + y*cos_t

        return x_r, y_r

    def deproject(self, inc, PA=0):
        """
        Deproject `(u,v)` coordinates.
        First, a rotation of a position angle `PA` is applied, and then a deprojection by inclination `inc`.

        Paramters
        ---------
        inc : float
            Inclination.
            **units**: rad
        PA : float, optional
            Position Angle (rad).
            **units**: rad

        Returns
        -------
        A copy of the UVTable object, with deprojected `(u,v)` coordinates.

        """
        cos_inc = np.cos(inc)

        # Rotation by -PA
        # Note: the right ascension is a reversed x axis, thus the anti-rotation
        # of a reversed PA Angle is the same of a direct rotation.
        u_deproj, v_deproj = self.rotate(self.u.copy(), self.v.copy(), PA)

        # Deprojection
        # Note u and v in the Fourier space, thus
        # instead of dividing by cos(), we must multiply
        u_deproj *= cos_inc
        return UVTable((u_deproj, v_deproj, self.re, self.im, self.weights))


    def uvplot(self):
        raise NotImplementedError

