#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

try:
    from json.decoder import JSONDecodeError
except ImportError:
    # Python 2.7 does not have JSONDecodeError (which is a subclass of ValueError)
    JSONDecodeError = ValueError

from .constants import clight

__all__ = ["UVTable"]

COLUMNS_V0 = ['u', 'v', 'Re', 'Im', 'weights']
COLUMNS_V1 = ['u', 'v', 'Re', 'Im', 'weights', 'freqs', 'spws']
COLUMNS_V2 = ['u', 'v', 'V', 'weights', 'freqs', 'spws']
COLUMNS_FORMATS = [COLUMNS_V0, COLUMNS_V1, COLUMNS_V2]

UNITS_V0 = ['lambda', 'lambda', 'Jy', 'Jy', 'None']
UNITS_V1 = ['lambda', 'lambda', 'Jy', 'Jy', 'None', 'Hz', 'None']
UNITS_V2 = ['lambda', 'lambda', 'Jy', 'None', 'Hz', 'None']
UNITS_FORMATS = [UNITS_V0, UNITS_V1, UNITS_V2]


def parse_columns(columns):
    return "# Columns\t{}".format(' '.join(columns))


COLUMNS_FORMATS_TEXT = "FORMAT\t\tCOLUMNS\t\t\t\t\t\t\tCOLUMNS_LINE " \
                       "(copy-paste as 2nd line in the ASCII file)\n" \
                       "COLUMNS_V0\t{0}\t\t\t'{1}'\n" \
                       "COLUMNS_V1\t{2}\t'{3}'\n" \
                       "COLUMNS_V2\t{4}\t\t'{5}'\n".format(COLUMNS_V0, parse_columns(COLUMNS_V0),
                                                           COLUMNS_V1, parse_columns(COLUMNS_V1),
                                                           COLUMNS_V2, parse_columns(COLUMNS_V2))


class UVTable(object):
    """
    UV table.

    """

    def __init__(self, uvtable=None, filename=None, format='ascii', columns=None, wle=1., **kwargs):
        """
        Create a UVTable object by importing the visibilities from ASCII file or from a uvtable.
        Provide either a `uvtable` or a `filename`, but not both of them.

        Parameters
        ----------
        uvtable : array-like, float
            A list of 1d arrays (or a 2d array) with `u, v, Re, Im, wegiths` columns.
        filename : str
            Name of the ASCII file containing the uv table.
        format : str
            If a filename is provided, the format of the uvtable: 'ascii' or 'binary' for .npz files.
        columns : list
            Columns mapping. Choose one among COLUMNS_V0, COLUMNS_V1, and COLUMNS_V2 (see Notes).
        wle : float, optional
            Observing wavelength. Default is 1, i.e. the `(u,v)` coordinates are
            assumed to be expressed in units of wavelength.
            **units**: same as the u, v coordinates in the uv table.

        Notes
        -----
        Inside the UVTable object the (u, v) coordinates are stored in units of `wle`.
        
        Columns formats:
        COLUMNS_V0 = ['u', 'v', 'Re', 'Im', 'weights']
        COLUMNS_V1 = ['u', 'v', 'Re', 'Im', 'weights', 'freqs', 'spws']
        COLUMNS_V2 = ['u', 'v', 'V', 'weights', 'freqs', 'spws']

        """
        if filename and uvtable:
            raise ValueError("Cannot provide both filename and uvtable.")

        if not filename and not uvtable:
            raise ValueError("Provide at least a filename or a uvtable.")

        self.wle = wle
        self._u = None
        self._v = None
        self._re = None
        self._im = None
        self._weights = None
        self._freqs = None
        self._spws = None
        self.header = None
        self.columns = None

        if filename:
            format = format.upper()

            if format == 'ASCII':
                self.read_ascii_uvtable(filename, columns=columns)

            if format == 'BINARY':
                self.read_binary_uvtable(filename, columns=columns)

        if uvtable:
            self.import_uvtable(uvtable, columns=columns)
            self.header = kwargs.get('header', None)

        self.ndat = len(self.u)
        self._uvdist = None
        self.bin_uvdist = None

        if self._freqs is not None:
            self._freqs_avg = np.mean(self._freqs)
            self._freqs_wrt_avg = self._freqs / self._freqs_avg

    def import_uvtable(self, uvtable, columns):

        self.columns = columns

        Ncols = len(uvtable)

        assert Ncols == len(columns), "Expect {} columns ({}), but the uvtable list " \
                                      "contains {} columns.".format(len(columns), columns, Ncols)

        if columns == COLUMNS_V0:
            self.u = uvtable[0]
            self.v = uvtable[1]
            self.re = uvtable[2]
            self.im = uvtable[3]
            self.weights = uvtable[4]

        elif columns == COLUMNS_V2:
            self.u = uvtable[0]
            self.v = uvtable[1]
            self.re = uvtable[2].real
            self.im = uvtable[2].imag
            self.weights = uvtable[3]
            self.freqs = uvtable[4]
            self.spws = uvtable[5]

        elif columns == COLUMNS_V1:
            self.u = uvtable[0]
            self.v = uvtable[1]
            self.re = uvtable[2]
            self.im = uvtable[3]
            self.weights = uvtable[4]
            self.freqs = uvtable[5]
            self.spws = uvtable[6]

        self.u = self.u / self.wle
        self.v = self.v / self.wle

    def read_ascii_uvtable(self, filename, columns=None, **kwargs):
        """ Read an ASCII uvtable """
        verbose = kwargs.get('verbose', True)

        if verbose:
            print("Reading uvtable from {} ...".format(filename))

        uvtable = np.require(np.loadtxt(filename, unpack=True), requirements='C')

        Ncols, Nuv_tot = uvtable.shape

        if not columns:
            # try to detect the columns line
            with open(filename, 'r') as f:
                line1 = f.readline()
                line2 = f.readline()

            try:
                header = json.loads(line1)
                columns = header['columns']

            except JSONDecodeError:
                if verbose:
                    print("Line 1 is not a dictionary containing the header.")
                    print("Trying to parse Columns format from line 2 of the ASCII file.")

                try:
                    if len(line2.split(':\t')) < 2:
                        print("Header not found in the ASCII file, please provide the "
                              "columns parameter when creating the UVTable object.")
                        exit(1)

                    head, cols = line2.split(':\t')

                    if line2[0] != '#' or 'COLUMNS' not in line2.upper():
                        print("Unable to find 'columns' line in ASCII file")
                        print("Trying to assume COLUMNS_V0: {}".format(COLUMNS_V0))
                        columns = COLUMNS_V0
                    else:
                        columns = cols.split('\n')[0].split(' ')

                    assert Ncols == len(columns), "Expect {} columns ({}), but the ASCII file " \
                                                  "contains {} columns.".format(len(columns),
                                                                                columns, Ncols)

                except AssertionError:
                    print("Expect the second line of the ASCII file to contain the columns format")
                    print("Alternatively, provide the columns parameter choosing from:")
                    print(COLUMNS_FORMATS_TEXT)

            if columns not in COLUMNS_FORMATS:
                raise AssertionError(
                    "Detected columns {} format is not among the valid formats"
                    "Pease provide `columns` format choosing one of\n"
                    " COLUMNS_V0: {}\n"
                    " COLUMNS_V1: {}\n"
                    " COLUMNS_V2: {}\n".format(columns, COLUMNS_V0, COLUMNS_V1, COLUMNS_V2))

        print("Assuming column format: {}".format(columns))

        self.import_uvtable(uvtable, columns)

        if verbose:
            print("Reading uvtable from {} ...done".format(filename))

    def read_binary_uvtable(self, filename, columns=None, **kwargs):
        """ Read binary uvtable """
        verbose = kwargs.get('verbose', True)

        if verbose:
            print("Reading uvtable from {} ...".format(filename), end='')

        loaded = np.load(filename)

        self.header = loaded['header'].item()

        if columns:
            assert columns == self.header['columns']
        else:
            columns = self.header['columns']

        self.columns = columns

        if columns == COLUMNS_V0:
            self.u = loaded['u']
            self.v = loaded['v']
            self.re = loaded['Re']
            self.im = loaded['Im']
            self.weights = loaded['weights']
        elif columns == COLUMNS_V2:
            self.u = loaded['u']
            self.v = loaded['v']
            self.re = loaded['V'].real
            self.im = loaded['V'].imag
            self.weights = loaded['weights']
            self.freqs = loaded['freqs']
            self.spws = loaded['spws']

        if verbose:
            print("done")

    def save_ascii_uvtable(self, filename, **kwargs):
        """ Save ascii uvtable """
        ascii_fmt = kwargs.get('ascii_fmt', '%10.6e')

        self.header.update(columns=COLUMNS_V1, units=UNITS_V2)

        np.savetxt(filename,
                   np.column_stack([self.u,
                                    self.v,
                                    self.V.real,
                                    self.V.imag,
                                    self.weights,
                                    self.freqs,
                                    self.spws]),
                   fmt=ascii_fmt.encode(), delimiter='\t',
                   header=json.dumps(self.header))

    def save_binary_uvtable(self, filename, **kwargs):
        """ Save binary uvtable """
        compressed = kwargs.get('compressed', True)

        self.header.update(columns=COLUMNS_V2, units=UNITS_V1)

        if compressed:
            np_save_func = np.savez_compressed
        else:
            np_save_func = np.savez

        np_save_func(filename,
                     u=self.u, v=self.v, V=self.V, weights=self.weights,
                     freqs=self.freqs, spws=self.spws,
                     header=self.header)

    def save(self, filename, export_fmt, **kwargs):
        """ Save uvtable """
        verbose = kwargs.get('verbose', True)

        if verbose:
            print("Write uvtable to file {} ...".format(filename), end='')
            sys.stdout.flush()

        export_fmt = export_fmt.upper()

        if export_fmt == "ASCII":
            self.save_ascii_uvtable(filename, **kwargs)

        elif export_fmt == "BINARY":
            self.save_binary_uvtable(filename, **kwargs)

        if verbose:
            print("done")

    @property
    def u(self):
        """ u coordinate (units: observing wavelength) """
        return self._u

    @u.setter
    def u(self, value):
        self._u = np.ascontiguousarray(value)

    @property
    def v(self):
        """ v coordinate (units: observing wavelength) """
        return self._v

    @v.setter
    def v(self, value):
        self._v = np.ascontiguousarray(value)

    @property
    def u_m(self):
        """ u coordinate (units: m) """
        return self.u / self.freqs * clight

    @property
    def v_m(self):
        """ v coordinate (units: m) """
        return self.v / self.freqs * clight

    @property
    def re(self):
        """ Real part of the Visibilities (units: Jy) """
        return self._re

    @re.setter
    def re(self, value):
        self._re = np.ascontiguousarray(value)

    @property
    def im(self):
        """ Imaginary part of the Visibilities (units: Jy) """
        return self._im

    @im.setter
    def im(self, value):
        self._im = np.ascontiguousarray(value)

    @property
    def weights(self):
        """ Weights of the Visibilities """
        return self._weights

    @weights.setter
    def weights(self, value):
        self._weights = np.ascontiguousarray(value)

    @property
    def freqs(self):
        """ Frequency of the Visibilities (units: Hz) """
        return self._freqs

    @freqs.setter
    def freqs(self, value):
        self._freqs = np.ascontiguousarray(value)

    @property
    def freqs_avg(self):
        """ Average frequency of the Visibilities (units: Hz) """
        return self._freqs_avg

    @property
    def freqs_wrt_avg(self):
        """ Frequency of the Visibilities normalised to the average frequency. (units: pure no.) """
        return self._freqs_wrt_avg

    @property
    def spws(self):
        """ Spectral window ID of the visibilities """
        return self._spws

    @spws.setter
    def spws(self, value):
        self._spws = np.ascontiguousarray(value)

    @property
    def V(self):
        """ Complex Visibilities (units: Jy) """
        return self.re + 1j * self.im

    @property
    def uvdist(self):
        """ Uv-distance (units: observing wavelength) """
        self._uvdist = np.hypot(self._u, self._v)

        return self._uvdist

    def uvbin(self, uvbin_size, **kwargs):
        """
        Compute the intervals (bins) of the uv-distances given the size of the bin (bin_size_wle).

        Parameters
        ----------
        bin_size_wle : float
            Bin size in units of the wavelength.

        Note
        ----
        To compute the weights, we do not need to divide by the weight_corr factor since it cancels out when we compute

        """
        self.nbins = np.ceil(self.uvdist.max() / uvbin_size).astype('int')
        self.bin_uvdist = np.zeros(self.nbins)
        self.bin_weights = np.zeros(self.nbins)
        self.bin_count = np.zeros(self.nbins, dtype='int')
        self.uv_intervals = []

        self.uv_bin_edges = np.arange(self.nbins + 1, dtype='float64') * uvbin_size

        for i in range(self.nbins):
            uv_interval = np.where((self.uvdist >= self.uv_bin_edges[i]) &
                                   (self.uvdist < self.uv_bin_edges[i + 1]))
            self.bin_count[i] = len(uv_interval[0])

            if self.bin_count[i] != 0:
                self.bin_uvdist[i] = self.uvdist[uv_interval].sum() / self.bin_count[i]
                self.bin_weights[i] = np.sum(self.weights[uv_interval])
            else:
                self.bin_uvdist[i] = self.uv_bin_edges[i] + 0.5 * uvbin_size

            self.uv_intervals.append(uv_interval)

        self.bin_re, self.bin_re_err = self.bin_quantity(self.re, **kwargs)
        self.bin_im, self.bin_im_err = self.bin_quantity(self.im, **kwargs)

    def bin_quantity(self, x, use_std=False):
        """
        Compute bins of the quantity x based on the intervals of the uv-distances of the current Uvtable.
        To compute the uv-distances use Uvtable.compute_uvdist() and to compute the intervals use Uvtable.compute_uv_intervals().

        Parameters
        ----------
        x : array-like
            Quantity to be binned.
        use_std : bool, optional
            If provided, the error on each bin is computed as the standard deviation divided by sqrt(npoints_per_bin).

        Returns
        -------
        bin_x, bin_x_err : array-like
            Respectively the bins and the bins error of the quantity x.

        """
        bin_x, bin_x_err = np.zeros(self.nbins), np.zeros(self.nbins)

        for i in range(self.nbins):

            if self.bin_count[i] != 0:
                bin_x[i] = np.sum(x[self.uv_intervals[i]] * self.weights[self.uv_intervals[i]]) / \
                           self.bin_weights[i]

                if use_std is True:
                    bin_x_err[i] = np.std(x[self.uv_intervals[i]])
                else:
                    bin_x_err[i] = 1. / np.sqrt(self.bin_weights[i])

        return bin_x, bin_x_err

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
        if dRA == 0 and dDec == 0:
            return

        dRA *= 2. * np.pi
        dDec *= 2. * np.pi

        phi = self.u * dRA + self.v * dDec
        vis = (self._re + 1j * self._im) * (np.cos(phi) + 1j * np.sin(phi))

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

        x_r = x * cos_t - y * sin_t
        y_r = x * sin_t + y * cos_t

        return x_r, y_r

    def deproject(self, inc, PA=0, inplace=True):
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
        inplace : bool, optional
            By default, the `(u,v)` coordinates stored in the UVTable object are overwritten.
            If False, `deproject` returns a copy of the UVTable object with deprojected `(u,v)` coordinates.

        Returns
        -------
        If inplace is True, a copy of the UVTable object, with deprojected `(u,v)` coordinates.

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

        if inplace:
            self.u = u_deproj
            self.v = v_deproj

        else:
            if self.columns == COLUMNS_V0:
                return UVTable(uvtable=[u_deproj, v_deproj, self.re, self.im, self.weights],
                               columns=COLUMNS_V0, header=self.header)
            elif self.columns == COLUMNS_V1:
                return UVTable(uvtable=[u_deproj, v_deproj, self.re, self.im, self.weights,
                                        self.freqs, self.spws],
                               columns=COLUMNS_V1, header=self.header)
            elif self.columns == COLUMNS_V2:
                return UVTable(uvtable=[u_deproj, v_deproj, self.V, self.weights, self.freqs,
                                        self.spws],
                               columns=COLUMNS_V2, header=self.header)

    def uvcut(self, maxuv, verbose=False):
        """
       Apply uv cut to a table. Consider only baselines shorter than maxuv.

       Parameters
       ----------
       maxuv : float
           Maximum baseline to be considered.
           **units**: observing wavelength.
       verbose : bool, optional
           If true, print an informative message.

       Returns:
       u, v, w, re, im, w : ndarrays
           Visibilities.

       """
        uvcut = self.uvdist <= maxuv

        if verbose:
            print("Consider only baselines up to {} klambda ({} out of {} uv-points)".format(
                maxuv / 1e3, np.count_nonzero(uvcut), self.ndat))

        if self.columns == COLUMNS_V0:
            return UVTable(uvtable=[a[uvcut] for a in [self.u, self.v, self.re, self.im,
                                                       self.weights]],
                           columns=COLUMNS_V0, header=self.header)
        elif self.columns == COLUMNS_V1:
            return UVTable(uvtable=[a[uvcut] for a in [self.u, self.v, self.re, self.im,
                                                       self.weights, self.freqs, self.spws]],
                           columns=COLUMNS_V1, header=self.header)
        elif self.columns == COLUMNS_V2:
            return UVTable(uvtable=[a[uvcut] for a in [self.u, self.v, self.V,
                                                       self.weights, self.freqs, self.spws]],
                           columns=COLUMNS_V2, header=self.header)

    def plot(self, fig_filename=None, color='k', linestyle='.', label='',
             fontsize=18, linewidth=2.5, alpha=1., yerr=True, caption=None, axes=None,
             uvbin_size=0, vis_filename=None, verbose=True):
        """
        Produce a uv plot.

        Parameters
        ----------
        fig_filename : str, optional
            File name of the output plot.
        color : str, optional
            Line color.
        linestyle : str, optional
            Line style.
        label : str, optional
            Legend label.
        fontsize : int, optional
            Font size to be used in the text of the plot.
        linewidth : float, optional
            Line width of the errobar.
        alpha : float, optional
            Transparency of the errorbar.
        yerr : bool, optional
            If True, the y errors are shown. Default is True.
        caption : dict, optional
            Caption for the whole plot. Must be a dictionary with following keys:
            'x' and 'y' (the coordinates of the caption in units of the plotted quantities),
            'fontsize' (the fontsize to be used), 'text' (the caption text).
        axes : matplotlib.Axes, optional
            If provided, the plots are done in axes, otherwise a new figure is created.
        uvbin_size : float, optional
            Size of the uv-distance bin.
            **units**: observing wavelength.
        vis_filename : str, optional
            File name where to store the binned visibiltiies, e.g. "vis_binned.txt".
        verbose : bool, optional
            If True, print informative messages.

        Returns
        -------
        axes : matplotlib.Axes
            The Real and Imaginary axes on which the uv plot is done.

        """
        if not axes:
            fig = plt.figure(figsize=(6, 6))
            gs = GridSpec(2, 1, height_ratios=[4, 1])
            axes = plt.subplot(gs[0]), plt.subplot(gs[1])

        assert len(axes) == 2
        ax_Re, ax_Im = axes

        if uvbin_size != 0:
            self.uvbin(uvbin_size)

        if self.bin_uvdist is None:
            raise AttributeError("Expected uv table with already binned data, or "
                                 "an input parameter uvbin_size != 0")

        uvbins = self.bin_uvdist / 1.e3

        mask = self.bin_count != 0

        if verbose:
            print("Masking {} uv bins".format(len(uvbins[~mask])))

        ax_Re.errorbar(uvbins[mask], self.bin_re[mask],
                       yerr=self.bin_re_err[mask] if yerr is True else None,
                       fmt=linestyle,
                       color=color, linewidth=linewidth, capsize=2.5,
                       markersize=13, elinewidth=0., label=label, alpha=alpha)

        ax_Im.errorbar(uvbins[mask], self.bin_im[mask],
                       yerr=self.bin_im_err[mask] if yerr is True else None,
                       fmt=linestyle,
                       color=color, linewidth=linewidth, capsize=2.5,
                       markersize=13, elinewidth=0., label=label, alpha=alpha)

        ax_Im.set_xlabel(r'uv-distance (k$\lambda$)', fontweight='bold',
                         fontsize=fontsize)
        ax_Re.set_ylabel('Re(V) (Jy)', fontweight='bold', fontsize=fontsize)
        ax_Im.set_ylabel('Im(V) (Jy)', fontweight='bold', fontsize=fontsize)
        ax_Re.yaxis.set_label_coords(-0.25, 0.5)
        ax_Im.yaxis.set_label_coords(-0.25, 0.5)

        ax_Re.tick_params(labelbottom=False)

        if vis_filename:
            np.savetxt(vis_filename,
                       np.stack([uvbins[mask],
                                 self.bin_re[mask], self.bin_re_err[mask],
                                 self.bin_im[mask], self.bin_im_err[mask]], axis=-1),
                       header="uv-distance\tRe(V)\te_Re(V)\tIm(V)\te_Im(V)\n"
                              "(klambda)\t(Jy)\t(Jy)\t(Jy)\t(Jy)")
            if verbose:
                print("Visibilities written to file {}".format(vis_filename))

        # ax[0].set_xlim(uvlim)
        # ax[1].set_xlim(uvlim)
        # ax[0].set_ylim(Jylims[0])
        # ax[1].set_ylim(Jylims[1])
        # ax_Re.set_xticklabels(["")
        # ax[1].set_yticklabels(Jyticks_im)

        if label:
            ax_Re.legend(fontsize=fontsize, frameon=False)

        if caption:
            caption_keys = ['x', 'y', 'text', 'fontsize']
            if all(k in caption for k in caption_keys):
                pass
            else:
                raise KeyError(
                    "Expected caption dict with following keys: {}, but got {}".format(
                        caption_keys, caption.keys()))

            ax_Re.text(caption['x'], caption['y'], caption['text'],
                       fontsize=caption['fontsize'], fontweight='bold')

        ax_Re.figure.subplots_adjust(left=0.25, right=0.97, hspace=0., bottom=0.15, top=0.98)

        if fig_filename:
            plt.savefig(fig_filename)
        else:
            return axes
