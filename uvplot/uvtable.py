#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from .constants import clight

__all__ = ["UVTable", "export_uvtable"]


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

            uvdata = self.read_uvtable(self.filename, kwargs.get('format', 'ascii'))

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
        self.bin_uvdist = None

    @staticmethod
    def read_uvtable(filename, format):
        """ Read uvtable from file, given a specific format. """
        if format == 'ascii':
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
        self.nbins = np.ceil(self.uvdist.max()/uvbin_size).astype('int')
        self.bin_uvdist = np.zeros(self.nbins)
        self.bin_weights = np.zeros(self.nbins)
        self.bin_count = np.zeros(self.nbins, dtype='int')
        self.uv_intervals = []

        self.uv_bin_edges = np.arange(self.nbins+1, dtype='float64')*uvbin_size

        for i in range(self.nbins):
            uv_interval = np.where((self.uvdist >= self.uv_bin_edges[i]) &
                                   (self.uvdist < self.uv_bin_edges[i+1]))
            self.bin_count[i] = len(uv_interval[0])

            if self.bin_count[i] != 0:
                self.bin_uvdist[i] = self.uvdist[uv_interval].sum()/self.bin_count[i]
                self.bin_weights[i] = np.sum(self.weights[uv_interval])
            else:
                self.bin_uvdist[i] = self.uv_bin_edges[i]+0.5*uvbin_size

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
                bin_x[i] = np.sum(x[self.uv_intervals[i]]*self.weights[self.uv_intervals[i]])/self.bin_weights[i]

                if use_std is True:
                    bin_x_err[i] = np.std(x[self.uv_intervals[i]])
                else:
                    bin_x_err[i] = 1./np.sqrt(self.bin_weights[i])

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
            return UVTable((u_deproj, v_deproj, self.re, self.im, self.weights))

    def plot(self, fig_filename=None, color='k', linestyle='.', label='',
             fontsize=18, yerr=True, caption=None, axes=None,
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
                       color=color, linewidth=2.5, capsize=2.5,
                       markersize=13, elinewidth=0., label=label)

        ax_Im.errorbar(uvbins[mask], self.bin_im[mask],
                       yerr=self.bin_im_err[mask] if yerr is True else None,
                       fmt=linestyle,
                       color=color, linewidth=2.5, capsize=2.5,
                       markersize=13, elinewidth=0., label=label)

        ax_Im.set_xlabel(r'uv-distance (k$\lambda$)', fontweight='bold',
                         fontsize=fontsize)
        ax_Re.set_ylabel('Re(V) (Jy)', fontweight='bold', fontsize=fontsize)
        ax_Im.set_ylabel('Im(V) (Jy)', fontweight='bold', fontsize=fontsize)
        ax_Re.yaxis.set_label_coords(-0.25, 0.5)
        ax_Im.yaxis.set_label_coords(-0.25, 0.5)

        ax_Re.set_xticklabels("")

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
        # ax[1].set_yticks(Jyticks_im)
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


def export_uvtable(uvtable_filename, tb, vis="", split_args=None, split=None,
                   dualpol=True, fmt='%10.6e', datacolumn="CORRECTED_DATA", keep_tmp_ms=False, verbose=False):
    """
    Export visibilities from an MS Table to a uvtable. Requires execution inside CASA.
    Currently the only uvtable format supported is ASCII.

    Typicall call signature::

    export_uvtable('uvtable_new.txt', tb, split=split,
                   split_args={'vis': 'sample.ms', 'datacolumn': 'data', spw='0,1'}, verbose=True)"

    Parameters
    ----------
    uvtable_filename : str
        Filename of the output uvtable, e.g. "uvtable.txt"
    tb : CASA `tb` object
        As tb parameter you **must** pass the tb object that is defined in the CASA shell.
        Since tb cannot be accessed outside CASA, export_uvtable('uvtable.txt', tb, ...) 
        can be executed only inside CASA.
    vis : str, optional
        MS Table filename, e.g. mstable.ms
    split_args : dict, optional
        Default is None. If provided, perform a split before exporting the uvtable.
        The split_args dictionary is passed to the CASA::split task.
        The CASA::split task must be provided in input as split.
    split : optional
        CASA split task
    dualpol : bool, optional
        If the MS Table contains dual polarisation data. Default is True.
    fmt : str, optional
        Format of the output ASCII uvtable.
    datacolumn: str, optional
        Data column to be extracted, e.g. "DATA", "CORRECTED_DATA", "MODEL_DATA".
    keep_tmp_ms : bool, optional
        If True, keeps the temporary outputvis created by the split command.
    verbose: bool, optional
        If True, print informative messages.

    Note
    ----
    By default, only the 1st spectral window (spw) is exported.
    To export all the spws in an MS table provide split_args, e.g.::

        split_args = {'vis': 'input.ms', 'outputvis': 'input_tmp.ms', spw: '*'}

    Example
    -------
    From within CASA, to extract all the visibilities from an MS table::

        export_uvtable('uvtable.txt', tb, vis='sample.ms')

    where `tb` is the CASA tb object (to inspect it type `tb` in the CASA shell).
    For more information on `tb` see `<https://casa.nrao.edu/docs/CasaRef/table-Module.html>`_
    
    From within CASA, to extract the visibilities in spectral windows 0 and 2 use
    the `split_args` parameter and the CASA `split` task::

        export_uvtable('uvtable.txt', tb, split=split,
         split_args={'vis': 'sample.ms' , 'datacolumn': 'data', 'spw':'0,2'})

    To perform these operations without running CASA interactively::

        casa --nologger --nogui -c "from uvplot import export_uvtable; export_uvtable(...)"

    ensuring that the strings are inside the '..' characters and not the "..." one.

    """
    if vis != "":
        MStb_name = vis

    if split_args:
        if split is None:
            raise RuntimeError("Missing split parameter: provide the CASA split object in input. See typical call signature.")
        if vis != "" and vis != split_args['vis']:
            # raise RuntimeError("extract_uvtable: either provide `vis` or `split_args` as input parameters, not both.")
            raise RuntimeError(
                "extract_uvtable: vis={} input parameter doesn't match with split_args['vis']={}".format(
                    vis, split_args['vis']))

        if not 'outputvis' in split_args.keys():
            split_args.update(outputvis='mstable_tmp.ms'.encode())

        MStb_name = split_args['outputvis']

        if verbose:
            print("Applying split. Creating temporary MS table {} from original MS table {}".format(MStb_name, split_args['vis']))

        split(**split_args)
    else:
        if vis == "":
            raise RuntimeError("Missing vis parameter: provide a valid MS table filename.")

    if verbose:
        print("Reading {}".format(MStb_name))

    tb.open(MStb_name)

    # get coordinates
    uvw = tb.getcol("UVW".encode())
    u, v, w = [uvw[i, :] for i in range(3)]

    # get weights
    weights_orig = tb.getcol("WEIGHT".encode())

    # get visibilities
    tb_columns = tb.colnames()

    if datacolumn.upper() in tb_columns:
        data = tb.getcol(datacolumn.encode())
    else:
        raise KeyError("datacolumn {} is not available.".format(datacolumn))

    spw = tb.getcol("DATA_DESC_ID".encode())
    nspw = len(np.unique(spw))
    if nspw > 1:
        if split_args is None or (split_args is not None and 'spw' not in split_args.keys()):
            print(
            "Warning: the MS table {} has {} spectral windows. By default all of them are exported."
            " To choose which spws to export, provide split_args with the spw parameter.".format(
                MStb_name, nspw))

    ispw = 0

    if dualpol:
        # dual polarisation: extract the polarised visibilities and weights
        V_XX = data[0, ispw, :]
        V_YY = data[1, ispw, :]
        weights_XX = weights_orig[0, :]
        weights_YY = weights_orig[1, :]

        # compute weighted average of the visibilities and weights
        V = (V_XX * weights_XX + V_YY * weights_YY) / (weights_XX + weights_YY)
        weights = 0.5 * (weights_XX + weights_YY)

    else:
        # single polarisation
        V = data[0, ispw, :]
        weigths = weights_orig

    spw_path = tb.getkeyword('SPECTRAL_WINDOW'.encode()).split()[-1]
    tb.close()

    # get the mean observing frequency
    tb.open(spw_path)
    freqs = tb.getcol('CHAN_FREQ'.encode())  # [GHz]
    tb.close()

    wle = clight / freqs.mean()  # [m]

    # export to file as ascii
    if verbose: print(
        "Exporting visibilities to {}".format(uvtable_filename))
    np.savetxt(uvtable_filename,
               np.column_stack([u, v, w, V.real, V.imag, weights]), fmt=fmt.encode(),
               delimiter='\t',
               header='Extracted from {}.\nwavelength[m] = {}\nColumns:\tu[m]\tv[m]\tw[m]\tRe(V)[Jy]\tIm(V)[Jy]\tweight'.format(
                   MStb_name, wle))

    if split_args:
        if not keep_tmp_ms:
            from subprocess import call
            if verbose:
                print("Removing temporary MS table {}".format(split_args['outputvis']))
            call("rm -rf {}".format(split_args['outputvis']), shell=True)
