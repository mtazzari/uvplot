#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np

from .constants import clight

__all__ = ["export_uvtable"]


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
            raise RuntimeError \
                ("Missing split parameter: provide the CASA split object in input. See typical call signature.")
        if vis != "" and vis != split_args['vis']:
            # raise RuntimeError("extract_uvtable: either provide `vis` or `split_args` as input parameters, not both.")
            raise RuntimeError(
                "extract_uvtable: vis={} input parameter doesn't match with split_args['vis']={}".format(
                    vis, split_args['vis']))

        if not 'outputvis' in split_args.keys():
            split_args.update(outputvis='mstable_tmp.ms'.encode())

        MStb_name = split_args['outputvis']

        if verbose:
            print \
                ("Applying split. Creating temporary MS table {} from original MS table {}".format
                    (MStb_name, split_args['vis']))

        split(**split_args)
    else:
        if vis == "":
            raise RuntimeError \
                ("Missing vis parameter: provide a valid MS table filename.")

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
        if split_args is None or \
                (split_args is not None and 'spw' not in split_args.keys()):
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
    if verbose: print("Exporting visibilities to {}".format(uvtable_filename))
    np.savetxt(uvtable_filename,
               np.column_stack([u, v, w, V.real, V.imag, weights]),
               fmt=fmt.encode(), delimiter='\t',
               header='Extracted from {}.\nwavelength[m] = {}\nColumns:\tu[m]\tv[m]\tw[m]\tRe(V)[Jy]\tIm(V)[Jy]\tweight'.format(
                   MStb_name, wle))

    if split_args:
        if not keep_tmp_ms:
            from subprocess import call
            if verbose:
                print("Removing temporary MS table {}".format(
                    split_args['outputvis']))
            call("rm -rf {}".format(split_args['outputvis']), shell=True)
