import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# def plot_uvplot(uvtables, uvbinsize, wavelength, inc, PA, Jylims, Jyticks_im, linestyle, colors, apply_shift=None, normalized=False,
#                 **kwargs):
#     """
#     Plots a uv-plot.
#
#     Parameters
#     ----------
#     uvtables:
#     :return:
#     """
#     nuvtables = len(uvtables)
#     uvlim = kwargs.get('uvlim', [0, 1000.])
#     inset_caption = kwargs.get('inset_caption', "")
#     weight_corr = kwargs.get('weight_corr', np.ones(nuvtables))
#     disk_ffcontr = kwargs.get('disk_ffcontr', np.zeros(nuvtables))
#     plot_filename = kwargs.get('plot_filename', 'uv_plot.pdf')
#     fontsize = kwargs.get('fontsize', 18)
#     ismodel = kwargs.get('ismodel', [False for i in range(nuvtables)])
#     labels = kwargs.get('labels', ["" for i in range(nuvtables)])
#     fig = plt.figure(figsize=(8, 8))
#     gs = GridSpec(2, 1, height_ratios=[4, 1])
#     ax = plt.subplot(gs[0]), plt.subplot(gs[1])
#
#     # create obs dataset
#     from ..observations import ObsData
#
#     par_obs = {'data_filenames': uvtables,
#                'wle_mm': wavelength,
#                'weight_corr': weight_corr,
#                'disk_ffcontr': disk_ffcontr}
#
#     def add_uvtable():
#
#     obs = ObsData('', par_obs=par_obs)
#     for i, uvtable in enumerate(obs.uvtables):
#         # shift the data
#         if apply_shift:
#             vis_shifted = apply_phase_vis(apply_shift['delta_alpha'][i], apply_shift['delta_delta'][i],
#                                                                  uvtable.u/uvtable.wle, uvtable.v/uvtable.wle, uvtable.re + 1j*uvtable.im)
#             uvtable.re, uvtable.im = vis_shifted.real, vis_shifted.imag
#
#         uvtable.bin_table(PA[i], inc[i], uvbinsize[i], **kwargs)
#
#         print("Initializing the observational datasets and binning the data.")
#
#         uvbins = uvtable.uvdist_bin/1.e3/uvtable.wle
#
#         if normalized == True:
#             uvtable.re_bin /= uvtable.re_bin[0]
#             uvtable.re_bin_err /= uvtable.re_bin[0]
#
#         ax[0].errorbar(uvbins, uvtable.re_bin,
#                        yerr=None if ismodel[i] is True else uvtable.re_bin_err,
#                        fmt=linestyle[i], color=colors[i],
#                        linewidth=2.5, capsize=2.5, markersize=13, elinewidth=0., label=labels[i])
#
#         ax[1].errorbar(uvbins, uvtable.im_bin,
#                        yerr=None if ismodel[i] is True else uvtable.im_bin_err,
#                        fmt=linestyle[i], color=colors[i],
#                        linewidth=2.5, markersize=13, capsize=2.5)
#
#     ax[1].set_xlabel(r'uv-distance (k$\lambda$)', fontweight='bold', fontsize=fontsize)
#     ax[0].set_ylabel('Re(V) (Jy)', fontweight='bold', fontsize=fontsize)
#     ax[1].set_ylabel('Im(V) (Jy)', fontweight='bold', fontsize=fontsize)
#     ax[0].set_xlim(uvlim)
#     ax[1].set_xlim(uvlim)
#     # ax[0].set_ylim(Jylims[0])
#     # ax[1].set_ylim(Jylims[1])
#     ax[0].set_xticklabels("")
#     ax[1].set_yticks(Jyticks_im)
#     ax[1].set_yticklabels(Jyticks_im)
#     ax[0].legend(fontsize=20)
#     ax[0].text(80., Jylims[0][1]*0.85, inset_caption, fontsize=20, fontweight='bold')
#     fig.subplots_adjust(left=0.2, right=0.95, hspace=0.)
#     fig.subplots_adjust(hspace=0)
#
#     return fig
#
#     # fig.savefig(plot_filename)
#     # print("Plot saved in {0}".format(plot_filename))

def uvplot(bin_uvdist, bin_re, bin_im, bin_re_err=None, bin_im_err=None,
           color='k', linestyle='.', label='', fontsize=18, caption=None):

    fig = plt.figure(figsize=(6, 6))
    gs = GridSpec(2, 1, height_ratios=[4, 1])
    ax = plt.subplot(gs[0]), plt.subplot(gs[1])

    print("Initializing the observational datasets and binning the data.")

    uvbins = bin_uvdist / 1.e3
    ax[0].errorbar(uvbins, bin_re,
                   yerr=bin_re_err,
                   fmt=linestyle, color=color, linewidth=2.5,
                   capsize=2.5, markersize=13, elinewidth=0., label=label)

    ax[1].errorbar(uvbins, bin_im,
                   yerr=bin_im_err,
                   fmt=linestyle, color=color, linewidth=2.5,
                   capsize=2.5, markersize=13, elinewidth=0., label=label)


    ax[1].set_xlabel(r'uv-distance (k$\lambda$)', fontweight='bold',
                     fontsize=fontsize)
    ax[0].set_ylabel('Re(V) (Jy)', fontweight='bold', fontsize=fontsize)
    ax[1].set_ylabel('Im(V) (Jy)', fontweight='bold', fontsize=fontsize)

    ax[0].set_xticklabels("")
    # ax[0].set_xlim(uvlim)
    # ax[1].set_xlim(uvlim)
    # ax[0].set_ylim(Jylims[0])
    # ax[1].set_ylim(Jylims[1])
    # ax[1].set_yticks(Jyticks_im)
    # ax[1].set_yticklabels(Jyticks_im)

    if label:
        ax[0].legend(fontsize=fontsize, frameon=False)

    if caption:
        caption_keys = ['x', 'y', 'text', 'fontsize']
        if all (k in caption for k in caption_keys):
            pass
        else:
            raise KeyError("Expected caption dict with following keys: {}, but got {}".format(caption_keys, caption.keys()))

        ax[0].text(caption['x'], caption['y'], caption['text'],
                   fontsize=caption['fontsize'], fontweight='bold')

    fig.subplots_adjust(left=0.2, right=0.95, hspace=0., top=0.95)
    fig.subplots_adjust(hspace=0)
    # gs.update(left=0.15, right=0.95, top=0.9, bottom=0.1,hspace=0) #, vspace=0.15)

    plt.savefig("uvplot.png")
