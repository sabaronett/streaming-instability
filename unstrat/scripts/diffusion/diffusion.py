#!/usr/bin/env python3
#==============================================================================
# diffusion.py
#
# Analyze the diffusion of each dust species[ and make relevant plots].
#
# Author: Stanley A. Baronett, Chao-Chin Yang
# Created: 2021-05-10
# Last Modified: 2022-06-28
#==============================================================================
import sys
sys.path.insert(0, '/home6/sbaronet/athena-dust/vis/python')
import athena_read
from dust import Diffusion, MultiSpecies
import math
# import matplotlib as mpl
# mpl.use("PDF")
# from matplotlib import cm
# from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib.colors import LogNorm
# from matplotlib.offsetbox import AnchoredText
# import matplotlib.pyplot as plt
import numpy as np
# import PencilCode as pc
import random
# from utilities import compose_title

# Get the stopping times, solid-to-gas ratios, and no. of outputs to process.
athinput = athena_read.athinput('athinput.si')
# tausmin, tausmax = 10**par.logtausmin, 10**par.logtausmax
# ms = MultiSpecies(par=par)
taus = athinput['problem']['omega'] * athinput['particles']['taus0']
epsilon = athinput['problem']['epsilon']
nout = 3 # Minimum number needed to calculate the diffusion coefficient.
if sys.argv[1]: nout = sys.argv[1]

# Find the displacement and the diffusion coefficient.
diff = Diffusion(athinput=athinput, nout=nout)
dpx, dpy, dpz = diff.coefficient()
np.savez(f'output/dcoeff-{nout}', taus=taus, dpx=dpx, dpy=dpy, dpz=dpz)

# # Select several snapshots.
# ntimes = len(diff.time)
# itimes = range(ntimes - 1, 0, -(ntimes - 1) // 4)
# tlabel = ["$t = {:.3G}P$".format(diff.time[i]) for i in itimes]

# # Compose a title.
# title = compose_title(par, nsp=ms.nsp, two_lines=True)
# dim = pc.read.dimensions()
# title += r" --- ${}\times{}$ grid".format(dim.nxgrid, dim.nzgrid)

# # Adjust Matplotlib style.
# plt.style.use("mine")
# mpl.rcParams["axes.autolimit_mode"] = "round_numbers"

# # Define function to initiate a figure.
# def initiate_figure():
#     fig = plt.figure(constrained_layout=True)
#     return fig, fig.gca()

# # Define function to finalize a figure.
# def finalize_figure(pdf, fig, ax, xlabel, ylabel, s=None):
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     ax.grid()
#     ax.set_title(title, size="small", pad=12)
#     if s is not None:
#         text = r"$\tau_{{\mathrm{{s}},{}}} = {:.3g}$".format(s+1, taus[s])
#         text += '\n' + r"$\epsilon_{{{}}} = {:.3g}$".format(s+1, epsilon[s])
#         at = AnchoredText(text, loc="upper left", prop=dict(linespacing=2))
#         ax.add_artist(at)
#     pdf.savefig(fig)
#     plt.close(fig)

# # Define function to make a colorbar.
# def colorbar(fig):
#     norm = LogNorm(tausmin, tausmax)
#     cmap = cm.get_cmap(lut=ms.nsp)
#     sm = cm.ScalarMappable(norm=norm, cmap=cmap)
#     cb = fig.colorbar(sm)
#     cb.set_label(r"$\tau_\mathrm{s}$")
#     return cmap.colors

# # Define function to process one of the dimensions.
# tf = np.linspace(0, diff.time[-1], 5 * (len(diff.time) - 1) + 1)
# def process_dim(dxp, mean, std, dpx, dim, pdf):

#     # Plot mean displacement for each species.
#     fig, ax = initiate_figure()
#     color = colorbar(fig)
#     for m, c in zip(mean, color):
#         ax.plot(diff.time, m, color=c, linewidth=0.5)
#     ax.set_xlim(0, diff.time[-1])
#     dlabel = rf"$\overline{{\Delta {dim}_\mathrm{{p}}}} / H_\mathrm{{g}}$"
#     finalize_figure(pdf, fig, ax, r"$t / P$", dlabel)

#     # Plot standard deviation.
#     fig, ax = initiate_figure()
#     color = colorbar(fig)
#     for s, d, c in zip(std, dpx, color):
#         ax.plot(tf, np.sqrt(2 * d * tf), color=c, linewidth=0.5)
#         ax.plot(diff.time, s, '.', color=c, markersize=1)
#     ax.set_xlim(0, diff.time[-1])
#     dlabel = rf"$\mathrm{{std}}(\Delta {dim}_\mathrm{{p}}) / H_\mathrm{{g}}$"
#     finalize_figure(pdf, fig, ax, r"$t / P$", dlabel)

#     # Plot the diffusion coefficient.
#     fig, ax = initiate_figure()
#     ax.semilogx(taus, dpx / (2 * math.pi), ".:")
#     ax.set_xlim(tausmin, tausmax)
#     tauslabel = r"\tau_\mathrm{s}"
#     dlabel = rf"$D_{{\mathrm{{p}},{dim}}}(" + tauslabel + ')'
#     dlabel += r"/ c_\mathrm{s} H_\mathrm{g}$"
#     tauslabel = '$' + tauslabel + '$'
#     finalize_figure(pdf, fig, ax, tauslabel, dlabel)

#     dlabel = rf"$\Delta {dim}_\mathrm{{p}} / H_\mathrm{{g}}$"
#     for s in range(ms.nsp):
#         indices = ms.get_indices_of_species(s)

#         # Plot trajectories of randomly chosen particles.
#         fig, ax = initiate_figure()
#         ax.plot(diff.time, dxp[:,random.choices(indices, k=8)], ".-")
#         ax.set_xlim(0, diff.time[-1])
#         finalize_figure(pdf, fig, ax, r"$t / P$", dlabel, s)

#         # Make the histograms of the displacement.
#         fig, ax = initiate_figure()
#         dxps = [dxp[i][indices] for i in itimes]
#         n, bins, patches = ax.hist(dxps, bins=50, density=True,
#                 histtype="step", label=tlabel)
#         x = np.linspace(bins[0], bins[-1], 200)
#         for k, itime in enumerate(itimes):
#             mu, sigma = mean[s][itime], std[s][itime]
#             y = np.exp(-0.5 * ((x - mu) / sigma)**2)
#             y /= math.sqrt(2 * math.pi) * sigma
#             ax.plot(x, y, ':', color=patches[k][0].get_edgecolor())
#         ax.legend(loc="upper right")
#         finalize_figure(pdf, fig, ax, dlabel, "Probability Density", s)

# # Make the analyses and plots.
# with PdfPages("diffx.pdf") as pdfx, PdfPages("diffz.pdf") as pdfz:
#     process_dim(diff.dxp, diff.dxp_mean, diff.dxp_std, dpx, 'x', pdfx)
#     process_dim(diff.dzp, diff.dzp_mean, diff.dzp_std, dpz, 'z', pdfz)
