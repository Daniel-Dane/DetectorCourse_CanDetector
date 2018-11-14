#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  31 08:19:00 2018

@author: dnielsen
"""

import numpy as np
import scipy.interpolate as interpolate
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from common import mca_to_hist, show_title, show_text
#import sys
#import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm

#from os import path
#sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

#from rootpy.plotting import Hist
#from ROOT import TF1, kRed
#from matplotlib import gridspec
#import matplotlib
#matplotlib.use('Agg')

"""
all of below need to be loaded after `matplotlib` and `matplotlib.use('Agg')`
"""
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt

np.random.seed(42)


######################################
# Fit spectra 
######################################

def get_draw_spline( fname, smoothing_strength=0.002, color_hist='k', color_spline='y', label="unlabeled", ax=None, axins=None ) :
    (hist, r_min, r_max, time) = mca_to_hist(fname, False)
    hist.Scale(1/time)
    hist.Rebin(4)
#    rplt.hist(hist, stacked=False, fill=False, axes=ax)
#    print(time)
    x = [hist.GetBinCenter(x) for x in range(0,hist.GetNbinsX())]
    y = [hist.GetBinContent(x) for x in range(0,len(x))]
    t, c, k = interpolate.splrep(x, y, s=smoothing_strength, k=3)
    xx = np.linspace(x[0], x[-1], 200)
    spline = interpolate.BSpline(t, c, k, extrapolate=False)
    ax.plot(xx, spline(xx), color_spline, label=label,zorder=10, linestyle="--")
    if axins is not None:
        axins.plot(xx, spline(xx), color_spline, label=label,zorder=10, linestyle="--")
#    plt.grid()
    ax.legend(loc='best')
    return [hist,spline]

def subtract_bkg( h_sig, sp_sig, sp_bkg, color='k', ax=None, axins=None ) :
    x = [ h_sig.GetBinCenter(x) for x in range(1, h_sig.GetNbinsX()) ]
    h_new = h_sig.Clone()
    for x in range( 0, h_sig.GetNbinsX() ) :
        h_new.SetBinContent( x, max(0,sp_sig(h_sig.GetBinCenter(x)) - sp_bkg(h_sig.GetBinCenter(x))) )
    rplt.hist(h_new, color=color, axes=ax)
    if axins is not None:
        rplt.hist(h_new, color=color, axes=axins)
    return h_new

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

# make zoomed in sub-figure
axins = zoomed_inset_axes(ax, 3.5, loc=5)

# Get histograms and make+draw splines
[h_fe,  spline_fe]  = get_draw_spline("../data/mca/fe_4_1937_spectrum.mca",  0.02,  'b', 'b', "Fe-55",      ax, axins)
[h_am,  spline_am]  = get_draw_spline("../data/mca/am_4_1937_spectrum.mca",  0.02,  'r', 'r', "Am-241",     ax, axins)
[h_bkg, spline_bkg] = get_draw_spline("../data/mca/bkg_4_1937_spectrum.mca", 0.002, 'g', 'g', "Background", ax, axins)

# subtract backround and draw
h_fe_new = subtract_bkg( h_fe, spline_fe, spline_bkg, "b", ax, axins )
h_am_new = subtract_bkg( h_am, spline_am, spline_bkg, "r", ax, axins )

# finished zoomed in sub-figure
axins.set_xlim(0, 150)
axins.set_ylim(0, 1.5)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

# spice it up and show
x=0.2
show_title(ax,x=x)
show_text("heyyyyyyy", ax, y=0.85,x=x)
ax.set_ylabel("Counts per second per 4 channels [1/s/bit]")
ax.set_xlabel("Channel [bit]")
fig.show()

# if using a terminal
#input("ready...")