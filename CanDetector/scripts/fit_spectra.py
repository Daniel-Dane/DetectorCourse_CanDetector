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
import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt
np.random.seed(42)
from scipy.optimize import curve_fit
import ROOT
from array import array



#%%#####################################
# Show data with bkg subtracted 
######################################

def get_draw_spline( fname, smoothing_strength=0.002, color_hist='k', color_spline='y', label="unlabeled", ax=None, axins=None, do_norm=True, time_to_norm_to=None ) :
    (hist, r_min, r_max, time_hist) = mca_to_hist(fname, False)
    if do_norm:
        if time_to_norm_to is None:
            hist.Scale(1/time_hist)
        else:
            hist.Scale(time_to_norm_to/time_hist)
    hist.Rebin(4)
#    rplt.hist(hist, stacked=False, fill=False, axes=ax)
#    print(time_hist)
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
    return [hist,spline,time_hist]

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

# Get histograms and make+draw splines (normalizing to fe)
#[h_fe,  spline_fe, time_fe]  = get_draw_spline("../data/mca/fe_4_1937_spectrum.mca",  0.02,  'b', 'b', "Fe-55",      ax, axins, False)
#[h_am,  spline_am, _]  = get_draw_spline("../data/mca/am_4_1937_spectrum.mca",  0.02,  'r', 'r', "Am-241",     ax, axins, True, time_fe)
#[h_bkg, spline_bkg, _] = get_draw_spline("../data/mca/bkg_4_1937_spectrum.mca", 0.002, 'g', 'g', "Background", ax, axins, True, time_fe)

# Get histograms and make+draw splines (normalizing to 1)
[h_fe,  spline_fe, time_fe]  = get_draw_spline("../data/mca/fe_4_1937_spectrum.mca",  0.02,  'b', 'b', "Fe-55",      ax, axins)
[h_am,  spline_am, _]  = get_draw_spline("../data/mca/am_4_1937_spectrum.mca",  0.02,  'r', 'r', "Am-241",     ax, axins)
[h_bkg, spline_bkg, _] = get_draw_spline("../data/mca/bkg_4_1937_spectrum.mca", 0.002, 'g', 'g', "Background", ax, axins)

# subtract backround and draw
h_fe_new = subtract_bkg( h_fe, spline_fe, spline_bkg, "b", ax, axins )
h_am_new = subtract_bkg( h_am, spline_am, spline_bkg, "r", ax, axins )

# finished zoomed in sub-figure
axins.set_xlim(0, 149)
#axins.set_ylim(0, 600)
axins.set_ylim(0, 1.5)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

# spice it up and show
x=0.2
show_title(ax,x=x)
show_text("heyyyyyyy", ax, y=0.85,x=x)
#ax.set_ylabel("Counts per {:.0f} seconds per 4 channels [1/s/bit]".format(time_fe))
ax.set_ylabel("Counts per second per 4 channels [1/s/bit]")
ax.set_xlabel("Channel [bit]")
fig.show()



#%%#####################################
# Fit spectra 
######################################

def gauss_single(x, c0, m0, s0):
    return c0*np.exp(-(x-m0)**2/(2*s0**2))
def gauss_double(x, c0, m0, s0, c1, m1, s1):
    return c0*np.exp(-(x-m0)**2/(2*s0**2)) + c1*np.exp(-(x-m1)**2/(2*s1**2))
def gauss_triple(x, c0, m0, s0, c1, m1, s1, c2, m2, s2):
    return c0*np.exp(-(x-m0)**2/(2*s0**2)) + c1*np.exp(-(x-m1)**2/(2*s1**2)) + c2*np.exp(-(x-m2)**2/(2*s2**2))
def gauss_single_root(x, par):
    return gauss_single(x[0], *par)
def gauss_double_root(x, par):
    return gauss_double(x[0], *par)
def gauss_triple_root(x, par):
    return gauss_triple(x[0], *par)

# energy of peaks in keV
fe_escape_energy = 2.96 # 60/76*2.958+16/76*2.956
fe_main_energy = 5.89 # 60/76*5.89875 + 16/76*5.88765
am_main_energy = 59.5409



######################################
# Fitting with ROOT (doesn't work) 
######################################
#import ROOT
#fitgausses = ROOT.TF1("fitgausses", "gaus+gaus(3)", 0, 1024, 6);
##fitgausses.SetParameters(2873.44120691, 90.35997591, 7.04385417)
#fitgausses.SetParameters(187.79413975, 46.42939315, 6.16071696, 2873.44120691, 90.35997591, 7.04385417)
#h_fe_new.Fit(fitgausses,"RLL")
##res = h_fe_new.Fit(fitgausses,"RS")
#print("Fit prob. = {:.1f}%".format(fitgausses.GetProb()*100))
##fig = plt.figure()
##ax = plt.subplot()
##rplt.hist(h_fe_new, color='r', axes=ax)
##fig.show()



######################################
# Fitting with RooFit (doesn't work) 
######################################

#from rootpy.stats import Workspace
#w = Workspace()
#w.factory('Gaussian::g(x[0,1024],mu[85,95],sigma[5,10])')
#w.factory("ExtendPdf:model(g,nevt[10,0,100000])")
##w.factory("x[0,200]")
#h = ROOT.TH1D("h","h",1024, 0, 1024)
#for i in range(0,1024):
#    h.SetBinContent(i+1,h_fe_new.GetBinContent(i+1))
#    h.SetBinError(i+1,h_fe_new.GetBinError(i+1))
##    if h.GetBinContent(i+1) > 0.0001 :
##        print(h.GetBinContent(i+1)-h_fe_new.GetBinContent(i+1))
#data = ROOT.RooDataHist("data","data",ROOT.RooArgList(w.var("x")),h)
#pdf = w.pdf('model')
##fitResult = pdf.fitTo(data,ROOT.RooFit.Save(),ROOT.RooFit.PrintLevel(-1),ROOT.RooFit.SumW2Error(ROOT.kTRUE))
#m = ROOT.RooMinuit(pdf.createNLL(data,ROOT.RooFit.NumCPU(20)))
#m.migrad()
##fitResult.Print()



######################################
# Fitting with scipy (works but no chisquare value)
######################################

def fit(hist, func, startval, ax, xrange=None, dont_plot_hist=False):
    x = np.array([hist.GetBinCenter(x) for x in range(0,hist.GetNbinsX())])
    y = np.array([hist.GetBinContent(x) for x in range(0,len(x))])
    if xrange is not None:
        y = y[(x>xrange[0])&(x<xrange[1])]
        x = x[(x>xrange[0])&(x<xrange[1])]
    popt, pcov = curve_fit(func, x, y, p0 = startval)
    if not dont_plot_hist:
        rplt.hist(hist, color='r', axes=ax)
    ax.plot(x,func(x, *popt),'b-', zorder=10)
    return popt

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

#x = [h_fe_new.GetBinCenter(x) for x in range(0,h_fe_new.GetNbinsX())]
#y = [h_fe_new.GetBinContent(x) for x in range(0,len(x))]
##popt, pcov = curve_fit(gauss_double, x, y, p0 = [0.62503264, 46.45198925,  6.1881488 ,  9.56303187, 90.35905582, 7.04286884])
#popt, pcov = curve_fit(gauss_triple, x, y, p0 = [12.94639985,   12.72625493,    3.11847216,  187.79332257, 46.4294929 ,    6.16076703, 2873.44121312,   90.35997589, 7.04385413])
#fig = plt.figure()
#ax = plt.subplot()
#rplt.hist(h_fe_new, color='r', axes=ax)
#plt.plot(x,gauss_triple(x, *popt))
#fig.show()
#c0, m0, s0, c1, fe_esc_mean, fe_esc_sigma, c2, fe_mean, fe_sigma = popt

c0, m0, s0, c1, fe_esc_mean, fe_esc_sigma, c2, fe_mean, fe_sigma = fit(h_fe_new, gauss_triple,[12.9, 12.7, 3.11, 187, 46.4, 6.16, 2873, 90.3, 7.04], ax)

fig.show()

#x = np.array([h_am_new.GetBinCenter(x) for x in range(0,h_am_new.GetNbinsX())])
#y = np.array([h_am_new.GetBinContent(x) for x in range(0,len(x))])
#x2 = x[(x>850)&(x<940)]
#y2 = y[(x>850)&(x<940)]
#popt, pcov = curve_fit(gauss_single, x2, y2, p0 = [1.03901811, 869.30521183,  38.45829851])
#fig = plt.figure()
#ax = plt.subplot()
#rplt.hist(h_am_new, color='r', axes=ax)
#plt.plot(x2,gauss_single(x2, *popt))
#fig.show()
#c0, am_mean, am_sigma = popt

# make figure and axes
#fig = plt.figure()
#ax = plt.subplot()

c0, am_mean, am_sigma = fit(h_am_new, gauss_single, [1.03, 869,  38.4], ax, [850,940])

fig.show()



######################################
# Make energy-channel calibration (linear fit with ROOT)
######################################

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

# do fit
y = [ fe_escape_energy, fe_main_energy, am_main_energy ]
x = [ fe_esc_mean, fe_mean, am_mean ]
yerr = [ 0, 0, 0 ]
xerr = [ fe_esc_sigma, fe_sigma, am_sigma ]
gr = ROOT.TGraphErrors( len(x), array('d',x), array('d',y), array('d',xerr), array('d',yerr) )
fit1 = ROOT.TF1("fit1","pol1", min(x), max(x));
fit1.SetParameters(-0.2232671292611002, 0.06796181642128599)
res = gr.Fit(fit1, "RS")
print("Fit prob. = {:.1f}%".format(fit1.GetProb()*100))

# plot points and fit result
plt.errorbar(x=x, xerr=xerr, y=y, yerr=yerr, fmt="none", color='r')
x = np.linspace(min(x), max(x), 1000)
y = [fit1.Eval(x) for x in x]
plt.plot(x, y, 'b-')

# spice it up and show
show_title(ax)
show_text("p( X², ndof ) = p( {:.2f}, {:d} ) = {:.1f}%".format(fit1.GetChisquare(), fit1.GetNDF(), fit1.GetProb()*100), ax, y=0.85)
show_text("y = {:.3f} + {:.3f}*x".format(fit1.GetParameter(0),fit1.GetParameter(1)), ax, y=0.80)
ax.set_ylabel("Energy [keV]")
ax.set_xlabel("Channel [bit]")
fig.show()



######################################
# Find additional peaks in Am
######################################

# transform channel to energy, neglecting uncertainties
#newbins = [ fit1.Eval(h_am_new.GetBinCenter(x)) for x in range(2,h_am_new.GetNbinsX()) ]
#h_am_new = ROOT.TH1D("bla","bla", len(newbins), array('d',newbins))

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

_, am1_mean, am1_sigma = fit(h_am_new, gauss_single, [1, 260, 10], ax, [245,275])
_, am2_mean, am2_sigma = fit(h_am_new, gauss_single, [1, 320, 10], ax, [300,335], True)
_, am3_mean, am3_sigma = fit(h_am_new, gauss_single, [1, 390, 10], ax, [370,415], True)

def evaluncertainty(val):
    x = array('d',[val])
    err = array('d',[0])
    res.GetConfidenceIntervals(1, 1, 1, x, err, 0.683, False)
    return err[0]

def energywithuncertainty(val):
    a = fit1.GetParameter(1)
    b = fit1.GetParameter(0)
    da = fit1.GetParError(1)
    db = fit1.GetParError(0)
    e = a*val+b
    de = np.sqrt((a*val*np.sqrt((da/a)**2+(db/b)**2))**2+db**2)
    return [e, de]

am1_energy = fit1.Eval(am1_mean)
am1_energy_unc = evaluncertainty(am1_mean)

am1_energy_sigma = fit1.Eval(am1_sigma)
am1_energy_sigma_unc = evaluncertainty(am1_sigma)

am1_energy = energywithuncertainty(am1_mean)
am2_energy = energywithuncertainty(am2_mean)
am3_energy = energywithuncertainty(am3_mean)

x=0.3
show_title(ax,x=x)
show_text("Channel {:.0f}: energy = {:.3f} +- {:.3f} keV".format(am1_mean,*am1_energy), ax, y=0.85, x=x)
#show_text("Peak at channel {:.0f}:\nenergy = {:.3f} +- {:.3f} keV,\nwidth = {:.3f} +- {:.3f} keV".format(am1_mean,am1_energy,am1_energy_unc,am1_energy_sigma,am1_energy_sigma_unc), ax, y=0.75, x=x)
show_text("Channel {:.0f}: energy = {:.3f} +- {:.3f} keV".format(am2_mean,*am2_energy), ax, y=0.80, x=x)
show_text("Channel {:.0f}: energy = {:.3f} +- {:.3f} keV".format(am3_mean,*am3_energy), ax, y=0.75, x=x)

fig.show()

# if using a terminal
#input("ready...")