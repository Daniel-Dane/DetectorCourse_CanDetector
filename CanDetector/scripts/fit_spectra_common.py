#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  31 08:19:00 2018

@author: dnielsen
"""

import numpy as np
import scipy.interpolate as interpolate
from common import mca_to_hist
import rootpy.plotting.root2matplotlib as rplt
np.random.seed(42)
#from scipy.optimize import curve_fit
import ROOT
from inspect import signature



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
#    hist.Rebin(4)
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
#        rplt.hist(hist, stacked=False, fill=False, axes=axins)
#    plt.grid()
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



#%%#####################################
# Fit spectra 
######################################

def gauss_single(x, c0, m0, s0):
    return c0/(np.sqrt(2*np.pi)*s0)*np.exp(-(x-m0)**2/(2*s0**2))
def gauss_double(x, c0, m0, s0, c1, m1, s1):
    return c0/(np.sqrt(2*np.pi)*s0)*np.exp(-(x-m0)**2/(2*s0**2)) + \
           c1/(np.sqrt(2*np.pi)*s1)*np.exp(-(x-m1)**2/(2*s1**2))
def gauss_triple(x, c0, m0, s0, c1, m1, s1, c2, m2, s2):
    return c0/(np.sqrt(2*np.pi)*s0)*np.exp(-(x-m0)**2/(2*s0**2)) + \
           c1/(np.sqrt(2*np.pi)*s1)*np.exp(-(x-m1)**2/(2*s1**2)) + \
           c2/(np.sqrt(2*np.pi)*s2)*np.exp(-(x-m2)**2/(2*s2**2))
def gauss_quad(x, c0, m0, s0, c1, m1, s1, c2, m2, s2, c3, m3, s3):
    return c0/(np.sqrt(2*np.pi)*s0)*np.exp(-(x-m0)**2/(2*s0**2)) + \
           c1/(np.sqrt(2*np.pi)*s1)*np.exp(-(x-m1)**2/(2*s1**2)) + \
           c2/(np.sqrt(2*np.pi)*s2)*np.exp(-(x-m2)**2/(2*s2**2)) + \
           c3/(np.sqrt(2*np.pi)*s3)*np.exp(-(x-m3)**2/(2*s3**2))

def gauss_p0(x, c0, m0, s0, p0):
    return gauss_single(x, c0, m0, s0) + p0

def gauss_plus_exp(x, c0, m0, s0, p0, t0):
    return gauss_single(x, c0, m0, s0) + p0*np.exp(-x*t0)

def gauss_p1(x, c0, m0, s0, p0, p1):
    return gauss_single(x, c0, m0, s0) + p0 + p1*x

def gauss_p2(x, c0, m0, s0, p0, p1, p2):
    return gauss_single(x, c0, m0, s0) + p0 + p1*x + p2*x**2

def gauss_triple_p1(x, c0, m0, s0, c1, m1, s1, c2, m2, s2, p0, p1):
    return gauss_triple(x, c0, m0, s0, c1, m1, s1, c2, m2, s2) + p0 + p1*x

def gauss_quad_p0(x, c0, m0, s0, c1, m1, s1, c2, m2, s2, c3, m3, s3, p0):
    return gauss_quad(x, c0, m0, s0, c1, m1, s1, c2, m2, s2, c3, m3, s3) + p0

def gauss_quad_p1(x, c0, m0, s0, c1, m1, s1, c2, m2, s2, c3, m3, s3, p0, p1):
    return gauss_quad(x, c0, m0, s0, c1, m1, s1, c2, m2, s2, c3, m3, s3) + p0 + p1*x

def gauss_quad_p2(x, c0, m0, s0, c1, m1, s1, c2, m2, s2, c3, m3, s3, p0, p1, p2):
    return gauss_quad(x, c0, m0, s0, c1, m1, s1, c2, m2, s2, c3, m3, s3) + p0 + p1*x + p2*x**2

def gauss_double_uncorr(x, N, r, m0, s0, m1, s1):
    return N*(r/(np.sqrt(2*np.pi)*s0)*np.exp(-(x-m0)**2/(2*s0**2)) + (1-r)/(np.sqrt(2*np.pi)*s1)*np.exp(-(x-m1)**2/(2*s1**2)))

# energy of peaks in keV
#fe_escape_energy = 2.96 # 60/76*2.958+16/76*2.956
fe_escape_energy = 3.19
fe_main_energy = 5.89
fe_sec_energy = 6.49045
am_main_energy = 59.5409

fe_escape_energy_unc = 0.01
fe_main_energy_unc = 0.00001
fe_sec_energy_unc = 0.01 # higher uncertainty due to being two peaks
am_main_energy_unc = 0.0001



######################################
# Fitting with ROOT (FINALLY WORKS!) 
######################################

fitWithROOT_counter = 0
def fit_with_ROOT(hist_orig, fitobject, fitoptions="RS"):
    global fitWithROOT_counter
    fitWithROOT_counter += 1
    hist = ROOT.TH1D("h"+str(fitWithROOT_counter), "h"+str(fitWithROOT_counter), hist_orig.GetNbinsX(), hist_orig.GetXaxis().GetXmin(), hist_orig.GetXaxis().GetXmax() )
    for i in range(1,hist_orig.GetNbinsX()+1):
        hist.SetBinContent(i,hist_orig.GetBinContent(i))
        hist.SetBinError(i,np.sqrt(hist_orig.GetBinContent(i)))
    hist.Fit(fitobject, fitoptions)
    return fitobject.GetChisquare(), fitobject.GetNDF(), fitobject.GetProb()*100

makeFitObject_counter = 0
def make_fit_object(funcOrExpr, xmin, xmax):
    global makeFitObject_counter
    makeFitObject_counter += 1
    if callable(funcOrExpr) :
        f = lambda x, pars : funcOrExpr(x[0], *pars)
        return ROOT.TF1("fit"+str(makeFitObject_counter), f, xmin, xmax, len(signature(funcOrExpr).parameters)-1)
    else:
        return ROOT.TF1("fit"+str(makeFitObject_counter), funcOrExpr, xmin, xmax, 4)

def fit_and_draw_ROOT(hist, func, startval, ax, xrange=None, dont_plot_hist=False, ax2=None, return_pcov=False, draw_individually=False, bounds=None, col='b-', dont_draw_fit=False, fitoptions="RS", label=None):
    if xrange is None:
        xrange = [0, 1024]
    fitobject = make_fit_object(func, xrange[0], xrange[1])
    if len(startval) < 10:
        fitobject.SetParameters(*startval)
    else:
        for i, val in enumerate(startval):
            fitobject.SetParameter(i, val)
    if bounds is not None:
        for i in range(0, len(bounds[0])):
            fitobject.SetParLimits(i, bounds[0][i], bounds[1][i])
    chi2, ndof, prob = fit_with_ROOT(hist, fitobject, fitoptions)
    if callable(func) :
        funcrange = range(0,len(signature(func).parameters)-1)
    else:
        funcrange = range(4)
    pars = [fitobject.GetParameter(i) for i in funcrange]
    errs = [fitobject.GetParError(i) for i in funcrange]
    
    x = np.linspace(xrange[0], xrange[1], 1000)
    if not dont_plot_hist:
        if ax is not None:
            rplt.hist(hist, color='r', axes=ax, label=label)
        if ax2 is not None:
            rplt.hist(hist, color='r', axes=ax2)
    if not dont_draw_fit:
        if ax is not None:
            ax.plot(x,func(x, *pars), col, zorder=10)
#           print(pars)
        if ax2 is not None:
            ax2.plot(x,func(x, *pars), col, zorder=10)
    if draw_individually:
        if func.__name__ == "gauss_double_uncorr":
            ax.plot(x, gauss_single(x, pars[0]*pars[1], pars[2], pars[3]), 'r--')
            ax.plot(x, gauss_single(x, pars[0]*(1-pars[1]), pars[4], pars[5]), 'g--')
        else:
            # Assumes an arbitrary number of gauss followed by a possible n-dim polynomial
            n = 0
            a = dict(signature(func).parameters)
            while True:
                if 's'+str(n) in a:
                    n+=1
                else:
                    break
            n*=3
            cols=['m:','g:','b:','k:']
            for i in range(0,n,3):
                ax.plot(x, gauss_single(x, pars[i], pars[i+1], pars[i+2]), cols[i//3])
            n = -1
            while True:
                if 'p'+str(n+1) in a:
                    n+=1
                else:
                    break
            if n > 0:
                y = 0
                for i in range(n+1):
                    y += pars[-1-n+i]*x**i
                ax.plot(x, y, 'g--')
    return pars, errs, chi2, ndof, prob



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

#def fit(hist, func, startval, ax, xrange=None, dont_plot_hist=False, ax2=None, return_pcov=False, draw_individually=False, bounds=(-np.inf,np.inf)):
#    x = np.array([hist.GetBinCenter(x) for x in range(1,hist.GetNbinsX())])
#    y = np.array([hist.GetBinContent(x) for x in range(1,len(x)+1)])
#    if xrange is not None:
#        y = y[(x>xrange[0])&(x<xrange[1])]
#        x = x[(x>xrange[0])&(x<xrange[1])]
#    popt, pcov = curve_fit(func, x, y, p0 = startval, method='trf', bounds=bounds)
#    if not dont_plot_hist:
#        if ax is not None:
#            rplt.hist(hist, color='r', axes=ax)
#        if ax2 is not None:
#            rplt.hist(hist, color='r', axes=ax2)
#    if ax is not None:
#        ax.plot(x,func(x, *popt), 'b-', zorder=10)
##        print(popt)
#    if ax2 is not None:
#        ax2.plot(x,func(x, *popt), 'b-')
#    if draw_individually:
#        ax.plot(x,gauss_single(x, popt[0]*popt[1], popt[2], popt[3]), 'r--')
#        ax.plot(x,gauss_single(x, popt[0]*(1-popt[1]), popt[4], popt[5]), 'g--')
##    print(pcov)
#    if return_pcov:
#        return [popt,np.sqrt(pcov.diagonal())]
#    else:
#        return popt



#%%#####################################
# Find additional peaks in Am
######################################

# E = p0 + p1*channel
# does error propagation
def energywithuncertainty(fitobj, val, dval, caluncoff=False):
    a = fitobj.GetParameter(1)
    b = fitobj.GetParameter(0)
    da = fitobj.GetParError(1)
    db = fitobj.GetParError(0)
    if caluncoff:
        da=0
        db=0
    e = a*val+b
    de = np.sqrt( (a*val)**2*((da/a)**2+(dval/val)**2) + db**2 )
    return [e, de]

def energyall(fitobj, val, dval):
    return [energywithuncertainty(fitobj, val, dval)[0],
            energywithuncertainty(fitobj, val, dval, True)[1],
            energywithuncertainty(fitobj, val, 0)[1]
            ]
