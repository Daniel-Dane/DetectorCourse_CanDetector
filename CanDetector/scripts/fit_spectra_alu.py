#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  31 08:19:00 2018

@author: dnielsen
"""

import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from common import show_title, show_text
import matplotlib.pyplot as plt
np.random.seed(42)
import ROOT
from array import array
from fit_spectra_common import get_draw_spline, subtract_bkg, fit_and_draw_ROOT, energyall, \
                               gauss_single, gauss_double_uncorr, gauss_p1, \
                               fe_escape_energy, fe_main_energy, fe_sec_energy, am_main_energy, \
                               fe_escape_energy_unc, fe_main_energy_unc, fe_sec_energy_unc, am_main_energy_unc
import common
common.show_title_subtitle = "(Al pipe)"



#%%#####################################
# Show data with bkg subtracted 
######################################

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

# make zoomed in sub-figure
axins = zoomed_inset_axes(ax, 5.0, loc=5)

# Get histograms and make+draw splines (normalizing to fe)
#[h_fe,  spline_fe, time_fe]  = get_draw_spline("../data/mcapipe/fe_10_1224_spec_cop.mca",  0.02,  'b', 'b', "Fe-55",      ax, axins, False)
#[h_am,  spline_am, _]  = get_draw_spline("../data/mcapipe/am_10_1224_spec_cop.mca",  0.02,  'r', 'r', "Am-241",     ax, axins, True, time_fe)
#[h_bkg, spline_bkg, _] = get_draw_spline("../data/mcapipe/bkg_10_1224_spec_cop.mca", 0.002, 'g', 'g', "Background", ax, axins, True, time_fe)
[h_fe,  spline_fe, time_fe]  = get_draw_spline("../data/mcapipe/fe_10_1217_spec_alu.mca",  0.02,  'b', 'b', "Fe-55",      ax, axins, False)
[h_am,  spline_am, _]  = get_draw_spline("../data/mcapipe/am_10_1217_spec_alu.mca",  0.02,  'r', 'r', "Am-241",     ax, axins, True, time_fe)
[h_bkg, spline_bkg, _] = get_draw_spline("../data/mcapipe/bkg_10_1217_spec_alu.mca", 0.002, 'g', 'g', "Background", ax, axins, True, time_fe)

# subtract backround and draw
h_fe_new = subtract_bkg( h_fe, spline_fe, spline_bkg, "b", ax, axins )
h_am_new = subtract_bkg( h_am, spline_am, spline_bkg, "r", ax, axins )

# finished zoomed in sub-figure
axins.set_xlim(0, 149)
axins.set_ylim(0, 99)
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

# spice it up and show
x=0.05
ax.set_ylim(top=1.2*ax.get_ylim()[1])
show_title(ax, x=x)
show_text("Dashed lines: Bsplines of original histograms", ax, y=0.85, x=x)
show_text("Full lines: Bsplined background subtracted", ax, y=0.80, x=x)
ax.set_ylabel("Counts for {:.0f} seconds per channel [1/s/bit]".format(time_fe))
#ax.set_ylabel("Counts per second per 4 channels [1/s/bit]")
ax.set_xlabel("Channel [bit]")
ax.legend(loc='best')
fig.show()
plt.savefig("../graphics/alubkgsubtraction.pdf", format='pdf')

# used later on for calculating the number of events in gauss
binwidth = h_am_new.GetBinWidth(1)



#%%#####################################
# Fit spectra 
######################################

# make figure and axes
f, (ax, ax2) = plt.subplots(1, 2)
f.set_figwidth(7)

# plot points and fit result
[fe_esc_c, fe_esc_mean, fe_esc_sigma], [_, fe_esc_unc, _], fe_esc_chi2, fe_esc_ndof, fe_esc_prob = fit_and_draw_ROOT( h_fe_new, gauss_single, [961, 53.39, 6.53], ax, [50,70])#, bounds=([0,45,2],[10_000,57,14]) )

[N, r, fe_mean, fe_sigma, fe_sec_mean, fe_sec_sigma], [_, _, fe_unc, _, fe_sec_unc, _], fe_chi2, fe_ndof, fe_prob = fit_and_draw_ROOT( h_fe_new, gauss_double_uncorr, [12755, 0.89, 95.64, 6.83, 107.2, 4.723], ax, [85,115], draw_individually=True, dont_plot_hist=True)# good bounds: , bounds=([5_000,0.5,85,2,102,1],[20_000,1,105,11,113,8]) )

#[_, fe_mean, fe_sigma], [_, fe_unc, _], fe_chi2, fe_ndof, fe_prob = fit_and_draw_ROOT( h_fe_new, gauss_single, [49022, 89, 10], ax, [80,95], bounds=([10,85,1],[100000,95,15]), fitoptions="RSLL" )

[am_c, am_mean, am_sigma, _, _], [_, am_unc, _, _, _], am_chi2, am_ndof, am_prob = fit_and_draw_ROOT( h_am_new, gauss_p1, [214.2, 939.23, 11.05, 8.437, -0.00133], None, [925, 960], ax2=ax2)#, bounds=([10,930,2,0,-100],[10_000,950,15,100,100]) )

#def crystalball_function(x, alpha, n, mu, sigma):
#    return ROOT.Math.crystalball_function(alpha, n, sigma, mu)
#[alpha, _, am_mean, _], [alpha_unc, _, am_unc, _], am_chi2, am_ndof, am_prob = fit_and_draw_ROOT( h_am_new, crystalball_function, [2, 4896, 883, 22], None, [850, 950], ax2=ax2, bounds=([0,10,875,20],[100,10000,885,28]) )

ax.set_xlim(0,149)
ax2.set_xlim(901,1024)

ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.yaxis.tick_right()
ax2.tick_params(labelleft='off')
ax2.tick_params(labelright='on')

d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1-d,1+d), (-d,+d), **kwargs)
ax.plot((1-d,1+d),(1-d,1+d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d), (1-d,1+d), **kwargs)
ax2.plot((-d,+d), (-d,+d), **kwargs)

# spice it up and show
x=0.9
ax2.set_ylim(top=30)
show_title(ax, subtitle_x_offset=0.4)
show_text("Fe-55", ax, y=0.1)
show_text("Am-241", ax2, y=0.1, x=0.7)
show_text("Esc. peak:  {: 5.1f} ± {:.1f} bit ({:.0f}/{:d}=  {:.0f}%)".format(fe_esc_mean,fe_esc_unc,fe_esc_chi2,fe_esc_ndof,fe_esc_prob), ax2, y=0.90, x=x, ha="right")
show_text("Fe K-α peak:  {: 5.1f} ± {:.1f} bit ({:.0f}/{:d}=  {:.0f}%)".format(fe_mean,fe_unc,fe_chi2,fe_ndof,fe_prob), ax2, y=0.85, x=x, ha="right")
show_text("Fe K-β peak:  {: 5.1f} ± {:.1f} bit ({:.0f}/{:d}=  {:.0f}%)".format(fe_sec_mean,fe_sec_unc,fe_chi2,fe_ndof,fe_prob), ax2, y=0.80, x=x, ha="right")
show_text("Am peak: {:05.1f} ± {:.1f} bit ({:.0f}/{:d}=  {:.0f}%)".format(am_mean,am_unc,am_chi2,am_ndof,am_prob), ax2, y=0.75, x=x, ha="right")
#ax.set_ylabel("Counts per second per 4 channels [1/s/bit]")
ax.set_ylabel("Counts for {:.0f} seconds per channel [1/s/bit]".format(time_fe))
ax2.yaxis.set_label_position("right")
ax2.yaxis.labelpad = 10
#ax2.set_ylabel("Counts per second per 4 channels [1/s/bit]")
ax2.set_ylabel("Counts for {:.0f} seconds per channel [1/s/bit]".format(time_fe))
ax.set_xlabel("Channel [bit]")
ax2.set_xlabel("Channel [bit]")
f.show()
plt.savefig("../graphics/aluchannelfits.pdf", format='pdf')



#%%#####################################
# Make energy-channel calibration (linear fit with ROOT)
######################################

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

#fe_esc_unc = fe_esc_sigma/np.sqrt(fe_esc_c/binwidth)
#fe_unc = fe_sigma/np.sqrt(N*r/binwidth)
#fe_sec_unc = fe_sec_sigma/np.sqrt(N*(1-r)/binwidth)
#am_unc = am_sigma/np.sqrt(am_c/binwidth)

# do fit
y = [ fe_escape_energy, fe_main_energy, fe_sec_energy, am_main_energy ]
x = [ fe_esc_mean, fe_mean, fe_sec_mean, am_mean ]
yerr = [ fe_escape_energy_unc, fe_main_energy_unc, fe_sec_energy_unc, am_main_energy_unc ]
xerr = np.array([ fe_esc_unc, fe_unc, fe_sec_unc, am_unc ])
gr = ROOT.TGraphErrors( len(x), array('d',x), array('d',y), array('d',xerr), array('d',yerr) )
fit1 = ROOT.TF1("fit1","pol1", min(x), max(x));
fit1.SetParameters(-0.2232671292611002, 0.06796181642128599)
res = gr.Fit(fit1, "RS")
print("Fit prob. = {:.1f}%".format(fit1.GetProb()*100))

# plot points and fit result
plt.errorbar(x=x, xerr=10*xerr, y=y, yerr=yerr, fmt="none", color='r',elinewidth=3, zorder=10)
x = np.linspace(min(x)-10, max(x)+10, 1000)
y = [fit1.Eval(x) for x in x]
plt.plot(x, y, 'b-')

# spice it up and show
show_title(ax)
show_text("Note: Channel uncertainties scaled by 10", ax, y=0.85)
show_text("Fit: y = b + ax", ax, y=0.80)
show_text("b = {:5.3f}     ± {:.3f} ({:3.1f}σ from 0)".format(fit1.GetParameter(0),fit1.GetParError(0),abs(fit1.GetParameter(0))/fit1.GetParError(0)), ax, y=0.75)
show_text("a = {: 7.5f} ± {:.5f}".format(fit1.GetParameter(1),fit1.GetParError(1)), ax, y=0.70)
show_text("p( X², ndof ) = p( {:.1f}, {:d} ) = {:.1f}%".format(fit1.GetChisquare(), fit1.GetNDF(), fit1.GetProb()*100), ax, y=0.65)
ax.set_ylabel("Energy [keV]")
ax.set_xlabel("Channel [bit]")
#ax.arrow(10,10,50,-2.7,width=0.5,head_length=15)
fig.show()
plt.savefig("../graphics/aluenergychannelcalib.pdf", format='pdf')



#%%#####################################
# Find additional peaks in Am
######################################

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

#[ am1_c, am1_mean, am1_sigma, am1_p0, am1_p1 ], [ am1_c_unc, am1_mean_unc, am1_sigma_unc, am1_p0_unc, am1_p1_unc ], \
#am1_chi2, am1_ndof, am1_prob = fit_and_draw_ROOT(h_am_new, gauss_p1, [161, 261, 11, 20, -0.03], ax, [240,285], col='g-', dont_draw_fit=True)
#
#[ am2_c, am2_mean, am2_sigma, am2_p0, am2_p1 ], [ am2_c_unc, am2_mean_unc, am2_sigma_unc, am2_p0_unc, am2_p1_unc ], \
#am2_chi2, am2_ndof, am2_prob = fit_and_draw_ROOT(h_am_new, gauss_p1, [58, 318, 7, 20, -0.03], ax, [300,335], True, \
#                                                 bounds=([10, 310, 2, 10, -0.2],[100, 325, 11, 10000, -0.001]), col='g-'
#                                                 )

[ am3_c, am3_mean, am3_sigma, am3_p0, am3_p1 ], [ am3_c_unc, am3_mean_unc, am3_sigma_unc, am3_p0_unc, am3_p1_unc ], \
am3_chi3, am3_ndof, am3_prob = fit_and_draw_ROOT(h_am_new, gauss_p1, [447, 418, 10, 20, -0.03], ax, [385,445], False, col='g-', label="Am-241")

#[ am4_c, am4_mean, am4_sigma, am4_p0, am4_p1 ], [ am4_c_unc, am4_mean_unc, am4_sigma_unc, am4_p0_unc, am4_p1_unc ], \
#am4_chi4, am4_ndof, am4_prob = fit_and_draw_ROOT(h_am_new, gauss_p1, [2000, 203, 10, 31, -0.03], ax, [190,210], True, \
#                                                 bounds=([1000,200,3,1,-2],[10000,210,25,10000,-0.001]), col='g-'
#                                                 )

#lower_bound=0.8
#upper_bound=1.2
#
#p0=20
#p1=-0.001
#
#[ am1_c, am1_mean, am1_sigma, am2_c, am2_mean, am2_sigma, \
#  am3_c, am3_mean, am3_sigma, am4_c, am4_mean, am4_sigma, \
#  am_final_p0, am_final_p1, am_final_p2 ], \
#[ am1_c_unc, am1_mean_unc, am1_sigma_unc, \
#  am2_c_unc, am2_mean_unc, am2_sigma_unc, \
#  am3_c_unc, am3_mean_unc, am3_sigma_unc, \
#  am4_c_unc, am4_mean_unc, am4_sigma_unc, \
#  am_final_p0_unc, am_final_p1_unc, am_final_p2_unc ], \
#am_final_chi4, am_final_ndof, am_final_prob = \
#    fit_and_draw_ROOT( h_am_new, gauss_quad_p2, \
#                       [ am1_c, am1_mean, am1_sigma, am2_c, am2_mean, am2_sigma, \
#                         am3_c, am3_mean, am3_sigma, am4_c, am4_mean, am4_sigma, \
#                         20, -0.001, 0
#                         ], \
#                       ax, [175,450], True, \
#                       draw_individually=True, \
#                       bounds=( [ 60,                am1_mean*lower_bound, am1_sigma*lower_bound, \
#                                  am2_c*lower_bound, am2_mean*lower_bound, am2_sigma*lower_bound, \
#                                  am3_c*lower_bound, am3_mean*lower_bound, am3_sigma*lower_bound, \
#                                  am4_c*lower_bound, am4_mean*lower_bound, am4_sigma*lower_bound, \
#                                  0.0, -1.0, 0.000000001 ],
#                                [ 1000,              am1_mean*upper_bound, am1_sigma*upper_bound, \
#                                  am2_c*upper_bound*100, am2_mean*upper_bound, am2_sigma*upper_bound, \
#                                  am3_c*upper_bound, am3_mean*upper_bound, am3_sigma*upper_bound, \
#                                  am4_c*upper_bound, am4_mean*upper_bound, am4_sigma*upper_bound, \
#                                  100, -0.001, 0.000000001 ]
#                               )
#                      )

# we fit with a gauss on top of a flat background
# the uncertainty on the mean is then the width divided by the square root of the number of entries
# the normalization constant divided by the binwidth gives us exactly the number of entries
# this elaborate exercise gives us the actual uncertainty of the mean for just the signal/gauss
#am1_energy = energyall(fit1, am1_mean, am1_sigma/np.sqrt(am1_c/binwidth))
#am2_energy = energyall(fit1, am2_mean, am2_sigma/np.sqrt(am2_c/binwidth))
am3_energy = energyall(fit1, am3_mean, am3_sigma/np.sqrt(am3_c/binwidth))
#am4_energy = energyall(fit1, am4_mean, am4_sigma/np.sqrt(am4_c/binwidth))

# spice it up and show
x=0.011
ax.set_ylim(top=1.3*ax.get_ylim()[1])
show_title(ax, x=x, y=0.92)
#show_text("Peak 1: E = {:.3f} ± {:.3f} (stat.) ± {:.3f} (cal.) ± {:.3f} (syst.) keV".format(*am4_energy, 0), ax, y=0.87, x=x)
#show_text("Peak 2: E = {:.3f} ± {:.3f} (stat.) ± {:.3f} (cal.) ± {:.3f} (syst.) keV".format(*am1_energy, 0), ax, y=0.82, x=x)
#show_text("Peak 3: E = {:.3f} ± {:.3f} (stat.) ± {:.3f} (cal.) ± {:.3f} (syst.) keV".format(*am2_energy, 0), ax, y=0.77, x=x)
show_text("Peak 4: E = {:.3f} ± {:.3f} (stat.) ± {:.3f} (cal.) ± {:.3f} (syst.) keV".format(*am3_energy, 0), ax, y=0.72, x=x)
ax.set_ylabel("Counts for {:.0f} seconds per channel [1/s/bit]".format(time_fe))
ax.set_xlabel("Channel [bit]")
ax.legend(loc='best')
fig.show()
plt.savefig("../graphics/alupeaksearch.pdf", format='pdf')

# if using a terminal
#input("ready...")