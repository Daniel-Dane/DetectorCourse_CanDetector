#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  31 08:19:00 2018

@author: dnielsen
"""

import numpy as np
import ROOT
from array import array
from common import show_title, show_text
import matplotlib.pyplot as plt
np.random.seed(42)



######################################
# Channel to voltage 
######################################

# load calibration from handwritten notes
#calib = np.loadtxt("../data/cidercan_calibration.txt", skiprows=1)
#mean = array('d', list(calib[:,2]))
#FWHM = array('d', list(calib[:,3]))
#volt = array('d', list(calib[:,0]))
#verr = array('d', list(calib[:,1]/1000))

# calculate calibration values from fits from files (commented out to increase SPEEED)
#from plot_HVscan import plot_confs
#cal_confs = ["10_12", "10_126", "10_144", "10_152", "10_32", "10_50", "10_79", "10_100"]
#plot_confs(cal_confs, "Calibration", "cal")

# calibration values from fits from files
volt = np.array([ 12.2, 32.4, 49.8, 78.8, 99.98, 125.6, 144.0])#, 151.6 ]
verr = np.array([ 0.220, 0.250, 0.260, 0.280, 0.300, 0.330, 0.200])#, 0.200 ]
mean = np.array([ 72.78174428570775, 210.0788135419234, 326.80617438591594, 515.2406913592274, 657.0662425249498, 832.8401720174578, 954.6987048468875])#, 994.3395288514907 ]
#FWHM = np.array([ 1.640102628564026, 1.6065840660143793, 1.551179962733933, 1.539019813215115, 1.517681317408192, 1.5699381520492788, 1.543115523427546])#, 1.5527226864220074 ]
mean_unc = np.array([ 0.00511876991278432, 0.00500321982192352, 0.00483271239072621, 0.00479870783938811, 0.00471302303328386, 0.00488808702421034, 0.00481756197898153])#, 0.00484989791657979 ]

# fit calibration values to find channel to volt
gr = ROOT.TGraphErrors( len(volt), array('d',volt), array('d',mean), array('d',verr), array('d',mean_unc) )
fit1 = ROOT.TF1("fit1","pol1", 0, 1024);
res = gr.Fit(fit1, "RS")
#res.Print()
print("Fit prob. = {:.1f}%".format(fit1.GetProb()*100))

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

# plot points and fit result
plt.errorbar(y=mean, yerr=5000*mean_unc, x=volt, xerr=10*verr, fmt="none", color='r')
x0 = volt[0]-10
xlast = volt[-1]+10
x = np.linspace(x0, xlast, 1000)
y = [fit1.Eval(x) for x in x]
plt.plot(x, y, 'b-')

# spice it up and show
show_title(ax)
show_text("p( X², ndof ) = p( {:.1f}, {:d} ) = {:.1f}%".format(fit1.GetChisquare(), fit1.GetNDF(), fit1.GetProb()*100), ax, y=0.85)
show_text("y = {:.2f} + {:.2f}*x".format(fit1.GetParameter(0),fit1.GetParameter(1)), ax, y=0.80)
show_text("Note: Channel uncertainties scaled by 5000", ax, y=0.75)
show_text("          Voltage uncertainties scaled by 10", ax, y=0.70)
ax.set_xlabel("Voltage [V]")
ax.set_ylabel("Channel [bit]")
fig.show()
plt.savefig("../graphics/voltchannelcalib.pdf", format='pdf')

# if using a terminal
#input("ready...")