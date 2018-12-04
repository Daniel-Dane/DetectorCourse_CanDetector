#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  31 08:19:00 2018

@author: 
"""

import numpy as np
from math import sqrt, log, exp
import ROOT
from array import array
from common import show_title, show_text, font
import matplotlib.pyplot as plt
np.random.seed(42)

######################################
# From perform_calibration.py (later the two scripts can be merged)
######################################

# calibration values from fits from files
volt = np.array([ 12.2, 32.4, 49.8, 78.8, 99.98, 125.6, 144.0])#, 151.6 ] # This is in mV
verr = np.array([ 0.220, 0.250, 0.260, 0.280, 0.300, 0.330, 0.200])#, 0.200 ] # This is in \mu V
mean = np.array([ 72.78174428570775, 210.0788135419234, 326.80617438591594, 515.2406913592274, 657.0662425249498, 832.8401720174578, 954.6987048468875])#, 994.3395288514907 ]
#FWHM = np.array([ 1.640102628564026, 1.6065840660143793, 1.551179962733933, 1.539019813215115, 1.517681317408192, 1.5699381520492788, 1.543115523427546])#, 1.5527226864220074 ]
mean_unc = np.array([ 0.00511876991278432, 0.00500321982192352, 0.00483271239072621, 0.00479870783938811, 0.00471302303328386, 0.00488808702421034, 0.00481756197898153])#, 0.00484989791657979 ]

C = 1.0 # pF

Q = C * volt * 1E-3#pC
Qerr = C * verr * 1E-3# pC

mean_sys_unc = 0#0.1*mean #[x*0.2 for x in mean]
mean_tot_unc = np.sqrt( mean_unc**2 + mean_sys_unc**2 ) #[sqrt( mean_unc[x]**2 + mean_sys_unc[x]**2 ) for x in range(len(mean))]

# fit calibration values to find channel to volt
gr_Q = ROOT.TGraphErrors( len(Q), array('d',mean), array('d',Q), array('d',mean_tot_unc), array('d',Qerr) )
fit_Q = ROOT.TF1("fit_Q","pol1", 70., 1000.);
res_Q = gr_Q.Fit(fit_Q, "RS")
#res.Print()
print("Fit prob. = {:.1f}%".format(fit_Q.GetProb()*100))

# make figure and axes
fig = plt.figure()
ax = plt.subplot()

# plot points and fit result
ax.errorbar(y=Q, yerr=10*Qerr, x=mean, xerr=5000*mean_tot_unc, color='r', linestyle='None')
x0 = mean[0]-10
xlast = mean[-1]+10
xx = np.linspace(x0, xlast, 1000)
yy = [fit_Q.Eval(x) for x in xx]
ax.plot(xx, yy, 'b-')

# spice it up and show
show_title(ax)
show_text("Note: Channel uncertainties scaled by 5000", ax, y=0.85)
show_text("          Charge uncertainties scaled by 10", ax, y=0.80)
show_text("Fit: y = b + ax", ax, y=0.75)
show_text("b = {:.5f}     ± {:.5f} ({:3.1f}σ from 0)".format(fit_Q.GetParameter(0),fit_Q.GetParError(0),abs(fit_Q.GetParameter(0))/fit_Q.GetParError(0)), ax, y=0.70)
show_text("a = {:.7f} ± {:.7f}".format(fit_Q.GetParameter(1),fit_Q.GetParError(1)), ax, y=0.65)
show_text("p( X², ndof ) = p( {:.1f}, {:d} ) = {:.1f}%".format(fit_Q.GetChisquare(), fit_Q.GetNDF(), fit_Q.GetProb()*100), ax, y=0.60)
ax.set_xlabel("Channel [bit]")
ax.set_ylabel("Charge [pC]")

ax2 = ax.twinx()
ax2.errorbar(y=volt, yerr=10*verr, x=mean, xerr=5000*mean_tot_unc, color='r', linestyle='None')
ax2.set_ylabel("Voltage [V]")

fig.show()
fig.savefig('../graphics/chargevschannel.pdf')

# if using a terminal
#input("ready...")


# calculate calibration values from fits from files (commented out to increase SPEEED)
#from plot_HVscan import plot_confs
#fe_confs = ["100_1420", "100_1470", "100_1523", "100_1572", "100_1617", "100_1717", "10_1978", "10_2000", "10_2042", "10_2080", "20_1901", "20_1978", "2_2201", "2_2254", "2_2303", "40_1717", "40_1801", "40_1853", "40_1901", "4_2080", "4_2124", "4_2201"]
#plot_confs(fe_confs, "Iron", "fe")
#am_confs = ["100_1136", "100_1191", "100_1244", "100_1297", "100_1351", "100_1399", "10_1665", "10_1712", "10_1758", "20_1559", "20_1603", "20_1666", "2_1900", "2_1951", "2_2001", "40_1397", "40_1455", "40_1502", "40_1559", "4_1757", "4_1808", "4_1858", "4_1899"]
#plot_confs(am_confs, "Americium", "am")

volt_Fe     = np.array([1901, 2080, 1523, 1717, 1978, 1470, 1901, 2254, 2201, 1801, 1572, 1717, 2124, 2080, 1617, 1978, 1420, 2303, 2042, 2000, 1853, 2201])
mean_Fe     = np.array([733.000238894175, 287.3984121692716, 128.25781370505598, 477.4041824598613, 319.19575942632196, 91.61897409206911, 361.783492611437, 544.9752538768655, 730.8688865143666, 350.97609048419804, 175.2700149709277, 194.383649040269, 398.8114766454168, 701.3327541796942,  238.737276971451, 644.0874131923914, 67.5509766439823, 800.2437373440258, 523.3445694884902, 379.94073634377907, 516.4059395041165, 359.34178590511846])
mean_unc_Fe = np.array([0.8094318219191396, 0.2665817419245241, 0.1323214953147828, 0.44532396703254734, 0.3106330019657088, 0.1091804108561152, 0.325044174096051, 0.517016213614744, 0.7392414225509276, 0.3207267761833333, 0.17191059551466836, 0.17865397608890807, 0.35836368013832987, 0.6616140317806388, 0.22162734028774378, 0.6152156370324625, 0.0914466961203270, 0.7733527351170573, 0.490041192006111, 0.354816410101314, 0.47682297677610685,  0.33031643094205165])
gain_Fe     = np.array([40, 4, 100, 100, 10, 100, 20, 2, 4, 40, 100, 40, 4, 10, 100, 20, 100, 2, 10, 10, 40, 2])
mean_syst_unc_Fe = 0#[x*0.2 for x in mean_Fe]
mean_tot_unc_Fe = np.sqrt(mean_unc_Fe**2+mean_syst_unc_Fe**2)#[sqrt( mean_unc_Fe[x]**2  +  mean_syst_unc_Fe[x]**2  ) for x in range(len(mean_Fe))]

volt_Am = np.array([1399, 1900, 1136, 1758, 1191, 1899, 1559, 1502, 1665, 1858, 1712, 1808, 1666, 1351, 1455, 1757, 1559, 1244, 1603, 1397, 2001, 1951, 1297])
mean_Am = np.array([544.5050345089651, 326.54317404616205, 127.43071791508439, 583.1179426444957, 168.17293382088152, 664.6309157208256, 616.6382473661497, 424.6967760119166, 307.2646013317212, 491.6720944076552, 425.2699471739359, 345.64279002057714, 619.5692060677644, 412.3728612037973, 315.6736531410237, 238.1303995063521, 304.003946203382, 225.20774106290764, 407.2465678168178, 221.10818867536287, 682.42716036647, 471.56052001424354, 301.853614468452])
mean_unc_Am = np.array([0.506039369118879, 0.3519355211309751, 0.1939561101144484,  0.5668291396588969, 0.2108980208118094, 0.5813120503493101, 0.4314897281094914, 0.29656483037448206, 0.2569726706757164, 0.366757556247128, 0.40212867545908054, 0.3087515122151757, 0.52404887890526, 0.355285435847614, 0.223769775890340, 0.2011415760430970, 0.2370662570751042, 0.272520173643901, 0.2950501496973, 0.174656632615955, 0.769557610293755, 0.53759025583648, 0.335624971474936])
gain_Am = np.array([100, 2, 100, 10, 100, 4, 40, 40, 10, 4, 10, 4, 20, 100, 40, 4, 20, 100, 20, 40, 2, 2, 100])
mean_syst_unc_Am = 0#[x*0.2 for x in mean_Am]
mean_tot_unc_Am = np.sqrt(mean_unc_Am**2+mean_syst_unc_Am**2)#[sqrt( mean_unc_Am[x]**2  +  mean_syst_unc_Am[x]**2  ) for x in range(len(mean_Am))]

num_of_electrons_Fe = np.zeros(len(volt_Fe))
num_of_electrons_error_Fe = np.zeros(len(volt_Fe))
for j in range(len(volt_Fe)):
	num_of_electrons_Fe[j] = ( (fit_Q.Eval(mean_Fe[j]) * (10./gain_Fe[j]) * 1E-12) / (1.602 * 1E-19))
	num_of_electrons_error_Fe[j] = sqrt( (1**2 * fit_Q.GetParError(0)**2)  +  (fit_Q.GetParError(1)**2 * (mean_Fe[j]**2)) + ((fit_Q.GetParameter(1)**2) * (mean_tot_unc_Fe[j]**2))   ) * ((10./gain_Fe[j]) * 1E-12) / (1.602 * 1E-19)

num_of_electrons_Am = np.zeros(len(volt_Am))
num_of_electrons_error_Am = np.zeros(len(volt_Am))
for j in range(len(volt_Am)):
	num_of_electrons_Am[j] = ( (fit_Q.Eval(mean_Am[j]) * (10./gain_Am[j]) * 1E-12) / (1.602 * 1E-19))
	num_of_electrons_error_Am[j] = sqrt(  (1**2 * fit_Q.GetParError(0)**2)  +  (fit_Q.GetParError(1)**2 * (mean_Am[j]**2)) + ((fit_Q.GetParameter(1)**2)*(mean_tot_unc_Am[j]**2)) ) * ((10./gain_Am[j]) * 1E-12) / (1.602 * 1E-19)

fig1, ax1 = plt.subplots()
ax1.set_xlabel("Voltage [V]")
ax1.set_ylabel("Number of electrons")
ax1.set_yscale("log", nonposy="clip")

ax1.errorbar(volt_Fe, num_of_electrons_Fe, num_of_electrons_error_Fe, color = "r", label = r"$\gamma = 5.9$ keV (Fe-55)", linestyle='None', marker='o')
ax1.errorbar(volt_Am, num_of_electrons_Am, num_of_electrons_error_Am, color = "b", label = r"$\gamma = 59.5$ keV (Am-241)", linestyle='None', marker='d')
ax1.text(0.5,0.9, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
            fontproperties=font, transform=ax.transAxes)
ax1.legend(loc="upper left", numpoints=1)
ax1.grid()
fig1.show()
fig1.savefig('../graphics/numbervsvoltage.pdf')
#plt.close(fig1)
#plt.clf()


M_Fe =  num_of_electrons_Fe * ( 1./227.)
M_unc_Fe = num_of_electrons_error_Fe * ( 1./227.)
M_Am = num_of_electrons_Am * (1./2290.)
M_unc_Am = num_of_electrons_error_Am * (1./2290.)

xx = np.linspace(1000., 2500., 1000)
#yy = np.zeros(1000)
def M_thoe(V):
	b = 6.58234/2 #cm
	a =  0.005/2 #cm 50 microns (from our measurement)
	DeltaU = 23.6 #V
	p =  1.03# 1 atm
	K = 4.8 * 1E+4 # 10^4 V/cm*atm 

	first_term = (V)/(log(b/a))
	second_term = (log(2.))/(DeltaU)
	third_term = log(   (V)  /  (K* p * a * log(b/a))    )  #- log(K)   
	#lnM = first_term * second_term * third_term
	lnM = ((V)/(log(b/a))) * ((log(2.))/(DeltaU)) * (log(   (V)  /  (K* p * a * log(b/a))    ))
	return lnM

def M_theo_unc(V):
	# assume DeltaU, K and V constant
	# add factos for a, b, p
	b = 6.58234/2 #cm
	a =  0.005/2 #cm 50 microns (from our measurement)
	DeltaU = 23.6 #V
	p =  1.03# 1 atm
	K = 4.8 * 1E+4 # 10^4 V/cm*atm 

	sigma_a = 0.001
	sigma_b = 0.009
	sigma_p = 0.001
	a_factor = ((log(2.)*V)/DeltaU) * ((   1 - log(K) - log(b/a)   + log((V)/(p*a*log(b/a)))  )/(  a * log(b/a) *log(b/a)  ))
	b_factor = ((log(2.)*V)/DeltaU) * ( (  1 - log(K) + log( (V)/(a*p*log(b/a)) ) )/( b*log(b/a) *log(b/a)  )  )
	p_factor = ((log(2.)*V)/DeltaU) * (-1./p) * (1./log(b/a))

	sigma_lnM = sqrt( sigma_a * sigma_a * a_factor*a_factor  + sigma_b*sigma_b *b_factor*b_factor )#+ sigma_p*sigma_p * p_factor*p_factor     )
	return sigma_lnM

yy = [exp(M_thoe(x)) for x in xx]
yy2 = [exp(M_thoe(x) + M_theo_unc(x)) for x in xx]
yy3 = [exp(M_thoe(x) - M_theo_unc(x)) for x in xx]

fig2, ax2 = plt.subplots()
ax2.set_xlabel("Voltage [V]")
ax2.set_ylabel("M")
ax2.set_yscale("log", nonposy="clip")
ax2.errorbar(volt_Fe, M_Fe, M_unc_Fe, color = "r", label = r"$\gamma = 5.9$ keV (Fe-55)", marker='o', linestyle='None')
ax2.errorbar(volt_Am, M_Am, M_unc_Am, color = "b", label = r"$\gamma = 59.5$ keV (Am-241)", marker='d', linestyle='None')
ax2.plot(xx,yy,color = "g", label=r"Prediction (± 1$\sigma$)")
ax2.fill_between(xx,yy3,yy2,facecolor="none", linewidth=0.0, edgecolor = "green", hatch="\\\\")#, label=r"Prediction ± 1$\sigma$")
#ax2.errorbar(volt_Fe, num_of_electrons_Fe, num_of_electrons_error_Fe, color = "k", label = r"$\gamma = 5.9 keV$ Fe-55", linestyle='None', marker='o')
#ax2.errorbar(volt_Am, num_of_electrons_Am, num_of_electrons_error_Am, color = "k", label = r"$\gamma = 59.5 keV$ Am-241", linestyle='None', marker='d')
ax2.text(0.5,0.9, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
            fontproperties=font, transform=ax2.transAxes)
ax2.legend(loc="upper left", numpoints=1)
ax2.grid()
fig2.show()
fig2.savefig('../graphics/Mvsvoltage.pdf')
#plt.close(fig2)
#plt.clf()



import sys
if sys.version_info[0]==2:
    raw_input("Press enter to finish")











