#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.optimize as optimization

def func(x, a, b):
    return a*x + b 

# gives the MCA-to-preamp voltage relationship
calib = np.genfromtxt("../data/cidercan_calibration.txt",skip_header=1)
input_v = calib[:,0]/1000.
outputv = calib[:,5]
Dy_preamp = calib[:,6]/1000.

plt.errorbar(input_v,outputv,yerr=Dy_preamp,fmt='s')
plt.title('Preamplifier gain')
plt.xlabel('input voltage (V)',fontsize=13)
plt.ylabel('output after preamplification (V)',fontsize=13) # As measured with the oscilloscope

# Fit a linear function. Print out the slope of the preamplifier gain and its uncertainty
Fit,covariance  = optimization.curve_fit(func,input_v ,outputv, sigma=Dy_preamp,
                                         absolute_sigma=True)

a = Fit[0]
b = Fit[1]
print covariance
perr = np.sqrt(np.diag(covariance)) # 1-sigma error on the parameter
print perr[0]

plt.plot(input_v,func(input_v,a,b),label='linear fit')
plt.text(0.08,4.,r'Preamplifier gain: %.02f $\pm$ %.02f'%(a,perr[0]))
plt.show()


# Gives coarse gain to MCA relationship
Data_fe = np.genfromtxt("../data/cidercan_fe_HVscan.txt",skip_header=1)
Data_am = np.genfromtxt("../data/cidercan_am_HVscan.txt",skip_header=1)

# columns are the following for the HV scans
#[Coarse gain]	[Voltage kV]	[MCA mean]	[MCA FWHM]	[FWHM/mean]


# Step 1: Estimate the uncertainty on the coarse gain amplifier
#
# For this we have 4 points where the same voltage is probed in two different setting

point_1 = Data_fe[5:7,:]
point_2 = Data_fe[9:11,:]
point_3 = Data_fe[11:13,:]
point_4 = Data_fe[15:17,:]
point_5 = Data_fe[18:20,:]

point_6 = Data_am[5:7,:]
point_7 = Data_am[9:11,:]
point_8 = Data_am[12:14,:]
point_9 = Data_am[15:17,:]
point_10 = Data_am[19:21,:]


alll_fe = [point_1,point_2,point_3,point_4, point_5]
alll_am = [point_6,point_7,point_8,point_9, point_10]

y_fe  = np.zeros(len(alll_fe))
dy_fe = np.zeros(len(alll_fe))
x_fe=[]

y_am  = np.zeros(len(alll_am))
dy_am = np.zeros(len(alll_am))
x_am=[]


x_values = np.linspace(0,6,5)

for p,i in zip(alll_fe,range(0,len(alll_fe))):
    
    pair = [p[0,0],p[1,0]]
    x_fe.append('%i->%i'%(min(pair),max(pair)))
    y_fe[i]  = p[0,2]/float(p[1,2]) # Ratio of output MCA
    dy_fe[i] = y_fe[i]*np.sqrt((p[0,3]/float(p[0,2]))**2.+(p[1,3]/p[1,2])**2.) # propagated uncertainties
    
plt.errorbar(x_values,y_fe,yerr=dy_fe,fmt='o',color='k',label=r'$Fe^{55}$')

for p,i in zip(alll_am,range(0,len(alll_am))):
    
    pair = [p[0,0],p[1,0]]
    x_am.append('%i->%i'%(min(pair),max(pair)))
    y_am[i]  = p[0,2]/float(p[1,2]) # Ratio of output MCA
    dy_am[i] = y_am[i]*np.sqrt((p[0,3]/float(p[0,2]))**2.+(p[1,3]/p[1,2])**2.) # propagated uncertainties
   
plt.errorbar(x_values+0.2,y_am,yerr=dy_am,fmt='o',color='b',label=r'$Am^{241}$')

# Get a combined measurement of Am and Fe data
# weighted average
coarse_gain = (y_fe/dy_fe**2.+y_am/dy_am**2.)/(1./dy_fe**2+1./dy_am**2)
coarse_gain_unc = np.sqrt(1./((1./dy_fe**2+1./dy_am**2)))
print x_fe
print coarse_gain
print coarse_gain_unc

plt.errorbar(x_values+0.4,coarse_gain,yerr=coarse_gain_unc,fmt='o',color='r',label=r'Combined')

# Change the xticks so that they reflect the labels
ax = plt.gca()
ax.set_xticks(x_values+0.1)
ax.set_xticklabels(x_fe)
plt.title("Coarse Gain uncertainty")
plt.xlabel('Coarse gain step',fontsize=13)
plt.ylabel('Measured Gain',fontsize=13)
plt.legend()
plt.show()

#******************************************************************

# Part where you plot the theoretical joker bands

voltages = X
plt.fill_between(voltages,one_pc,ninety9,label=r'$2\sigma$',color='m')
plt.fill_between(voltages,ten_pc,ninety,label=r'$1\sigma$',color='g')
plt.plot(voltages,one_pc,'m')
plt.plot(voltages,ninety9,'m')
plt.plot(voltages,ten_pc,'g')
plt.plot(voltages,ninety,'g')
plt.plot(voltages,median,'k',label='median')

plt.xlabel('Operating voltage (kV)')
plt.ylabel('M')
plt.legend()
plt.ylim([1,6.])
plt.show()

            
