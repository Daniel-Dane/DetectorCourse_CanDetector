#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.optimize as optimization

def func(x, a, b):
    return a*x + b 

#****************************************************************************
# Preamplifier business

# A. Extract the dimensionless preamp gain

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
perr = np.sqrt(np.diag(covariance)) # 1-sigma error on the parameter
plt.plot(input_v,func(input_v,a,b),label='linear fit')
plt.text(0.08,4.,r'Preamplifier gain: %.02f $\pm$ %.02f'%(a/10.,perr[0]/10.)) # divide by 10 to remove effect of the coarse gain setting
plt.show()


# B. Extract gain, but including the MCA-to-volt conversion factor
output_mca = calib[:,2]
Dy_mca = calib[:,4]

plt.errorbar(input_v,output_mca,yerr=Dy_mca,fmt='v',color='r')
plt.title('Preamplifier gain')
plt.xlabel('input voltage (V)',fontsize=13)
plt.ylabel('Digitized Output (MCA)',fontsize=13) # As measured with the MCA

# Fit a linear function. Print out the slope of the preamplifier gain and its uncertainty
Fit_mca,covariance_mca  = optimization.curve_fit(func,input_v ,output_mca, sigma=Dy_mca,
                                         absolute_sigma=True)

mca_preampgain = Fit_mca[0]
mca_b = Fit_mca[1]
mca_perr = np.sqrt(np.diag(covariance_mca)) # 1-sigma error on the parameter
plt.plot(input_v,func(input_v,mca_preampgain,mca_b),label='linear fit',color='k')
plt.text(0.06,200.,r'Preamplifier gain: %.02f $\pm$ %.02f MCA\V'%(mca_preampgain/10.,mca_perr[0]/10.))
plt.show()


# Gives coarse gain to MCA relationship
Data_fe = np.genfromtxt("../data/cidercan_fe_HVscan.txt",skip_header=1)
Data_am = np.genfromtxt("../data/cidercan_am_HVscan.txt",skip_header=1)




#*************************************************************************************
# Coarse gain business
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
   
plt.errorbar(x_values+0.2,y_am,yerr=dy_am,fmt='P',color='b',label=r'$Am^{241}$')

# Get a combined measurement of Am and Fe data
# weighted average
coarse_gain = (y_fe/dy_fe**2.+y_am/dy_am**2.)/(1./dy_fe**2+1./dy_am**2)
coarse_gain_unc = np.sqrt(1./((1./dy_fe**2+1./dy_am**2)))
print("Coarse gains: ",x_fe)
print("Coarse gain values: ",coarse_gain)
print("Coarse gain uncertainties: ",coarse_gain_unc)

plt.errorbar(x_values+0.4,coarse_gain,yerr=coarse_gain_unc,fmt='s',color='r',label=r'Combined')

# Change the xticks so that they reflect the labels
ax = plt.gca()
ax.set_xticks(x_values+0.1)
ax.set_xticklabels(x_fe)
plt.title("Coarse Gain Ratios")
plt.xlabel('Ratio between steps',fontsize=13)
plt.ylabel('Measured Gain',fontsize=13)
plt.legend()
plt.show()



#******************************************************************

# Part where you plot the theoretical joker bands
import pickle
import seaborn as sns

dataname = "./uncertainties_100000trials.p"
data = pickle.load(open(dataname))

ten_pc = data[0]
one_pc = data[1]
ninety = data[2]
ninety9= data[3]
median = data[4]
voltages=data[5]

plt.fill_between(voltages,one_pc,ninety9,label=r'$2\sigma$',color='#deebf7')
plt.fill_between(voltages,ten_pc,ninety,label=r'$1\sigma$',color='#9ecae1')
plt.plot(voltages,one_pc,'#deebf7')
plt.plot(voltages,ninety9,'#deebf7')
plt.plot(voltages,ten_pc,'#9ecae1')
plt.plot(voltages,ninety,'#9ecae1')
plt.plot(voltages,median,'k',label='median')

#**********************************************************************
# Part where you plot the data, divided by the proper factor
M_coarse_Fe = Data_fe[:,0] # coarse gain setting
M_X_Fe = Data_fe[:,1]      # operating voltage
M_MCA_Fe = Data_fe[:,2]    # Signal in MCA counts
M_MCA_Fe_unc = Data_fe[:,3]


# Instead of the raw files, I'll use the ones computed by Daniel for americium
# calibration values from fits from files
M_X_am = np.array([1399, 1900, 1136, 1758, 1191, 1899, 1559, 1502, 1665, 1858, 1712, 1808, 1666, 1351, 1455, 1757, 1559, 1244, 1603, 1397, 2001, 1951, 1297])/1000.
M_MCA_am = np.array([544.5050345089651, 326.54317404616205, 127.43071791508439, 583.1179426444957, 168.17293382088152, 664.6309157208256, 616.6382473661497, 424.6967760119166, 307.2646013317212, 491.6720944076552, 425.2699471739359, 345.64279002057714, 619.5692060677644, 412.3728612037973, 315.6736531410237, 238.1303995063521, 304.003946203382, 225.20774106290764, 407.2465678168178, 221.10818867536287, 682.42716036647, 471.56052001424354, 301.853614468452])
M_MCA_am_unc = np.array([0.506039369118879, 0.3519355211309751, 0.1939561101144484,  0.5668291396588969, 0.2108980208118094, 0.5813120503493101, 0.4314897281094914, 0.29656483037448206, 0.2569726706757164, 0.366757556247128, 0.40212867545908054, 0.3087515122151757, 0.52404887890526, 0.355285435847614, 0.223769775890340, 0.2011415760430970, 0.2370662570751042, 0.272520173643901, 0.2950501496973, 0.174656632615955, 0.769557610293755, 0.53759025583648, 0.335624971474936])
M_coarse_am = np.array([100, 2, 100, 10, 100, 4, 40, 40, 10, 4, 10, 4, 20, 100, 40, 4, 20, 100, 20, 40, 2, 2, 100])


"""
# Part where you plot the data, divided by the proper factor
M_coarse_am = Data_am[:,0] # coarse gain setting
M_X_am = Data_am[:,1]/1000.      # operating voltage
M_MCA_am = Data_am[:,2]    # Signal in MCA counts
M_MCA_am_unc = Data_am[:,3]
"""


preamp_gain = mca_preampgain/10. # preamp data was taken at a coarse gain of 10. Removing it form the equation
preamp_unc  = mca_perr[0]/10.

# Write out the ratios explicitely:
G100_G40=coarse_gain[0]
G40_G20=coarse_gain[1]
G20_G10=coarse_gain[2]
G10_G4=coarse_gain[3]
G4_G2=coarse_gain[4]

G100_G40_e=coarse_gain_unc[0]
G40_G20_e =coarse_gain_unc[1]
G20_G10_e =coarse_gain_unc[2]
G10_G4_e  =coarse_gain_unc[3]
G4_G2_e   =coarse_gain_unc[4]

#Define the correct coarse gain factor to use
M_Gratio  = np.zeros(len(M_coarse_Fe))
M_Gratio_unc = np.zeros(len(M_coarse_Fe))

for g,i in zip(M_coarse_Fe,range(0,len(M_coarse_Fe))):

    if g==2:
        M_Gratio[i] = 10./G10_G4/G4_G2
        M_Gratio_unc[i] = M_Gratio[i]*np.sqrt((G10_G4_e/G10_G4)**2.+(G4_G2_e/G4_G2)**2.)
    elif g==4:
        M_Gratio[i] = 10./G10_G4
        M_Gratio_unc[i] = M_Gratio[i]*(G10_G4_e/G10_G4)
    elif g==10:
        M_Gratio[i] = 10.
        M_Gratio_unc[i] = 0.0
    elif g==20:
        M_Gratio[i] = 10.*G20_G10
        M_Gratio_unc[i] = M_Gratio[i]*(G20_G10_e/G20_G10)
    elif g==40:
        M_Gratio[i] = 10.*G20_G10*G40_G20
        M_Gratio_unc[i] = M_Gratio[i]*np.sqrt((G20_G10_e/G20_G10)**2.+(G40_G20_e/G40_G20)**2.)
    elif g==100:
        M_Gratio[i] = 10.*G20_G10*G40_G20*G100_G40
        M_Gratio_unc[i] = M_Gratio[i]*np.sqrt((G20_G10_e/G20_G10)**2.+(G40_G20_e/G40_G20)**2.+(G100_G40_e/G100_G40)**2.)
        
Q_detector_fe = M_MCA_Fe/preamp_gain*(1e-12)/1.602e-19/M_Gratio
Q_detector_fe_unc = Q_detector_fe*np.sqrt((preamp_unc/preamp_gain)**2.+(M_Gratio_unc/M_Gratio)**2.+(M_MCA_Fe_unc/M_MCA_Fe)**2.)
M_Fe = Q_detector_fe/227.
M_Fe_unc = Q_detector_fe_unc/227.
ln_M_Fe = np.log(M_Fe)
ln_M_Fe_unc = M_Fe_unc/M_Fe


#***************************************************************************
# Doing the Americium stuff
M_Gratio_am  = np.zeros(len(M_coarse_am))
M_Gratio_am_unc = np.zeros(len(M_coarse_am))
         
for g,i in zip(M_coarse_am,range(0,len(M_coarse_am))):

    if g==2:
        M_Gratio_am[i] = 10./G10_G4/G4_G2
        M_Gratio_am_unc[i] = M_Gratio_am[i]*np.sqrt((G10_G4_e/G10_G4)**2.+(G4_G2_e/G4_G2)**2.)
    elif g==4:
        M_Gratio_am[i] = 10./G10_G4
        M_Gratio_am_unc[i] = M_Gratio_am[i]*(G10_G4_e/G10_G4)
    elif g==10:
        M_Gratio_am[i] = 10.
        M_Gratio_am_unc[i] = 0.0
    elif g==20:
        M_Gratio_am[i] = 10.*G20_G10
        M_Gratio_am_unc[i] = M_Gratio_am[i]*(G20_G10_e/G20_G10)
    elif g==40:
        M_Gratio_am[i] = 10.*G20_G10*G40_G20
        M_Gratio_am_unc[i] = M_Gratio_am[i]*np.sqrt((G20_G10_e/G20_G10)**2.+(G40_G20_e/G40_G20)**2.)
    elif g==100:
        M_Gratio_am[i] = 10.*G20_G10*G40_G20*G100_G40
        M_Gratio_am_unc[i] = M_Gratio_am[i]*np.sqrt((G20_G10_e/G20_G10)**2.+(G40_G20_e/G40_G20)**2.+(G100_G40_e/G100_G40)**2.)
             

Q_detector_am = M_MCA_am/preamp_gain*(1e-12)/1.602e-19/M_Gratio_am
Q_detector_am_unc = Q_detector_am*np.sqrt((preamp_unc/preamp_gain)**2.+(M_MCA_am_unc/M_MCA_am)**2.+(M_Gratio_am_unc/M_Gratio_am)**2.)#+

M_am = Q_detector_am/(2290.)
M_am_unc = Q_detector_am_unc/(2290.)

ln_M_am = np.log(M_am)
ln_M_am_unc = M_am_unc/M_am

plt.errorbar(M_X_am,ln_M_am,yerr=ln_M_am_unc,fmt='sr',label=r'$Am^{241}$',color='#de2d26')
plt.errorbar(M_X_Fe,ln_M_Fe,yerr=ln_M_Fe_unc,fmt='o',label=r'$Fe^{55}$',color='#31a354')

plt.xlabel('Operating voltage (kV)',fontsize=13)
plt.ylabel(r'$\ln(M)$',fontsize=13)
plt.legend()
plt.xlim([0.8,2.5])
#plt.yscale('log')
plt.ylim([-2.,14.5])
plt.show()

            
