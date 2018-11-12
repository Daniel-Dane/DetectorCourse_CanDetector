# Remember to set up ROOT in your environment before running this script, as several ROOT functions are used
import numpy as np
from ROOT import TCanvas, TGraphErrors, TGraph, TF1, kBlack, gStyle, TRandom3, TH1D, TLatex
from array import array
import sys
from math import sqrt, exp
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages

# Factor to match the impedances of the two oscilloscopes used to measure

# The amplitude of the calibration curve was acquired with one oscilloscope while the amplitude of the signal was 
# acquired with another oscilloscope and those two have two be matched with a factor, since the impedances are different
impedance_factor = 0.7
rc_factor_10_90  = 2.197

def LoadCalibrationFile(calibration_filename):
    calibration_diode = np.loadtxt(calibration_filename,skiprows=1)
    C           = calibration_diode[:,0]
    e_C         = np.zeros([len(C)])
    noise_rms   = calibration_diode[:,1]
    d_noise_rms = calibration_diode[:,2]
    fall_time   = calibration_diode[:,3]
    d_fall_time = calibration_diode[:,4]
    amplitude   = calibration_diode[:,5]
    d_amplitude = calibration_diode[:,6]
    return C, e_C, noise_rms, d_noise_rms, fall_time, d_fall_time, amplitude, d_amplitude


def OtherPlots(calibration_filename, IV_filename, CV_filename_openneedle, CV_filename):
    # Read data fron the calibration curve dataset
    C, e_C, noise_rms, d_noise_rms, fall_time, d_fall_time, amplitude, d_amplitude = LoadCalibrationFile(calibration_filename)

    font0 = FontProperties()
    font = font0.copy()
    font.set_style('italic')
    font.set_weight('bold')
    font.set_size('x-large')
    font_par = font0.copy()
    font_par.set_style('italic')
    font_par.set_weight('medium')

    # Plot of the Noise RMS vs capacitance
    fig, ax = plt.subplots()
    ax.set_xlabel('C [pF]')
    ax.set_xlim(0.0, 120)
    ax.set_ylabel('Noise RMS [\mu V]')
    ax.errorbar(C, noise_rms, d_noise_rms, marker='o', color='b', linestyle='None', label='Data')
    ax.text(0.5,0.9, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
    ax.legend(loc='upper left',numpoints=1)
    plt.grid()
    #plt.show()
    fig.savefig('../graphics/noise_vs_capacitance.pdf')
    plt.close(fig)
    plt.clf()

    # Plot of thw Amplitude vs the capacitance
    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('C [pF]')
    ax1.set_xlim(0.0, 120)
    ax1.set_ylabel('Amplitude [mV]')
    ax1.errorbar(C, amplitude, d_amplitude, marker='o', color='b', linestyle='None', label='Data')
    ax1.text(0.5,0.8, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
    ax1.legend(loc='upper right', numpoints=1)
    plt.grid()
    #plt.show()
    fig1.savefig('../graphics/amplitude_vs_capacitance.pdf')
    plt.close(fig1)
    plt.clf()

    IV = np.loadtxt(IV_filename,skiprows=1)
    V   = IV[:,0]
    pc1 = IV[:,1]
    tc1 = IV[:,2]
    pc2 = IV[:,3]
    tc2 = IV[:,4]
    pc3 = IV[:,5]
    tc3 = IV[:,6]
    pc4 = IV[:,7]
    tc4 = IV[:,8] 

    pc_all = np.array([pc1,pc2,pc3,pc4])
    tc_all = np.array([tc1,tc2,tc3,tc4])
    pc_mean = np.mean(pc_all, axis = 0)
    pc_std  = np.std(pc_all, axis = 0)
    tc_mean = np.mean(tc_all, axis = 0)
    tc_std  = np.std(tc_all, axis = 0) 

    # Plot of the IV curve
    fig3, ax3 = plt.subplots(2, sharex= True)
    ax3[1].set_xlabel('Voltage [v]')
    ax3[0].set_ylabel('Pad current [A]')
    ax3[1].set_ylabel('Total current [A]')
    ax3[0].set_ylim(0.0, -4.2E-8)
    ax3[1].set_ylim(0.0, -4.2E-8)
    ax3[0].set_xlim(0.0, -160)
    ax3[1].set_xlim(0.0, -160)
    ax3[0].errorbar(V, pc_mean, pc_std, marker='o', color='k', label='Pad current')
    ax3[1].errorbar(V, tc_mean, tc_std, marker='o', color='k', label='Total current')
    ax3[0].text(0.5,0.8, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
    ax3[1].text(0.5,0.8, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
    ax3[1].annotate('Depletion voltage', xy=((V[6] - V[5])/2 + V[5], tc_mean[6]), xytext=(-75, -2.8E-8), verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font_par, transform=ax.transAxes,arrowprops=dict(arrowstyle="->",facecolor='black'))
    ax3[1].legend(loc='upper left', numpoints=1)
    ax3[0].legend(loc='upper left', numpoints=1)
    ax3[0].grid()
    ax3[1].grid()
    #plt.show()
    fig3.savefig('../graphics/IV_V_vs_A.pdf')
    plt.close(fig3)
    plt.clf()

    CV_openneedle = np.loadtxt(CV_filename_openneedle,skiprows=1)
    V_openneedle  = CV_openneedle[:,0]
    C1_openneedle = CV_openneedle[:,1]
    D1_openneedle = CV_openneedle[:,2]
    C2_openneedle = CV_openneedle[:,3]
    D2_openneedle = CV_openneedle[:,4]
    C3_openneedle = CV_openneedle[:,5]
    D3_openneedle = CV_openneedle[:,6]
    C4_openneedle = CV_openneedle[:,7]
    D4_openneedle = CV_openneedle[:,8] 

    CV = np.loadtxt(CV_filename,skiprows=1)
    V  = CV[:,0]
    C1 = CV[:,1]
    D1 = CV[:,2]
    C2 = CV[:,3]
    D2 = CV[:,4]
    C3 = CV[:,5]
    D3 = CV[:,6]
    C4 = CV[:,7]
    D4 = CV[:,8] 

    # Plot of the CV curve
    C_openneedle_all = np.array([C1_openneedle, C2_openneedle, C3_openneedle, C4_openneedle])
    C_all = np.array([C1, C2, C3, C4])
    C_final_all = C_all - C_openneedle_all
    C_openneedle_mean = np.mean(C_openneedle_all, axis = 0)
    C_openneedle_std  = np.std(C_openneedle_all, axis = 0) 
    C_mean = np.mean(C_all, axis = 0)
    C_std  = np.std(C_all, axis = 0)
    C_final_mean = np.mean(C_final_all, axis = 0)
    C_final_std  = np.std(C_final_all, axis = 0)

    fig3, ax3 = plt.subplots(1)
    ax3.set_xlabel('Voltage [v]')
    ax3.set_ylabel('Capacitance [A]')
    ax3.set_ylim(1E-12, 1E-9)
    ax3.set_xlim(0.0, -160)
    ax3.set_xlim(0.0, -160)
    ax3.set_yscale("log", nonposy='clip')
    ax3.errorbar(V, C_openneedle_mean, C_openneedle_std, marker='o', color='b', label='Open needle ')
    ax3.errorbar(V, C_mean, C_std, marker='o', color='g', label='Diode')
    ax3.errorbar(V, C_final_mean, C_final_std,  marker='o', color='r', label='Diode - open needle pad 1')
    ax3.text(0.1,0.8, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
    #ax3.text(0.4,0.6, 'C = {0:.2f} +/- {1:.2f} pF'.format(C_final_mean[7]*1E12, C_final_std[7]*1E12), verticalalignment='bottom', horizontalalignment='left',
    #            fontproperties=font_par, transform=ax.transAxes)
    ax3.annotate('C = {0:.2f} +/- {1:.2f} pF'.format(C_final_mean[7]*1E12, C_final_std[7]*1E12), xy=(V[7], C_final_mean[7]), xytext=(-75, 1.1E-10), verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font_par, transform=ax.transAxes,arrowprops=dict(arrowstyle="->",facecolor='black'))
    ax3.legend(loc='upper right', numpoints=1)
    ax3.grid()
    #plt.show()
    fig3.savefig('../graphics/V_vs_C.pdf')
    plt.close(fig3)
    plt.clf()

    return

def CalibrationRiseTime(calibration_filename):
    # Read data fron the calibration curve dataset
    C, e_C, noise_rms, d_noise_rms, fall_time, d_fall_time, amplitude, d_amplitude = LoadCalibrationFile(calibration_filename)

    gStyle.SetOptStat(1111)

    # Construct the  curve
    calib_curve = TGraphErrors(len(C), array('d',C.tolist()), array('d', fall_time.tolist()), array('d', e_C.tolist()), array('d',d_fall_time.tolist()) )

    # Construct the fit with a 2nd order polynomial
    fit_curve = TF1('fit_curve','[0]*x*x + [1]*x + [2]',0.,110.)
    fit_curve.SetParameter(0, 0.0)
    fit_curve.SetParameter(1, 0.0)
    fit_curve.SetParameter(2, 0.0)
    fit_curve.SetLineColor(2)

    # FIt the calibration curve
    calib_curve.Fit(fit_curve,'R')

    # Extract the results fron the fit
    par0   = fit_curve.GetParameter(0)
    e_par0 = fit_curve.GetParError(0)
    par1   = fit_curve.GetParameter(1)
    e_par1 = fit_curve.GetParError(1)
    par2   = fit_curve.GetParameter(2)
    e_par2 = fit_curve.GetParError(2)

    xx = np.linspace(0., 105, 1000)
    yy = [fit_curve.Eval(x) for x in xx]
    par = [par2, par1, par0]
    l = 'fit pol2 $p_0$: {:.3g} $p_1$: {:.3g} $p_2$: {:.2f}'.format(*par)

    # Plot the calibration curve with matplotlib
    font0 = FontProperties()
    font = font0.copy()
    font.set_style('italic')
    font.set_weight('bold')
    font.set_size('x-large')

    fig, ax = plt.subplots()
    ax.set_xlabel('C [pF]')
    ax.set_xlim(0.0, 120)
    ax.set_ylabel('Rise time [ns]')
    ax.errorbar(C, fall_time, d_fall_time, marker='o', color='k', linestyle='None', label='Data')
    ax.plot(xx, yy, label=l, color='r')
    #plt.set_title('Diode calibration curve')
    ax.text(0.5,0.9, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
    ax.legend(loc='upper left', numpoints=1)
    plt.grid()
    #plt.show()
    fig.savefig('../graphics/calibration_diode.pdf')
    plt.close(fig)
    plt.clf()

    return par0, e_par0, par1, e_par1, par2, e_par2

def CalibrationENC(calibration_filename):
    # Read data fron the calibration curve dataset
    C, e_C, noise_rms, d_noise_rms, fall_time, d_fall_time, amplitude, d_amplitude = LoadCalibrationFile(calibration_filename)

    noise_rms = noise_rms / 1000. # convert form \muV to mV
    d_noise_rms = d_noise_rms / 1000.

    ENC = (noise_rms/amplitude) * (C * 20.0) * pow(10,-15) / (1.602 * pow(10,-19) ) #pF * mV = 10^-12 * 10^-3 = 10^-15 C / 10^-19 C = 10^4 electrons
    d_ENC = np.zeros([len(ENC)])
    for i in range(len(ENC)):
        d_ENC[i] = sqrt( d_noise_rms[i]*d_noise_rms[i]*(C[i] * 20.0/amplitude[i])*(C[i] * 20.0/amplitude[i])  +  d_amplitude[i]*d_amplitude[i]*((noise_rms[i]/(amplitude[i]*amplitude[i])) * (C[i] * 20.0))*((noise_rms[i]/(amplitude[i]*amplitude[i])) * (C[i] * 20.0))    )
    
    d_ENC = d_ENC * pow(10,-15) / (1.602 * pow(10,-19) ) # conversted to # of electrons

    calib_curve = TGraphErrors(len(ENC), array('d', C[:-1].tolist()), array('d', ENC.tolist()), array('d', e_C.tolist()), array('d', d_ENC) )

    # Construct the fit with a 1st and 2nd order polynomial
    fit_curve = TF1('fit_curve','[1]*x + [0]',0.,110.)
    fit_curve2 = TF1('fit_curve2','[0]*x*x + [1]*x + [2]',0.,110.)

    # Fit the calibration curve
    calib_curve.Fit(fit_curve,'R')
    calib_curve.Fit(fit_curve2, 'R')

    lin_par0   = fit_curve.GetParameter(0)
    lin_e_par0 = fit_curve.GetParError(0)
    lin_par1   = fit_curve.GetParameter(1)
    lin_e_par1 = fit_curve.GetParError(1)

    xx = np.linspace(0., 105, 1000)
    yy = [fit_curve.Eval(x) for x in xx]
    yy2 = [fit_curve2.Eval(x) for x in xx]

    par = [lin_par0, lin_par1]
    l = 'fit pol2 $p_0$: {:.3g} $p_1$: {:.3g}'.format(*par)


    # Plot the calibration curve with matplotlib
    font0 = FontProperties()
    font = font0.copy()
    font.set_style('italic')
    font.set_weight('bold')
    font.set_size('x-large')

    fig, ax = plt.subplots()
    ax.set_xlabel('C [pF]')
    ax.set_ylabel('ENC [# of electrons]')
    ax.errorbar(C, ENC, d_ENC, marker='o', color='k', linestyle='None', label='Data')
    ax.plot(xx, yy, color='r', label=l)
    ax.plot(xx, yy2, color='g', label='Fit2')
    ax.text(0.5,0.9, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
    ax.legend(loc='upper left', numpoints=1)
    plt.grid()
    #plt.show()
    fig.savefig('../graphics/calibrationENC_diode.pdf')
    plt.close(fig)
    plt.clf()

    return lin_par0, lin_e_par0, lin_par1, lin_e_par1

# Extract the noise and the rise time from the data collected
def SignalCurves(datasets, path):
    decay_data             = []
    decay_curve            = []
    rise_time              = []
    rise_time_error        = []
    raise_fit              = []
    noise_std_holesystem   = []
    noise_mean_holesystem  = []
    noise_fit              = []

    for i in datasets:
        q = np.loadtxt(path + i, skiprows=6, delimiter=',')
        decay_data.append(q)
        print "MESSAGE INFO: The sample %s has been loaded "%(i)

    j = -1
    for data in decay_data:
        j += 1

        data[:,1] = data[:,1]

        curve = TGraph(len(data[:,0]), array('d', data[:,0]), array('d', data[:,1]) )
        decay_curve.append(curve)

        # The limits for the rise time fits need to be hard coded, since not all the 
        # curved start and finish at the same point.
        # I tried to develop a more independent method, by calculating the derivarived of the curve and hence finding
        # the limits of the rising slope, but since the data have significant oscillations the only way was to read this by eye
        limits_down = [-0.12E-6, -0.06E-6, -0.16E-6, -0.03E-6, 
                        -0.1E-6, -0.05E-6, -0.02E-6, -0.095E-6, 
                        -0.03E-6, -0.06E-6, -0.075E-6, -0.085E-6]
        limits_up = [0.0, 0.025E-6, -0.1E-6, 0.05E-6, 
                    0.02E-6, 0.05E-6, 0.08E-6, -0.01E-6, 
                    0.05E-6, 0.03E-6, 0.01E-6, 0.0]

        # Rise time fit
        r_fit = TF1('r_fit{:}'.format(str(j)), '[0]*(1.0 - exp(-(x)/[1])) + [2]' , limits_down[j] , limits_up[j])
        r_fit.SetParameter(0, 0.04)
        r_fit.SetParameter(1, 0.0000003)
        r_fit.SetParameter(2, 0.002)
        curve.Fit(r_fit, 'R')
        raise_fit.append(r_fit)
        rise_time.append(r_fit.GetParameter(1))
        rise_time_error.append(r_fit.GetParError(1))
        
        # Noise fit
        n_fit = TF1('n_fit_{:}'.format(str(j)), '[0]', -2.5E-6, -0.1E-6)
        n_fit.SetParameter(0, 0.0)
        curve.Fit(n_fit, 'R')
        noise_fit.append(n_fit)
        #noise_std_holesystem.append(n_fit.GetParError(0))
        noise_mean_holesystem.append(n_fit.GetParameter(0))
        noise_std_holesystem.append(np.std(data[data[:,0] < -0.5E-6,1]))
        print len(data[data[:,0] < -0.1E-6])
        print np.std(data[data[:,0] < -0.1E-6,1])

        def RiseFunct(x, par0, par1, par2):
            return par0*(1.0 - exp(-x/par1)) + par2


        def NoiseFunct(x, parN):
            return np.ones(len(x)) * parN

        xx_noise = np.linspace(-2.5E-6, -0.1E-6, 1000)
        xx_rise = np.linspace(limits_down[j], limits_up[j], 1000)
        yy_rise  = [r_fit.Eval(x) for x in xx_rise]
        yy_noise = [n_fit.Eval(x) for x in xx_noise]


        pdfP = PdfPages('../graphics/data_{:}.pdf'.format(str(j)))
        font0 = FontProperties()
        font = font0.copy()
        font.set_style('italic')
        font.set_weight('bold')
        font.set_size('x-large')
        font_par = font0.copy()
        font_par.set_style('italic')
        font_par.set_weight('medium')

        fig= plt.figure(dpi=100)
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel('Time [s]')
        ax.set_xlim(-0.7E-6, 0.5E-6)
        ax.set_ylim(-0.005, 0.05)
        ax.set_ylabel('Amplitude [V]')
        ax.plot(data[:,0], data[:,1], marker='o', color='k', label='Data', linestyle='None')
        #ax.plot(xx_noise, yy_noise, color='g', label='Noise Fit')
        ax.plot(xx_rise, yy_rise, color = 'r', label='Rise time fit', linewidth=3)
        ax.text(0.45,0.9, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                    fontproperties=font, transform=ax.transAxes)
        #ax.text(0.05,0.7, 'Noise mean = {0:.3f} [\mu V] '.format(noise_mean_holesystem[j]*1000000.), verticalalignment='bottom', horizontalalignment='left',
        #            fontproperties=font_par, transform=ax.transAxes)
        ax.text(0.05,0.65, 'Noise rms = {0:.3f} [\mu V]'.format(noise_std_holesystem[j]*1000000.), verticalalignment='bottom', horizontalalignment='left',
                    fontproperties=font_par, transform=ax.transAxes)
        ax.text(0.05,0.6, 'Rise time = {0:.3f} +/- {1:.3f} [ns]'.format(rise_time[j]*(1E+9),rise_time_error[j]*(1E+9) ) , verticalalignment='bottom', horizontalalignment='left',
                    fontproperties=font_par, transform=ax.transAxes)
        ax.legend(loc='upper left',numpoints=1)
        plt.grid()
        plt.savefig(pdfP, format='pdf')
        plt.savefig('../graphics/data_{:}.png'.format(str(j)))
        pdfP.close()
        plt.close(fig)
        plt.clf()

    # Make to time constant
    #rise_time = [rise_time[i]*2.197 for i in range(len(rise_time))]
    return rise_time, rise_time_error, noise_std_holesystem

def ExtractCapacitance(rise_time, rise_time_error, par0, e_par0, par1, e_par1, par2, e_par2):
    # Make to time constant
    rise_time = [rise_time[i]*2.197 for i in range(len(rise_time))]

    # Mean and error on the mean of the results
    mean_rise_time = np.mean(rise_time) 
    mean_error_rise_time = np.std(rise_time)/sqrt(len(rise_time))

    # Weighted mwan and weighted error on the mean of the results
    #weights = np.array(rise_time_error)**(-2)
    #weighted_rise_time = np.average(rise_time, weights=weights)
    #weighted_error_rise_time = sqrt(1.0/(np.sum(weights)))
    #print "MEan rise time: ", mean_rise_time

    # Find the value of the capacitance from the extracted parameters and fall time average
    def Capacitance(par0, par1, par2, rise_time):
        if (par1**2 -4.0*par0*(par2 - rise_time*pow(10,9)) < 0): 
            print "Negative value!!!"
            return False, 0.0
        C = (-par1 + sqrt(par1**2 -4.0*par0*(par2 - rise_time*pow(10,+9))))/(2*par0)
        return True, C

    b, C = Capacitance(par0, par1, par2, mean_rise_time)
    

    # Draw random numbers from gaussian histograms for each of the variables and then calculate for each a capacitance
    # and take the standard deviation of this
    
    C_gaus = []
    par0_gaus = []
    par1_gaus = []
    par2_gaus = []
    rise_time_gaus = []

    r = TRandom3()

    C_hist    = TH1D('C_hist', 'C_hist', 100, C*0.1, C*3 )
    par2_hist = TH1D('par2_hist', 'par2_hist', 25, par2*0.001, par2*1.999)
    par1_hist = TH1D('par1_hist', 'par1_hist', 25, par1*0.3, par1*1.7)
    par0_hist = TH1D('par0_hist', 'par0_hist', 25, par0*0.8, par0*1.2)
    risetime_hist = TH1D('risetime_hist', 'risetime_hist', 25, mean_rise_time*0.85, mean_rise_time*1.15)

    neg = 0
    for i in range(100000):
        if (i % 1000 == 0): print "MESSAGE INFO: pseudoexperiment: ", i
        par0_tmp = r.Gaus(par0, e_par0)
        par1_tmp = r.Gaus(par1, e_par1)
        par2_tmp = r.Gaus(par2, e_par2)
        rise_time_tmp = r.Gaus(mean_rise_time, mean_error_rise_time)
        b, C_tmp = Capacitance(par0_tmp, par1_tmp, par2_tmp, rise_time_tmp)

        par0_gaus.append(par0_tmp)
        par1_gaus.append(par1_tmp)
        par2_gaus.append(par2_tmp)
        rise_time_gaus.append(rise_time_tmp)
        if b == True: 
            C_hist.Fill(C_tmp)
            C_gaus.append(C_tmp)
        if b == False: neg += 1
        par0_hist.Fill(par0_tmp)
        par1_hist.Fill(par1_tmp)
        par2_hist.Fill(par2_tmp)
        risetime_hist.Fill(rise_time_tmp)

    print "MESSAGE INFO: The negative happens: ", neg

    g = TF1('fit_C_stat_uncert', 'gaus', C + 0.2*C, C -0.2*C)
    g.SetLineColor(2)
    gStyle.SetOptStat(1111)
    gStyle.SetOptFit(1111)
    C_hist.Fit(g, 'R')

    xx = np.linspace(C + 0.2*C, C -0.2*C, 1000)
    yy = [g.Eval(x) for x in xx]

    font0 = FontProperties()
    font = font0.copy()
    font.set_style('italic')
    font.set_weight('bold')
    font.set_size('x-large')
    font_par = font0.copy()
    font_par.set_style('italic')
    font_par.set_weight('medium')

    fig, ax = plt.subplots()
    ax.set_xlabel('C [pF]')
    #ax.set_xlim(0.0, 120)
    ax.set_ylabel('Pseudoexperiment')
    ax.hist(C_gaus, 100, range=(C*0.1, C*3), label='C pseudoexperiments')
    ax.text(0.7,0.7, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
    ax.plot(xx,yy, label='Gaussian Fit',color='r',linewidth=2)
    ax.text(0.7,0.65, 'C std = {0:.3f} [pF]'.format(g.GetParameter(2)), verticalalignment='bottom', horizontalalignment='left',
                    fontproperties=font_par, transform=ax.transAxes)

    ax.legend(loc='upper right',numpoints=1)
    plt.grid()
    fig.savefig('../graphics/stat_error_on_C_pseudoexperiments.pdf')
    plt.close(fig)
    plt.clf()

    e_C = g.GetParameter(2)

    return C, e_C

def ExtractENC(C, e_C, lin_par0, lin_e_par0, lin_par1, lin_e_par1):
    # Insert the ENC in the linear equation
    ENC = lin_par0 + lin_par1 * C
    # Use error propagation to find the error on the ENC
    e_ENC = sqrt( (lin_e_par0*lin_e_par0)*(1.0)  +  (lin_e_par1*lin_e_par1)*(C*C)  +  (e_C*e_C)*(lin_par1*lin_par1) )

    return ENC, e_ENC

if __name__ == "__main__":

    ######################################
    # Plot other curves measured 
    ######################################

    OtherPlots('../data/diodetest/diode_calibration.txt', '../data/diodetest/sample_name_IV_31_10_2018_14_17.txt', '../data/diodetest/diode2_CV_31_10_2018_14_23OpenNeedle.txt', '../data/diodetest/diode2_CV_31_10_2018_14_23.txt')
    
    ######################################
    # Find the capacitance of the diode from
    # the calibration curve
    ######################################
    # Find the parameters from the Calibration curve for the rise time and ENC
    lin_par0, lin_e_par0, lin_par1, lin_e_par1 = CalibrationENC('../data/diodetest/diode_calibration.txt')
    
    par0, e_par0, par1, e_par1, par2, e_par2   = CalibrationRiseTime('../data/diodetest/diode_calibration.txt')

    path = '../data/measurements/'
    # Read in the data from the Strontiom decay and fit the decay curve to an exponential to get the decay time
    datasets = ['C1_goodones00000.txt', 'C1_goodones00001.txt', 'C1_goodones00002.txt', 'C1_goodones00003.txt',
            'C1_goodones00004.txt', 'C1_goodones00005.txt', 'C1_goodones00006.txt', 'C1_goodones00007.txt',
            'C1_goodones00008.txt', 'C1_goodones00009.txt', 'C1_goodones00010.txt', 'C1_goodones00011.txt']
    rise_time, rise_time_error, noise_std_holesystem = SignalCurves(datasets, path)

    # Find the capacitance from the rise time and estimate the error
    C, e_C     = ExtractCapacitance(rise_time, rise_time_error, par0, e_par0, par1, e_par1, par2, e_par2)
    ENC, e_ENC = ExtractENC(C, e_C, lin_par0, lin_e_par0, lin_par1, lin_e_par1) 

    print "##########################"
    print "Mean rise time\t: {:} s ".format(np.mean(rise_time))
    print "Std rise time\t: {:} s".format(np.std(rise_time))
    print "##########################"
    print "##########################"
    print "Capacitance \t\t\t: {:} pF".format(C)
    print "Stat error on capacitance\t: {:} pF".format(e_C)
    print "##########################"
    print "##########################"
    print "Mean Noise RMS from diode + circuit\t: {:} V".format(np.mean(noise_std_holesystem))
    print "Std Noise RMS from diode + circuit\t: {:} V".format(np.std(noise_std_holesystem))
    print "##########################"
    print "##########################"
    print "Mean ENC in diode \t: {:} # of electrons".format(ENC)
    print "Std ENC in diode \t: {:} # of electrons".format(e_ENC)
    print "##########################"
    
    raw_input("Press enter to finish")
