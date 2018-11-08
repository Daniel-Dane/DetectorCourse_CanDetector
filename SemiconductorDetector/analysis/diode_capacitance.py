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

def Calibration(calibration_filename, plot_Calibration_curve):
    # Read data fron the calibration curve dataset
    calibration_diode = np.loadtxt(calibration_filename,skiprows=1)
    C           = calibration_diode[:,0]
    e_C         = np.zeros([len(C)])
    noise_rms   = calibration_diode[:,1]
    d_noise_rms = calibration_diode[:,2]
    fall_time   = calibration_diode[:,3]
    d_fall_time = calibration_diode[:,4]
    amplitude   = calibration_diode[:,5]
    d_amplitude = calibration_diode[:,6]

    gStyle.SetOptStat(1111)

    # Construct the  curve
    calib_curve = TGraphErrors(len(C), array('d',C.tolist()), array('d', fall_time.tolist()), array('d', e_C.tolist()), array('d',d_fall_time.tolist()) )
    calib_curve.SetMarkerColor(kBlack)
    calib_curve.SetMarkerSize(1)
    calib_curve.SetMarkerStyle(20)
    calib_curve.GetXaxis().SetTitle('C [pF]')
    calib_curve.GetYaxis().SetTitle('Rise time [ns]')

    # Construct the fit with a 2nd order polynomial
    fit_curve = TF1('fit_curve','[0]*x*x + [1]*x + [2]',0.,110.)
    fit_curve.SetParameter(0, 0.0)
    fit_curve.SetParameter(1, 0.0)
    fit_curve.SetParameter(2, 0.0)
    fit_curve.SetLineColor(2)

    # FIt the calibration curve
    calib_curve.Fit(fit_curve,'R')

    # Plot the curve and the fit with ROOT
    #c = TCanvas('c', 'c', 10, 10, 600, 600)
    #calib_curve.Draw('APE')
    #fit_curve.Draw('SAME')
    #c.Update()
    #c.SaveAs('calib.pdf')

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
    if plot_Calibration_curve:
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
        ax.legend(loc='upper left')
        plt.grid()
        #plt.show()
        fig.savefig('../graphics/calibration_diode.pdf')
        plt.close(fig)
        plt.clf()

    return par0, e_par0, par1, e_par1, par2, e_par2

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
        curve = TGraph(len(data[:,0]), array('d', data[:,0]), array('d', data[:,1]) )
        #curve.SetMarkerColor(kBlack)
        #curve.SetMarkerSize(0.5)
        #curve.SetMarkerStyle(20)
        #curve.SetTitle('')
        #curve.GetXaxis().SetTitle('Time')
        #curve.GetYaxis().SetTitle('Amplitude')
        #curve.GetXaxis().SetRangeUser(-0.15, 0.15)
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
        noise_std_holesystem.append(np.std(data[data[:,0] < -0.1E-6,1]))
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
        ax.plot(xx_rise, yy_rise, color = 'r', label='Rise time fit')
        ax.text(0.45,0.9, 'Group 1', verticalalignment='bottom', horizontalalignment='left',
                    fontproperties=font, transform=ax.transAxes)
        #ax.text(0.05,0.7, 'Noise mean = {0:.3f} [\mu V] '.format(noise_mean_holesystem[j]*1000000.), verticalalignment='bottom', horizontalalignment='left',
        #            fontproperties=font_par, transform=ax.transAxes)
        ax.text(0.05,0.65, 'Noise rms = {0:.3f} [\mu V]'.format(noise_std_holesystem[j]*1000000.), verticalalignment='bottom', horizontalalignment='left',
                    fontproperties=font_par, transform=ax.transAxes)
        ax.text(0.05,0.6, 'Rise time = {0:.3f} +/- {1:.3f} [ns]'.format(rise_time[j]*(1E+9),rise_time_error[j]*(1E+9) ) , verticalalignment='bottom', horizontalalignment='left',
                    fontproperties=font_par, transform=ax.transAxes)
        ax.legend(loc='upper left')
        plt.grid()
        plt.savefig(pdfP, format='pdf')
        #plt.savefig('data_{:}.png'.format(str(j)))
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
    #c5 = TCanvas('c5', 'c5', 10, 10, 600, 600)
    #C_hist.GetXaxis().SetTitle('C [pF]')
    #C_hist.Draw('HIST')
    #g.Draw('SAME')
    #c5.Update()

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

    ax.legend(loc='upper right')
    plt.grid()
    #plt.show()
    fig.savefig('../graphics/stat_error_on_C_pseudoexperiments.pdf')
    #plt.show()
    plt.close(fig)
    plt.clf()

    """
    c6 = TCanvas('c6', 'c6', 10, 10, 600, 600)
    par0_hist.GetXaxis().SetTitle('par0')
    par0_hist.Draw('HIST')
    c6.Update()

    c7 = TCanvas('c7', 'c7', 10, 10, 600, 600)
    par1_hist.GetXaxis().SetTitle('par1')
    par1_hist.Draw('HIST')
    c7.Update()

    c8 = TCanvas('c8', 'c8', 10, 10, 600, 600)
    par2_hist.GetXaxis().SetTitle('par2')
    par2_hist.Draw('HIST')
    c8.Update()

    c9 = TCanvas('c9', 'c9', 10, 10, 600, 600)
    risetime_hist.GetXaxis().SetTitle('rise time')
    risetime_hist.Draw('HIST')
    c9.Update()
    """

    e_C = g.GetParameter(2)

    return C, e_C

if __name__ == "__main__":
    # Find the parameters from the Calibration curve
    #sys.append("../data/diodetest")
    par0, e_par0, par1, e_par1, par2, e_par2 = Calibration('../data/diodetest/diode_calibration.txt', 1)

    path = '../data/measurements/'
    # Read in the data from the Strontiom decay and fit the decay curve to an exponential to get the decay time
    datasets = ['C1_goodones00000.txt', 'C1_goodones00001.txt', 'C1_goodones00002.txt', 'C1_goodones00003.txt',
            'C1_goodones00004.txt', 'C1_goodones00005.txt', 'C1_goodones00006.txt', 'C1_goodones00007.txt',
            'C1_goodones00008.txt', 'C1_goodones00009.txt', 'C1_goodones00010.txt', 'C1_goodones00011.txt']
    rise_time, rise_time_error, noise_std_holesystem = SignalCurves(datasets, path)


    # Find the capacitance from the rise time and estimate the error
    C, e_C = ExtractCapacitance(rise_time, rise_time_error, par0, e_par0, par1, e_par1, par2, e_par2)

    print "##########################"
    print "Capacitance: {:} pF".format(C)
    print "Stat error on capacitance: {:} pF".format(e_C)
    print "##########################"

    raw_input("Press enter to finish")

