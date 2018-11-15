#!/usr/bin/env python2

import logging
import argparse
import sys
import math
from numpy import genfromtxt
from rootpy.plotting import Hist
from ROOT import TF1, kRed
from matplotlib import gridspec
import matplotlib
matplotlib.use('Agg')

"""
all of below need to be loaded after `matplotlib` and `matplotlib.use('Agg')`
"""
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt

from common import mca_to_hist, font, font_wip

#######################################################################
# Package information
#######################################################################
__author__ = "Fabian A.J. Thiele"
__copyright__ = "2018, Thiele"
__credits__ = ["Fabian Thiele"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Fabian Thiele"
__email__ = "fabian.thiele@cern.ch"
__status__ = "Development"


#######################################################################
# configure logging
#######################################################################
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(name)-12s '
                           '%(levelname)-8s %(message)s',
                    datefmt='%d-%m-%y %H:%M')
console = logging.StreamHandler()
log = logging.getLogger('plotHVscan')


def parseArguments(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sample', nargs='+', help='Choose `fe` and/or `am`', default=['am', 'fe'])
    args = parser.parse_args()

    return args


def plotSpec(name, h, r_min=0, r_max=0):
    pp = PdfPages(name)

    plt.figure(figsize=(15, 8), dpi=300)
    plt.xlabel('xtick', fontsize=5)
    plt.ylabel('ytick', fontsize=5)
    gs = gridspec.GridSpec(1, 1)

    gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0],
                                           width_ratios=[1],
                                           height_ratios=[1],
                                           wspace=0.1,
                                           hspace=0.2)

    ax1 = plt.subplot(gs1[:, :])
    axes = [ax1]

    for ax in axes:

        # name axes and set maximum of y-axis to accomodate
        # all histograms inside histstack
        ax.set_ylabel("Counts")
        ax.set_xlabel("Channel]")
        ax.set_title(name.replace(".pdf", ""))

        #ax.errorbar(d_volt, d_y, d_y_unc, d_volt_unc, linestyle='None', color='red')
        f = TF1("tf_{:}".format(name), "gaus", float(r_min), float(r_max))
        h.Fit("tf_{:}".format(name), "LLR") # , "R")
        rplt.hist(h, stacked=False, fill=False, axes=ax)

        print("constant = {:}".format(f.GetParameter(0)))
        print("mean = {:}".format(f.GetParameter(1)))
        print("mean unc = {:}".format(f.GetParError(1)))
        print("std = {:}".format(f.GetParameter(2)))
        fwhm = f.GetParameter(2)*2.*math.sqrt(2*math.log(2))
        print("FWHM = {:}".format(fwhm))
        print("eval for r_min = {:}: {:}".format(r_min, f.Eval(float(r_min))))
        xx = range(int(r_min), int(r_max))
        yy = [f.Eval(x) for x in xx]
        plt.plot(xx, yy, 'r')
        ax.axvline(x=r_min)
        ax.axvline(x=r_max)


        # make ATLAS fonts in plot
        ax.text(0.2, 0.9, '1-Group',
                verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
        ax.text(0.2, 0.87, 'Data',
                verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font_wip, transform=ax.transAxes)
        ax.text(0.8, 0.9, "Mean = {:.2f}".format(f.GetParameter(1)),
                verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font_wip, transform=ax.transAxes)
        ax.text(0.8, 0.87, "FWHM = {:.2f}".format(fwhm),
                verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font_wip, transform=ax.transAxes)

        ax.legend(fontsize='small')
        ax.grid()

    plt.show()
    plt.savefig(pp, format='pdf')
    pp.close()

    return (f.GetParameter(1), f.GetParameter(2), fwhm, f.GetParError(1), f.GetParError(2)*2.*math.sqrt(2*math.log(2)), name.split("_")[1],f.GetParError(1))

def plot_confs(confs, title, abbrev):

    fit_pars = {}
    for v in confs:
        name = "../data/mca/"+abbrev+"_"+v+".mca"
        path = "../graphics/"+abbrev+"_"+v+".pdf"
        (h, r_min, r_max, time) = mca_to_hist(name)
        print(v.split("_"))
        fit_pars[int(v.split("_")[1])] = plotSpec(path, h, r_min, r_max)

    gain_style = {100: 'blue',
                  10: 'red',
                  20: 'pink',
                  2: 'purple',
                  40: 'orange',
                  4: 'green'}

    d_y = {}
    d_y_unc = {}
    d_volt = {}

    for gain in gain_style:
        d_volt[gain] = []
        d_y_unc[gain] = []
        d_y[gain] = []

    for key, value in fit_pars.items():
        # mean, std, fwhm, mean_err, fwhm_err
        gain = int(value[5])
        d_y[gain].append(value[2]/value[0])
        rel_unc_mean = value[3]/value[0]
        d_volt[gain].append(key)
        print("volt",key," mean",value[0]," FWHM",value[2], " mean unc",value[6])

        rel_unc_fwhm = value[4]/value[2]
        rel_unc_fwhm = math.sqrt(math.pow(rel_unc_fwhm, 2) + math.pow(0.1, 2))
        d_y_unc[gain].append(math.sqrt(math.pow(rel_unc_mean, 2)+
                                       math.pow(rel_unc_fwhm, 2))*value[2]/value[0])

    pp = PdfPages("../graphics/"+title.lower()+"_scan.pdf")


    plt.figure(figsize=(15, 8), dpi=300)
    plt.xlabel('xtick', fontsize=5)
    plt.ylabel('ytick', fontsize=5)
    gs = gridspec.GridSpec(1, 1)

    gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0],
                                           width_ratios=[1],
                                           height_ratios=[1],
                                           wspace=0.1,
                                           hspace=0.2)

    ax1 = plt.subplot(gs1[:, :])
    axes = [ax1]

    for ax in axes:

        # name axes and set maximum of y-axis to accomodate
        # all histograms inside histstack
        ax.set_ylabel("FWHM/Mean")
        ax.set_xlabel("Voltage [V]")
        ax.set_title(title)

        for key, value in gain_style.items():
            # mean, std, fwhm, mean_err, fwhm_err
            ax.errorbar(d_volt[key], d_y[key], d_y_unc[key], marker='o', color=gain_style[key], linestyle='None')

        # define a legend
        l = []
        for key in sorted(gain_style.keys(), reverse=True):
             l.append("Gain: {:}".format(key))
        ax.legend(l, numpoints=1)


        # make ATLAS fonts in plot
        ax.text(0.2, 0.9, '1-Group',
                verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font, transform=ax.transAxes)
        ax.text(0.2, 0.87, 'Data',
                verticalalignment='bottom', horizontalalignment='left',
                fontproperties=font_wip, transform=ax.transAxes)

        ax.legend(fontsize='small')
        ax.grid()

    plt.show()
    plt.savefig(pp, format='pdf')
    pp.close()



def main(argv):
    args = parseArguments(argv)
    
    print("am_40_1502 fit is not exactly right")

    am_confs = ["100_1136", "100_1191", "100_1244", "100_1297", "100_1351", "100_1399", "10_1665", "10_1712", "10_1758", "20_1559", "20_1603", "20_1666", "2_1900", "2_1951", "2_2001", "40_1397", "40_1455", "40_1502", "40_1559", "4_1757", "4_1808", "4_1858", "4_1899"]

    fe_confs = ["100_1420", "100_1470", "100_1523", "100_1572", "100_1617", "100_1717", "10_1978", "10_2000", "10_2042", "10_2080", "20_1901", "20_1978", "2_2201", "2_2254", "2_2303", "40_1717", "40_1801", "40_1853", "40_1901", "4_2080", "4_2124", "4_2201"]

    if 'am' in args.sample:
        plot_confs(am_confs, "Americium", "am")
    if 'fe' in args.sample:
        plot_confs(fe_confs, "Iron", "fe")


if __name__ == "__main__":
    main(sys.argv[1:])
