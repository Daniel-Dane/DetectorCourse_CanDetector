#!/usr/bin/env python3

import logging
import argparse
import sys
import math
from numpy import genfromtxt
from rootpy.plotting import Hist
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
import matplotlib
matplotlib.use('Agg')

"""
all of below need to be loaded after `matplotlib` and `matplotlib.use('Agg')`
"""
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import rootpy.plotting.root2matplotlib as rplt

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
    parser.add_argument("--file", type=str,
                        help="file with CSV data separated by whitespace",
                        required=True)
    args = parser.parse_args()

    return args


def mca_to_hist(filename):
    roi_on = False
    data_on = False
    data_line = 1
    bins = 1024
    r_min = 0
    r_max = 0

    h = Hist(1024, 0, 1024)

    for line in open(filename, 'r'):
        if roi_on:
            val = line.split()
            r_min = val[0]
            r_max = val[1]
            roi_on = False

        if data_on:
            h.SetBinContent(data_line, line.split()[0])
            data_line += 1

        if line.find("ROI") > -1:
            roi_on = True

        if line.find("DATA") > -1:
            data_on = True

    print("region of interest: {:} to {:}".format(r_min, r_max))

    return h


def main(argv):
    args = parseArguments(argv)

    data = genfromtxt(args.file)  # By default, any consecutive whitespaces
                                  # act as delimiter.

    d_volt = data[:, 1]
    d_volt_unc = [2 for x in d_volt]
    d_mean = data[:, 2]
    d_mean_unc = data[:, 4]
    d_fwhm = data[:, 3]
    d_fwhm_unc = [0.2*x for x in d_fwhm]

    am_confs = ["100_1136", "100_1191", "100_1244", "100_1297", "100_1351", "100_1399", "10_1665", "10_1712", "10_1758", "20_1559", "20_1603", "20_1666", "2_1900", "2_1951", "2_2001", "40_1397", "40_1455", "40_1502", "40_1559", "4_1757", "4_1808", "4_1858", "4_1899"]

    print(mca_to_hist("../data/mca/am_"+am_confs[5]+".mca"))


    d_y = []
    d_y_unc = []

    for idx in range(len(d_volt)):
        d_y.append(d_fwhm[idx]/d_mean[idx])
        rel_unc_mean = d_mean_unc[idx]/d_mean[idx]
        rel_unc_fwhm = d_fwhm_unc[idx]/d_fwhm[idx]
        d_y_unc.append(math.sqrt(math.pow(rel_unc_mean, 2)+
                                 math.pow(rel_unc_fwhm, 2))*d_y[idx])

    pp = PdfPages(args.file.replace(".txt", "")+'.pdf')
    font0 = FontProperties()
    font = font0.copy()
    font.set_style('italic')
    font.set_weight('bold')
    font.set_size('x-large')
    font_wip = font0.copy()
    font_wip.set_style('italic')
    font_wip.set_weight('medium')

    plt.figure(figsize=(10, 8), dpi=300)
    plt.xlabel('xtick', fontsize=5)
    plt.title("test")
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

        ax.errorbar(d_volt, d_y, d_y_unc, d_volt_unc, linestyle='None', color='red')

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



if __name__ == "__main__":
    main(sys.argv[1:])
