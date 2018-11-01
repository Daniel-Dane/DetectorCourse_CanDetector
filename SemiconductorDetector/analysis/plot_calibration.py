#!/usr/bin/env python

capacitance = []

signal = []
unc_signal = []

noise = []
unc_noise = []

rise_time = []
unc_rise_time = []


with open("../data/diodetest/diode_calibration.txt") as fin:

    fin.readline()
    for line in fin:
        line = line.split()
        c = float(line[0])
        n = float(line[1])/1000.
        un= float(line[2])/1000.
        rt= float(line[3])
        rtu=float(line[4])
        s = float(line[5])
        us= float(line[6])

        capacitance.append(c)
        signal.append(s)
        noise.append(n/1000.)
        unc_noise.append(un/1000.)
        unc_signal.append(us)
        rise_time.append(rt)
        unc_rise_time.append(rtu)
        

import matplotlib.pyplot as plt
import numpy as np
plt.errorbar(capacitance,signal,yerr=unc_signal,fmt='o')
plt.xlabel("Capacitance (pF)")
plt.ylabel("Signal (mV)")
plt.show()

plt.errorbar(capacitance,noise,yerr=unc_noise,fmt='rs')
plt.xlabel("Capacitance")
plt.ylabel("RMS noise (mV)")
plt.show()


plt.plot(np.array(capacitance),np.array(signal)/np.array(noise),'gs')
plt.xlabel("capacitance (pF)")
plt.ylabel("ENC")
plt.show()

plt.errorbar(capacitance,rise_time,yerr=unc_rise_time,fmt='^m')
plt.xlabel("capacitance (pF)")
plt.ylabel("rise time (ns)")
plt.show()

    
