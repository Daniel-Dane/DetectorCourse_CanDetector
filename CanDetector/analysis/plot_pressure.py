#!/usr/bin/env python



import numpy as np
import datetime

data = open('../data/pressure_data_helsinki_oct31_nov2_2018.csv')
data.readline()

X = []
P = []

for line in data:
    d = line.split(',')

    X.append(datetime.datetime(int(d[0]),
                               int(d[1]),
                               int(d[2]),
                               int(d[3].split(':')[0]),
                               int(d[3].split(':')[1])))
    P.append(float(d[-1])/1000.)


import matplotlib.pyplot as plt

plt.plot(X,P)
plt.xlabel('Time')
plt.ylabel('Pressure (w.r.t 1 bar)')
plt.show()
