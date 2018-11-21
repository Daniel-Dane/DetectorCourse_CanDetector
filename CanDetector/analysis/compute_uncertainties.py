#!/usr/bin/env python

#####################################################
# Code that calculates the gas multiplication factor
# and propagates uncertainties due to the can geometry
# and environmental condition
#
#####################################################
import numpy as np
import sys


# Store the gas properties
# Source: Blum, Riegler, Rolandi. Particle Detection with drift chambers.
#         page 136, table 4.2

class gas_properties:

    def __init__(self):
        self.name = 'P10'
        self.Emin = 48.   # units: kV/cm
        self.Emin_unc =3.
        self.dv = 23.6 /1000.# units: kV
        self.dv_unc = 5.4/1000.
        self.pressure = 1.0 # Units = fraction of standard pressure (p/[1 bar])
        self.pressure_unc = 0.001 # Based on monitoring over a day
        self.temperature = 1.0# Units = fraction of standard temperature (T/[273.15 K])

    def __getitem__(self,key):
        return self.__dict__[key]
    
        
    def items(self):
        return self.__dict__.items()
    
    def __setitem__(self, key, item):
        self.__dict__[key] = item
        
    def wigle(self,list_to_wigle = None):
        #wigle the geometry according to
        #the quoted uncertainties

        # if no list of element to wigle is provided,
        # will wigle all elements that are not uncertainties
        if list_to_wigle is None:
            list_to_wigle = []
            for element in self.__dict__.keys():
                if '_unc' not in element: 
                    list_to_wigle.append(element)

        for k in list_to_wigle:
            if k+'_unc' in self.__dict__.keys():
                # normal-distributed uncertainty
                wig = np.random.normal(loc = self[k],
                                       scale = self[k+'_unc'])
                self[k] = wig


    def compute_density_ratio(self):
        self.rho = self.pressure/self.temperature


    

class cangeo:

    def __init__(self):
        #initialize with the values in
        # our lab report (can values)
        #outer can diameter)
        self.out_d = 8.     #units: cm
        self.out_d_unc = 0.001
        #thickness
        self.thick = 0.0105 # units cm
        self.thick_unc = 0.001
        #Wire diameter      #units cm
        self.wire = 0.005
        self.wire_unc = 0.001

    def __getitem__(self,key):
        return self.__dict__[key]
    
    def items(self):
        return self.__dict__.items()

    def __setitem__(self, key, item):
        self.__dict__[key] = item
        

    def wigle(self,list_to_wigle = None):
        #wigle the geometry according to
        #the quoted uncertainties

        # if no list of element to wigle is provided,
        # will wigle all elements that are not uncertainties
        if list_to_wigle is None:
            list_to_wigle = []
            for element in self.__dict__.keys():
                if '_unc' not in element: 
                    list_to_wigle.append(element)
            
        for k in list_to_wigle:
            if k+'_unc' in self.__dict__.keys():
                # normal-distributed uncertainty
                wig = np.random.normal(loc = self[k],
                                       scale = self[k+'_unc'])
                self[k] = wig
        

    # uncertainties shall be propagated using the wigle method
    # prior to calling this function
    
    def compute_useful_geo(self):
        #anode radius
        self.ra = self.wire/2.  #units cm
        #cathode radius
        self.rc = self.out_d/.2-self.thick  #units cm             
    
        
        

def lnm(gas,geo,V):
    '''
    Compute the natural log of
    the multiplication factor.
    
    inputs are:
    gas: gas_property object (see definition above)
    geo: detector_geometry object (see above)
    V: operating voltage (in kV)
    '''
    a = geo['ra']
    b = geo['rc']

    rho_ratio = gas['rho']
    Emin = gas['Emin']
    dV = gas['dv']
    
    return np.log(V)/np.log(b/a)*np.log(2)/dV*(np.log(V/(a*np.log(b/a)))-np.log(Emin*rho_ratio))


if __name__=='__main__':

    ntrials = 1000000

    lnm_distributions = {}
    voltages = np.linspace(1.3,4.,101)
    
    for OV in voltages: #scan accross operating voltages

        lnm_distributions[OV] = []
   
        for i in range(0,ntrials):

            gas = gas_properties()
            gas.wigle()
            gas.compute_density_ratio()

            geo = cangeo()
            geo.wigle()
            geo.compute_useful_geo()

            lm = lnm(gas,geo,OV)


            lnm_distributions[OV].append(lm)

    # plot everything

    import matplotlib.pyplot as plt
    
    ten_pc = []
    one_pc = []
    ninety = []
    ninety9= []
    median = []
    
    X = []
    for k,v in sorted(lnm_distributions.items()):
        X.append(k)
        v = np.array(sorted(v))

        ten_pc.append(np.percentile(v,10))
        one_pc.append(np.percentile(v,1))
        ninety.append(np.percentile(v,90))
        ninety9.append(np.percentile(v,99))
        median.append(np.median(v))

    import pickle
    pickle.dump([ten_pc,
                 one_pc,
                 ninety,
                 ninety9,
                 median],open("uncertainties_lnm_sub_p0.001_p_uncertainty.p","wb"))

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

            

