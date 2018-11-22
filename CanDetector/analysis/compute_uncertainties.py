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
        self.K = 4.8e4   #Units of V/cm/atm
        self.K_unc = 0.3e4   #Units of V/cm/atm
        self.dU = 23.6   #Units of eV
        self.dU_unc = 5.4   #Units of eV
        
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
        # theoretical_calc.py
        # this file has measurements in mm
        sys.path.append("../scripts/")
        from theoretical_calc import cider_diameter,std_caliper
        from theoretical_calc import cider_wall,std_micrometerscrew
        from theoretical_calc import anodewire_diameter, std_micrometerscrew
        
        #outer can diameter)
        self.out_d = np.mean(cider_diameter)/10.  #units: cm
        self.out_d_unc = std_caliper/10. #units: cm
        
        #thickness
        self.thick = cider_wall/10.              # units cm
        self.thick_unc = std_micrometerscrew/10. #units: cm
        
        #Wire diameter      #units cm
        self.wire = anodewire_diameter/10.  #
        self.wire_unc = std_micrometerscrew/10.

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
        self.ra = self['wire']/2.  #units cm
        #cathode radius
        self.rc = self['out_d']/2.-self['thick']  #units cm             
    
        
        

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
    
    return V/np.log(b/a)*np.log(2)/dV*(np.log(V/(a*np.log(b/a)))-np.log(Emin*rho_ratio))


def lnm_beer(gas,geo,V):
    '''
    Compute the berrcan version of the lnm formula
    
    inputs are:
    gas: gas_property object (see definition above)
    geo: detector_geometry object (see above)
    V: operating voltage (in kV)
    '''
    a = geo['ra']
    b = geo['rc']
    p = gas['pressure']
    K = gas['K']
    dU = gas['dU']

    V = 1000.*V # voltage is in volt in this formula

    A =  V/np.log(b/a)

    B = np.log(2.)/dU

    C = np.log(V/(p*a*np.log(b/a)))-np.log(K)

    return A*B*C


if __name__=='__main__':

    import argparse

    parser = argparse.ArgumentParser("generate MC data of Gain")
    parser.add_argument('--ntrials',
                        type=int,
                        help='# of trials to run per voltage value',
                        default=1000)
    parser.add_argument('--output',
                        help='output pickle file that stores the sigma contours.',
                        default='uncertainties_lnm.p')
    args = parser.parse_args()


    gas = gas_properties()
    gas.compute_density_ratio()

    geo = cangeo()
    geo.compute_useful_geo()

    lnm_distributions = {}
    voltages = np.linspace(0.8,4.,101)

    distrib_ra = []
    distrib_rc = []
    
    for OV in voltages: #scan accross operating voltages

        lnm_distributions[OV] = []
   
        for i in range(0,args.ntrials):

            gas = gas_properties()
            gas.wigle()
            gas.compute_density_ratio()

            geo = cangeo()
            geo.wigle()
            geo.compute_useful_geo()

            distrib_ra.append(geo['ra'])
            distrib_rc.append(geo['rc'])
            
            lm = lnm(gas,geo,OV)


            lnm_distributions[OV].append(lm)

    print np.mean(distrib_ra),np.std(distrib_ra)
    print np.mean(distrib_rc),np.std(distrib_rc)

    
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
                 median,
                 voltages],open(args.output,"wb"))
    print "saved output to ",args.output

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
    plt.ylim([0,60.])
    plt.show()

            

