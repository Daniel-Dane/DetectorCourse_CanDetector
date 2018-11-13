# milimeters
std_caliper = 0.05 #?
std_micrometerscrew = 0.005 # micrometer screw gauge

cider_diameter        = np.array( [66.04, 65.50, 65.70, 65.90, 65.98] )
cider_wall            = np.array( [0.105] ) # micrometerscrew
plastictube_diameter  = np.array( [5.89, 5.95, 6.00, 6.01, 6.01] )
brasstube_diameter    = np.array( [1.0] )
nylonscrew_diameter   = np.array( [7.8, 7.75] )
HV_connector_diameter = np.array( [9.37] )
anodewire_diameter    = np.array( [50] )

coppertube_length     = np.array( [149.06, 149.10, 149.06, 149.09] )
coppertube_wall       = np.array( [1.04, 1.00, 1.02, 1.02, 1.02, 1.14, 1.05, 1.04, 1.04] )
coppertube_inner_diameter = np.array( [19.97, 19.98, 19.99, 19.88, 19.70, 19.80, 19.96, 20.00] )
coppertube_radiationhole = np.array( [5.03, 5.04, 5.08] )
# all the endcaps for the copper tube go here <--
brasstube_diameter = np.array( [1.0] )
anodewire_coppertube_diameter = np.array( [0.025] )

# to get mean and std for a variable
#coppertube_inner_diameter.mean()
#coppertube_inner_diameter.std(ddof=1)
