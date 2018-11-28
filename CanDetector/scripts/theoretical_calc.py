import numpy as np

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

coppertube_frontend_length_tot = np.array( [36.82, 36.76] )
coppertube_frontend_length_top = np.array( [15.91, 15.94, 15.86] )
coppertube_frontend_tube_radius = np.array( [19.84, 19.85, 19.86] )

coppertube_backend_length_tot = np.array( [15.04, 14.97, 14.89, 14.93, 14.95, 15.06] )
coppertube_backend_length_top = np.array( [5.05, 4.98, 4.95] )
coppertube_backend_tube_radius = np.array( [19.88, 19.86] )

brasstube_diameter = np.array( [1.0] )
anodewire_coppertube_diameter = np.array( [0.025] )

# to get mean and std for a variable
#coppertube_inner_diameter.mean()
#coppertube_inner_diameter.std(ddof=1)


print("can detector")
print( "${:.2f} \pm {:.2f}$".format(cider_diameter.mean(),cider_diameter.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(cider_wall.mean(),cider_wall.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(plastictube_diameter.mean(),plastictube_diameter.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(brasstube_diameter.mean(),brasstube_diameter.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(nylonscrew_diameter.mean(),nylonscrew_diameter.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(HV_connector_diameter.mean(),HV_connector_diameter.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(anodewire_diameter.mean(),anodewire_diameter.std(ddof=1)) )







print("copper tube")
print( "${:.2f} \pm {:.2f}$".format(coppertube_length.mean(),coppertube_length.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(coppertube_wall.mean(),coppertube_wall.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(coppertube_inner_diameter.mean(),coppertube_inner_diameter.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(coppertube_radiationhole.mean(),coppertube_radiationhole.std(ddof=1)) )

print( "${:.2f} \pm {:.2f}$".format(coppertube_frontend_length_tot.mean(),coppertube_frontend_length_tot.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(coppertube_frontend_length_top.mean(),coppertube_frontend_length_top.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(coppertube_frontend_tube_radius.mean(),coppertube_frontend_tube_radius.std(ddof=1)) )

print( "${:.2f} \pm {:.2f}$".format(coppertube_backend_length_tot.mean(),coppertube_backend_length_tot.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(coppertube_backend_length_top.mean(),coppertube_backend_length_top.std(ddof=1)) )
print( "${:.2f} \pm {:.2f}$".format(coppertube_backend_tube_radius.mean(),coppertube_backend_tube_radius.std(ddof=1)) )

print( "${:.2f} \pm {:.2f}$".format(brasstube_diameter.mean(),brasstube_diameter.std(ddof=0)) )
print( "${:.2f} \pm {:.2f}$".format(anodewire_coppertube_diameter.mean(),anodewire_coppertube_diameter.std(ddof=0)) )