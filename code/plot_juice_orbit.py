import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice

import sys
sys.path.append('../lib')
import juice_lib

# SPICE test
spice.tkvrsn('TOOLKIT')

# load SPICE ketnels
spice_dir = 'C:/share/Linux/doc/spice\kernels\\'
kaguya_lib.spice_ini(kaguya_spice_dir)

# set date/time
utctim='2008-07-02T07:00:00'    # start date/time
et_ex=spice.str2et(utctim)      # seconds
print('ET:', et_ex)
nd = 1440 # number of data
dt= 60   # time step [second]
et = et_ex + dt * np.arange(0, nd)

# calculate Kaguya position with spice
x, y, z, r, lat, lon = kaguya_lib.get_kaguya_pos(et)

# plot
plt.scatter(lon,lat,s=10)
plt.show()
