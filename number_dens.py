import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

G=6.67e-11
msol = 2e30
kpc = 3.086e19
km = 1000
bins = 15
mi = np.pi/180/60

#============Ontly Variable that Needs Changing================
#============Set between 0 and 3, inclusive====================

#========================User input ======================================================|
#-------------------------------------------------Select Index--------------------------------------------------------------------------------------------------------|
#   0               1               2               3
#   Sextans     Sculptor    Fornax       Carina
gal_numb = 3
#=====================User input Complete==================================================|


cols = ['r','g','b','k']
diri ='data_files/'
filis = ['Sextans.xlsx', 'Sculptor.xlsx', 'Fornax.xlsx', 'Carina.xlsx']

gs_dist = np.array([ 86, 79, 138, 101])

rl = np.array([0.768,0.282,0.714,0.254])

dist_gal = [90,90,140,110]

num_file = diri + 'nom_dens.xlsx'

df = pd.read_excel(num_file,skiprows=1)

ri = df.iloc[:,6-2*gal_numb]
ro = df.iloc[:,7-2*gal_numb]

mi = mi*gs_dist[gal_numb]

ro = np.array(ro)
ri = np.array(ri)
ro=ro/mi/mi
ri=ri*mi/rl[gal_numb]

plt.plot(ri,ro,'.')
plt.xlabel('R/$R_e$')
plt.ylabel('Projected Number Density ($kpc^{-2}$)')
plt.title(filis[gal_numb][0:-5])

plt.xscale('log')

plt.show()