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
mi = mi*101


cols = ['r','g','b','k']
diri ='data_files/'
filis = ['Sextans.xlsx', 'Sculptor.xlsx', 'Fornax.xlsx', 'Carina.xlsx']
#filis =['Sculptor.xlsx', 'Fornax.xlsx', 'Carina.xlsx']
maxis = [245,150,90,245]
minis = [205,95,20,205]

gs_dist = np.array([ 86, 79, 138, 101])

rl = np.array([0.768,0.282,0.714,0.254])

dist_gal = [90,90,140,110]

num_file = diri + 'nom_dens.xlsx'

df = pd.read_excel(num_file,skiprows=1)

ri = df.iloc[:,6]
ro = df.iloc[:,7]

ro = np.array(ro)
ri = np.array(ri)
ro=ro/mi/mi
ri=ri*mi
plt.plot(ri,ro,'.')
plt.xlabel('R (kpc)')
plt.ylabel('Projected Number Density ($kpc^{-2}$)')

plt.show()