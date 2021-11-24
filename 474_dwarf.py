import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

G=6.67e-11
msol = 2e30
kpc = 3.086e19
km = 1000
bins = 25

diri ='data_files/'
filis = ['Sextans.xlsx', 'Sculptor.xlsx', 'Fornax.xlsx', 'Carina.xlsx']
#filis =['Sculptor.xlsx', 'Fornax.xlsx', 'Carina.xlsx']
maxis = [245,130,80,250]
minis = [205,90,30,200]

dist_gal = [90,90,140,110]

#maxis = [155,100,260]
#minis = [60,0,190]

for bb in range(0,4):
    df = pd.read_excel(diri +filis[bb])
    
    dist = np.array(df.iloc[:,-1],dtype=float)
    vel = np.array(df.iloc[:,-2],dtype=float)
    dist = dist[np.isnan(vel) == False]
    vel = vel[np.isnan(vel) == False]
     
    dist = dist[vel > minis[bb]]
    vel = vel[vel > minis[bb]]
    
    dist = dist[vel < maxis[bb]]
    d_dist = 10*dist/dist_gal[bb]
    
    vel = vel[vel < maxis[bb]]
    S_sq = ((vel-vel.mean())**2).sum()/len(vel)
    
    m = (((vel-vel.mean())*km)**2)*dist*kpc/(G*msol)
    #m=abs(vel)
    mt =[]
    sigo = []
    dsi = []
    yt = []
    rt = []
    dr = max(dist)/bins
    dt = []
    for ii in range(0,bins):
        count = 0
        sum_of = 0
        rt.insert(ii,dr*(ii+0.5))
        for jj in range(0,len(dist)):
            if ii*dr < dist[jj] and dist[jj] < (ii+1)*dr:
                count+=1
        if count >0:
            sigo.insert(ii,np.sqrt(count*S_sq/(count-1)))
            dsi.insert(ii,sigo[ii]/np.sqrt(2*(count-1)))
        else:
            sigo.insert(ii,0)
            dsi.insert(ii,0)
    sigo = np.array(sigo)
    rt = np.array(rt)
    dsi = np.array(dsi)
    '''for ii in range(0,25):
        count = 0
        sum_of = 0
        for jj in range(0,len(dist)):
            if ii*dr < dist[jj] and dist[jj] < (ii+1)*dr:
                count+=1
                sum_of += m[jj]
        if count >0:
            mt.insert(ii,sum_of/count)
        else:
            mt.insert(ii,0)
    for ii in range(0,25):
        count = 0
        sum_of = 0
        for jj in range(0,len(dist)):
            if ii*dr < dist[jj] and dist[jj] < (ii+1)*dr:
                count+=1
                sum_of += (m[jj]-mt[ii])**2
        if count >0:
            yt.insert(ii,np.sqrt(sum_of)/count)
        else:
            yt.insert(ii,0)'''
    mt = (sigo*sigo)*rt*km*km*kpc/(G*msol)
    yt = 2*sigo*dsi*rt*km*km*kpc/(G*msol)
    for ii in range(1,bins):
        dt.insert(ii-1,3*(mt[ii]-mt[ii-1])/(4*np.pi*(pow(rt[ii],3)-pow(rt[ii-1],3))))
    dt.insert(0,dt[0])
    plt.errorbar(rt,sigo,yerr=dsi,label = filis[bb][0:-5],linestyle=' ')
    plt.xlabel('Distance (kpc)')
    plt.ylabel('$\sigma$ (km/s)')
    plt.legend()
    plt.show()
    '''plt.plot(dist,vel,'.',label = filis[bb][0:-5])
    plt.xlabel('Distance (kpc)')
    plt.ylabel('$v_{gs}$ (km/s)')
    plt.legend()
    plt.show()'''
    plt.errorbar(rt,mt, yerr =yt,label = filis[bb][0:-5],linestyle = ' ')
    plt.xlabel('Distance (kpc)')
    plt.ylabel('Mass ($M_{sol}$)')
    plt.legend()
    plt.show()
    '''plt.plot(rt,dt,'.',label = filis[bb][0:-5])
    plt.xlabel('Distance (kpc)')
    plt.ylabel('Density ($M_{sol}$/$kpc^3$)')
    plt.legend()
    plt.show()'''
    
