import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd


#========================User input ======================================================|
#-------------------------------------------------Select Number of Bins------------------------------------------------------------------------------------------|
bins = 15
#---------------------------------Select Desired Variabe---------------------------------------------------------------------------------------------------------| 
#   1   =   Mass
#   2   =   Velocity
#   3   =   Density
des_var = 1
#---------------------------------Include Pisces and Pegasus (True or False)--------------------------------------------------------------------------|
include_pis = False
#---------------------------------Log Y Scale-------------------------------------------------------------------------------------------------------------------------|
log_y = False
#=====================User input Complete==================================================|

#=======================Calculate Distances if Needed=========================================
def eri_proto(df_j,distj,vcol):
    
    vel_star = df_j.iloc[:,vcol]
    dec_star = df_j.iloc[:,3]
    ra_star = df_j.iloc[:,2]
    for iij in range(0,len(vel_star)):
        if '+' in str(vel_star[iij]):
            vel_star[iij] = str(vel_star[iij]).split("+")[0]
    raj = np.array(ra_star,dtype=float)/180*np.pi
    raj = raj - np.mean(raj)
    decj = np.array(dec_star,dtype=float)/180*np.pi
    decj = decj - np.mean(decj)
    rj = np.sqrt(decj**2 + raj**2)*distj
    
    vel_j = np.array(vel_star,dtype=float)
    return rj, vel_j

#=================Function for Getting Data from Excel==========================================    
def get_draco(filisi):
    df_i = pd.read_excel(filisi)
    if 'erida' in filisi:
        radi, veli = eri_proto(df_i,366,7)
    elif 'Draco' in filisi:
        disti = 80
        radi = np.array(df_i.iloc[:,17],dtype = float)
        veli = -np.array(df_i.iloc[:,9],dtype = float)
        radi = radi/180*np.pi*disti
    elif 'Antli' in filisi:
        radi,veli = eri_proto(df_i,130,4)
    elif 'Crater' in filisi:
        radi, veli = eri_proto(df_i,118,7)
    else:
        disti = 210
        radi = np.array(df_i.iloc[:,7],dtype = float)
        veli = np.array(df_i.iloc[:,8],dtype = float)
        radi = radi/180/60*np.pi*disti
    return radi, veli

#========================Some Constants==================================================
G=6.67e-11
msol = 2e30
kpc = 3.086e19
km = 1000

#========================Data File Locations=============================================
diri ='data_files/'
filis = ['Sextans.xlsx', 'Sculptor.xlsx', 'Fornax.xlsx', 'Carina.xlsx','Draco.xlsx', 'LeoII.xlsx','AntliaII.xlsx','eridanas.xlsx','CraterII.xlsx','LeoI.xlsx']

#=================Initial Velocity Locations from Visual Inspection==================================
maxis = [245,150,90,245,1000,105,310,100,110,325]
minis = [205,95,20,205,0,55,270,40,80,240]

#================Distance to Galaxy and 2D Half Light Radius (in KPC)=========================
gs_dist = np.array([ 86, 79, 138, 101,80,210,130,366,118,250])
rl = np.array([0.768,0.282,0.714,0.254,0.221,0.176,2.8,0.28,1.07,0.251])

#==========Lower Res Galaxy Distances, Used in Some of the Data files=============================
#===========These are corrected to the more accurate values later on and are no longer used===============
dist_gal = [90,90,140,110]

#================Iterate Over all Galaxies=================================================
for bb in range(0,10):
    #=========================Get Data==========================================
    if bb > 3 and bb < 9:
        dist, vel = get_draco(diri+filis[bb])
    elif bb < 4:
        df = pd.read_excel(diri +filis[bb])
    
        dist = np.array(df.iloc[:,-1],dtype=float)
        dist = dist/dist_gal[bb]*gs_dist[bb]
    
        vel = np.array(df.iloc[:,-2],dtype=float)
    else:
        df = pd.read_excel('data_files/mateo_table.xls')
        dist = np.array(df.iloc[:,1],dtype = float)
        vel = np.array(df.iloc[:,3],dtype = float)
        dist = dist[np.isnan(vel) == False]
        vel = vel[np.isnan(vel) == False]
        dist = dist*250*np.pi/(180*3600)
        
    vel=vel[dist<1.5]
    dist=dist[dist<1.5]
    dist = dist[np.isnan(vel) == False]
    vel = vel[np.isnan(vel) == False]
    
    dist = dist[vel > minis[bb]]
    vel = vel[vel > minis[bb]]
    
    dist = dist[vel < maxis[bb]]
    d_dist = 10*dist/gs_dist[bb]
    
    vel = vel[vel < maxis[bb]]
    
    S_sq = 0
    
    #===========Eliminate bad values for velocity====================================
    for ii in range(0,5):
        main_sig = sum(((vel - vel.mean())*km)**2)/len(vel)
        dmain_sig = main_sig/np.sqrt(2*(len(vel)-1))
        vel_all = np.sqrt(3)*np.sqrt(main_sig)/km
        dvel_all = np.sqrt(3)*dmain_sig/(2*np.sqrt(main_sig))/km
        dist = dist[abs(vel-vel.mean())<np.sqrt(2)*vel_all]
        vel = vel[abs(vel-vel.mean())<np.sqrt(2)*vel_all]
    
    
    #==============Initialize arrays======================================
    mt =[]
    sigo = []
    dsi = []
    yt = []
    rt = []
    dr = max(dist)/bins
    dt = []
    ddt = []
    
    #=====================Get Velocity Dispersion and Error===============================
    for ii in range(0,bins):
        count = 0
        sum_of = 0
        rt.insert(ii,dr*(ii+0.5))
        vel_m = vel[ii*dr < dist]
        dist_m = dist[ii*dr < dist]
        vel_m = vel_m[dist_m < (ii+1)*dr]
        count = len(vel_m)
        S_sq = sum((vel_m-vel_m.mean())**2)
        if count >0:
            S_sq = S_sq/count
            sigo.insert(ii,np.sqrt(count*S_sq/(count-1)))
            dsi.insert(ii,sigo[ii]/np.sqrt(2*(count-1)))
        else:
            sigo.insert(ii,0)
            dsi.insert(ii,0)
    sigo = np.array(sigo)
    rt = np.array(rt)
    dsi = np.array(dsi)
    vel2 = np.sqrt(3)*sigo
    dvel = np.sqrt(3)*dsi

    #=============================Get Masses========================================
    mt = 4*(sigo*sigo)*rt*km*km*kpc/(G*msol)
    yt = 8*sigo*dsi*rt*km*km*kpc/(G*msol)
    
    #===================================Get Density=================================
    for ii in range(1,bins):
        dt.insert(ii-1,(mt[ii]-mt[ii-1])/(np.pi*(pow(rt[ii],2)-pow(rt[ii-1],2))))
        ddt.insert(ii-1,np.sqrt(yt[ii]**2+yt[ii-1]**2)/(np.pi*(pow(rt[ii],2)-pow(rt[ii-1],2))))
    dt.insert(0,dt[0])
    ddt.insert(0,ddt[0])
    dt = np.array(dt)
    ddt = np.array(ddt)
    
    #================Get mass from total velocity dispersion===================
    mass_o = 4*rl[bb]*main_sig*kpc/G/msol
    
    
    #=========Get Unitless Density=========================================================
    d_half = mass_o/(np.pi*rl[bb]**2)
    
    dt = dt/d_half
    ddt = ddt/d_half
    
    #==================Plot mass, velocity, or density============================================
    if bb != 0: #Doesn't include Sextans, change to arbitrarily high number to get Sextans
        
        #=============Mass
        if des_var ==1:
            plt.errorbar(rt[0:int(bins)],mt[0:int(bins)], yerr =yt[0:int(bins)],label = filis[bb][0:-5],linestyle = ' ',color='k',marker='.')
            mass_check = mass_o*(rt/rl[bb])**2
            plt.plot(rt,mass_check,'r--')
            plt.plot([0.35,0.35],[10**6,10**8],'b--')
            plt.plot([1,1],[10**6,10**8],'b--')
        
        #=============Velocity
        elif des_var ==2:
            plt.errorbar(rt[0:int(bins)],vel2[0:int(bins)], yerr =dvel[0:int(bins)]/vel_all,marker='.',label = filis[bb][0:-5],linestyle = ' ',color='k')
        
        #============Density
        elif des_var ==3:
            plt.errorbar(rt[3:int(bins)]/rl[bb],dt[3:int(bins)], yerr =ddt[3:int(bins)],label = filis[bb][0:-5],linestyle = ' ',marker='.',color='k')
    
    #======================Label Axes===============================================
    if des_var<3:
        plt.xlabel('Distance (kpc)')
    else:
        plt.xlabel('r/$r_e$')
    
    if des_var ==1:
        plt.ylabel('Mass ($M_{sol}$)')
    elif des_var ==2:
        plt.ylabel('$V_c$ (km/s)')
    elif des_var ==3:
        plt.ylabel(r'$\Sigma/<\Sigma_{1/2}>$')
    plt.legend()
    
    print(filis[bb][0:-5])
    print("Vc_{c1/2}: "+str(vel_all)+ " +- "+ str(dvel_all) +" km/s")
    print("M_{1/2}: "+str(mass_o/1e6) + "x 10^6 solar M")
    
    #===================Half Light Mass and Velocity =========================================
    if des_var == 1:
        plt.plot(rl[bb],mass_o,'.',label="r_e,M_e",color='k',marker='p')
    elif des_var == 2:
        plt.plot(rl[bb],vel_all,'.',label="r_e,V_c",color='k',marker='p')
    
    #==========================Set log scales=============================================
    plt.xscale('log')
    if log_y:
        plt.yscale('log')
    plt.legend()
    plt.show()
#plt.show()
#========================Make plots of Re vs Vc ====================================================
#==========================Include PiscesII and PegasusIII=====================================
if include_pis:
    names_final =  ['Sextans', 'Sculptor', 'Fornax', 'Carina','Leo_I','Draco','Leo_II','Eridanas_II','AntliaII','CraterII','PiscesII','PegasusIII','TriangulumII','SegueI']
    v_circ = np.array([12.3,14.5,20.9,12.7,16.3,18.0,13.3,19,11.1,9.1,9,9,9,7])
    dv_circ = np.array([0.2,0.1,0.1,0.1,0.3,0.5,0.4,1,0.5,0.3,5,5,4,2])
    rad_ef = np.array([768,282,714,254,251,221,176,280,2800,1066,58,53,34,30])
#=========================Don't Include PiscesII and PegasusIII===============================
else:
    names_final =  ['Sextans', 'Sculptor', 'Fornax', 'Carina','Leo_I','Draco','Leo_II','Eridanas_II','AntliaII','CraterII']
    v_circ = np.array([12.3,14.5,20.9,12.7,16.3,18.0,13.3,19,11.1,9.1])
    dv_circ = np.array([0.2,0.1,0.1,0.1,0.3,0.5,0.4,1,0.5,0.3])
    rad_ef = np.array([768,282,714,254,251,221,176,280,2800,1066])


drad_ef = np.array([44,45,77,39,27,19,42,10,100])
markers = ['d','_','1','.','8','p','<','>','^','v','*','s','$a$','$b$']

lin_x = np.array([0.07,0.3])
lin_y = 200*lin_x**(3/2)
lin_x2 = np.array([0.03,0.3])
lin_y2 = 400*lin_x2**(1.45)

lin_x3 = np.array([0.07,0.3])
lin_y3 = 100*lin_x3**(3/2)
lin_x4 = np.array([0.03,0.3])
lin_y4 = 250*lin_x4**(1.45)

for ii in range(0,len(names_final)):
    plt.errorbar(rad_ef[ii]/1000,v_circ[ii],yerr=dv_circ[ii],linestyle='',marker=markers[ii],label=names_final[ii],color='k')

plt.plot(lin_x,lin_y,'k--')
plt.plot(lin_x2,lin_y2,'b--')
plt.plot(lin_x3,lin_y3,'k--')
plt.plot(lin_x4,lin_y4,'b--')

plt.xlabel('$R_e$ (kpc)')
plt.ylabel('$V_{ce}$ (km/s)')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.show()