import urllib
import os
import glob
import subprocess
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.time import Time
import csv
import json
from pprint import pprint
import os.path
import sys
from scipy.interpolate import UnivariateSpline
from valenti_ni56 import *
from ni56_mass_print import *
#I = r - 1.2444*(r - i) - 0.3820 lupton 2005
sys.path.insert(0, '/home/afsari/')
from SNAP5.Analysis import *

with open("Data/sn_names_tail1.txt") as f:
    file_names = f.readlines()

file_names = [line.rstrip('\n') for line in file_names]
files_count=len(file_names)
#s&d 2011
coef = {'B': 3.626, 'V': 2.742, 'I': 1.505, "i'": 1.698, "r'":2.285, "R":2.169}
add=[]
for i,sn_name in enumerate(['SN2011dh']):
    print(sn_name)
    print "name",sn_name
    mag = np.zeros(shape=(0, 6))
    location='./Data/SN_json/'
    print os.path.isfile(location+sn_name+'.json')
    if os.path.isfile(location+sn_name+'.json')==False:
        url='https://sne.space/astrocats/astrocats/supernovae/output/json/'+sn_name+'.json'
        urllib.urlretrieve(url,location+sn_name+'.json')
    with open(location+sn_name+'.json') as data_file:
        print data_file
        data = json.load(data_file)
    if sn_name=='SN2011dh':
        ebv=0.07
        d=7.8e6
        de=1e6
        t0=55712.5
        DM=5*np.log10(d)-5
        DMe=np.abs(5*de/(d*np.log(10)))
        Band=["R","V"]

    for dat in data[sn_name]["photometry"]:
        if  np.any(np.in1d(dat["band"],Band)):
            if "e_magnitude" in dat:
                error=float(dat["e_magnitude"])
            else:
                error=0
            if '17' in dat["source"].split(','):
                add = np.concatenate(([float(dat["time"]), dat["magnitude"], error],[deredMag(float(dat["magnitude"]), float(ebv), coef[dat["band"]])-DM, np.sqrt(error**2+DMe**2),dat["band"]]))
                add = np.reshape(add, (1, 6))
                mag = np.concatenate((mag, add), axis=0)
    np.save("Data/SN_photometry/"+sn_name + ".npy", mag)

    plt.scatter(mag[:,0][mag[:,5]=='V'].astype(float)-t0,mag[:,3][mag[:,5]=='V'].astype(float))
    plt.scatter(mag[:,0][mag[:,5]=="R"].astype(float)-t0,mag[:,3][mag[:,5]=="R"].astype(float))
    #plt.scatter(mag[:,0][mag[:,5]=="B"].astype(float)-t0,mag[:,3][mag[:,5]=="B"].astype(float))

    t_u = np.arange(10, 130, 0.1)
    M_u=np.zeros(shape=(len(Band)+2,np.shape(t_u)[0]))
    for i,band in enumerate(Band):
        t=mag[:, 0][mag[:, 5] == band].astype(float)
        M=mag[:,3][mag[:,5]==band].astype(float)
        M = M[np.argsort(t)]
        t = t[np.argsort(t)]-t0
        fit = UnivariateSpline(t,M, s=0.2)
        M_u[i,:]=fit(t_u)
        plt.scatter(t_u,M_u[i,:])

    # M_u[3, :]=M_u[1, :]- 1.2444*(M_u[1, :] - M_u[0, :]) - 0.3820 #I band filter
    # M_u[4, :] = M_u[1, :] - 0.2936*(M_u[1, :] - M_u[0, :]) - 0.1439 #R band filter
    # BC_V_I=0.213-0.203*(M_u[2,:]-M_u[3,:])-0.079*(M_u[2,:]-M_u[3,:])**2
    BC_V_R = 0.197 - 0.183 * (M_u[1, :] - M_u[0, :]) - 0.419 * ((M_u[1, :] - M_u[0, :])**2)
    Mbol=BC_V_R+M_u[1,:]-0
    #plt.scatter(t_u, Mbol, label='Mbol')
    plt.legend()
    plt.gca().invert_yaxis()
    plt.show()
    Msun = 4.74
    Lsun = 3.84e33
    lbol = Lsun * np.power(10, ((Msun - Mbol) / 2.5))
    print "lbol_max:", np.max(lbol), t_u[np.argmax(lbol)]
    print "M_ni(beta)=", nickel_mass(t_u[np.argmax(lbol)], np.max(lbol), 9/8)
    beta = np.arange(0, 5, 0.005)
    m56 = np.zeros(shape=beta.shape)
    for i, b in enumerate(beta):
        m56[i] = nickel_mass(t_u[np.argmax(lbol)], np.max(lbol), b).value
    E=np.arange(2.7, 3.5, 0.1)
    Mni56=np.zeros(shape=E.shape)
    Mni56_std=np.zeros(shape=E.shape)
    for i,e in enumerate(E):
        M_ni56=valenti_ni56(t_u[t_u >= 60],lbol[t_u >=60],1,e)
        Mni56_std[i]=np.std(M_ni56)
        Mni56[i]=np.mean(M_ni56)
    print Mni56[np.argmin(Mni56_std)],np.min(Mni56_std),E[np.argmin(Mni56_std)]
    for i, b in enumerate(beta):
        m56[i] = nickel_mass(t_u[np.argmax(lbol)], np.max(lbol), b).value
    print "m56, beta", m56[np.argmin(np.abs(m56-Mni56[np.argmin(Mni56_std)]))],beta[np.argmin(np.abs(m56-Mni56[np.argmin(Mni56_std)]))]
    alpha_co = -1 / 111.26
    beta_co = 43
    co = np.zeros(shape=(2,))
    co[1] = beta_co
    co[0] = alpha_co
    fit_3 = np.poly1d(co)

    alpha_ni = -1 / 6.1
    beta_ni = 44
    ni = np.zeros(shape=(2,))
    ni[1] = beta_ni
    ni[0] = alpha_ni
    fit_4 = np.poly1d(ni)
    #print M_u[2,:]-M_u[4,:]
    plt.plot(np.linspace(100, 140), fit_3(np.linspace(100, 140)), '--', color='blue', label=r'$^{56}$Co', lw=2)
    plt.plot(np.linspace(10, 15), fit_4(np.linspace(10, 15)), '-.', color='green', label=r'$^{56}$Ni', lw=2)
    plt.xlabel('Time [days]')
    plt.ylabel(r'L_bol [erg/s]')
    ax2=plt.subplot(111)
    plt.legend()
    plt.scatter(t_u, np.log10(lbol))
    #plt.gca().invert_yaxis()
    plt.show()

