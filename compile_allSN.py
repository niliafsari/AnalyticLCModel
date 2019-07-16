import csv
import requests
import numpy as np
import urllib
import os
import glob
import subprocess
import commands
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
from error_interpolation import *
from astropy.constants import M_sun
#from ni56_mass_print import *
#I = r - 1.2444*(r - i) - 0.3820 lupton 2005
sys.path.insert(0, '/home/afsari/')
from SNAP5.Analysis import *
from lyman_BC import  *

def nickel_mass(t_peak,L_peak,beta):
    return (L_peak*(beta**2)*(t_peak/8.8)**2)/(2*3.9e10*((0.83*(1-beta*t_peak/8.8)*np.exp(-beta*t_peak/8.8))+(26.56*(1-((1+beta*t_peak/111.3)*np.exp(-beta*t_peak/111.3))))))/M_sun.to("g")


fig = plt.figure()
url='https://docs.google.com/spreadsheets/d/e/2PACX-1vQsvlrp1xYmZcWSy11LxxPg_X5qNPxlzjRUyQ_iGGSXcYESRxub6loaMWdD4J9jtqc8B4dKoKW_Fs0i/pub?gid=0&single=true&output=csv'
plt.close("all")
with requests.Session() as s:
    download = s.get(url)
    decoded_content = download.content.decode('utf-8')
    cr = csv.reader(decoded_content.splitlines(), delimiter=',')
    my_list = list(cr)
sn_fit_max={'SN2016gkg':1}
sn_s={'SN2016gkg':0.05}
info=np.array(my_list)
coef = {'B': 3.626, 'V': 2.742, 'I': 1.505, "i'": 1.698, "r'": 2.285, "R": 2.169}
add = []
sug_beta={'Ib':1.125, 'Ic':1.125, 'IIb':0.82, 'Ia':1.6}
f = plt.figure(1)
f.set_figheight(7)
f.set_figwidth(12)
f = plt.figure(2)
f.set_figheight(7)
f.set_figwidth(12)
for i,sn_name in enumerate(info[1:,0]):
    print "SN name:",sn_name
    location='./Data/SN_json/'
    if os.path.isfile(location+sn_name+'.json')==False:
        url='https://sne.space/astrocats/astrocats/supernovae/output/json/'+sn_name+'.json'
        urllib.urlretrieve(url,location+sn_name+'.json')
    with open(location+sn_name+'.json') as data_file:
        data = json.load(data_file)
    ebv=0
    ebv_e=0
    d=0
    de=0
    t0=0
    t0_e=0
    ebv = float(info[i+1,3].split('(')[0])
    if (len(info[i+1,3].split('('))>1):
        ebv_e=float(info[i+1,3].split('(')[1].replace(")", ""))
    d = float(info[i+1,5].split('(')[0])*1.0e6
    if (len(info[i+1,5].split('('))>1):
        de=float(info[i+1,5].split('(')[1].replace(")", ""))
    t0 = float(info[i+1,7].split('(')[0])
    if (len(info[i+1,7].split('('))>1):
        t0_e=float(info[i+1,7].split('(')[1].replace(")", ""))
    print t0
    t0 = t0-2400000.5
    print t0
    if (d>0):
        DM = 5 * np.log10(d) - 5
    else:
        DM=0
    if (de>0):
        DMe = np.abs(5 * de / (d * np.log(10)))
    else:
        DMe=0
    Band=info[i+1,10].split(',')[0].split("-")
    print ebv,ebv_e,d,de,t0,t0_e,Band
    if sn_name=='SN2008D':
        continue
    mag = np.zeros(shape=(0, 6))
    for dat in data[sn_name]["photometry"]:
       if "band" in dat.keys():
            if  np.any(np.in1d(dat["band"],Band)):
                if "e_magnitude" in dat:
                    error=float(dat["e_magnitude"])
                else:
                    error=0
                if info[i+1,11] in dat["source"].split(','):
                    add = np.concatenate(([float(dat["time"]), dat["magnitude"], error],[deredMag(float(dat["magnitude"]), float(ebv), coef[dat["band"]])-DM, np.sqrt(error**2+DMe**2),dat["band"]]))
                    add = np.reshape(add, (1, 6))
                    mag = np.concatenate((mag, add), axis=0)
    np.save("Data/SN_photometry/"+sn_name + ".npy", mag)

    if sn_name == 'SN2016gkg':
        t_u = np.arange(10, 100, 0.1)
        s=0.05
    elif sn_name == 'iPTF13bvn':
        t_u = np.arange(10, 230, 0.1)
        s=0.05
    elif sn_name == 'SN2009jf':
        t_u = np.arange(10, 130, 0.1)
        s=0.05
    elif sn_name == 'SN2013df':
        t_u = np.arange(10, 250, 0.1)
        s=0.25
    else:
        t_u = np.arange(10, 130, 0.1)
        s=0.3
    M_u=np.zeros(shape=(len(Band)+2,np.shape(t_u)[0]))
    M_e=np.zeros(shape=(len(Band)+2,np.shape(t_u)[0]))
    for j,band in enumerate(Band):
        t=mag[:, 0][mag[:, 5] == band].astype(float)
        M=mag[:,3][mag[:,5]==band].astype(float)
        Me = mag[:, 4][mag[:, 5] == band].astype(float)
        M = M[np.argsort(t)]
        Me = Me[np.argsort(t)]
        t = t[np.argsort(t)] - t0
        t,u=np.unique(t,return_index=True)
        fit = ErrorPropagationSpline(t,M[u],Me[u], s=s)
        M_u[j,:],M_e[j,:]=fit(t_u)
        if (sn_name=='iPTF13bvn' )| (sn_name=='SN2016gkg') | (sn_name=='SN2013df') :
            fit = UnivariateSpline(t, M[u], s=s)
            M_u[j, :]= fit(t_u)
    if 1:
        plt.figure(1)
        plt.subplot(2, 4, i)
        plt.scatter(mag[:,0][mag[:,5]==Band[0]].astype(float)-t0,mag[:,3][mag[:,5]==Band[0]].astype(float),label=Band[0])
        plt.scatter(mag[:,0][mag[:,5]==Band[1]].astype(float)-t0,mag[:,3][mag[:,5]==Band[1]].astype(float),label=Band[1])
        plt.scatter(t_u, M_u[0, :],s=0.5)
        plt.scatter(t_u, M_u[1, :],s=0.5)
        plt.legend()
        plt.title(sn_name)
        plt.gca().invert_yaxis()
        lbol,le=lyman_BC(np.reshape(M_u[0, :],(1,np.shape(M_u[0, :])[0])), np.reshape(M_u[1, :],(1,np.shape(M_u[1, :])[0])),Band[0],Band[1])
        plt.figure(2)
        plt.subplot(2, 4, i)
        plt.errorbar(np.reshape(t_u, (1,t_u.shape[0])),lbol,yerr=le)
        if sn_name=='SN2013ge':
            plt.figure(3)
            plt.errorbar(np.reshape(t_u-t_u[np.argmax(lbol)], (1, t_u.shape[0])), lbol, yerr=le)
            data_ge=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SN2013ge_PseudoBol.dat', delimiter=',')
            plt.errorbar(data_ge[:,1], data_ge[:,3], yerr=data_ge[:,6])
            lbol=data_ge[:,3]
            t_u=data_ge[:,0]-t0
            fig.subplots_adjust(hspace=0, wspace=0)
            plt.tight_layout()
            plt.show()
        plt.legend()
        plt.title(sn_name)

    print "lbol_max:", np.max(lbol), t_u[np.argmax(lbol)]
    print "M_ni(beta)=", nickel_mass(t_u[np.argmax(lbol)], np.max(lbol),sug_beta[info[i+1,2]])
    beta = np.arange(0, 5, 0.005)
    m56 = np.zeros(shape=beta.shape)
    for uu, b in enumerate(beta):
        m56[uu] = nickel_mass(t_u[np.argmax(lbol)], np.max(lbol), b).value
    E = np.arange(2.7, 5, 0.1)
    Mni56 = np.zeros(shape=E.shape)
    Mni56_std = np.zeros(shape=E.shape)
    lbol=np.reshape(lbol,t_u.shape)
    for i, e in enumerate(E):
        M_ni56 = valenti_ni56(t_u[(t_u >= 60) & (t_u < 120)], lbol[(t_u >= 60)& (t_u <120) ], 1, e)
        Mni56_std[i] = np.std(M_ni56)
        Mni56[i] = np.mean(M_ni56)
    print Mni56[np.argmin(Mni56_std)], np.min(Mni56_std), E[np.argmin(Mni56_std)]
    for i, b in enumerate(beta):
        m56[i] = nickel_mass(t_u[np.argmax(lbol)], np.max(lbol), b).value
    print "m56, beta", m56[np.argmin(np.abs(m56 - Mni56[np.argmin(Mni56_std)]))], beta[
        np.argmin(np.abs(m56 - Mni56[np.argmin(Mni56_std)]))]
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

#plt.xlabel('Time (days)')
#plt.ylabel('Magnitude (mag)')
plt.figure(1)
fig.subplots_adjust(hspace=0, wspace=0)
plt.tight_layout()
plt.show()

plt.figure(2)
fig.subplots_adjust(hspace=0, wspace=0)
plt.tight_layout()
plt.show()