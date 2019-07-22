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

def nickel_mass(t_peak,L_peak,L_peak_e,beta):
    n=1000
    L= np.multiply(np.tile(L_peak_e,n), np.random.randn(n)+np.tile(L_peak,n))
    t_peak_t=np.tile(t_peak, n)
    MNi=np.divide(np.multiply(L_peak*(beta**2),(t_peak_t/8.8)**2),(2*3.9e10*((0.83*np.multiply((1-beta*t_peak_t/8.8),np.exp(-beta*t_peak_t/8.8)))+(26.56*(1-(np.multiply((1+beta*t_peak_t/111.3),np.exp(-beta*t_peak_t/111.3))))))))/M_sun.to("g").value
    return np.mean(MNi,axis=0), np.std(MNi,axis=0)

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
results = np.zeros(shape=(info[1:,0].shape[0]+1, 7), dtype=object)
results[0,:]=['Lbol_peak', 't_peak', '56Ni from Suggested beta', '56Ni tail','M_ej/E', '56Ni tail error', 'beta required']
print results[0,:]
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
        s=0.15
    elif sn_name == 'iPTF13bvn':
        t_u = np.arange(10, 90, 0.1)
        s=0.2
    elif sn_name == 'SN1993J':
        t_u = np.arange(10, 100, 0.1)
        s=0.25
    elif sn_name == 'SN2013df':
        t_u = np.arange(10, 200, 0.1)
        s=0.35
    else:
        t_u = np.arange(13, 130, 0.1)
        s=0.2
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
        if (sn_name=='SN2016gkg')  :
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
        plt.errorbar(np.reshape(t_u, (1,t_u.shape[0])),lbol,yerr=le, fmt='b', ecolor='yellow')
        # if sn_name=='SN1993J':
        #     plt.figure(3)
        #     plt.scatter(mag[:, 0][mag[:, 5] == Band[0]].astype(float) - t0,
        #                 mag[:, 3][mag[:, 5] == Band[0]].astype(float), label=Band[0])
        #     plt.scatter(mag[:, 0][mag[:, 5] == Band[1]].astype(float) - t0,
        #                 mag[:, 3][mag[:, 5] == Band[1]].astype(float), label=Band[1])
        #     plt.errorbar(np.reshape(t_u-t_u[np.argmax(lbol)], (1, t_u.shape[0])), lbol, yerr=le, label='Lyman16 BC')
        #     #data_ge=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SN2013ge_PseudoBol.dat', delimiter=',')
        #     #plt.errorbar(data_ge[:,1], data_ge[:,3], yerr=data_ge[:,6], label='drout16')
        #     #lbol=data_ge[:,3]
        #     #t_u=data_ge[:,0]-t0
        #     fig.subplots_adjust(hspace=0, wspace=0)
        #     plt.tight_layout()
        #     plt.xlabel('Time')
        #     plt.ylabel('Luminosity [erg/s]')
        #     plt.legend()
        #     plt.show()
        plt.legend()
        plt.title(sn_name)
    results[i + 1, 1]=t_u[np.argmax(lbol)]
    results[i + 1, 0]=np.max(lbol)
    print "lbol_max:", np.max(lbol), t_u[np.argmax(lbol)]
    print i,le.shape, lbol.shape
    results[i + 1, 2] =nickel_mass(t_u[np.argmax(lbol)], np.max(lbol),le[0,np.argmax(lbol)],sug_beta[info[i+1,2]])[0]
    print "M_ni(beta)=", results[i+1,2]
    beta = np.arange(0, 5, 0.005)
    m56 = np.zeros(shape=beta.shape)
    m56_e = np.zeros(shape=beta.shape)
    # for uu, b in enumerate(beta):
    #     m56[uu],m56_e[uu]  = nickel_mass(t_u[np.argmax(lbol)], np.max(lbol),le[np.argmax(lbol)], b)
    E = np.arange(2.7, 5, 0.1)
    Mni56 = np.zeros(shape=E.shape)
    Mni56_e = np.zeros(shape=E.shape)
    Mni56_std = np.zeros(shape=E.shape)
    lbol=np.reshape(lbol,t_u.shape)
    le = np.reshape(le, t_u.shape)
    for k, e in enumerate(E):
        M_ni56, Mni56_e = valenti_ni56(t_u[(t_u >= 60) & (t_u < 120)], lbol[(t_u >= 60)& (t_u <120) ], 1, e,le[(t_u >= 60)& (t_u <120) ])
        Mni56_std[k] = np.sqrt(np.std(M_ni56)**2+np.mean(Mni56_e)**2)
        Mni56[k] = np.mean(M_ni56)
    print Mni56[np.argmin(Mni56_std)], np.min(Mni56_std), E[np.argmin(Mni56_std)]
    results[i + 1, 3] =Mni56[np.argmin(Mni56_std)]
    results[i + 1,4] =np.min(Mni56_std)
    results[i + 1, 5] =E[np.argmin(Mni56_std)]
    for k, b in enumerate(beta):
        m56[k], m56_e[k]= nickel_mass(t_u[np.argmax(lbol)], np.max(lbol),le[np.argmax(lbol)], b)
    print "m56, beta", m56[np.argmin(np.abs(m56 - Mni56[np.argmin(Mni56_std)]))], beta[
        np.argmin(np.abs(m56 - Mni56[np.argmin(Mni56_std)]))]
    results[i + 1, 6] = beta[np.argmin(np.abs(m56 - Mni56[np.argmin(Mni56_std)]))]

np.savetxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/Results.txt',results,fmt='%12s',delimiter=',')
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