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

from scipy.optimize import curve_fit
from scipy.optimize import leastsq

#essential files
from SNAP4.Analysis.LCRoutines import*
from SNAP4.Analysis.LCFitting import*
from scipy.optimize import leastsq


fig = plt.figure()
url='https://docs.google.com/spreadsheets/d/e/2PACX-1vQsvlrp1xYmZcWSy11LxxPg_X5qNPxlzjRUyQ_iGGSXcYESRxub6loaMWdD4J9jtqc8B4dKoKW_Fs0i/pub?gid=0&single=true&output=csv'
plt.close("all")
with requests.Session() as s:
    download = s.get(url)
    decoded_content = download.content.decode('utf-8')
    cr = csv.reader(decoded_content.splitlines(), delimiter=',')
    my_list = list(cr)

url='https://docs.google.com/spreadsheets/d/e/2PACX-1vTRiEGa3J6bnG5R4JhDz6Wx2gb_KK9TnRdkn5YdWIoar6_ugH1p42KYLkDMHiMjjOUxXPxaPH22wjRN/pub?gid=0&single=true&output=csv'
plt.close("all")
with requests.Session() as s:
    download = s.get(url)
    decoded_content = download.content.decode('utf-8')
    cr = csv.reader(decoded_content.splitlines(), delimiter=',')
    my_list1 = list(cr)

info=np.array(my_list)
prentice=np.array(my_list1)
prentice_name=[]
for j,sn in enumerate(prentice):
    prentice_name.append(sn[0])
coef = {'B': 3.626, 'V': 2.742, 'I': 1.505, "i'": 1.698, "r'": 2.285, "R": 2.169}
add = []
sug_beta={'Ib':1.125, 'Ic':1.125, 'IIb':0.82, 'Ia':1.6}
f = plt.figure(1)
f.set_figheight(10)
f.set_figwidth(12)
f = plt.figure(6)
f.set_figheight(10)
f.set_figwidth(12)
f = plt.figure(2)
f.set_figheight(10)
f.set_figwidth(12)
results = np.zeros(shape=(info[1:,0].shape[0]+1, 14), dtype=object)
results[0,:]=['Name', 'type', 'Lbol_peak', 't_peak', '56Ni from Suggested beta', '56Ni tail', '56Ni tail error','F','M_ejEk', 'beta required', 'Lp_prentice', 'tp_prentice', 'Ni_prentice','Arnett Ni']

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
    if sn_name=='SN2008D':
        Band = ["V","i'", "r'"]
    print info[i+1,11].split(',')
    print ebv,ebv_e,d,de,t0,t0_e,Band
    mag = np.zeros(shape=(0, 6))
    for dat in data[sn_name]["photometry"]:
       if "band" in dat.keys():
            if  np.any(np.in1d(dat["band"],Band)):
                if "e_magnitude" in dat:
                    error=float(dat["e_magnitude"])
                else:
                    error=0
                if np.any(np.in1d(dat["source"].split(','),info[i+1,11].split(','))):
                    add = np.concatenate(([float(dat["time"]), dat["magnitude"], error],[deredMag(float(dat["magnitude"]), float(ebv), coef[dat["band"]])-DM, np.sqrt(error**2+DMe**2),dat["band"]]))
                    add = np.reshape(add, (1, 6))
                    mag = np.concatenate((mag, add), axis=0)
    print info[i+1,11]
    np.save("Data/SN_photometry/"+sn_name + ".npy", mag)

    if sn_name == 'SN2016gkg':
        t_u = np.arange(10, 100, 0.1)
        s=0.15
    elif sn_name == 'iPTF13bvn':
        t_u = np.arange(10, 90, 0.1)
        s=0.2
    elif sn_name == 'SN2008D':
        t_u = np.arange(10, 120, 0.1)
        s=0.25
    elif sn_name == 'SN1993J':
        t_u = np.arange(10, 100, 0.1)
        s=0.25
    elif sn_name == 'SN2013df':
        t_u = np.arange(10, 200, 0.1)
        s=0.35
    elif sn_name == 'SN2006jc':
        t_u = np.arange(2, 65, 0.1)
        s=0.2
    elif sn_name == 'SN1994I':
        t_u = np.arange(5, 100, 0.1)
        s=0.25
    elif sn_name == 'SN2004aw':
        t_u = np.arange(5, 85, 0.1)
        s=0.25
    elif sn_name == 'SN2009bb':
        t_u = np.arange(3, 45, 0.1)
        s=0.35
    elif sn_name == 'SN2007ru':
        t_u = np.arange(1, 130, 0.1)
        s=0.1
    elif sn_name == 'SN2002ap':
        t_u = np.arange(1, 40, 0.1)
        s=0.2
    else:
        t_u = np.arange(5, 130, 0.1)
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
    t_u_t=np.clip(t_u,a_min=60, a_max=130)
    if sn_name == 'SN2009bb':
        t_u_t = np.arange(60, 75, 0.1)
    if sn_name == 'SN2002ap':
        t_u_t = np.arange(50, 230, 0.1)
    if sn_name == 'SN2006jc':
        t_u_t = np.arange(50, 100, 0.1)
    t_u_t = np.unique(t_u_t)
    M_u_t=np.zeros(shape=(len(Band)+2,np.shape(t_u_t)[0]))
    M_e_t=np.zeros(shape=(len(Band)+2,np.shape(t_u_t)[0]))
    # if sn_name=='SN2002ap':
    #     continue
    for j,band in enumerate(Band):
        t=mag[:, 0][mag[:, 5] == band].astype(float)
        M=mag[:,3][mag[:,5]==band].astype(float)
        Me = mag[:, 4][mag[:, 5] == band].astype(float)
        M = M[np.argsort(t)]
        Me = Me[np.argsort(t)]
        t = t[np.argsort(t)] - t0
        t,u=np.unique(t,return_index=True)
        M=M[u]
        Me=Me[u]
        fit = ErrorPropagationLinear(t[(t>50)],M[(t>50) ],Me[(t>50)])
        M_u_t[j,:],M_e_t[j,:]=fit(t_u_t)
    if 1:
        plt.figure(1)
        plt.subplot(3,6, i+1)
        plt.scatter(mag[:,0][mag[:,5]==Band[0]].astype(float)-t0,mag[:,3][mag[:,5]==Band[0]].astype(float),label=Band[0])
        plt.scatter(mag[:,0][mag[:,5]==Band[1]].astype(float)-t0,mag[:,3][mag[:,5]==Band[1]].astype(float),label=Band[1])
        plt.scatter(t_u_t, M_u_t[0, :],s=0.5)
        plt.scatter(t_u_t, M_u_t[1, :],s=0.5)
        plt.legend()
        plt.title(sn_name)
        plt.gca().invert_yaxis()
        plt.figure(6)
        plt.subplot(3,6, i+1)
        plt.scatter(mag[:,0][mag[:,5]==Band[0]].astype(float)-t0,mag[:,3][mag[:,5]==Band[0]].astype(float),label=Band[0])
        plt.scatter(mag[:,0][mag[:,5]==Band[1]].astype(float)-t0,mag[:,3][mag[:,5]==Band[1]].astype(float),label=Band[1])
        plt.scatter(t_u, M_u[0, :],s=0.5)
        plt.scatter(t_u, M_u[1, :],s=0.5)
        plt.legend()
        plt.title(sn_name)
        plt.gca().invert_yaxis()
        if sn_name=='SN2008D':
            M_u[1, :] = M_u[2, :] - 1.2444 * (M_u[2, :] - M_u[1, :]) - 0.3820
            Band[1]='I'
        lbol,le=lyman_BC(np.reshape(M_u[0, :],(1,np.shape(M_u[0, :])[0])), np.reshape(M_u[1, :],(1,np.shape(M_u[1, :])[0])),Band[0],Band[1])
        lbol_t,le_t=lyman_BC(np.reshape(M_u_t[0, :],(1,np.shape(M_u_t[0, :])[0])), np.reshape(M_u_t[1, :],(1,np.shape(M_u_t[1, :])[0])),Band[0],Band[1])
        plt.figure(2)
        plt.subplot(3, 6, i+1)
        plt.gca().annotate(sn_name,xy=(0.7, 0.9), xycoords='axes fraction')
        plt.errorbar(np.reshape(t_u, (1,t_u.shape[0])),lbol,yerr=le, fmt='b', ecolor='yellow')
        xx = 1
        if (sn_name=='SN2013ge') | (sn_name=='iPTF13bvn') | (sn_name=='SN2009jf'):
            #plt.figure(3)
            #plt.errorbar(np.reshape(t_u-t_u[np.argmax(lbol)], (1, t_u.shape[0])), lbol, yerr=le, label='Lyman16 BC')
            if sn_name == 'SN2013ge':
                plt.figure(3)
                plt.errorbar(np.reshape(t_u - t_u[np.argmax(lbol)], (1, t_u.shape[0])), lbol, yerr=le,
                             label='Lyman16 BC')
                data_ge = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/' + sn_name + '_PseudoBol.dat',
                                     delimiter=',')
                plt.errorbar(data_ge[:,1], data_ge[:,3], yerr=data_ge[:,6], label=sn_name+' drout16')
                lbol = data_ge[:, 3]
                t_u = data_ge[:, 0] - t0
                plt.xlabel('Time')
                plt.ylabel('Luminosity [erg/s]')
                plt.legend()
                plt.title(sn_name)
            elif sn_name=='iPTF13bvn' :
                plt.figure(4)
                plt.errorbar(np.reshape(t_u, (1, t_u.shape[0])), lbol, yerr=le, label='Lyman16 BC')
                data_ge = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/' + sn_name + '_PseudoBol.dat',
                        delimiter=',')
                plt.scatter(data_ge[:, 0], 10**data_ge[:, 1],s=2, label=sn_name+' drout16')
                lbol = 10**data_ge[:, 1]
                t_u = data_ge[:, 0]
                plt.xlabel('Time')
                plt.ylabel('Luminosity [erg/s]')
                plt.legend()
                plt.title(sn_name)
            elif sn_name=='SN2009jf' :
                plt.figure(5)
                plt.errorbar(np.reshape(t_u, (1, t_u.shape[0])), lbol, yerr=le, label='Lyman16 BC')
                data_ge = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/' + sn_name + '_PseudoBol.dat',
                        delimiter=',')
                plt.scatter(data_ge[:, 0], 10**((43-data_ge[:, 1])+41),s=2, label=sn_name+' drout16')
                lbol = 10**((43-data_ge[:, 1])+41)
                t_u = data_ge[:, 0]
                plt.xlabel('Time')
                plt.ylabel('Luminosity [erg/s]')
                plt.legend()
                plt.title(sn_name)
            xx=0
            # fig.subplots_adjust(hspace=0, wspace=0)
            # plt.tight_layout()
    results[i + 1, 3]=t_u[np.argmax(lbol)]
    results[i + 1,2]=np.log10(np.max(lbol))
    print "lbol_max:", np.max(lbol), t_u[np.argmax(lbol)]
    print i,le.shape, lbol.shape
    results[i + 1, 4] =nickel_mass_khatami(t_u[np.argmax(lbol)], np.max(lbol),le[0,np.argmax(lbol)],sug_beta[info[i+1,2].split(' ')[0]])[0]
    print "M_ni(beta)=", results[i+1,2]
    lbol=np.reshape(lbol,t_u.shape)
    lbol_t=np.reshape(lbol_t,t_u_t.shape)
    if xx:
        le_t = np.reshape(le_t, t_u_t.shape)
        le = np.reshape(le, t_u.shape)
    else:
        le_t = np.zeros(t_u_t.shape)
        le = np.zeros(t_u.shape)

    x=[0.1, 3]

    p, p_err = fit_bootstrap(x, t_u_t[(t_u_t >= 60) & (t_u_t < 120)], lbol_t[(t_u_t >= 60)& (t_u_t <120) ],np.clip(le_t[(t_u_t >= 60)& (t_u_t <120) ],a_min=0.0001e42,a_max=1e51), wheelerErr, errfunc=True, perturb=True, n=1000, nproc=4)

    F =32*p[1]
    print p[0], p_err[0], p[1], "F:",F, 32*p_err[1]
    results[i + 1, 5] =p[0]
    results[i + 1,6] =p_err[0]
    results[i + 1, 7] = p[1]*32
    results[i + 1, 8] =p[1]
    plt.figure(2)
    plt.subplot(3, 6, i + 1)
    plt.plot(t_u_t[(t_u_t >= 60) & (t_u_t < 120)],wheeler_bol(t_u_t[(t_u_t >= 60) & (t_u_t < 120)],p[0],p[1]),color="orange")
    beta0=0.3
    betafit = leastsq(nickel_mass_khatami_err, beta0, args=(t_u[np.argmax(lbol)], np.max(lbol), le[np.argmax(lbol)],p[0],p_err[0]), full_output=0, maxfev=10000000)
    results[i + 1, 9] = betafit[0][0]
    try:
        ind = prentice_name.index(sn_name.strip('SN'))
        results[i+1,10]=prentice[ind,2]
        results[i + 1,11] = prentice[ind, 3]
        results[i + 1,12] = prentice[ind, 4]
    except:
        print 'no data'
    results[i + 1, 13]= valenti_ni56(results[i + 1, 3],10**results[i + 1,2])[0]
    results[i + 1, 0] =sn_name
    results[i+1,1]=info[i+1,2]
    # if i==17:
    #     break

np.savetxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/Results_linear.txt',results,fmt='%12s',delimiter=',')
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

plt.figure(3)
plt.show()

plt.figure(4)
plt.show()

plt.figure(5)
plt.show()