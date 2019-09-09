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

SN_verified_t0=['SN2008D','']

coef = {'B': 3.626, 'V': 2.742, 'I': 1.505, "i'": 1.698, "r'": 2.285, "R": 2.169}

class supernova:
    def __init__(self, name,  host=None,sn_type=None,distance=[None, None], band_extinction=None, band_bc=None , source=None, ebv_gal=[None, None],t0=[None, None] ,ebv_host=[None, None],source_bc=None,smoothing_bc=None):
        self.name=name
        self.sn_type=sn_type
        self.host=host
        self.distance=distance
        self.band_extinction=band_extinction
        self.band_bc=band_bc
        self.source=source
        self.t0=t0
        self.ebv_gal=ebv_gal
        self.ebv_host=ebv_host
        self.verbose=0
        self.source_bc=source_bc
        self.smoothing_bc=smoothing_bc
        self.tail_ni_mass=[None, None]
        self.tail_meje=[None, None]
        self.peakL=None
        self.peakt=None
        self.arnett_ni_mass=[None, None]
        self.peakL=[None, None]
    def fill(self,info):
        location = './Data/SN_json/'
        if os.path.isfile(location + self.name + '.json') == False:
            url = 'https://sne.space/astrocats/astrocats/supernovae/output/json/' + self.name  + '.json'
            urllib.urlretrieve(url, location + self.name  + '.json')
        with open(location + self.name  + '.json') as data_file:
            data = json.load(data_file)
        index=np.where(info[:,0]==self.name)[0][0]
        self.sn_type = info[index , np.where(info[0,:]=='Type')[0][0]].split(' ')[0]
        self.host = info[index, np.where(info[0,:]=='Host')[0][0]]
        self.distance[0] = float(info[index , np.where(info[0,:]=='luminosity distance')[0][0]].split('(')[0])*1e6
        if (len(info[index , np.where(info[0,:]=='luminosity distance')[0][0]].split('(')) > 1):
            self.distance[1] = float(info[index, np.where(info[0,:]=='luminosity distance')[0][0]].split('(')[1].replace(")", ""))*1e6
        self.band_extinction = info[index, np.where(info[0,:]=='extinction color')[0][0]].split(" ")[0].split("-")
        self.band_bc = info[index ,  np.where(info[0,:]=='Color')[0][0]].split(',')[0].split("-")
        self.source = info[index, np.where(info[0,:]=='source')[0][0]].split(',')
        self.ebv_gal[0] = float(info[index , np.where(info[0,:]=='AV galactic(S&F2011)')[0][0]].split('(')[0]) / 3.1#Rv
        if (len(info[index , np.where(info[0,:]=='AV galactic(S&F2011)')[0][0]].split('(')) > 1):
            self.ebv_gal[1] = float(info[index , np.where(info[0,:]=='AV galactic(S&F2011)')[0][0]].split('(')[1].replace(")", ""))
        self.source_bc=info[index, np.where(info[0,:]=='source_number')[0][0]].split(',')
        self.smoothing_bc = float(info[index, np.where(info[0, :] == 'smoothing_bc')[0][0]])
        self.t0[0] = float(info[index, np.where(info[0, :] == 't0 (JD)')[0][0]].split('(')[0])-2400000.5
        if (len(info[index , np.where(info[0,:]=='t0 (JD)')[0][0]].split('(')) > 1):
            self.t0[1] = float(info[index , np.where(info[0,:]=='t0 (JD)')[0][0]].split('(')[1].replace(")", ""))
        if 1:
            print self.name,self.host,self.sn_type,self.distance,self.band_extinction,self.band_bc,self.source,self.ebv_gal,self.t0

    def stritzinger_color(self,s=0.2):
        add = []
        location='./Data/SN_json/'
        with open(location+self.name+'.json') as data_file:
            data = json.load(data_file)
        mag = np.zeros(shape=(0, 4))
        for dat in data[self.name]["photometry"]:
           if "band" in dat.keys():
                if  np.any(np.in1d(dat["band"],self.band_extinction)):
                    if "e_magnitude" in dat:
                        error=float(dat["e_magnitude"])
                    else:
                        error=0
                    if np.any(np.in1d(dat["source"].split(','),self.source)):
                        add = [float(dat["time"]), deredMag(float(dat["magnitude"]), float(self.ebv_gal[0]), coef[dat["band"]]), error, dat["band"]]
                        add = np.reshape(add, (1, 4))
                        mag = np.concatenate((mag, add), axis=0)
        t_u = np.arange(0, 10, 0.1)
        M_u = np.zeros(shape=(len(self.band_extinction) + 2, np.shape(t_u)[0]))
        M_e = np.zeros(shape=(len(self.band_extinction) + 2, np.shape(t_u)[0]))
        for j, band in enumerate(self.band_extinction):
            t = mag[:, 0][mag[:, 3] == band].astype(float)
            M = mag[:, 1][mag[:, 3] == band].astype(float)  # type: object
            Me = mag[:, 2][mag[:, 3] == band].astype(float)
            M = M[np.argsort(t)]
            Me = Me[np.argsort(t)]
            t = t[np.argsort(t)]
            t, u = np.unique(t, return_index=True)
            if j==0:
                t0=t[np.argmin(M[u])]
            if (self.name=='SN2016gkg'):
                t0=57668
            if (self.name=='SN2013df'):
                t0=56469
            if (self.name=='SN1993J'):
                t0=49093
            t=t-t0
            fit = ErrorPropagationSpline(t, M[u], Me[u], s=s)
            M_u[j, :], M_e[j, :] = fit(t_u)
            if (self.name == 'SN2016gkg'):
                fit = UnivariateSpline(t, M[u], s=0.1)
                M_u[j, :] = fit(t_u)
                M_e[j, :]=0
        folder='BmX_color_templates'
        scale=1
        if (self.band_extinction[0]=='V' ) & (self.band_extinction[1]=='R'):
            Mr, Mr_e=convertRtor(M_u[0, :], M_e[0, :],M_u[1, :], M_e[1, :])
            Mg, Mg_e = convertVtog(M_u[0, :], M_e[0, :], M_u[1, :], M_e[1, :],Mr, Mr_e)
            M_u[0, :], M_e[0, :]=Mg, Mg_e
            M_u[1, :], M_e[1, :]=Mr, Mr_e
            self.band_extinction[0]='g'
            folder='GmX_color_templates'
            scale=0.94
        filename='{}m{}_template_{}_kcorr.mat.txt'.format(self.band_extinction[0].lower(),self.band_extinction[1].lower(),self.sn_type)
        temp=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/templates/{}/'.format(folder)+filename)
        if self.verbose:
            plt.figure(1)
            plt.plot(t_u, M_u[0, :],'b')
            plt.plot(t_u, M_u[1, :], 'r')
            plt.scatter(t,M[u],s=5,color='g')
            plt.scatter(mag[:, 0][mag[:, 3] == self.band_extinction[0]].astype(float)-t0, mag[:, 1][mag[:, 3] == self.band_extinction[0]].astype(float), s=5)
            plt.title(self.name)
            plt.gca().invert_yaxis()
            plt.figure(2)
            plt.plot(t_u, M_u[0, :] - M_u[1, :], 'r')
            plt.plot(temp[:,0],temp[:,1])
            plt.title(self.name)
            plt.show()

        color_obs=M_u[0, :] - M_u[1, :]
        color_obs_e=np.sqrt(M_e[0, :]**2+M_e[1, :]**2)
        fit_extinction = ErrorPropagationSpline(temp[:,0], temp[:,1], temp[:,2], s=0.2)
        color_in, color_in_e = fit_extinction(t_u)
        self.ebv_host[0]=np.mean(color_obs[(t_u < 10) & (t_u > 5)]-color_in[(t_u < 10) & (t_u > 5)])/scale
        if self.ebv_host[0]<0:
            self.ebv_host[0]=0
        self.ebv_host[1]= np.mean(np.sqrt(color_obs_e[(t_u < 10) & (t_u > 5)]**2+ color_in_e[(t_u < 10) & (t_u > 5)]**2))
    def tailNickel(self):
        add = []
        location='./Data/SN_json/'
        with open(location+self.name+'.json') as data_file:
            data = json.load(data_file)
        mag = np.zeros(shape=(0, 6))
        DM = 5 * np.log10(self.distance[0]) - 5
        DMe = np.abs(5 * self.distance[1] / ( self.distance[0]* np.log(10)))
        for dat in data[self.name]["photometry"]:
            if "band" in dat.keys():
                if np.any(np.in1d(dat["band"], self.band_bc)):
                    if "e_magnitude" in dat:
                        error = float(dat["e_magnitude"])
                    else:
                        error = 0
                    if np.any(np.in1d(dat["source"].split(','), self.source_bc)):
                        add = np.concatenate(([float(dat["time"]), dat["magnitude"], error],
                                              [deredMag(float(dat["magnitude"]), self.ebv_host[0]+self.ebv_gal[0], coef[dat["band"]]) - DM,
                                               np.sqrt(error ** 2 + DMe ** 2), dat["band"]]))
                        add = np.reshape(add, (1, 6))
                        mag = np.concatenate((mag, add), axis=0)
        s=self.smoothing_bc
        s=1
        if self.name == 'SN2016gkg':
            t_u = np.arange(10, 100, 0.1)
        elif self.name== 'iPTF13bvn':
            t_u = np.arange(10, 90, 0.1)
        elif self.name== 'SN2008D':
            t_u = np.arange(10, 120, 0.1)
        elif self.name == 'SN1993J':
            t_u = np.arange(10, 100, 0.1)
        elif self.name == 'SN2013df':
            t_u = np.arange(10, 200, 0.1)
        elif self.name == 'SN2006jc':
            t_u = np.arange(2, 65, 0.1)
        elif self.name == 'SN1994I':
            t_u = np.arange(5, 100, 0.1)
        elif self.name == 'SN2004aw':
            t_u = np.arange(5, 80, 0.1)
        elif self.name == 'SN2009bb':
            t_u = np.arange(3, 45, 0.1)
        elif self.name == 'SN2007ru':
            t_u = np.arange(1, 130, 0.1)
        elif self.name == 'SN2002ap':
            t_u = np.arange(1, 40, 0.1)
        elif self.name=='SN2007gr':
            t_u = np.arange(10, 130, 0.1)
        else:
            t_u = np.arange(5, 130, 0.1)
        t_u_t = np.clip(t_u, a_min=60, a_max=130)
        if self.name == 'SN2009bb':
            t_u_t = np.arange(60, 75, 0.1)
        if self.name  == 'SN2002ap':
            t_u_t = np.arange(50, 230, 0.1)
        if self.name  == 'SN2006jc':
            t_u_t = np.arange(50, 100, 0.1)
        t_u_t = np.unique(t_u_t)
        M_u_t = np.zeros(shape=(len(self.band_bc), np.shape(t_u_t)[0]))
        M_e_t = np.zeros(shape=(len(self.band_bc), np.shape(t_u_t)[0]))
        t_u=t_u[t_u<30]
        M_u = np.zeros(shape=(len(self.band_bc), np.shape(t_u)[0]))
        M_e = np.zeros(shape=(len(self.band_bc), np.shape(t_u)[0]))
        for j, band in enumerate(self.band_bc):
            t = mag[:, 0][mag[:, 5] == band].astype(float)
            M = mag[:, 3][mag[:, 5] == band].astype(float)
            Me = mag[:, 4][mag[:, 5] == band].astype(float)
            M = M[np.argsort(t)]
            Me = Me[np.argsort(t)]
            t = t[np.argsort(t)] - self.t0[0]
            t, u = np.unique(t, return_index=True)
            M=M[u]
            Me=Me[u]
            fit = ErrorPropagationSpline(t[t<30], M[t<30], Me[t<30], s=s)
            M_u[j, :], M_e[j, :] = fit(t_u)
            fit_linear = ErrorPropagationLinear(t[(t > 50)], M[(t > 50)], Me[(t > 50)])
            M_u_t[j, :], M_e_t[j, :] = fit_linear(t_u_t)
            if (self.name == 'SN2016gkg'):
                fit = UnivariateSpline(t, M[u], s=0.1)
                M_u[j, :] = fit(t_u)
        if self.name=='SN2008D':
            M_u[1, :] = M_u[2, :] - 1.2444 * (M_u[2, :] - M_u[1, :]) - 0.3820
            self.band_bc[1]='I'
        lbol,le=lyman_BC(np.reshape(M_u[0, :],(1,np.shape(M_u[0, :])[0])), np.reshape(M_u[1, :],(1,np.shape(M_u[1, :])[0])),self.band_bc[0],self.band_bc[1])
        lbol_t,le_t=lyman_BC(np.reshape(M_u_t[0, :],(1,np.shape(M_u_t[0, :])[0])), np.reshape(M_u_t[1, :],(1,np.shape(M_u_t[1, :])[0])),self.band_bc[0],self.band_bc[1])
        lbol = np.reshape(lbol, t_u.shape)
        lbol_t = np.reshape(lbol_t, t_u_t.shape)
        le = np.reshape(le, t_u.shape)
        le_t = np.reshape(le_t, t_u_t.shape)
        if self.verbose:
            plt.figure(1)
            plt.scatter(mag[:, 0][mag[:, 5] == self.band_bc[0]].astype(float) - self.t0[0],
                        mag[:, 3][mag[:, 5] == self.band_bc[0]].astype(float), label=self.band_bc[0])
            plt.scatter(mag[:, 0][mag[:, 5] == self.band_bc[1]].astype(float) - self.t0[0],
                        mag[:, 3][mag[:, 5] == self.band_bc[1]].astype(float), label=self.band_bc[1])
            plt.scatter(t_u_t, M_u_t[0, :], s=0.5)
            plt.scatter(t_u_t, M_u_t[1, :], s=0.5)
            plt.legend()
            plt.title(self.name)
            plt.gca().invert_yaxis()
            plt.figure(2)
            plt.scatter(mag[:, 0][mag[:, 5] == self.band_bc[0]].astype(float) - self.t0[0],
                        mag[:, 3][mag[:, 5] == self.band_bc[0]].astype(float), label=self.band_bc[0])
            plt.scatter(mag[:, 0][mag[:, 5] == self.band_bc[1]].astype(float) - self.t0[0],
                        mag[:, 3][mag[:, 5] == self.band_bc[1]].astype(float), label=self.band_bc[1])
            plt.scatter(t_u, M_u[0, :], s=0.5)
            plt.scatter(t_u, M_u[1, :], s=0.5)
            plt.legend()
            plt.title(self.name)
            plt.gca().invert_yaxis()
            plt.figure(3)
            plt.plot(t_u_t[(t_u_t >= 60) & (t_u_t < 120)],lbol_t[(t_u_t >= 60) & (t_u_t < 120)], color="orange")
            plt.plot(t_u,lbol, color="red")
            plt.show()

        self.peakt = t_u[np.argmax(lbol)]
        self.peakL[0]= np.log10(np.max(lbol))
        self.peakL[1]=np.log10(le[np.argmax(lbol)])
        p, p_err = fit_bootstrap([0.1, 3], t_u_t[(t_u_t >= 60) & (t_u_t < 120)], lbol_t[(t_u_t >= 60) & (t_u_t < 120)],
                                 np.clip(le_t[(t_u_t >= 60) & (t_u_t < 120)], a_min=0.0001e42, a_max=1e51), valentiErr,
                                 errfunc=True, perturb=True, n=1000, nproc=4)
        self.tail_ni_mass=[p[0], p_err[0]]
        self.tail_meje=[p[1], p_err[1]]
        self.arnett_ni_mass=list(valenti_ni56(self.peakt,10**self.peakL[0],0,0,10**self.peakL[1]))
    def khatami_model(self):
        beta0 = 0.3
        betafit = fit_bootstrap([0.3], [self.peakt, self.tail_ni_mass[0]],[10**self.peakL[0]],[10**self.peakL[1]], khatami_err,
                                 errfunc=True, perturb=True, n=1000, nproc=4)
        self.beta_req= [betafit[0][0],betafit[1][0]]
    def prentice_data(self):
        url='https://docs.google.com/spreadsheets/d/e/2PACX-1vTRiEGa3J6bnG5R4JhDz6Wx2gb_KK9TnRdkn5YdWIoar6_ugH1p42KYLkDMHiMjjOUxXPxaPH22wjRN/pub?gid=0&single=true&output=csv'
        with requests.Session() as s:
            download = s.get(url)
            decoded_content = download.content.decode('utf-8')
            cr = csv.reader(decoded_content.splitlines(), delimiter=',')
            my_list= list(cr)
        prentice = np.array(my_list)
        prentice_names = []
        for j, sn in enumerate(prentice):
            prentice_names.append(sn[0])
        try:
            ind = prentice_names.index(self.name.strip('SN'))
            self.prentice_peakL = ' '.join(prentice[ind, 2].split()).split(" ")
            self.prentice_peakt= ' '.join(prentice[ind, 3].split()).split(" ")
            try:
                self.prentice_ni= ' '.join(prentice[ind,4].split()).split(" ")
            except:
                self.prentice_ni = None
        except:
            self.prentice_peakL = None
            self.prentice_peakt = None
            self.prentice_ni = None


plt.close("all")
url='https://docs.google.com/spreadsheets/d/e/2PACX-1vS3r49T74bVSHHYNJrkSwcc7O7-nXdmd9ERmfdO8lNKAocdrhjT4H7znvILr8nN8BERXZpM2_cx3ge0/pub?gid=0&single=true&output=csv'
with requests.Session() as s:
    download = s.get(url)
    decoded_content = download.content.decode('utf-8')
    cr = csv.reader(decoded_content.splitlines(), delimiter=',')
    my_list = list(cr)
info=np.array(my_list)

sn = supernova('SN2008D')
sn.fill(info)
index = np.where(info[:, 0] == sn.name)[0][0]
sn.stritzinger_color(float(info[index, 6]))
print "host E(B-V)", sn.ebv_host
sn.tailNickel()
print "Tail Ni", sn.tail_ni_mass, "Tail Mej/E", sn.tail_meje, "Tail Arnett", sn.arnett_ni_mass, "Lp", sn.peakL, "tp", sn.peakt
sn.khatami_model()
sn.prentice_data()
print "beta", sn.beta_req
print "t0", sn.t0

with open('Data/out.csv','w') as f:
    w = csv.DictWriter(f,fieldnames=sorted(vars(sn)))
    w.writeheader()
    w.writerow({k: getattr(sn, k) for k in vars(sn)})
    for i,sn_name in enumerate(info[1:,0]):
        if  (sn_name =="SN2006jc") |(sn_name=='SN2002ap') | (sn_name=='SN2008D'):
             continue
        # if (sn_name!='SN1998bw'):
        #     continue
        print sn_name
        sn = supernova(sn_name)
        sn.fill(info)
        index = np.where(info[:, 0] == sn.name)[0][0]
        sn.stritzinger_color(float(info[index, 6]))
        print "host E(B-V)",sn.ebv_host
        sn.tailNickel()
        print "Tail Ni",sn.tail_ni_mass, "Tail Mej/E",sn.tail_meje, "Tail Arnett",sn.arnett_ni_mass,"Lp", sn.peakL, "tp",sn.peakt
        sn.khatami_model()
        sn.prentice_data()
        print "beta",sn.beta_req
        print "t0",sn.t0
        w.writerow({k: getattr(sn, k) for k in vars(sn)})



