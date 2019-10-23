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
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
matplotlib.rcParams['axes.linewidth'] = 1.5 #set the value globally
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['xtick.major.width'] = 2
matplotlib.rcParams['xtick.minor.size'] = 2
matplotlib.rcParams['xtick.minor.width'] = 1.5
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['ytick.major.width'] = 2
matplotlib.rcParams['ytick.minor.size'] = 2
matplotlib.rcParams['ytick.minor.width'] = 1.5
matplotlib.rcParams.update({'font.size': 14})
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'Courier New'})
SN_verified_t0=['SN2008D','SN2016gkg','SN2011dh','SN2013df','SN1993J','SN1998bw']
sug_beta={'Ib':1.125, 'Ic':1.125, 'IIb':0.82, 'Ia':1.6}
coef = {'U':4.334,'g': 3.303,'i':1.698, 'r':2.285,'B': 3.626, 'V': 2.742, 'I': 1.505, "i'": 1.698, "r'": 2.285, "R": 2.169, 'r':2.285}


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
        self.verbose=1
        self.source_bc=source_bc
        self.smoothing_bc=smoothing_bc
        self.tail_ni_mass=[None, None]
        self.tail_meje=[None, None]
        self.peakL=[None, None]
        self.peakt=[None, None]
        self.arnett_ni_mass=[None, None]
        self.band_t0=None
        self.lb=None
        self.ub=None
        self.ni_khatami=[None, None]
    def fill(self,info):
        location = './Data/SN_json/'
        if os.path.isfile(location + self.name + '.json') == False:
            url='https://sne.space/astrocats/astrocats/supernovae/output/json/'+self.name+'.json'
            urllib.urlretrieve(url, location + self.name  + '.json')
        with open(location + self.name  + '.json') as data_file:
            data = json.load(data_file)
        index=np.where(info[:,0]==self.name)[0][0]
        self.sn_type = info[index , np.where(info[0,:]=='Type')[0][0]].split(' ')[0]
        self.host = info[index, np.where(info[0,:]=='Host')[0][0]]
        self.distance[0] = float(info[index , np.where(info[0,:]=='luminosity distance')[0][0]].split('(')[0])*1e6
        if (len(info[index , np.where(info[0,:]=='luminosity distance')[0][0]].split('(')) > 1):
            self.distance[1] = float(info[index, np.where(info[0,:]=='luminosity distance')[0][0]].split('(')[1].replace(")", ""))*1e6
        else:
            self.distance[1]=0
        self.band_extinction = info[index, np.where(info[0,:]=='extinction color')[0][0]].split(" ")[0].split("-")
        self.band_bc = info[index ,  np.where(info[0,:]=='Color')[0][0]].split(',')[0].split("-")
        self.source = info[index, np.where(info[0,:]=='source')[0][0]].split(',')
        self.ebv_gal[0] = float(info[index , np.where(info[0,:]=='E(B-V) (S&F2011)')[0][0]].split('(')[0])
        if (len(info[index , np.where(info[0,:]=='E(B-V) (S&F2011)')[0][0]].split('(')) > 1):
            self.ebv_gal[1] = float(info[index , np.where(info[0,:]=='E(B-V) (S&F2011)')[0][0]].split('(')[1].replace(")", ""))
        else:
            self.ebv_gal[1]=0
        try:
            self.ebv_host[0] = float(info[index , np.where(info[0,:]=='E(B-V)_hostf')[0][0]].split('(')[0])
        except:
            self.ebv_host[0] =0.0
        if (len(info[index , np.where(info[0,:]=='E(B-V)_hostf')[0][0]].split('(')) > 1):
            self.ebv_host[1] = float(info[index , np.where(info[0,:]=='E(B-V)_hostf')[0][0]].split('(')[1].replace(")", ""))
        else:
            self.ebv_host[1]=0
        self.source_bc=info[index, np.where(info[0,:]=='source_number')[0][0]].split(',')
        self.smoothing_bc = float(info[index, np.where(info[0, :] == 'smoothing_bc')[0][0]])
        try:
            self.t0[0] = float(info[index, np.where(info[0, :] == 't0_fit')[0][0]].split('(')[0])-2400000.5
        except:
            self.t0[0] =0
        if self.name not in SN_verified_t0:
            self.band_t0 = info[index, np.where(info[0, :] == 't0_band')[0][0]]
        if (len(info[index , np.where(info[0,:]=='t0_fit')[0][0]].split('(')) > 1):
            self.t0[1] = float(info[index , np.where(info[0,:]=='t0_fit')[0][0]].split('(')[1].replace(")", ""))
        else:
            self.t0[1] =0.0
        try:
            self.lb = float(info[index, np.where(info[0, :] == 'lb,ub')[0][0]].split(',')[0])
            self.ub = float(info[index, np.where(info[0, :] == 'lb,ub')[0][0]].split(',')[1])
        except:
            self.lb =None
            self.ub =None
        if 1:
            print self.name,self.host,self.sn_type,self.distance,self.band_extinction,self.band_bc,self.source,self.ebv_gal,self.t0
    def fit_t0(self):
        from scipy.optimize import curve_fit
        add = []
        location='./Data/SN_json/'
        with open(location+self.name+'.json') as data_file:
            data = json.load(data_file)
        mag = np.zeros(shape=(0, 8))
        DM = 5 * np.log10(self.distance[0]) - 5
        DMe = np.abs(5 * self.distance[1] / ( self.distance[0]* np.log(10)))
        if (self.name!='SN2017gax') & (self.name!='SN2003jd'):
            for dat in data[self.name]["photometry"]:
                try:
                    if np.any(np.in1d(dat["band"], self.band_t0)):
                            if "e_magnitude" in dat:
                                error = float(dat["e_magnitude"])
                            else:
                                error = 0
                            if 1:#np.any(np.in1d(dat["source"].split(','), self.source_bc)):
                                add = np.concatenate(([float(dat["time"]), dat["magnitude"], error],
                                                      [deredMag(float(dat["magnitude"]), self.ebv_host[0]+self.ebv_gal[0], coef[dat["band"]]/0.86) - DM,
                                                       np.sqrt(error ** 2 + DMe ** 2 + (coef[dat["band"]] / 0.86) ** 2 * (
                                                                   self.ebv_host[1] ** 2 + self.ebv_gal[1] ** 2)), dat["band"],0,0]))
                                add  = np.reshape(add, (1, 8))
                                mag = np.concatenate((mag, add), axis=0)
                except:
                    continue
        else:
            data = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/{}_lcophot.csv'.format(self.name), delimiter=',',
                              usecols=(0, 1, 2))
            data=data.reshape((data.shape[0], 3))
            mag = np.zeros(shape=(data.shape[0], 8))
            filter = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/{}_lcophot.csv'.format(self.name),
                                delimiter=',', usecols=(3,), dtype='|S8').reshape((data.shape[0],1))
            x=np.zeros((data.shape[0],1))
            mags = np.zeros((data.shape[0], 2))
            for  p,lab in enumerate(filter):
                try:
                    x[p,0]=coef[lab[0]] / 0.86
                except:
                    continue
                mags[p,0]=deredMag(data[p, 1], self.ebv_host[0] + self.ebv_gal[0],x[p,0] / 0.86)-DM
                mags[p,1] =np.sqrt(data[p,2] ** 2 + DMe ** 2 + (x[p,0] / 0.86) ** 2 * (self.ebv_host[1] ** 2 + self.ebv_gal[1] ** 2))
            mag=np.concatenate((data,mags,filter),axis=1)
            mag=np.concatenate((mag,np.zeros((data.shape[0],2))),axis=1)
            mag=mag[mag[:,5]==self.band_t0,:]
        for j,f in enumerate(mag[:,0]):
            mag[j,6:8]=Mag_toFlux(self.band_t0,float(mag[j,3]),float(mag[j,4]))
        flux=np.zeros(shape=(mag.shape[0],3))
        flux[:,1:3]= mag[:,6:8].astype(float)
        flux[:, 1:3]=flux[:,1:3]/np.max(flux[:,1])
        flux[:,0]= mag[:,0].astype(float)
        t_mjd=np.min(flux[:, 0])
        flux[:,0]=flux[:,0]-np.min(flux[:,0])
        u=np.argsort(flux[:,0])
        flux[flux[:,2]>1,2]=0.1
        flux[:, 1:3]=flux[u, 1:3]
        flux[:, 0] = flux[u, 0]
        flux=flux[(flux[:,0]<self.ub) & (flux[:,0]>=self.lb),:]
        t, indices =np.unique(flux[:, 0],return_index=True)
        flux1 = np.zeros(shape=(t.shape[0], 3))
        flux1[:, 1:3] = flux[indices, 1:3]
        flux1[:,0]=flux[indices, 0]
        print flux1
        def power_law(t, t0, a, n):
            return a*np.sign(t - t0) * (np.abs(t - t0))**n
        p0=[-5 ,  0.0001 ,  2.5]
        lb = [-20, 0.00001, 0.01]
        ub = [0,1,5]
        #p, p_err = fit_bootstrap(p0, flux1[:,0], flux1[:,1], flux1[:,2], powerLawErr, errfunc=True, perturb=False, n=3, nproc=4)
        # # print p,p_err
        n=2000
        randomDelta = np.random.normal(0., flux1[:,2]/10, (n,len(flux1[:,1])))
        randomdataY = flux1[:,1] + randomDelta
        # p, pcov = curve_fit(power_law, flux[:,0], flux[:,1],p0=p0,maxfev=10000000)
        x=np.zeros(shape=(3,n))
        x_err = np.zeros(shape=(3, 3,n))
        for i,randY in enumerate(randomdataY):
            x[:,i],x_err[:,:,i] = curve_fit(power_law, flux1[:, 0], randY, p0=p0,  bounds=(lb,ub),maxfev=1000000)

        p0=np.nanmedian(x,1)
        p0_err=np.nanstd(x,1)
        print p0, p0_err
        if self.verbose:
            plt.errorbar(flux[:,0],flux[:,1],yerr=flux[:,2],fmt='o',ls=None)
            plt.xlabel('Time (days)')
            plt.ylabel('Scaled flux')
            plt.plot(flux[:,0],power_law(flux[:,0],p0[0],p0[1],p0[2]))
            plt.show()
        self.t0=[p0[0]+t_mjd, p0_err[0]]
        print (self.t0)
    def drout_color(self,s=0.2):
        add = []
        location='./Data/SN_json/'
        with open(location+self.name+'.json') as data_file:
            data = json.load(data_file)
        mag = np.zeros(shape=(0, 4))
        if self.name!='SN2017gax':
            for dat in data[self.name]["photometry"]:
               if "band" in dat.keys():
                    if  np.any(np.in1d(dat["band"],self.band_extinction)):
                        if "e_magnitude" in dat:
                            error=float(dat["e_magnitude"])
                        else:
                            error=0
                        if np.any(np.in1d(dat["source"].split(','),self.source)):
                            add = [float(dat["time"]), deredMag(float(dat["magnitude"]), float(self.ebv_gal[0]), coef[dat["band"]]/0.86), np.sqrt(error**2 + ((coef[dat["band"]]/0.86)**2*float(self.ebv_gal[1])**2)), dat["band"]]
                            add = np.reshape(add, (1, 4))
                            mag = np.concatenate((mag, add), axis=0)
        else:
            mag = np.zeros(shape=(770, 4))
            data = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SN2017gax_lcophot.csv', delimiter=',',
                              usecols=(0, 1, 2)).reshape((770, 3))
            filter = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SN2017gax_lcophot.csv',
                                delimiter=',', usecols=(3,), dtype='|S8').reshape((770, 1))
            x = np.zeros((770, 1))
            mags = np.zeros((770, 2))
            for p, lab in enumerate(filter):
                x[p, 0] = coef[lab[0]] / 0.86
                mags[p, 0] = deredMag(data[p, 1], self.ebv_host[0] + self.ebv_gal[0], x[p, 0] / 0.86) - DM
                mags[p, 1] = np.sqrt(
                    data[p, 2] ** 2 + DMe ** 2 + (x[p, 0] / 0.86) ** 2 * (self.ebv_host[1] ** 2 + self.ebv_gal[1] ** 2))
            mag = np.concatenate((data[0:770,0], mags, filter), axis=1)
            mag = mag[np.isin(mag[:, 3],self.band_extinction), :]
            print mag.shape
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
            M_u[j, :], M_e[j, :],dummy = fit(t_u)
            if (self.name == 'SN2016gkg'):
                fit = UnivariateSpline(t, M[u], s=0.02)
                M_u[j, :] = fit(t_u)
                M_e[j, :]=0
        scale=1
        if (self.band_extinction[0]=='V' ) & (self.band_extinction[1]=='R'):
            scale=(coef['V']-coef['R'])/0.86
        if self.verbose:
            plt.figure(1)
            plt.plot(t_u, M_u[0, :],'b')
            plt.plot(t_u, M_u[1, :], 'r')
            plt.scatter(t,M[u],s=5,color='g',label=self.band_extinction[1])
            plt.scatter(mag[:, 0][mag[:, 3] == self.band_extinction[0]].astype(float)-t0, mag[:, 1][mag[:, 3] == self.band_extinction[0]].astype(float), s=5,label=self.band_extinction[0])
            plt.legend()
            plt.title(self.name)
            plt.gca().invert_yaxis()
            plt.figure(2)
            plt.plot(t_u, M_u[0, :] - M_u[1, :], 'r')
            plt.title(self.name)
            plt.show()
        color_obs=M_u[0, :] - M_u[1, :]
        color_obs_e=np.sqrt(M_e[0, :]**2+M_e[1, :]**2)
        print scale,color_obs[np.argmin(np.abs(t_u -10))]
        self.ebv_host[0]=(color_obs[np.argmin(np.abs(t_u -10))]-0.26)/scale
        if self.ebv_host[0]<0:
            self.ebv_host[0]=0
        self.ebv_host[1]= np.sqrt(color_obs_e[np.argmin(np.abs(t_u -10))]**2+ 0.06**2)

    def stritzinger_color(self,s=0.2):
        add = []
        location='./Data/SN_json/'
        with open(location+self.name+'.json') as data_file:
            data = json.load(data_file)
        mag = np.zeros(shape=(0, 4))
        if (self.name!='SN2017gax') & (self.name!='SN2003jd'):
            for dat in data[self.name]["photometry"]:
               if "band" in dat.keys():
                    if  np.any(np.in1d(dat["band"],self.band_extinction)):
                        if "e_magnitude" in dat:
                            error=float(dat["e_magnitude"])
                        else:
                            error=0
                        if np.any(np.in1d(dat["source"].split(','),self.source)):
                            add = [float(dat["time"]), deredMag(float(dat["magnitude"]), float(self.ebv_gal[0]), coef[dat["band"]]/0.86), np.sqrt(error**2 + ((coef[dat["band"]]/0.86)**2*float(self.ebv_gal[1])**2)), dat["band"]]
                            add = np.reshape(add, (1, 4))
                            mag = np.concatenate((mag, add), axis=0)
        else:
            data = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/{}_lcophot.csv'.format(self.name), delimiter=',',
                              usecols=(0, 1, 2))
            data=data.reshape((data.shape[0], 3))
            mag = np.zeros(shape=(data.shape[0], 4))
            filter = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/{}_lcophot.csv'.format(self.name),
                                delimiter=',', usecols=(3,), dtype='|S8').reshape((data.shape[0], 1))
            x = np.zeros((data.shape[0], 1))
            mags = np.zeros((data.shape[0], 2))
            for p, lab in enumerate(filter):
                try:
                    x[p, 0] = coef[lab[0]] / 0.86
                except:
                    continue
                mags[p, 0] = deredMag(data[p, 1], self.ebv_gal[0], x[p, 0] / 0.86)
                mags[p, 1] = np.sqrt(
                    data[p, 2] ** 2 + (x[p, 0] / 0.86) ** 2 * (self.ebv_gal[1] ** 2))
            mag = np.concatenate((data[0:data.shape[0],0:1], mags, filter), axis=1)
            mag = mag[np.in1d(mag[:, 3],self.band_extinction), :]
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
            M_u[j, :], M_e[j, :],dumm = fit(t_u)
            if (self.name == 'SN2016gkg'):
                fit = UnivariateSpline(t, M[u], s=0.02)
                M_u[j, :] = fit(t_u)
                M_e[j, :]=0
        folder='{}mX_color_templates'.format(self.band_extinction[0])
        scale=1
        if (self.band_extinction[0]=='V' ) & (self.band_extinction[1]=='R'):
            Mr, Mr_e=convertRtor(M_u[0, :], M_e[0, :],M_u[1, :], M_e[1, :])
            M_u[1, :], M_e[1, :]=Mr, Mr_e
            self.band_extinction[1]='r'
            folder='VmX_color_templates'
            scale=(coef['V']-coef['R'])/0.86
        # if (self.band_extinction[0]=='V' ) & (self.band_extinction[1].replace("'","")=="r"):
        #     scale = (coef[self.band_extinction[0]] - coef[self.band_extinction[1]])/0.86
        scale = (coef[self.band_extinction[0]] - coef[self.band_extinction[1]]) / 0.86
        filename='{}m{}_template_{}_kcorr.mat.txt'.format(self.band_extinction[0].lower(),self.band_extinction[1].lower().strip("'"),self.sn_type)
        print filename,folder
        temp=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/templates/{}/'.format(folder)+filename)
        if self.verbose:
            plt.figure(1)
            plt.plot(t_u, M_u[0, :],'b')
            plt.plot(t_u, M_u[1, :], 'r')
            plt.scatter(t,M[u],s=5,color='g',label=self.band_extinction[1])
            plt.scatter(mag[:, 0][mag[:, 3] == self.band_extinction[0]].astype(float)-t0, mag[:, 1][mag[:, 3] == self.band_extinction[0]].astype(float), s=5,label=self.band_extinction[0])
            plt.xlim(-10, 30)
            plt.legend()
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
        color_in, color_in_e,dumm = fit_extinction(t_u)
        #print color_in_e
        self.ebv_host[0]=np.mean(color_obs[(t_u < 10) & (t_u > 5)]-color_in[(t_u < 10) & (t_u > 5)])/scale
        if self.ebv_host[0]<0:
            self.ebv_host[0]=0
        self.ebv_host[1]=np.sqrt(np.mean(np.sqrt(color_obs_e[(t_u < 10) & (t_u > 5)]**2+ color_in_e[(t_u < 10) & (t_u > 5)]**2))**2+np.std(color_obs[(t_u < 10) & (t_u > 5)]-color_in[(t_u < 10) & (t_u > 5)])**2)
           # np.mean(np.sqrt(color_obs_e[(t_u < 10) & (t_u > 5)]**2+ color_in_e[(t_u < 10) & (t_u > 5)]**2))
    def tailNickel(self):
        add = []
        location='./Data/SN_json/'
        with open(location+self.name+'.json') as data_file:
            data = json.load(data_file)
        mag = np.zeros(shape=(0, 6))
        DM = 5 * np.log10(self.distance[0]) - 5
        DMe = np.abs(5 * self.distance[1] / ( self.distance[0]* np.log(10)))
        if (self.name != 'SN2017gax') & (self.name != 'SN2003jd'):
            for dat in data[self.name]["photometry"]:
                if "band" in dat.keys():
                    if np.any(np.in1d(dat["band"], self.band_bc)):
                        if "e_magnitude" in dat:
                            error = float(dat["e_magnitude"])
                        else:
                            error = 0
                        if np.any(np.in1d(dat["source"].split(','), self.source_bc)):
                            add = np.concatenate(([float(dat["time"]), dat["magnitude"], error],
                                                  [deredMag(float(dat["magnitude"]), self.ebv_host[0]+self.ebv_gal[0], coef[dat["band"]]/0.86) - DM,
                                                   np.sqrt(error ** 2 + DMe ** 2+(coef[dat["band"]]/0.86)**2*(self.ebv_host[1]**2+self.ebv_gal[1]**2)), dat["band"]]))
                            add = np.reshape(add, (1, 6))
                            mag = np.concatenate((mag, add), axis=0)
        else:
            data = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/{}_lcophot.csv'.format(self.name), delimiter=',',
                              usecols=(0, 1, 2))
            num_datapnts = data.shape[0]
            data=data.reshape((num_datapnts, 3))
            mag = np.zeros(shape=(num_datapnts, 6))
            filter = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/{}_lcophot.csv'.format(self.name),
                                delimiter=',', usecols=(3,), dtype='|S8').reshape((num_datapnts, 1))
            x = np.zeros((num_datapnts, 1))
            mags = np.zeros((num_datapnts, 2))
            data[:,0]=data[:,0]-2400000.5

            for p, lab in enumerate(filter):
                try:
                    x[p, 0] = coef[lab[0]] / 0.86
                except:
                    continue
                mags[p, 0] = deredMag(data[p, 1], self.ebv_host[0] + self.ebv_gal[0], x[p, 0] / 0.86) - DM
                mags[p, 1] = np.sqrt(
                    data[p, 2] ** 2 + DMe ** 2 + (x[p, 0] / 0.86) ** 2 * (self.ebv_host[1] ** 2 + self.ebv_gal[1] ** 2))
            mag = np.concatenate((data[0:num_datapnts,0:3], mags, filter), axis=1)
            mag = mag[np.in1d(mag[:, 5],self.band_bc), :]
        s=self.smoothing_bc
        if self.name == 'SN2016gkg':
            t_u = np.arange(10, 100, 0.1)
        if self.name == 'SN2008ax':
            t_u = np.arange(10, 130, 0.1)
        elif self.name== 'iPTF13bvn':
            t_u = np.arange(10, 90, 0.1)
        elif self.name== 'SN2008D':
            t_u = np.arange(10, 120, 0.1)
        elif self.name == 'SN1993J':
            t_u = np.arange(10, 100, 0.1)
        elif self.name == 'SN2013df':
            t_u = np.arange(20, 200, 0.1)
        elif self.name == 'SN2006jc':
            t_u = np.arange(2, 65, 0.1)
        elif self.name == 'SN1994I':
            t_u = np.arange(3, 100, 0.1)
        elif self.name == 'SN2004aw':
            t_u = np.arange(5, 80, 0.1)
        elif self.name == 'SN2009bb':
            t_u = np.arange(3, 45, 0.1)
        elif self.name == 'SN2007ru':
            t_u = np.arange(3, 190, 0.1)
        elif self.name == 'SN2002ap':
            t_u = np.arange(1, 40, 0.1)
        elif self.name=='SN2007gr':
            t_u = np.arange(3, 130, 0.1)
        elif self.name=='SN1996cb':
            t_u = np.arange(3, 100, 0.1)
        elif self.name=='SN2006el':
            t_u = np.arange(12, 130, 0.1)
        elif self.name=='SN2011bm':
            t_u = np.arange(10, 130, 0.1)
        elif self.name=='SN2005bf':
            t_u = np.arange(30, 140, 0.1)
        elif self.name=='SN2009jf':
            t_u = np.arange(10, 140, 0.1)
        elif self.name=='SN2003jd':
            t_u = np.arange(7, 130, 0.1)
        else:
            t_u = np.arange(5, 130, 0.1)
        if (self.name=='SN2016gkg') | (self.name=='SN1994I'):
            t_u_t = np.clip(t_u, a_min=50, a_max=130)
        elif self.name=='SN2013df':
            t_u_t = np.clip(t_u, a_min=40, a_max=200)
        elif self.name=='SN2006ep':
            t_u_t = np.clip(t_u, a_min=40, a_max=200)
        elif self.name=='SN2007ru':
            t_u_t = np.clip(t_u, a_min=60, a_max=200)
        else:
            t_u_t = np.clip(t_u, a_min=60, a_max=130)
        if self.name == 'SN2009bb':
            t_u_t = np.arange(50, 100, 0.1)
        if self.name  == 'SN2002ap':
            t_u_t = np.arange(38, 230, 0.1)
        if self.name  == 'SN2006jc':
            t_u_t = np.arange(50, 100, 0.1)
        if self.name=='SN2005bf':
            t_u_t = np.arange(75, 130, 0.1)
        t_u_t = np.unique(t_u_t)
        M_u_t = np.zeros(shape=(len(self.band_bc), np.shape(t_u_t)[0]))
        M_e_t = np.zeros(shape=(len(self.band_bc), np.shape(t_u_t)[0]))
        if (self.name=='SN2011bm'):
            ub=50
            t_u=t_u[(t_u<(ub-10))]
        elif (self.name=='SN2005bf'):
            ub=58
            t_u=t_u[(t_u<(ub-10))]
        elif (self.name=='SN1994I'):
            ub=30
            t_u=t_u[(t_u<(ub-10))]
        elif (self.name=='SN2009jf'):
            ub=50
            t_u=t_u[(t_u<(ub-10))]
        elif (self.name=='SN2003jd'):
            ub=40
            t_u=t_u[(t_u<(ub-10))]
        else:
            ub=50
            t_u=t_u[(t_u<(ub-10))]
        M_u = np.zeros(shape=(len(self.band_bc), np.shape(t_u)[0]))
        M_e = np.zeros(shape=(len(self.band_bc), np.shape(t_u)[0]))
        tpeak_err = np.zeros(shape=(len(self.band_bc), np.shape(t_u)[0]))
        for j, band in enumerate(self.band_bc):
            t = mag[:, 0][mag[:, 5] == band].astype(float)
            M = mag[:, 3][mag[:, 5] == band].astype(float)
            Me = mag[:, 4][mag[:, 5] == band].astype(float)
            M = M[np.argsort(t)]
            Me = Me[np.argsort(t)]
            #print t,self.t0, t-self.t0[0]
            t = t[np.argsort(t)] - self.t0[0]
            t, u = np.unique(t, return_index=True)
            M=M[u]
            Me=Me[u]
            if (self.name == 'SN2016gkg'):
                fit = UnivariateSpline(t, M[u], s=0.02)
                M_u[j, :] = fit(t_u)
                tpeak_err[j, :]=4
            else:
                if (self.name=='SN2006el' )| (self.name=='SN2011bm') | (self.name=='SN2003jd'):
                    Me=Me/5
                #print M,t
                fit = ErrorPropagationSpline(t[(t<ub)], M[(t<ub)], Me[(t<ub)], s=1)
                M_u[j, :], M_e[j, :],tpeak_err[j, :] = fit(t_u)
            fit_linear = ErrorPropagationLinear(t[(t > np.min(t_u_t)) & (t<np.max(t_u_t))], M[(t > np.min(t_u_t))& (t<np.max(t_u_t))], Me[(t > np.min(t_u_t))& (t<np.max(t_u_t))])
            M_u_t[j, :], M_e_t[j, :] = fit_linear(t_u_t)
        if (self.name=='SN2008D') | (self.name=='SN2006el'):
            M_u[1, :] = M_u[2, :] - 1.2444 * (M_u[2, :] - M_u[1, :]) - 0.3820 #i->I
            self.band_bc[1]='I'
        lbol,le=lyman_BC(np.reshape(M_u[0, :],(1,np.shape(M_u[0, :])[0])), np.reshape(M_u[1, :],(1,np.shape(M_u[1, :])[0])),self.band_bc[0],self.band_bc[1],np.reshape(M_e[0, :],(1,np.shape(M_e[0, :])[0])), np.reshape(M_e[1, :],(1,np.shape(M_e[1, :])[0])))
        lbol_t,le_t=lyman_BC(np.reshape(M_u_t[0, :],(1,np.shape(M_u_t[0, :])[0])), np.reshape(M_u_t[1, :],(1,np.shape(M_u_t[1, :])[0])),self.band_bc[0],self.band_bc[1],np.reshape(M_e_t[0, :],(1,np.shape(M_u_t[0, :])[0])), np.reshape(M_e_t[1, :],(1,np.shape(M_u_t[1, :])[0])))
        lbol = np.reshape(lbol, t_u.shape)
        lbol_t = np.reshape(lbol_t, t_u_t.shape)
        le = np.reshape(le, t_u.shape)
        le_t = np.reshape(le_t, t_u_t.shape)
        tpeak_e=np.mean(tpeak_err)
        if self.verbose:
            plt.figure(1)
            plt.errorbar(mag[:, 0][mag[:, 5] == self.band_bc[0]].astype(float) - self.t0[0],
                        mag[:, 3][mag[:, 5] == self.band_bc[0]].astype(float), yerr=mag[:, 4][mag[:, 5] == self.band_bc[0]].astype(float),fmt='o', lw=0.2, label=self.band_bc[0],color='r', ls='none',markersize=2.5)
            plt.errorbar(mag[:, 0][mag[:, 5] == self.band_bc[1]].astype(float) - self.t0[0],
                        mag[:, 3][mag[:, 5] == self.band_bc[1]].astype(float),yerr=mag[:, 4][mag[:, 5] == self.band_bc[1]].astype(float),fmt='o', lw=0.2, label=self.band_bc[1],color='b', ls='none',markersize=2.5)
            plt.plot(t_u, M_u[0, :], lw=1,color='maroon', label='fit '+self.band_bc[0])
            #plt.fill_between(t_u,M_u[0, :]-M_e[0, :],M_u[0, :]+M_e[0, :],alpha=0.3,facecolor='r',edgecolor='tomato')
            plt.plot(t_u, M_u[1, :], lw=1,color='navy', label='fit '+self.band_bc[1])
            #plt.fill_between(t_u, M_u[1, :] - M_e[1, :], M_u[1, :] + M_e[1, :], alpha=0.3, facecolor='b',edgecolor='lightblue')
            plt.plot(t_u_t, M_u_t[0, :], lw=1,color='maroon')
            plt.plot(t_u_t, M_u_t[1, :], lw=1,color='navy')
            plt.xlabel(r'Time (days)', fontname='Sans')
            plt.ylabel(r'$\rm M_{\rm bol}(\rm mag)$')
            plt.xlim([-5,150])
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
            plt.legend()
            plt.title(self.name)
            plt.gca().invert_yaxis()
            plt.savefig('./Plots/'+self.name+'_mag.pdf')
            plt.figure(2)
            plt.plot(t_u_t[(t_u_t >= 60) & (t_u_t < 120)],lbol_t[(t_u_t >= 60) & (t_u_t < 120)], color="orange")
            plt.fill_between(t_u_t[(t_u_t >= 60) & (t_u_t < 120)],lbol_t[(t_u_t >= 60) & (t_u_t < 120)]-le_t[(t_u_t >= 60) & (t_u_t < 120)],
            lbol_t[(t_u_t >= 60) & (t_u_t < 120)] + le_t[(t_u_t >= 60) & (t_u_t < 120)],alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
            plt.plot(t_u,lbol, color="orange")
            plt.fill_between(t_u,lbol - le, lbol + le, alpha=0.5,  edgecolor='#CC4F1B', facecolor='#FF9848')
            plt.xlabel(r'Time (days)',fontname='Sans')
            plt.ylabel(r'$\rm L_{\rm bol}(\rm erg \ s^{-1})$')
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
            plt.title(self.name)
            plt.savefig('./Plots/'+self.name+'_Lbol.pdf')
            plt.show()


        self.peakt[0] = t_u[np.argmax(lbol)]
        self.peakL[0]= np.log10(np.max(lbol))
        self.peakL[1]=np.log10(le[np.argmax(lbol)])
        if (self.name in SN_verified_t0):
            n=1
            s = [0]
        else:
            n=100
            s = np.random.normal(0, self.t0[1], n)
        p=np.zeros(shape=(2,n))
        p_err = np.zeros(shape=(2, n))
        for k,offset in enumerate(s):
            x, x_err = fit_bootstrap([0.1, 3], t_u_t[(t_u_t >= 60) & (t_u_t < 120)]+offset, lbol_t[(t_u_t >= 60) & (t_u_t < 120)],
                                     np.clip(le_t[(t_u_t >= 60) & (t_u_t < 120)], a_min=0.0001e42, a_max=1e51), valentiErr,
                                     errfunc=True, perturb=True, n=1000, nproc=4)
            p[:,k]=x
            p_err[:,k]=x_err
        self.tail_ni_mass=[np.mean(p[0,:]), np.sqrt(np.mean(p_err[0,:])**2+np.std(p[0,:])**2)]
        self.tail_meje=[np.mean(p[1,:]), np.sqrt(np.mean(p_err[1,:])**2+np.std(p[1,:])**2)]
        self.peakt[1]=np.sqrt(self.t0[1]**2+tpeak_e**2)
        n=100
        s = np.random.normal(0, self.peakt[1], n)
        p=np.zeros(n)
        p_err = np.zeros(n)
        for k,offset in enumerate(s):
            p[ k], p_err[ k]=valenti_ni56((self.peakt[0]+offset),10**self.peakL[0],0,0,10**self.peakL[1])
        self.arnett_ni_mass=[np.mean(p), np.sqrt(np.mean(p_err)**2+np.std(p)**2)]
        if 0:
            import matplotlib.gridspec as gridspec
            fig = plt.figure(3,figsize=(6,8))
            gs1 = gridspec.GridSpec(2, 1)
            gs1.update(wspace=0.0, hspace=0.0)
            ax=plt.subplot(gs1[0])
            plt.errorbar(mag[:, 0][mag[:, 5] == self.band_bc[0]].astype(float) - self.t0[0],
                         mag[:, 3][mag[:, 5] == self.band_bc[0]].astype(float),
                         yerr=mag[:, 4][mag[:, 5] == self.band_bc[0]].astype(float), fmt='o', lw=1,
                         label=self.band_bc[0], color='limegreen', ls='none', markersize=7, mec='k',mew=0.8)
            plt.errorbar(mag[:, 0][mag[:, 5] == self.band_bc[1]].astype(float) - self.t0[0],
                         mag[:, 3][mag[:, 5] == self.band_bc[1]].astype(float),
                         yerr=mag[:, 4][mag[:, 5] == self.band_bc[1]].astype(float), fmt='o', lw=1,
                         label=self.band_bc[1], color='r', ls='none', markersize=7, mec='k',mew=0.8)
            plt.plot(t_u, M_u[0, :], lw=2, color='darkgreen', label=self.band_bc[0]+' fit')
            # plt.fill_between(t_u,M_u[0, :]-M_e[0, :],M_u[0, :]+M_e[0, :],alpha=0.3,facecolor='r',edgecolor='tomato')
            plt.plot(t_u, M_u[1, :], lw=2, color='maroon', label=self.band_bc[1]+' fit')
            plt.plot(t_u_t, M_u_t[0, :], lw=2,color='darkgreen')
            plt.plot(t_u_t, M_u_t[1, :], lw=2,color='maroon')
            #plt.xlabel(r'MJD $-$ \ {}'.format(np.round(self.t0[0],2)), fontname='Sans')
            plt.ylabel(r'$\rm M_{\rm abs} \ (\rm mag)$')
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))

            plt.gca().set_xlim([-5, 140])
            plt.gca().set_ylim([-18, -13.3])
            ax.set_xticklabels([])
            plt.legend()
            plt.title(self.name, position=(0.5, 0.9))
            plt.gca().invert_yaxis()
            plt.subplots_adjust(wspace=0.0, hspace=0.0)
            plt.subplot(gs1[1])
            plt.plot(t_u_t[(t_u_t >= 60) & (t_u_t < 120)],lbol_t[(t_u_t >= 60) & (t_u_t < 120)]/1e42, color="olive", lw=2)
            plt.fill_between(t_u_t[(t_u_t >= 60) & (t_u_t < 120)],lbol_t[(t_u_t >= 60) & (t_u_t < 120)]/1e42-le_t[(t_u_t >= 60) & (t_u_t < 120)]/1e42,
            lbol_t[(t_u_t >= 60) & (t_u_t < 120)] /1e42+ le_t[(t_u_t >= 60) & (t_u_t < 120)]/1e42,alpha=0.5, edgecolor="olive", facecolor='palegoldenrod',lw=2)
            plt.plot(t_u,lbol/1e42, color="olive", lw=2)
            plt.fill_between(t_u,lbol /1e42- le/1e42, lbol/1e42 + le/1e42, alpha=0.5,  edgecolor="olive", facecolor='palegoldenrod',lw=2)
            plt.plot(np.arange(-5,self.peakt[0],0.1), np.repeat(10**self.peakL[0]/1e42,np.arange(-5,self.peakt[0],0.1).shape), color="orange", lw=2)
            plt.fill_between(np.arange(-5,self.peakt[0],0.1),np.repeat(10**self.peakL[0],np.arange(-5,self.peakt[0],0.1).shape)/1e42-np.repeat(10**self.peakL[1],np.arange(-5,self.peakt[0],0.1).shape)/1e42, np.repeat(10**self.peakL[0],np.arange(-5,self.peakt[0],0.1).shape)/1e42+np.repeat(10**self.peakL[1]/1e42,np.arange(-5,self.peakt[0],0.1).shape), alpha=0.1,  edgecolor="orange", facecolor='yellow',lw=2)
            #plt.plot(np.arange(-5,self.peakt[0],0.1), np.repeat(10**self.peakL[0]/1e42,np.arange(-5,self.peakt[0],0.1).shape), color="orange", lw=2)
            plt.plot(np.repeat(self.peakt[0],np.arange(0,10**self.peakL[0]/1e42,0.01).shape),np.arange(0,10**self.peakL[0]/1e42,0.01),color='orange',lw=2)
            plt.fill_betweenx(np.arange(0,10**self.peakL[0]/1e42,0.1),np.repeat(self.peakt[0],np.arange(0,10**self.peakL[0]/1e42,0.1).shape)-np.repeat(self.peakt[1],np.arange(0,10**self.peakL[0]/1e42,0.1).shape), np.repeat(self.peakt[0],np.arange(0,10**self.peakL[0]/1e42,0.1).shape)+np.repeat(self.peakt[1],np.arange(0,10**self.peakL[0]/1e42,0.1).shape), alpha=0.1,  edgecolor="orange", facecolor='yellow',lw=2)
            plt.plot(np.repeat(0,np.arange(0,2.4,0.01).shape),np.arange(0,2.4,0.01), color="mediumturquoise", lw=2)
            plt.fill_betweenx(np.arange(0,2.4,0.01),np.repeat(0,np.arange(0,2.4,0.01).shape)-np.repeat(self.t0[1],np.arange(0,2.4,0.01).shape),np.repeat(0,np.arange(0,2.4,0.01).shape)+np.repeat(self.t0[1],np.arange(0,2.4,0.01).shape), alpha=0.4,  edgecolor="mediumturquoise", facecolor='lightcyan',lw=2.5)

            plt.xlabel(r'MJD $-$ \ {}'.format(np.round(self.t0[0],2)), fontname='Sans')
            plt.ylabel(r'$\rm L_{\rm bol} \ (\rm 10^{42} \ erg \ s^{-1})$')
            plt.annotate(r'L$_p$',
                         xy=(5, 2),
                         xycoords='data',
                         textcoords='offset points',fontsize=18)
            plt.gca().set_xlim([-5, 140])
            plt.annotate(r't$_p$',
                         xy=(17, 0.8),
                         xycoords='data',
                         textcoords='offset points',fontsize=18)
            plt.annotate(r't$_0$',
                         xy=(0.7, 0.8),
                         xycoords='data',
                         textcoords='offset points',fontsize=18)

            plt.gca().set_xlim([-5, 140])
            plt.gca().set_ylim([0, 2.4])
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))
            plt.subplots_adjust(wspace=0.0, hspace=0.0)
            #plt.tight_layout()
            plt.savefig('./Plots/'+self.name+'_lbolphot.pdf')
            plt.show()
    def khatami_model(self):
        from SNAP5.Analysis.LCFitting import fit_bootstrap
        beta0 = 0.8
        betafit = fit_bootstrap([beta0], [self.peakt[0], self.tail_ni_mass[0]],[10**self.peakL[0]],[10**self.peakL[1]], khatami_err,bounds=([0.0], [np.inf]),
                                 errfunc=True, perturb=True, n=1000, nproc=4)
        self.beta_req= [betafit[0][0],betafit[1][0]]
    def khatami_model_litbeta(self):
        self.ni_khatami=nickel_mass_khatami(sn.peakt[0], 10 ** sn.peakL[0], 10 ** sn.peakL[1], sug_beta[self.sn_type])
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

# sn = supernova('SN2008D')
# sn.fill(info)
# index = np.where(info[:, 0] == sn.name)[0][0]
# #sn.stritzinger_color(float(info[index, 6]))
# sn.tailNickel()
# print "Tail Ni", sn.tail_ni_mass, "Tail Mej/E", sn.tail_meje, "Tail Arnett", sn.arnett_ni_mass, "Lp", sn.peakL, "tp", sn.peakt
# sn.khatami_model()
# sn.khatami_model_litbeta()
# sn.prentice_data()

# print "beta", sn.beta_req
# print "t0", sn.t0
# ni=np.zeros((2,7))
# arnett=np.zeros((2,7))
# for i in np.arange(0,7):
#     print i
#     sn.t0[0]=sn.t0[0]-i
#     sn.tailNickel()
#     sn.band_bc[1]="r'"
#     ni[0,i],ni[1,i]=sn.tail_ni_mass
#     arnett[0, i], arnett[1, i] =sn.arnett_ni_mass
# plt.errorbar(np.arange(0,7),ni[0,:],yerr=ni[1,:],fmt='d',linestyle='-',label='Tail Ni')
# plt.errorbar(np.arange(0,7),arnett[0,:],yerr=arnett[1,:],fmt='s',linestyle='-',label='Arnett Ni')
# plt.legend()
# plt.xlabel('dt [days]')
# plt.ylabel('M_ni')
# plt.show()
#sn.fit_t0()


# print "host E(B-V)", sn.ebv_host
# sn.tailNickel()
# print "Tail Ni", sn.tail_ni_mass, "Tail Mej/E", sn.tail_meje, "Tail Arnett", sn.arnett_ni_mass, "Lp", sn.peakL, "tp", sn.peakt
# sn.khatami_model()
# sn.prentice_data()
# print "beta", sn.beta_req
# counter=0
# print  np.arange(0.05,1.5,0.05).shape
# ni=np.zeros((29,1))
# for i,beta in enumerate(np.arange(0.05,1.5,0.05)):
#     print beta
#     ni[counter],x=nickel_mass_khatami(11.4,10**42.6,10**41,beta)
#     counter+=1
#
# plt.plot(np.arange(0.05,1.5,0.05),ni)
# plt.show()

for i,sn_name in enumerate(['SN2003jd']):
    if  ((sn_name =="SN2006jc")) :#| (i<18) | (i>27):
         continue
    print sn_name
    sn = supernova(sn_name)
    sn.fill(info)
    index = np.where(info[:, 0] == sn.name)[0][0]
    #sn.fit_t0()
    #sn.stritzinger_color(float(info[index, 6]))
    print "E(B-V)",sn.ebv_host[0]+sn.ebv_gal[0], np.sqrt(sn.ebv_host[1]**2+sn.ebv_gal[1]**2),"host",sn.ebv_host, "t0", sn.t0
    sn.tailNickel()
    print "Tail Ni",sn.tail_ni_mass, "Tail Mej/E",sn.tail_meje, "Arnett",sn.arnett_ni_mass,"Lp", sn.peakL, "tp",sn.peakt
    sn.khatami_model()
    sn.khatami_model_litbeta()
    print sn.ni_khatami
    print sn.beta_req
    sn.prentice_data()
    with open('Data/'+sn.name+'_results.csv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=sorted(vars(sn)))
        w.writeheader()
        w.writerow({k: getattr(sn, k) for k in vars(sn)})


# with open('Data/out_new1.csv','w') as f:
#     w = csv.DictWriter(f,fieldnames=sorted(vars(sn)))
#     w.writeheader()
#     w.writerow({k: getattr(sn, k) for k in vars(sn)})
#     for i,sn_name in enumerate(info[1:,0]):
#         if  (sn_name =="SN2006jc")  | (sn_name=='SN2008D') | (i>17):
#              continue
#         print sn_name
#         sn = supernova(sn_name)
#         sn.fill(info)
#         # index = np.where(info[:, 0] == sn.name)[0][0]
#         # sn.stritzinger_color(float(info[index, 6]))
#         print "host E(B-V)",sn.ebv_host
#         sn.tailNickel()
#         print "Tail Ni",sn.tail_ni_mass, "Tail Mej/E",sn.tail_meje, "Tail Arnett",sn.arnett_ni_mass,"Lp", sn.peakL, "tp",sn.peakt
#         sn.khatami_model()
#         sn.khatami_model_litbeta()
#         sn.prentice_data()
#         print "beta",sn.beta_req
#         print "t0",sn.t0
#         w.writerow({k: getattr(sn, k) for k in vars(sn)})

