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



coef = {'B': 3.626, 'V': 2.742, 'I': 1.505, "i'": 1.698, "r'": 2.285, "R": 2.169}

class supernova:
    def __init__(self, name,  host=None,sn_type=None,distance=[None, None], band_extinction=None, band_bc=None , source=None, ebv_gal=[None, None],t0=[None, None] ,ebv_host=[None, None]):
        self.name=name
        self.type=sn_type
        self.host=host
        self.distance=distance
        self.band_extinction=band_extinction
        self.band_bc=band_bc
        self.source=source
        self.t0=t0
        self.ebv_gal=ebv_gal
        self.ebv_host=ebv_host
        self.verbose=1
    def fill(self,info):
        location = './Data/SN_json/'
        if os.path.isfile(location + self.name + '.json') == False:
            url = 'https://sne.space/astrocats/astrocats/supernovae/output/json/' + self.name  + '.json'
            urllib.urlretrieve(url, location + self.name  + '.json')
        with open(location + self.name  + '.json') as data_file:
            data = json.load(data_file)
        index=np.where(info[:,0]==self.name)[0][0]
        self.sn_type = info[index , 2].split(' ')[0]
        self.host = info[index, 1]
        self.distance[0] = float(info[index , 10].split('(')[0])
        if (len(info[index , 10].split('(')) > 1):
            self.distance[1] = float(info[index, 10].split('(')[1].replace(")", ""))
        self.band_extinction = info[index, 4].split(" ")[0].split("-")
        self.band_bc = info[index , 17].split(',')[0].split("-")
        print self.band_bc
        self.source = info[index, 5].split(',')
        self.ebv_gal[0] = float(info[index , 3].split('(')[0]) / 3.1
        if (len(info[index , 2].split('(')) > 1):
            self.ebv_gal[1] = float(info[index , 2].split('(')[1].replace(")", ""))
        if self.verbose:
            print self.name,self.host,self.sn_type,self.distance,self.band_extinction,self.band_bc,self.source,self.ebv_gal

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
            plt.scatter(mag[:, 0][mag[:, 3] == 'V'].astype(float)-t0, mag[:, 1][mag[:, 3] == 'V'].astype(float), s=5)
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
        self.ebv_host[1]= np.mean(np.sqrt(color_obs_e[(t_u < 10) & (t_u > 5)]**2+ color_in_e[(t_u < 10) & (t_u > 5)]**2))



plt.close("all")
url='https://docs.google.com/spreadsheets/d/e/2PACX-1vS3r49T74bVSHHYNJrkSwcc7O7-nXdmd9ERmfdO8lNKAocdrhjT4H7znvILr8nN8BERXZpM2_cx3ge0/pub?gid=0&single=true&output=csv'
with requests.Session() as s:
    download = s.get(url)
    decoded_content = download.content.decode('utf-8')
    cr = csv.reader(decoded_content.splitlines(), delimiter=',')
    my_list = list(cr)
info=np.array(my_list)

sn=supernova('SN2013ge')
sn.fill(info)
index=np.where(info[:,0]==sn.name)[0][0]
print float(info[index,6])
sn.stritzinger_color(float(info[index,6]))
print sn.ebv_host
# for i,sn_name in enumerate(info[1:,0]):
#     distance=[0,0]
#     ebv = [0, 0]
#     if  (sn_name =="SN2006jc") |(sn_name=='SN2002ap'):
#          continue
#     if  (sn_name !="SN2013ge"):
#          continue
#     #sn=supernova(sn_name,host,sn_type,distance,band_extinction,band_bc,source,ebv)
#     sn=supernova(sn_name)
#     sn.fill(info)
#     sn.stritzinger_color(float(info[i+1,6]))
#     print sn.ebv_host