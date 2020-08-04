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

# data = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SN2003jd_lcophot.csv', delimiter=',',
#                   usecols=(0, 1, 2))
# data=data.reshape((data.shape[0], 3))
# filter = np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SN2003jd_lcophot.csv',
#                     delimiter=',', usecols=(3,), dtype='|S8').reshape((data.shape[0], 1))
# # mags = [deredMag(data[0:770, 1], 0.08, coef[filter] / 0.86) - DM, np.sqrt(
# #     data[0:770, 2] ** 2 + DMe ** 2 + (coef[filter] / 0.86) ** 2 * (self.ebv_host[1] ** 2 + self.ebv_gal[1] ** 2))]
# data = np.concatenate((np.reshape(data, (data.shape[0], 3)), np.reshape(filter, (data.shape[0], 1))), axis=1)
# labels=['B','V','I','R']
# colors=['r','b','k','g','orange']
# f=plt.figure(1)
# for i,lab in enumerate(labels):
#     plt.scatter(data[data[:,3]==lab,0].astype(np.float),data[data[:,3]==lab,1].astype(np.float),label=lab,color=colors[i])
#
# plt.gca().invert_yaxis()

tp=15
beta0=1
beta=[]
Lp=np.arange(41.8,45,0.1)
#Lp=42.5
#tp=np.logspace(0,3,100)
mni=0.05
for i,lp in enumerate(Lp):
    print i,lp
    betafit = fit_bootstrap([beta0], [tp, mni], [10**lp], [1], khatami_err, bounds=([-10], [np.inf]),
                            errfunc=True, perturb=False, n=300, nproc=8)
    beta.append(betafit[0][0])
    print betafit[0][0]
    #mni[i]=nickel_mass_khatami(tp,10**Lp,0.1,b)[0]



plt.plot(Lp,beta)
#plt.legend()
plt.show()