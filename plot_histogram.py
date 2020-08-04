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
from astropy.constants import M_sun
sys.path.insert(0, '/home/afsari/')
from SNAP5.Analysis import *
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib
import itertools
from matplotlib.ticker import AutoMinorLocator
import pandas as pd


def sample_cdf(ni,ni_e,ax,color,n=1000):
    data=np.multiply(np.repeat(np.expand_dims(ni_e, axis=1),n,axis=1), np.random.randn(np.shape(ni)[0],n)) +np.repeat(np.expand_dims(ni, axis=1),n,axis=1)
    print data.shape
    for c in range(0,n):
        h, edges = np.histogram(data[:,c], density=True, bins=20000)
        h = np.cumsum(h) / np.cumsum(h).max()
        X = edges.repeat(2)[:-1]
        x_s = X
        y = np.zeros_like(X)
        y[1:] = h.repeat(2)
        ax.plot(X, y, color=color, lw=0.1,alpha=0.07)


df=pd.read_csv('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/khatamiMni.csv')
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
matplotlib.rcParams.update({'font.size': 18})
plt.rc('text', usetex=True)
marker = itertools.cycle(('8', 'p','>','^','v','<','s', 'o','h','<','*','P','X','+','H','1','2','3','4'))
colors = itertools.cycle(("black","salmon","yellow","green","m","sienna","gray","blue",
                          "darkkhaki","peru","gold","deepskyblue","olive"))

#data=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/Results_linear.txt')
sn_names=[]
sn_type=[]
arnett_ni=[]
ni=[]
types=[]
arnett_ni_e=[]
ni_e=[]
line=0

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


from scipy import stats
class deterministic_gen(stats.rv_continuous):
    def _cdf(self, x):
        dat=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/CDF_typeII.csv',delimiter=',')
        dat=dat[np.argsort(dat[:,0]),:]
        res=np.zeros((x.shape[0]))
        for i,value in enumerate(list(x)):
            res[i]=dat[find_nearest(dat[:,0],value),1]
        return res

deterministic = deterministic_gen(name="type_ii")

plt.figure(3)
plt.plot(np.arange(0, 0.7, 0.001), deterministic.cdf(np.arange(0, 0.7, 0.001)))



tpeaks=[]
lpeaks=[]

dict={'Ic':'red', 'Ib':'blue','Ic BL':'orange', 'IcGRB':'orange','IIb':'green','Ib pec':'blue', 'Ibn':'blue'}
# Load all rows from the file into a variable called rows
with open("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda_26dec.csv", "r") as f_input:
    csv_reader = csv.reader(f_input, delimiter=",")
    rows = list(csv_reader)
    for row in rows:
        if line==0:
            index_name=row.index('name')
            index_lpeak=row.index('peakL')
            index_tpeak=row.index('peakt')
            index_ni_khatami=row.index('ni_khatami')
            index_arnett_ni_mass=row.index('arnett_ni_mass')
            index_beta_req=row.index('beta_req')
            index_prentice_ni= row.index('prentice_ni')
            index_tail_ni_mass=row.index('tail_ni_mass')
            index_tail_meje=row.index('tail_meje')
            index_prentice_peakt=row.index('prentice_peakt')
            index_prentice_peakL=row.index('prentice_peakL')
            index_sn_type=row.index('sn_type_full')
            line=line+1
        else:
            # mark = marker.next()
            # col = colors.next()
            sn_names.append(row[index_name])
            sn_type.append(row[index_sn_type])
            tpeaks.append(row[index_tpeak])
            lpeaks.append(row[index_lpeak])
            line=line+1
            arnett_ni.append(float(row[index_arnett_ni_mass].split(";")[0]))
            ni.append(float(row[index_tail_ni_mass].split(";")[0]))
            arnett_ni_e.append(float(row[index_arnett_ni_mass].split(";")[1]))
            ni_e.append(float(row[index_tail_ni_mass].split(";")[1]))
            types.append(row[index_sn_type].split(";")[0])

# import matplotlib.gridspec as gridspec
# f = plt.figure(1,figsize=(5,13))
# gs1 = gridspec.GridSpec(4, 1)
# gs1.update(wspace=0.0, hspace=0.0)
# ax0 = plt.subplot(gs1[0])
# types1=np.array(types)
# arnett_ni=np.array(arnett_ni)
# arnett_ni_e=np.array(arnett_ni_e)
# tpeaks=np.array(tpeaks)
# lpeaks=np.array(lpeaks)
# ni=np.array(ni)
# ni_e=np.array(ni_e)
# khatami_ni=np.array(df['Mni'].tolist())
# khatami_ni_e=np.array(df['Mni_e'].tolist())
# khatami_type=np.array(df['type'].tolist())
# h, edges = np.histogram(arnett_ni, density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax0.plot(X,y,color='olive',lw=2,label='Arnett')
# y=np.arange(0,1.1,0.01)
# X=np.repeat(np.mean(ni),y.shape)
# ax0.plot(X,y,color='#1f77b4',lw=1.2)
# sample_cdf(arnett_ni,arnett_ni_e,ax0,color='olive',n=500)
h, edges = np.histogram(ni, density=True, bins=20000)
h = np.cumsum(h)/np.cumsum(h).max()

X = edges.repeat(2)[:-1]
x_ss=X
y = np.zeros_like(X)
y[1:] = h.repeat(2)
y_ss=y
# ax0.plot(X,y,color='#1f77b4',lw=2,label='Tail')
# sample_cdf(ni,ni_e,ax0,color='#1f77b4',n=500)
# y=np.arange(0,1.1,0.01)
# X=np.repeat(np.mean(arnett_ni),y.shape)
# print (np.mean(arnett_ni), np.mean(ni), np.std(arnett_ni), np.std(ni))
# ax0.plot(X,y,color='olive',lw=1.2)
#
# h, edges = np.histogram(khatami_ni, density=True, bins=20000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# x_s=X
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# y_s=y
# ax0.plot(X,y,color='salmon',lw=2,label='KK19')
# y=np.arange(0,1.1,0.01)
# X=np.repeat(np.mean(khatami_ni),y.shape)
# ax0.plot(X,y,color='salmon',lw=1.2)
# sample_cdf(khatami_ni,khatami_ni_e,ax0,color='salmon',n=500)
#
# plt.legend(frameon=False)
# ax0.set_xlim([0, 0.7])
# ax0.set_ylim([0.01, 1.01])
# ax0.set_xticklabels([])
# ax0.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax0.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax0.xaxis.set_tick_params(width=1.5)
# ax0.yaxis.set_tick_params(width=1.5)
# ax0.yaxis.set_ticks_position('both')
# ax0.tick_params(direction = 'in',which ='both')
# ax0.set_ylabel(r'Overall CDF',fontsize=18)
#
# ax1 = plt.subplot(gs1[1])
# h, edges = np.histogram(ni[types1=='Ic BL'], density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax1.plot(X,y,color='orange',lw=3,label='Ic-BL')
# sample_cdf(ni[types1=='Ic BL'],ni_e[types1=='Ic BL'],ax1,color='orange',n=500)
# h, edges = np.histogram(ni[types1=='Ic'], density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax1.plot(X,y,color='red',lw=3,label='Ic')
# sample_cdf(ni[types1=='Ic'],ni_e[types1=='Ic'],ax1,color='red',n=500)
#
# h, edges = np.histogram(ni[types1=='Ib'], density=True, bins=10000 )
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax1.plot(X,y,color='blue',lw=3,label='Ib')
# sample_cdf(ni[types1=='Ib'],ni_e[types1=='Ib'],ax1,color='blue',n=500)
#
# h, edges = np.histogram(ni[types1=='IIb'], density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
#
# ax1.plot(X,y,color='green',lw=3,label='IIb')
# sample_cdf(ni[types1=='IIb'],ni_e[types1=='IIb'],ax1,color='green',n=500)
#
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(ni[types1=='Ic BL']),y.shape)
# ax1.plot(X,y,color='orange',lw=1.5,ls=':',alpha=1)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(ni[types1=='IIb']),y.shape)
# ax1.plot(X,y,color='green',lw=1.5,ls=':',alpha=0.9)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(ni[types1=='Ic']),y.shape)
# ax1.plot(X,y,color='red',lw=1.5,ls=':',alpha=0.9)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(ni[types1=='Ib']),y.shape)
# ax1.plot(X,y,color='blue',lw=1.5,ls=':',alpha=0.9)
#
# print ("tail:", np.mean(ni[types1=='Ic BL']), np.mean(ni[types1=='IIb']), np.mean(ni[types1=='Ic']), np.mean(ni[types1=='Ib']))
#
# ax1.set_xticklabels([])
# ax1.set_xlim([0, 0.7])
# ax1.set_ylim([0.01, 1])
# ax1.set_ylabel(r'Tail CDF',fontsize=18)
# ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax1.xaxis.set_tick_params(width=1.5)
# ax1.yaxis.set_tick_params(width=1.5)
# ax1.yaxis.set_ticks_position('both')
# ax1.tick_params(direction = 'in',which ='both')
# plt.gca().legend(loc='lower right',frameon=False,fontsize=16)
#
#
# ax2 = plt.subplot(gs1[2])
#
#
# h, edges = np.histogram(arnett_ni[types1=='Ic BL'], density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax2.plot(X,y,color='orange',lw=3,label='Ic-BL',ls='-')
# sample_cdf(arnett_ni[types1=='Ic BL'],arnett_ni_e[types1=='Ic BL'],ax2,color='orange',n=500)
# h, edges = np.histogram(arnett_ni[types1=='Ic'], density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax2.plot(X,y,color='red',lw=3,label='Ic',ls='-')
# sample_cdf(arnett_ni[types1=='Ic'],arnett_ni_e[types1=='Ic'],ax2,color='red',n=500)
#
#
#
# h, edges = np.histogram(arnett_ni[types1=='Ib'], density=True, bins=10000 )
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax2.plot(X,y,color='blue',lw=3,label='Ib',ls='-')
# sample_cdf(arnett_ni[types1=='Ib'],arnett_ni_e[types1=='Ib'],ax2,color='blue',n=500)
#
#
# from scipy import stats
# h = np.histogram(ni[(types1=='Ib')| (types1=='Ic')], density=True, bins=10000)
# hist_dist = stats.rv_histogram(h)
# d ,p = stats.kstest(list(ni[types1=='Ic BL']),hist_dist.cdf)
# print d, 100-p*100
#
# h = np.histogram(ni[(types1=='Ib')], density=True, bins=10000)
# hist_dist = stats.rv_histogram(h)
# d ,p = stats.kstest(list(ni[types1=='Ic']),hist_dist.cdf)
# print d, p*100
#
#
# h = np.histogram(ni[(types1=='IIb')], density=True, bins=10000)
# hist_dist = stats.rv_histogram(h)
# d ,p = stats.kstest(list(ni[ (types1=='Ib') | (types1=='Ic')]),hist_dist.cdf)
# print d, p*100
#
#
#
# h = np.histogram(ni[(types1=='IIb')], density=True, bins=10000)
# hist_dist = stats.rv_histogram(h)
# d ,p = stats.kstest(list(ni[ (types1=='Ib')]),hist_dist.cdf)
# print d, p*100
#
# h, edges = np.histogram(arnett_ni[types1=='IIb'], density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
#
#
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax2.plot(X,y,color='green',lw=3,label='IIb',ls='-')
# sample_cdf(arnett_ni[types1=='IIb'],arnett_ni_e[types1=='IIb'],ax2,color='green',n=500)
#
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(arnett_ni[types1=='Ic BL']),y.shape)
# ax2.plot(X,y,color='orange',lw=1.5,ls=':',alpha=1)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(arnett_ni[types1=='IIb']),y.shape)
# ax2.plot(X,y,color='green',lw=1.5,ls=':',alpha=0.9)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(arnett_ni[types1=='Ic']),y.shape)
# ax2.plot(X,y,color='red',lw=1.5,ls=':',alpha=0.9)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(arnett_ni[types1=='Ib']),y.shape)
# ax2.plot(X,y,color='blue',lw=1.5,ls=':',alpha=0.9)
#
# ax2.set_xlim([0.001, 0.7])
# ax2.set_ylim([0.01, 1])
# ax2.set_xlabel(r'$M_{\rm Ni} (M_\odot$)',fontsize=18)
# ax2.set_ylabel(r'Arnett CDF',fontsize=18)
# ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax2.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax2.xaxis.set_tick_params(width=1.5)
# ax2.yaxis.set_tick_params(width=1.5)
# ax2.yaxis.set_ticks_position('both')
# ax2.set_xticklabels([])
# ax2.tick_params(direction = 'in',which ='both')
# plt.gca().legend(frameon=False,fontsize=16)
#
# print ("IIb & %10.2f & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f " )%(np.round(np.mean(ni[types1=='IIb']),2), np.round(np.median(ni[types1=='IIb']),2), np.std(ni[types1=='IIb']),np.mean(arnett_ni[types1=='IIb']), np.median(arnett_ni[types1=='IIb']), np.std(arnett_ni[types1=='IIb']))
# print ("Ib & %10.2f & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f " )%(np.round(np.mean(ni[types1=='Ib']),2), np.round(np.median(ni[types1=='Ib']),2), np.std(ni[types1=='Ib']),np.mean(arnett_ni[types1=='Ib']), np.median(arnett_ni[types1=='Ib']), np.std(arnett_ni[types1=='Ib']))
# print ("Ic & %10.2f & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f " )%(np.round(np.mean(ni[types1=='Ic']),2), np.round(np.median(ni[types1=='Ic']),2), np.std(ni[types1=='Ic']),np.mean(arnett_ni[types1=='Ic']), np.median(arnett_ni[types1=='Ic']), np.std(arnett_ni[types1=='Ic']))
# print ("Ic-BL & %10.2f & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f " )%(np.round(np.mean(ni[types1=='Ic BL']),2), np.round(np.median(ni[types1=='Ic BL']),2), np.std(ni[types1=='Ic BL']),np.mean(arnett_ni[types1=='Ic BL']), np.median(arnett_ni[types1=='Ic BL']), np.std(arnett_ni[types1=='Ic BL']))
# print ("all & %10.2f & %5.2f & %5.2f & %5.2f & %5.2f & %5.2f " )%(np.round(np.mean(ni),2), np.round(np.median(ni),2), np.std(ni),np.mean(arnett_ni), np.median(arnett_ni), np.std(arnett_ni))
#
#
# ax3 = plt.subplot(gs1[3])
#
# print "khatami",khatami_type=='Ic BL',khatami_type
# h, edges = np.histogram(khatami_ni[khatami_type=='Ic BL'], density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax3.plot(X,y,color='orange',lw=3,label='Ic-BL',ls='-')
# sample_cdf(khatami_ni[types1=='Ic BL'],khatami_ni_e[types1=='Ic BL'],ax3,color='orange',n=500)
#
# h, edges = np.histogram(khatami_ni[khatami_type=='Ic'], density=True, bins=10000)
# h = np.cumsum(h)/np.cumsum(h).max()
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax3.plot(X,y,color='red',lw=3,label='Ic',ls='-')
# sample_cdf(khatami_ni[types1=='Ic'],khatami_ni_e[types1=='Ic'],ax3,color='red',n=500)
#
#
#
# h, edges = np.histogram(khatami_ni[khatami_type=='Ib'], density=True, bins=10000 )
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax3.plot(X,y,color='blue',lw=3,label='Ib',ls='-')
# sample_cdf(khatami_ni[types1=='Ib'],khatami_ni_e[types1=='Ib'],ax3,color='blue',n=500)
#
#
# h, edges = np.histogram(khatami_ni[khatami_type=='IIb'], density=True, bins=10000 )
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# ax3.plot(X,y,color='green',lw=3,label='IIb',ls='-')
# sample_cdf(khatami_ni[types1=='IIb'],khatami_ni_e[types1=='IIb'],ax3,color='green',n=500)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(khatami_ni[khatami_type=='Ic BL']),y.shape)
# ax3.plot(X,y,color='orange',lw=1.5,ls=':',alpha=1)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(khatami_ni[khatami_type=='IIb']),y.shape)
# ax3.plot(X,y,color='green',lw=1.5,ls=':',alpha=0.9)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(khatami_ni[khatami_type=='Ic']),y.shape)
# ax3.plot(X,y,color='red',lw=1.5,ls=':',alpha=0.9)
#
# y=np.arange(0,1.2,0.01)
# X=np.repeat(np.mean(khatami_ni[khatami_type=='Ib']),y.shape)
# ax3.plot(X,y,color='blue',lw=1.5,ls=':',alpha=0.9)
#
#
# ax3.set_xlim([0.001, 0.7])
# ax3.set_ylim([0, 1])
# ax3.set_xlabel(r'$M_{\rm Ni} (M_\odot$)',fontsize=18)
# ax3.set_ylabel(r'KK19 CDF',fontsize=18)
# ax3.xaxis.set_minor_locator(AutoMinorLocator(2))
# ax3.yaxis.set_minor_locator(AutoMinorLocator(2))
# ax3.xaxis.set_tick_params(width=1.5)
# ax3.yaxis.set_tick_params(width=1.5)
# ax3.yaxis.set_ticks_position('both')
# ax3.tick_params(direction = 'in',which ='both')
#
# plt.gca().legend(frameon=False,fontsize=16)
#
# plt.gcf().savefig('/home/afsari/PycharmProjects/typeIbcAnalysis/Plots/cdfs.pdf', bbox_inches='tight')
#
f2 = plt.figure(2,figsize=(5.5,5))
ax3=plt.subplot(111)
d ,p = stats.kstest(ni,deterministic.cdf)
data=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/CDF_typeII.csv',delimiter=',')
# print "typeII",d, p


ax3.plot(data[:,0],data[:,1],color='green',lw=3,label='H-rich Type II (Anderson19)')
data=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/CDF_SESNe.csv',delimiter=',')
ax3.plot(data[:,0],data[:,1],color='salmon',lw=3,label='SESN (Anderson19)')
ax3.plot(x_ss,y_ss,color='orange',lw=3,label='SESN (this work, Tail)')

df=pd.read_csv('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/Meza2020.csv')
print np.array(df['Tail Nickel'].mask( df['Tail Nickel'].eq('None')).dropna()).astype(float)
h, edges = np.histogram(np.array(df['Tail Nickel'].mask( df['Tail Nickel'].eq('None')).dropna()).astype(float), density=True, bins=20000)
h = np.cumsum(h)/np.cumsum(h).max()

X = edges.repeat(2)[:-1]
x_s=X
y = np.zeros_like(X)
y[1:] = h.repeat(2)
y_s=y
#ax0.plot(X,y,color='#1f77b4',lw=2,label='Tail')

ax3.plot(x_s,y_s,color='b',lw=3,label='SESN (Meza20, Tail)',alpha=0.2)
#print np.array(df['Tail Nickel'].mask( df['Tail Nickel'].eq('None')).dropna()).astype(float)
# h, edges = np.histogram(np.array(df['Arnett Nickel']).astype(float), density=True, bins=20000)
# h = np.cumsum(h)/np.cumsum(h).max()
#
# X = edges.repeat(2)[:-1]
# x_s=X
# y = np.zeros_like(X)
# y[1:] = h.repeat(2)
# y_s=y
# #ax0.plot(X,y,color='#1f77b4',lw=2,label='Tail')
#
# ax3.plot(x_s,y_s,color='cyan',lw=3,label='SESN (Meza20, Arnett)',alpha=0.2)


ax3.set_xlim([0.001, 1])
ax3.set_ylim([0, 1])
ax3.set_xlabel(r'$M_{\rm Ni} (M_\odot$)',fontsize=18)
ax3.set_ylabel(r'CDF',fontsize=18)
ax3.xaxis.set_minor_locator(AutoMinorLocator(5))
ax3.yaxis.set_minor_locator(AutoMinorLocator(10))
ax3.xaxis.set_tick_params(width=1.5)
ax3.yaxis.set_tick_params(width=1.5)
ax3.yaxis.set_ticks_position('both')
ax3.tick_params(direction = 'in',which ='both')
plt.gca().legend(frameon=False,fontsize=14)
plt.gcf().savefig('/home/afsari/PycharmProjects/typeIbcAnalysis/Plots/cdf_compare.pdf', bbox_inches='tight')

plt.show()
