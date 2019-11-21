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
from valenti_ni56 import *
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
marker = itertools.cycle(( 'x','p','>','^','v','<','s', 'o','h','d','D','*','P','X','+','H'))
colors = itertools.cycle(("black","salmon","yellow","green","m","sienna","gray","blue",
                          "darkkhaki","peru","gold","deepskyblue","olive"))

#data=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/Results_linear.txt')
sn_names=[]
sn_type=[]
arnett_ni=[]
ni=[]
line=0
f = plt.figure(1,figsize=(6,6))
ax=plt.subplot(111)




dict={'Ic':'red', 'Ib':'blue','Ic BL':'orange', 'IcGRB':'orange','IIb':'green','Ib pec':'blue', 'Ibn':'blue'}
dict_marker={'Ic':'o', 'Ib':'v','Ic BL':'s', 'IIb':'d','Ib pec':'v', 'Ibn':'v'}
# Load all rows from the file into a variable called rows

line=0
plot_lines = []
orange_marker=[]
blue_marker=[]
red_marker=[]
green_marker=[]
types=[]
beta_req=[]
dict_list={'orange':orange_marker,'blue':blue_marker, 'red':red_marker, 'green':green_marker}
with open("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda.csv", "r") as f_input:
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
            mark = marker.next()
            col = colors.next()
            # while mark in dict_list[dict[row[index_sn_type]]]:
            #     mark = marker.next()
            sn_names.append(row[index_name])
            sn_type.append(row[index_sn_type])
            if (len(orange_marker)==0) & (dict[row[index_sn_type]]=='orange'):
                label = 'Ic BL'
            elif (len(green_marker)==0) & (dict[row[index_sn_type]]=='green'):
                label = 'IIb'
            elif (len(red_marker) == 0) & (dict[row[index_sn_type]] == 'red'):
                label = 'Ic'
            elif (len(blue_marker) == 0) & (dict[row[index_sn_type]] == 'blue'):
                label = 'Ib'
            else:
                label=None
            line=line+1

            arnett_ni.append(float(row[index_arnett_ni_mass].split(";")[0]))
            ni.append(float(row[index_tail_ni_mass].split(";")[0]))
            types.append(row[index_sn_type])
            beta_req.append(float(row[index_beta_req].split(";")[0]))
            ax.errorbar(float(row[index_tail_ni_mass].split(";")[0]), float(row[index_beta_req].split(";")[0]), yerr=float(row[index_beta_req].split(";")[1]),
                         xerr=float(row[index_tail_ni_mass].split(";")[1]),
                                                                             marker=dict_marker[row[index_sn_type]],
                                                                             color='deepskyblue',
                                                                             label=label, linewidth=0.4,
                                                                             elinewidth=0.5, mec='k',markersize=8)
            if dict[row[index_sn_type]]=='orange':
                orange_marker.append(mark)
            elif dict[row[index_sn_type]]=='blue':
                blue_marker.append(mark)
            elif dict[row[index_sn_type]]=='green':
                green_marker.append(mark)
            elif dict[row[index_sn_type]]=='red':
                red_marker.append(mark)



types=np.array(types)
beta_req=np.array(beta_req)
#ax.legend(loc=1,bbox_to_anchor=(1.37,1.1),
#          fancybox=True, ncol=1, fontsize =9)
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
# ax.annotate(r'Ic-BL',color=dict['Ic BL'],
#              xy=(0.1, 0.2),
#              xycoords='axes fraction', fontsize=16)
# ax.annotate(r'IIb',color=dict['IIb'],
#              xy=(0.1, 0.15),
#              xycoords='axes fraction', fontsize=16)
# ax.annotate(r'Ib',color=dict['Ib'],
#              xy=(0.1, 0.1),
#              xycoords='axes fraction', fontsize=16)
# ax.annotate(r'Ic',color=dict['Ic'],
#              xy=(0.1, 0.05),
#              xycoords='axes fraction', fontsize=16)

ax.set_ylabel(r'$\beta$',fontsize=18)
ax.set_xlabel(r'$\rm M_{Ni} \ (M_\odot$)',fontsize=18)
ax.set_xscale("log")
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction = 'in',which ='both')
ax.set_xlim(0.01, 1)
ax.set_ylim(0.0001, 2)

from SNAP5.Analysis.LCFitting import fit_bootstrap
mni=np.logspace(np.log10(0.02),np.log10(1),5, endpoint=False)
tpeaks=np.logspace(np.log10(5),np.log10(80),4, endpoint=False)
#print tpeaks
lpeak=np.zeros((mni.size,tpeaks.size))
betas=np.zeros((mni.size,tpeaks.size))
count=0
for i, ni in enumerate(mni):
    for j,tp in enumerate(tpeaks):
        beta0 = 1.5
        lpeak[i,j]=valenti_bol(tp,ni,0)
        betafit = fit_bootstrap([beta0], [tp, ni],[lpeak[i,j]],[10**0.1], khatami_err,bounds=([0.0], [np.inf]),
                                 errfunc=True, perturb=False, n=30, nproc=8)
        betas[i,j]= betafit[0][0]
        if i==0:
            ax.text( mni[-1]+0.05,betas[i,j]-0.03, '{}'.format(np.round(tp,2)), fontsize=10)
        #print mni[-1],lpeak[i,j],tp, betas[i,j]
        count=count+1
    ax.plot(np.repeat(ni,betas[i,:].shape),betas[i,:],marker='x',ls='-',color='k',lw=0.5)
    ax.text(ni, 1.47, '{}'.format(np.round(ni,2)),fontsize=10)
ax.plot(mni,betas,marker='x',ls='-',color='k',lw=0.5)
ax.text(0.1, 1.4, '{}'.format(r'M$_{\rm Ni}$'),fontsize=13)
ax.text(0.67, 1.75, '{}'.format(r't$_{\rm peak}$'),fontsize=13)
betafit = fit_bootstrap([beta0], [17.5, 0.24], [6e42], [10 ** 0.1], khatami_err, bounds=([0.0], [np.inf]),
                        errfunc=True, perturb=False, n=50, nproc=8)
print betafit[0][0]
ax.scatter(0.24,betafit[0][0], marker='*',s=140,color='salmon',edgecolor='k')
ax.text(0.24, betafit[0][0]-0.085, '{}'.format('Barnes19'),fontsize=13)
    # ax.annotate(r'Ic-BL', color='k',
    #              xy=(0.1, 0.9),
    #              xycoords='axes fraction', fontsize=16)
        #print tp, ni,  lpeak[i,j],betas[i,j]
        #ax.scatter(ni,betas[i,j],marker='s',s=40,color='lightgray')

mni=np.logspace(np.log10(0.01),np.log10(1),5, endpoint=True)
tpeaks=np.linspace(5,80,5, endpoint=True)
ax.fill_between(mni, np.min(np.min(betas)), np.max(np.max(betas)), color='lightgray')
index_name = 0
index_lpeak = 2
index_tpeak = 1
index_tail_ni_mass = 3
names=[]
lpeaks=[]
tpeaks=[]
nimasses=[]
with open("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/dessart_models.txt", "r") as f_input:
    csv_reader = csv.reader(f_input, delimiter=",")
    rows = list(csv_reader)
    for row in rows:
        names.append(row[index_name])
        lpeaks.append(float(row[index_lpeak].replace(")","").split("(")[0])*(10**float(row[index_lpeak].replace(")","").split("(")[1])))
        tpeaks.append(float(row[index_tpeak].replace(")", "").split("(")[0]) *(10** float(
            row[index_tpeak].replace(")", "").split("(")[1])))
        nimasses.append(float(row[index_tail_ni_mass].replace(")", "").split("(")[0]) * (10**float(
            row[index_tail_ni_mass].replace(")", "").split("(")[1])))
e_Ni=3.90e10
M_sun=1.9884e33
plt.figure(2)
ax2=plt.subplot(111)
for i,ni in enumerate(nimasses):
    #print ni
    if i == 0:
        label = 'Dessart16'
    else:
        label = None
    betafit = fit_bootstrap([beta0], [tpeaks[i], ni], [lpeaks[i]], [1], khatami_err, bounds=([0.0], [np.inf]),
                            errfunc=True, perturb=False, n=30, nproc=8)
    ax.scatter(ni, betafit[0][0], marker='^', s=120, color='yellow', edgecolor='k',label=label)
    #print tpeaks[i],lpeaks[i],ni
    ax2.scatter(tpeaks[i],lpeaks[i]/(M_sun*ni*e_Ni) , marker='^', s=120, color='yellow', edgecolor='k')
    mni=nickel_mass_khatami(tpeaks[i],lpeaks[i],10,9.0/8.0)
    ax2.scatter( tpeaks[i],lpeaks[i] / (M_sun * mni[0] * e_Ni), marker='.', s=5, color='yellow', edgecolor='k')

    #print arnett_param(tpeaks[i])
ax.legend(loc=3,fontsize=12,ncol=1,fancybox=True)
ax.set_ylim([0.005, 2])
for tpeak in range(10,50,1):
    ax2.set_xlim([10,50])
    ax2.set_ylim([0.1,0.4])
    ax2.scatter(tpeak,arnett_param(tpeak), marker='.', s=5, color='yellow', edgecolor='k')
    ax2.set_xlim([10,50])
    ax2.set_ylim([0.1,0.4])

print nickel_mass_khatami(17.5,6e42,10,1.3)

f.savefig('./Plots/models.pdf')

plt.show()