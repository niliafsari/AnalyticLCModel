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
import pandas as pd

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

f3=plt.figure(3)
ax3=plt.subplot(111)




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
lfac=[]
lpf=[]
lpeaks=[]
tpeaks=[]
tpeaks_e=[]
lpeaks_e=[]
dict_list={'orange':orange_marker,'blue':blue_marker, 'red':red_marker, 'green':green_marker}
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
            mark = marker.next()
            col = colors.next()
            while mark in dict_list[dict[row[index_sn_type]]]:
                mark = marker.next()
            sn_names.append(row[index_name])
            sn_type.append(row[index_sn_type])
            if (len(orange_marker)==0) & (dict[row[index_sn_type]]=='orange'):
                label = 'Ic-BL'
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
            lpeaks.append(float(row[index_lpeak].split(";")[0]))
            tpeaks.append(float(row[index_tpeak].split(";")[0]))
            lpeaks_e.append(float(row[index_lpeak].split(";")[1]))
            tpeaks_e.append(float(row[index_tpeak].split(";")[1]))
            types.append(row[index_sn_type])
            beta_req.append(float(row[index_beta_req].split(";")[0]))
            ax.errorbar(float(row[index_tail_ni_mass].split(";")[0]), float(row[index_beta_req].split(";")[0]), yerr=float(row[index_beta_req].split(";")[1]),
                         xerr=float(row[index_tail_ni_mass].split(";")[1]),
                                                                             marker=dict_marker[row[index_sn_type]],
                                                                             color='deepskyblue',
                                                                             label=label, linewidth=0.4,
                                                                             elinewidth=0.5, mec='k',markersize=8)
            lbol,lbol_err=lbol_khatami(float(row[index_tpeak].split(";")[0]), float(row[index_tail_ni_mass].split(";")[0]),
                         float(row[index_tail_ni_mass].split(";")[1]), 9.0 / 8.0)
            print "khatami", np.log10(lbol), "lpeak",float(row[index_lpeak].split(";")[0]),row[index_name]
            lfrac=(10**float(row[index_lpeak].split(";")[0]))/ lbol
            lfac.append(1-1/lfrac)
            lpf.append((10**float(row[index_lpeak].split(";")[0]))*(1-1/lfrac))
            lfrac_err=lfrac*(np.sqrt((lbol_err/lbol)**2+(10**float(row[index_lpeak].split(";")[1])/10**float(row[index_lpeak].split(";")[0]))**2))
            ax3.errorbar(float(row[index_beta_req].split(";")[0]), lfrac, xerr=float(row[index_beta_req].split(";")[1]), yerr=lfrac_err,markersize=8,  marker=mark, mec='k',color=dict[row[index_sn_type]], label=row[index_name], linewidth=0.4)
            if dict[row[index_sn_type]]=='orange':
                orange_marker.append(mark)
            elif dict[row[index_sn_type]]=='blue':
                blue_marker.append(mark)
            elif dict[row[index_sn_type]]=='green':
                green_marker.append(mark)
            elif dict[row[index_sn_type]]=='red':
                red_marker.append(mark)

print np.median(lpeaks),np.median(tpeaks),np.median(lfac)

columns=['tp','Lp','beta','Mni','model']
df_models= pd.DataFrame(columns=columns)
df_models['tp'] = tpeaks
df_models['Lp']=lpeaks
df_models['beta']=beta_req
df_models['Mni']=ni
df_models['model']=sn_names
print(df_models)

columns=['name','Lp','f','Lp|f|']
df2= pd.DataFrame(columns=columns)
df2['Lp']=lpeaks
df2['Lp|f|']=np.log10(np.array(lpf))
df2['f']=lfac
df2['name']=sn_names
df2.to_csv('Data/lpf.csv', index=False)



lpeaks=np.array(lpeaks)
sn_type=np.array(sn_type)
lfac=np.array(lfac)


print "IIb",np.log10(np.mean(10**(lpeaks[sn_type=='IIb']))),np.log10(np.median(10**(lpeaks[sn_type=='IIb']))), np.std(np.log10(10**(lpeaks[sn_type=='IIb'])))
print "Ib",np.log10(np.mean(10**(lpeaks[sn_type=='Ib']))),np.log10(np.median(10**(lpeaks[sn_type=='Ib']))), np.std(np.log10(10**(lpeaks[sn_type=='Ib'])))
print "Ic",np.log10(np.mean(10**(lpeaks[sn_type=='Ic']))),np.log10(np.median(10**(lpeaks[sn_type=='Ic']))), np.std(np.log10(10**(lpeaks[sn_type=='Ic'])))
print "Ic-BL:",np.log10(np.mean(10**(lpeaks[sn_type=='Ic BL']))),np.log10(np.median(10**(lpeaks[sn_type=='Ic BL']))), np.std(np.log10(10**(lpeaks[sn_type=='Ic BL'])))
print "all:",np.log10(np.mean(10**(lpeaks))),np.log10(np.median(10**(lpeaks))), np.std(np.log10(10**(lpeaks)))

# lpeaks=lpeaks[lfac>0]
# print lpeaks
# sn_type=sn_type[lfac>0]
# lfac=lfac[lfac>0]
# lpeaks=np.log10(np.multiply(10**lpeaks,lfac))
# print "IIb",np.log10(np.mean(10**(lpeaks[sn_type=='IIb']))),np.log10(np.median(10**(lpeaks[sn_type=='IIb']))), np.std(np.log10(10**(lpeaks[sn_type=='IIb'])))
# print "Ib",np.log10(np.mean(10**(lpeaks[sn_type=='Ib']))),np.log10(np.median(10**(lpeaks[sn_type=='Ib']))), np.std(np.log10(10**(lpeaks[sn_type=='Ib'])))
# print "Ic",np.log10(np.mean(10**(lpeaks[sn_type=='Ic']))),np.log10(np.median(10**(lpeaks[sn_type=='Ic']))), np.std(np.log10(10**(lpeaks[sn_type=='Ic'])))
# print "Ic-BL:",np.log10(np.mean(10**(lpeaks[sn_type=='Ic BL']))),np.log10(np.median(10**(lpeaks[sn_type=='Ic BL']))), np.std(np.log10(10**(lpeaks[sn_type=='Ic BL'])))
# print "all:",np.log10(np.mean(10**(lpeaks))),np.log10(np.median(10**(lpeaks))), np.std(np.log10(10**(lpeaks)))


types=np.array(types)
beta_req=np.array(beta_req)
ni_plot=np.array(ni)
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))


ax.set_ylabel(r'$\beta$',fontsize=18)
ax.set_xlabel(r'$ M_{\rm Ni} \ (M_\odot$)',fontsize=18)
ax.set_xscale("log")
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.tick_params(direction = 'in',which ='both')
ax.set_xlim(0.01, 1)
ax.set_ylim(0.0001, 2)


ax3.set_ylabel(r'$ L_{\rm p} / L_{\rm khatami} $',fontsize=18)
ax3.set_xlabel(r'Tuned $\beta$ ',fontsize=18)
ax3.legend(loc=1,bbox_to_anchor=(1.35,1.1),
          fancybox=True, ncol=1, fontsize =8)
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')
ax3.tick_params(direction = 'in',which ='both')

from SNAP5.Analysis.LCFitting import fit_bootstrap
mni=np.logspace(np.log10(0.02),np.log10(1),5, endpoint=False)
tpeaks_a=np.logspace(np.log10(5),np.log10(80),4, endpoint=False)
lpeak=np.zeros((mni.size,tpeaks_a.size))
betas=np.zeros((mni.size,tpeaks_a.size))
count=0
for i, ni in enumerate(mni):
    for j,tp in enumerate(tpeaks_a):
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
    #ax.text(ni, 1.49, '{}'.format(np.round(ni,2)),fontsize=10)
ax.plot(mni,betas,marker='x',ls='-',color='k',lw=0.5)
#ax.text(0.05, 1.47, '{}'.format(r'$M_{\rm Ni}$'),fontsize=13)
ax.text(0.67, 1.75, '{}'.format(r'$t_{\rm p}$'),fontsize=13)




mni=np.logspace(np.log10(0.01),np.log10(1),5, endpoint=True)
#tpeaks=np.linspace(5,80,5, endpoint=True)
ax.fill_between(mni, np.min(np.min(betas)), np.max(np.max(betas)), color='lightgray')
index_name = 0
index_lpeak = 2
index_tpeak = 1
index_tail_ni_mass = 3
names=[]
lpeaks_d=[]
tpeaks_d=[]
nimasses=[]
with open("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/dessart_models.txt", "r") as f_input:
    csv_reader = csv.reader(f_input, delimiter=",")
    rows = list(csv_reader)
    for row in rows:
        names.append(row[index_name])
        lpeaks_d.append(float(row[index_lpeak].replace(")","").split("(")[0])*(10**float(row[index_lpeak].replace(")","").split("(")[1])))
        tpeaks_d.append(float(row[index_tpeak].replace(")", "").split("(")[0]) *(10** float(
            row[index_tpeak].replace(")", "").split("(")[1])))
        nimasses.append(float(row[index_tail_ni_mass].replace(")", "").split("(")[0]) * (10**float(
            row[index_tail_ni_mass].replace(")", "").split("(")[1])))
e_Ni=3.90e10
M_sun=1.9884e33
# plt.figure(2)
# ax2=plt.subplot(111)
dessart_beta=[]
for i,ni in enumerate(nimasses):
    #print ni
    if i == 0:
        label = r'Dessart16 (IIb \& Ib/c)'
    else:
        label = None
    betafit = fit_bootstrap([beta0], [tpeaks_d[i], ni], [lpeaks_d[i]], [1], khatami_err, bounds=([0.0], [np.inf]),
                            errfunc=True, perturb=False, n=30, nproc=8)
    dessart_beta.append(betafit[0][0])
    df_models = df_models.append(pd.Series([tpeaks_d[i], np.log10(lpeaks_d[i]), betafit[0][0], ni, 'Dessart'], index=df_models.columns), ignore_index=True)
    ax.scatter(ni, betafit[0][0], marker='^', s=120, color='yellow', edgecolor='k',label=label)

print "dessart" , np.min(np.array (dessart_beta)),np.mean(np.array (dessart_beta)), np.max(np.array (dessart_beta))
print (df_models)
dat=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/Nilou.dat', delimiter=',')
print dat.shape[1]
ertl_beta=[]
for i in range(0,dat.shape[0]):
    if dat[i,0]<8:
        if i == 0:
            label = 'Ertl19 (Ib/c)'
        else:
            label = None
        #print dat[i,1], dat[i,3]
        betafit = fit_bootstrap([beta0], [dat[i,1], dat[i,3]], [10**dat[i,2]], [1], khatami_err, bounds=([0.0], [np.inf]),
                                errfunc=True, perturb=False, n=30, nproc=8)
        ertl_beta.append(betafit[0][0])
        df_models = df_models.append(
            pd.Series([dat[i,1], dat[i,2], betafit[0][0], dat[i,3], 'Ertl'], index=df_models.columns),
            ignore_index=True)

        ax.scatter(dat[i,3], betafit[0][0], marker='p', s=100, color='springgreen', edgecolor='k',label=label)

dat=np.loadtxt('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/Io_data.csv', delimiter=',')
print dat.shape[1]
io_beta=[]
Msun = 4.74
Lsun = 3.84e33

for i in range(0,dat.shape[0]):
    if i == 0:
        label = 'Kleiser18 (Ib/c)'
    else:
        label = None
    #print dat[i,1], dat[i,3]
    lbol = Lsun * np.power(10, ((Msun - dat[i,1]) / 2.5))
    betafit = fit_bootstrap([beta0], [dat[i,0], dat[i,2]], [lbol], [1], khatami_err, bounds=([0.0], [np.inf]),
                            errfunc=True, perturb=False, n=30, nproc=8)
    #ertl_beta.append(betafit[0][0])
    df_models = df_models.append(
        pd.Series([dat[i,0], np.log10(lbol), betafit[0][0], dat[i,2], 'Kleiser'], index=df_models.columns),
        ignore_index=True)

    ax.scatter(dat[i,2], betafit[0][0], marker='<', s=100, color='m', edgecolor='k',label=label)

print "ertl:", np.min(np.array(ertl_beta)),np.mean(np.array(ertl_beta)), np.max(np.array(ertl_beta))
b=glob.glob("Data/*_barnes.csv")
i=0
beta_barnes=[]
for lc in b:
    if i == 0:
        label = 'Barnes18 (Ic-BL)'
    else:
        label = None
    LC=np.loadtxt(lc, delimiter=',')
    LC_p=np.max(LC[:,1])
    t_p=np.max(LC[np.argmax(LC[:,1]),0])-np.min(LC[:,0])
    print LC_p
    nickelmass=0.24
    t_p=17.5
    betafit = fit_bootstrap([beta0], [t_p, nickelmass], [LC_p], [1], khatami_err, bounds=([0.0], [np.inf]),
                            errfunc=True, perturb=False, n=300, nproc=8)
    beta_barnes.append(betafit[0][0])
    df_models = df_models.append(
        pd.Series([t_p, np.log10(LC_p), betafit[0][0],nickelmass , 'Barnes'], index=df_models.columns),
        ignore_index=True)
    ax.scatter(nickelmass, betafit[0][0], marker='*', s=120, color='salmon', edgecolor='k', label=label)
    i=i+1
print(df_models)
print "barnes",np.min(np.array(beta_barnes)),np.max(np.array(beta_barnes))
df_models.to_csv('Data/all_models_data.csv', index=False)
ax.plot(np.arange(0.001,1.1,0.1),np.repeat(1.125,(11,)),ls='--',color='k')
ax.text(0.011, 1.125+0.03, r'$\beta=1.125$',fontsize=13)
ax.plot(np.arange(0.001,1.1,0.1),np.repeat(1.6,(11,)),ls='--',color='k')
ax.text(0.011, 1.6+0.03, r'$\beta=1.6$',fontsize=13)
ax.plot(np.arange(0.001,1.1,0.1),np.repeat(0.82,(11,)),ls='--',color='k')
ax.text(0.011, 0.82+0.03, r'$\beta=0.82$',fontsize=13)
handles, labels = ax.get_legend_handles_labels()
print labels
order=[4,7,5,6,1,2,3,0]
ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc=3,fontsize=11,ncol=1,fancybox=True, frameon=False)

ax3.legend(fontsize=12,ncol=1,fancybox=True, frameon=False)
ax.set_ylim([-0.005, 2])


f.savefig('./Plots/models.pdf')


tpeaks=np.array(tpeaks)
lpeaks=np.array(lpeaks)
lfac=np.array(lfac)

f4=plt.figure(4,figsize=(6,5))
plt.subplots_adjust(wspace=None, hspace=None)
ax4=plt.subplot(111)
# plt.subplots_adjust(wspace=None, hspace=None)
# plt.hist(tpeaks,bins=10, color='orange', ec='k')
# ax4.set_xlabel(r'$t_{\rm p}$ (days)',fontsize=15)
# ax4.set_ylabel(r'Count',fontsize=15)
# ax5=plt.subplot(212)
# plt.subplots_adjust(wspace=None, hspace=None)
# plt.hist(lpeaks,bins=10, color='yellow', ec='k')
# ax5.set_xlabel(r'Log $L_{\rm p} \ (\rm erg \ s^{-1})$',fontsize=15)
# ax5.set_ylabel(r'Count',fontsize=15)
# plt.subplots_adjust(wspace=None, hspace=None)

sc=ax4.scatter(tpeaks,10**np.array(lpeaks),zorder=100,c=ni_plot,s=35,marker='o',cmap='YlGnBu',edgecolors='k',linewidths=1,
                norm=matplotlib.colors.LogNorm(vmin=0.01,vmax=1))
ax4.errorbar(tpeaks,10**np.array(lpeaks),xerr=tpeaks_e,yerr=10**np.array(lpeaks_e),ls=None,fmt='.',markersize=0.1, color=None, capsize=3, capthick=0.5,elinewidth=0.5,zorder=0,ecolor='k')
ax4.set_yscale('log')
cb=plt.colorbar(sc)
ax4.set_ylabel(r'$L_{\rm p} \ (\rm erg \ s^{-1})$',fontsize=15)
ax4.set_xlabel(r'$t_{\rm p}$ (days)',fontsize=15)
cb.set_label(r'$M_{\rm Ni}$ ($M_\odot$)')
from matplotlib.ticker import LogLocator
# cb.ax.yaxis.set_major_locator(LogLocator())  # <- Why? See above.
# cb.set_ticks(cb.ax.yaxis.get_major_locator().tick_values(0.01, 1))
cb.set_ticks([0.01,0.1,1])
cb.set_ticklabels([0.01,0.1,1])
cb.update_ticks()
#ax4.yaxis.set_minor_locator(AutoMinorLocator(10))
ax4.xaxis.set_minor_locator(AutoMinorLocator(10))
ax4.xaxis.set_ticks_position('both')
ax4.yaxis.set_ticks_position('both')
ax4.tick_params(direction = 'in',which ='both')
plt.tight_layout()
types=np.array(types)
f4.savefig('./Plots/tp_vs_lp.pdf')

print np.mean(tpeaks[types=='IIb']), np.median(tpeaks[types=='IIb']), np.std(tpeaks[types=='IIb']), np.mean(lpeaks[types=='IIb']), np.median(lpeaks[types=='IIb']), np.std(lpeaks[types=='IIb']), np.mean(lfac[types=='IIb']), np.median(lfac[types=='IIb']), np.std(lfac[types=='IIb'])
print np.mean(tpeaks[types=='Ib']), np.median(tpeaks[types=='Ib']), np.std(tpeaks[types=='Ib']), np.mean(lpeaks[types=='Ib']), np.median(lpeaks[types=='Ib']), np.std(lpeaks[types=='Ib']), np.mean(lfac[types=='Ib']), np.median(lfac[types=='Ib']), np.std(lfac[types=='Ib'])
print np.mean(tpeaks[types=='Ic']), np.median(tpeaks[types=='Ic']), np.std(tpeaks[types=='Ic']), np.mean(lpeaks[types=='Ic']), np.median(lpeaks[types=='Ic']), np.std(lpeaks[types=='Ic']), np.mean(lfac[types=='Ic']), np.median(lfac[types=='Ic']), np.std(lfac[types=='Ic'])
print np.mean(tpeaks[types=='Ic BL']), np.median(tpeaks[types=='Ic BL']), np.std(tpeaks[types=='Ic BL']), np.mean(lpeaks[types=='Ic BL']), np.median(lpeaks[types=='Ic BL']), np.std(lpeaks[types=='Ic BL']), np.mean(lfac[types=='Ic BL']), np.median(lfac[types=='Ic BL']), np.std(lfac[types=='Ic BL'])
print np.mean(tpeaks), np.median(tpeaks), np.std(tpeaks), np.mean(lpeaks), np.median(lpeaks), np.std(lpeaks), np.mean(lfac), np.median(lfac), np.std(lfac)
print 10**np.mean(lpeaks)
print zip(sn_names, lfac)
print lfac
import pandas as pd
d = {'SN_name': sn_names, 'lfac': lfac}
df = pd.DataFrame(data=d)


res = np.zeros((5,), dtype=[('Beta mean', np.float64), ('Beta median', np.float64),('Beta std', np.float64),
                             ('tp mean', np.float64), ('tp median', np.float64), ('tp std', np.float64),
                             ('lp mean', np.float64), ('lp median', np.float64), ('lp std', np.float64),
                             ('f mean', np.float64), ('f median', np.float64), ('f std', np.float64)])
print np.round(np.mean(np.array(beta_req[types=='IIb'])),2)
res['Beta mean']=[np.round(np.mean(np.array(beta_req[types=='IIb'])),2),np.round(np.mean(np.array(beta_req[types=='Ib'])),2),
                    np.round(np.mean(np.array(beta_req[types=='Ic'])),2),np.round(np.mean(np.array(beta_req[types=='Ic BL'])),2),
                    np.round(np.mean(np.array(beta_req)),2)]

print  res['Beta mean']

res['Beta median']=[np.round(np.median(np.array(beta_req[types=='IIb'])),2),np.round(np.median(np.array(beta_req[types=='Ib'])),2)
    ,np.round(np.median(np.array(beta_req[types=='Ic'])),2),np.round(np.median(np.array(beta_req[types=='Ic BL'])),2),
                      np.round(np.median(np.array(beta_req)),2)]
res['Beta std']=[np.round(np.std(np.array(beta_req[types=='IIb'])),2),np.round(np.std(np.array(beta_req[types=='Ib'])),2),
                   np.round(np.std(np.array(beta_req[types=='Ic'])),2),np.round(np.std(np.array(beta_req[types=='Ic BL'])),2),
                   np.round(np.std(np.array(beta_req)),2)]
res['tp mean']=[np.round(np.mean(np.array(tpeaks[types=='IIb'])),1),np.round(np.mean(np.array(tpeaks[types=='Ib'])),1),
                  np.round(np.mean(np.array(tpeaks[types=='Ic'])),1),np.round(np.mean(np.array(tpeaks[types=='Ic BL'])),1),
                  np.round(np.mean(np.array(tpeaks)),1)]
res['tp median']=[np.round(np.median(np.array(tpeaks[types=='IIb'])),1),np.round(np.median(np.array(tpeaks[types=='Ib'])),1),
                  np.round(np.median(np.array(tpeaks[types=='Ic'])),1),np.round(np.median(np.array(tpeaks[types=='Ic BL'])),1),
                  np.round(np.median(np.array(tpeaks)),1)]
res['tp std']=[np.round(np.std(np.array(tpeaks[types=='IIb'])),1),np.round(np.std(np.array(tpeaks[types=='Ib'])),1),
               np.round(np.std(np.array(tpeaks[types=='Ic'])),1),np.round(np.std(np.array(tpeaks[types=='Ic BL'])),1)
    ,np.round(np.std(np.array(tpeaks)),1)]
res['lp mean']=[np.round(np.mean(np.array(lpeaks[types=='IIb'])),2),np.round(np.mean(np.array(lpeaks[types=='Ib'])),2),
                np.round(np.mean(np.array(lpeaks[types=='Ic'])),2),np.round(np.mean(np.array(lpeaks[types=='Ic BL'])),2),
                np.round(np.mean(np.array(lpeaks)),2)]
res['lp median']=[np.round(np.median(np.array(lpeaks[types=='IIb'])),2),np.round(np.median(np.array(lpeaks[types=='Ib'])),2),
                  np.round(np.median(np.array(lpeaks[types=='Ic'])),2),np.round(np.median(np.array(lpeaks[types=='Ic BL'])),2),
                  np.round(np.median(np.array(lpeaks)),2)]
res['lp std']=[np.round(np.std(np.array(lpeaks[types=='IIb'])),2),np.round(np.std(np.array(lpeaks[types=='Ib'])),2),
               np.round(np.std(np.array(lpeaks[types=='Ic'])),2),np.round(np.std(np.array(lpeaks[types=='Ic BL'])),2),
               np.round(np.std(np.array(lpeaks)),2)]
res['f mean']=[np.round(np.mean(np.array(lfac[types=='IIb'])),2),np.round(np.mean(np.array(lfac[types=='Ib'])),2),
               np.round(np.mean(np.array(lfac[types=='Ic'])),2),np.round(np.mean(np.array(lfac[types=='Ic BL'])),2),
               np.round(np.mean(np.array(lfac)),2)]
res['f median']=[np.round(np.median(np.array(lfac[types=='IIb'])),2),np.round(np.median(np.array(lfac[types=='Ib'])),2),
                 np.round(np.median(np.array(lfac[types=='Ic'])),2),np.round(np.median(np.array(lfac[types=='Ic BL'])),2),
                 np.round(np.median(np.array(lfac)),2)]
res['f std']=[np.round(np.std(np.array(lfac[types=='IIb'])),2),np.round(np.std(np.array(lfac[types=='Ib'])),2),
              np.round(np.std(np.array(lfac[types=='Ic'])),2),np.round(np.std(np.array(lfac[types=='Ic BL'])),2),
              np.round(np.std(np.array(lfac)),2)]

np.savetxt("Data/table4.csv", res, fmt="%10.2f & %10.2f  & %10.2f  & %10.1f & %10.1f & %10.1f & %10.2f & %10.2f & %10.2f & %10.2f & %10.2f & %10.2f ", newline=' \\\\\n')


df.to_csv('Data/lfrac.dat', index=False)


plt.show()