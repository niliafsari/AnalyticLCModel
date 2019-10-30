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
marker = itertools.cycle(('8', 'p','>','^','v','<','s', 'o','h','<','*','P','X','+','H','1','2','3','4'))
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
f2 = plt.figure(2)
ax2=plt.subplot(111)
f3 = plt.figure(3)
ax3=plt.subplot(111)
f4 = plt.figure(4)
ax4=plt.subplot(111)
# f5 = plt.figure(5)
# ax5=plt.subplot(111)
f6 = plt.figure(6)
ax6=plt.subplot(111)
f7 = plt.figure(7)
ax7=plt.subplot(111)
f8 = plt.figure(8)
ax8=plt.subplot(111)

f9 = plt.figure(9,figsize=(8,5))
ax9=plt.subplot(111)
f10 = plt.figure(10)
ax10=plt.subplot(111)



dict={'Ic':'red', 'Ib':'blue','Ic BL':'orange', 'IcGRB':'orange','IIb':'green','Ib pec':'blue', 'Ibn':'blue'}
# Load all rows from the file into a variable called rows
yflag=0
gflag=0
rflag=0
bflag=0
line=0
plot_lines = []
with open("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata.csv", "r") as f_input:
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
            sn_names.append(row[index_name])
            sn_type.append(row[index_sn_type])
            line=line+1
            arnett_ni.append(row[index_arnett_ni_mass].split(";")[0])
            ni.append(row[index_tail_ni_mass].split(";")[0])
            ax.errorbar(float(row[index_tail_ni_mass].split(";")[0]),float(row[index_arnett_ni_mass].split(";")[0]),
                        xerr=float(row[index_tail_ni_mass].split(";")[1])*5,
                        yerr=float(row[index_arnett_ni_mass].split(";")[1]),
                        marker=mark,color=dict[row[index_sn_type]],label=row[index_name], linewidth=0.4, elinewidth=0.5)
            ax10.errorbar(float(row[index_tail_ni_mass].split(";")[0]),float(row[index_arnett_ni_mass].split(";")[0]),
                        xerr=float(row[index_tail_ni_mass].split(";")[1])*5,
                        yerr=float(row[index_arnett_ni_mass].split(";")[1]),
                        marker=mark,color=dict[row[index_sn_type]],label=row[index_name], linewidth=0.4, elinewidth=0.5,mec='k')
            # ax9.errorbar(float(row[index_tail_ni_mass].split(";")[0]),float(row[index_arnett_ni_mass].split(";")[0]),
            #             xerr=float(row[index_tail_ni_mass].split(";")[1])*5,
            #             yerr=float(row[index_arnett_ni_mass].split(";")[1]),
            #             marker=mark,color=dict[row[index_sn_type]],label=row[index_name], linewidth=0.4, elinewidth=0.5)
            ax7.scatter(float(row[index_tail_ni_mass].split(";")[0]), float(row[index_beta_req].split(";")[0]), s=10, marker=mark, color=col, label=row[index_name], edgecolors='black',
                       linewidth=0.2)
            ax8.scatter(float(row[index_tail_ni_mass].split(";")[0]), float(row[index_ni_khatami].strip("(").strip(")").split(";")[0]), s=10, marker=mark, color=col, label=row[index_name], edgecolors='black',
                       linewidth=0.2)
            ax9.scatter(float(row[index_tpeak].split(";")[0]),float(row[index_lpeak].split(";")[0]),s=15, marker=mark, color=dict[row[index_sn_type]], label=row[index_name])
#             if (dict[row[1].replace(" ","")]=='yellow') and (yflag==0):
#                 ax5.scatter(float(row[2]), float(row[5]), s=25, marker='^', color=dict[row[1].replace(" ","")],label='Ic BL', edgecolors='black',  linewidth=0.2)
#                 yflag=1
#             elif (dict[row[1].replace(" ","")]=='red')and (rflag==0):
#                 ax5.scatter(float(row[2]), float(row[5]), s=25, marker='^', color=dict[row[1].replace(" ", "")],
#                             label='Ib', edgecolors='black', linewidth=0.2)
#                 rflag=1
#             elif (dict[row[1].replace(" ","")]=='blue')and (bflag==0):
#                 ax5.scatter(float(row[2]), float(row[5]), s=25, marker='^', color=dict[row[1].replace(" ", "")],
#                             label='Ic', edgecolors='black', linewidth=0.2)
#                 bflag=1
#             elif (dict[row[1].replace(" ","")]=='green')and (gflag==0):
#                 ax5.scatter(float(row[2]), float(row[5]), s=25, marker='^', color=dict[row[1].replace(" ", "")],
#                             label='IIb', edgecolors='black', linewidth=0.2)
#                 gflag=1
#             else:
#                 ax5.scatter(float(row[2]), float(row[5]), s=25, marker='^', color=dict[row[1].replace(" ", "")], edgecolors='black', linewidth=0.2)
            if (row[index_prentice_ni].split(";")[0]=='') | (row[index_prentice_ni].split(";")[0]=='NA'):
                continue
            ax2.errorbar(float(row[index_tail_ni_mass].split(";")[0]), float(row[index_prentice_ni].split(";")[0]),xerr=float(row[index_tail_ni_mass].split(";")[1]),yerr=float(row[index_prentice_ni].split(";")[1]), fmt='o', markersize=5, marker=mark, color=col, label=row[index_name], linewidth=0.2, capsize=2, capthick=2)
            ax3.scatter(float(row[index_lpeak].split(";")[0]),  float(row[index_prentice_peakL].split(";")[0]), s=10, marker=mark, color=col, label=row[index_name], edgecolors='black',
                       linewidth=0.2)
            ax4.scatter(float(row[index_tpeak].split(";")[0]),  float(row[index_prentice_peakt].split(";")[0]), s=10, marker=mark, color=col, label=row[index_name], edgecolors='black',
                       linewidth=0.2)
            ax6.scatter( float(row[index_prentice_ni].split(";")[0]), float(row[index_arnett_ni_mass].split(";")[0]), s=10, marker=mark, color=col, label=row[index_name], edgecolors='black',
                       linewidth=0.2)
#
#
#plt.figure(1,figsize=(10,10))
plt.figure(1)
ax.legend(loc=1,bbox_to_anchor=(1.3, 1.1),
          fancybox=True, ncol=1, fontsize =8)
plt.plot(np.arange(0, 1, 0.1),np.arange(0, 1, 0.1),'--',color='black',linewidth=1)
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
# ax.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax.set_ylabel(r'Arnett $\rm M_{Ni} \ (M_\odot$)')
ax.set_xlabel(r'Tail $\rm M_{Ni} \ (M_\odot$)')
plt.xlim(0, 0.9)
plt.ylim(0, 0.9)
ax.xaxis.set_ticks(np.arange(0, 1, 0.1))
ax.yaxis.set_ticks(np.arange(0, 1, 0.1))
plt.tight_layout()
ax.set_aspect('equal', adjustable='box')
f.savefig('/home/afsari/PycharmProjects/typeIbcAnalysis/Plots/Ni_Arnett_Tail1.pdf', bbox_inches='tight')

plt.figure(10)
#plt.figure(9,figsize=(10,10))
ax10.legend(loc=1,bbox_to_anchor=(1.35,1.1),
          fancybox=True, ncol=1, fontsize =8)
#legend1 = plt.legend(plot_lines[0], ["algo1", "algo2", "algo3"], loc=1)
lines = ax10.get_lines()
plt.annotate(r'Ic-BL',color=dict['Ic BL'],
             xy=(0.1, 0.9),
             xycoords='axes fraction', fontsize=16)
plt.annotate(r'IIb',color=dict['IIb'],
             xy=(0.1, 0.85),
             xycoords='axes fraction', fontsize=16)
plt.annotate(r'Ib',color=dict['Ib'],
             xy=(0.1, 0.75),
             xycoords='axes fraction', fontsize=16)
plt.annotate(r'Ic',color=dict['Ic'],
             xy=(0.1, 0.8),
             xycoords='axes fraction', fontsize=16)

ax10.set_xscale("log")
ax10.set_yscale("log")
plt.loglog(np.arange(0.001, 1.1, 0.01),np.arange(0.001, 1.1, 0.01),'--',color='black',linewidth=1)
plt.xlim(0.01, 1)
plt.ylim(0.01, 1)
#plt.loglog(np.arange(0, 1, 0.1),np.arange(0, 1, 0.1),'--',color='black',linewidth=1)
# ax9.yaxis.set_minor_locator(AutoMinorLocator(5))
# ax9.xaxis.set_minor_locator(AutoMinorLocator(5))
# ax.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax10.xaxis.set_ticks_position('both')
ax10.yaxis.set_ticks_position('both')
ax10.tick_params(direction = 'in',which ='both')
ax10.set_aspect('equal', adjustable='box')
ax10.set_ylabel(r'Arnett $\rm M_{Ni} \ (M_\odot$)')
ax10.set_xlabel(r'Tail $\rm M_{Ni} \ (M_\odot$)')
f10.savefig('/home/afsari/PycharmProjects/typeIbcAnalysis/Plots/Ni_Arnett_Tail_loglog.pdf', bbox_inches='tight')


plt.figure(9)
#plt.figure(9,figsize=(10,10))
ax9.legend(loc='best',bbox_to_anchor=(1.9, 1),
          fancybox=True, ncol=2, fontsize =6)
#plt.loglog(np.arange(0, 1, 0.1),np.arange(0, 1, 0.1),'--',color='black',linewidth=1)
# ax9.yaxis.set_minor_locator(AutoMinorLocator(5))
# ax9.xaxis.set_minor_locator(AutoMinorLocator(5))
# ax.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax9.set_ylabel(r'$\rm log( L_{p} \ (erg \ s^{-1}))$',fontsize=12)
ax9.set_xlabel(r't$_p$ (days)')
# ax9.set_xscale("log")
# ax9.set_yscale("log")
# plt.xlim(0, 0.9)
# plt.ylim(0, 0.9)
# ax9.xaxis.set_ticks(np.arange(0, 1, 0.1))
# ax9.yaxis.set_ticks(np.arange(0, 1, 0.1))
#plt.tight_layout()
#ax9.set_aspect('equal', adjustable='box')
#f9.savefig('/home/afsari/PycharmProjects/typeIbcAnalysis/Plots/Ni_Arnett_Tail_loglog.pdf', bbox_inches='tight')
#
#
ax2.legend(loc='best', bbox_to_anchor=(1.3, -0.2),
          fancybox=True, ncol=4, fontsize =10)
ax2.plot(np.arange(0, 0.8, 0.1),np.arange(0, 0.8, 0.1),'--',color='black',linewidth=1)

ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
# ax.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax2.set_ylabel(r'Prentice $\rm M_{Ni} \ (M_\odot$)')
ax2.set_xlabel(r'Tail $\rm M_{Ni} \ (M_\odot$)')
ax2.set_xlim(0, 0.6)
ax2.set_ylim(0, 0.6)
ax2.xaxis.set_ticks(np.arange(0, 0.7, 0.1))
ax2.yaxis.set_ticks(np.arange(0, 0.7, 0.1))
ax2.set_aspect('equal', adjustable='box')
#plt.tight_layout()
f2.savefig('/home/afsari/PycharmProjects/typeIbcAnalysis/Plots/Ni_prenticeComparison1.pdf', bbox_inches='tight')

ax3.legend(loc='best', bbox_to_anchor=(1.3, -0.3),
          fancybox=True, ncol=3, fontsize =10)
ax3.plot(np.arange(42.2, 43.6, 0.2),np.arange(42.2, 43.6, 0.2),'--',color='black',linewidth=1)

ax3.yaxis.set_minor_locator(AutoMinorLocator(4))
ax3.xaxis.set_minor_locator(AutoMinorLocator(4))
# ax.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax3.set_ylabel(r'Prentice $\rm log( L_{p} \ (erg \ s^{-1}))$',fontsize=12)
ax3.set_xlabel(r'$\rm log( L_{p} \ (erg \ s^{-1}))$',fontsize=12)
ax3.set_xlim(42.2, 43.4)
ax3.set_ylim(42.2, 43.4)
ax3.xaxis.set_ticks(np.arange(42.2, 43.6, 0.4))
ax3.yaxis.set_ticks(np.arange(42.2, 43.6, 0.4))
plt.tight_layout()
ax3.set_aspect('equal', adjustable='box')
f3.savefig('./Plots/Lpeak_prenticeComparison1.pdf', bbox_inches='tight')
#
ax4.legend(loc='best', bbox_to_anchor=(1.3, -0.3),
          fancybox=True, ncol=3, fontsize =10)
ax4.plot(np.arange(10, 30, 5),np.arange(10, 30, 5),'--',color='black',linewidth=1)

ax4.yaxis.set_minor_locator(AutoMinorLocator(4))
ax4.xaxis.set_minor_locator(AutoMinorLocator(4))
# ax.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax4.set_ylabel(r'Prentice $\rm t_{p} \ (days)$',fontsize=12)
ax4.set_xlabel(r'$\rm t_{p} \ (days)$',fontsize=12)
ax4.set_xlim(10, 25)
ax4.set_ylim(10, 25)
ax4.xaxis.set_ticks(np.arange(10, 30, 5))
ax4.yaxis.set_ticks(np.arange(10, 30, 5))
plt.tight_layout()
ax4.set_aspect('equal', adjustable='box')
f4.savefig('./Plots/tpeak_prenticeComparison1.pdf', bbox_inches='tight')
#
# ax5.legend(loc='best', bbox_to_anchor=(0.5, -0.3),
#           fancybox=True, ncol=3, fontsize =10)
# #ax5.plot(np.arange(10, 30, 5),np.arange(10, 30, 5),'--',color='black',linewidth=1)
#
# ax5.yaxis.set_minor_locator(AutoMinorLocator(4))
# ax5.xaxis.set_minor_locator(AutoMinorLocator(4))
# # ax.xaxis.set_tick_params(width=1.5)
# # ax.yaxis.set_tick_params(width=1.5)
# ax5.set_xlabel(r'$\rm log( L_{p} \ (erg \ s^{-1}))$',fontsize=12)
# ax5.set_ylabel(r'Tail $\rm M_{Ni} \ (M_\odot$)')
# ax5.set_xlim(42.2, 43.2)
# ax5.set_ylim(0, 0.3)
# ax5.xaxis.set_ticks(np.arange(42.2, 43.6, 0.4))
# ax5.yaxis.set_ticks(np.arange(0, 0.3+0.05, 0.05))
# f5.savefig('./Plots/Ni_Tail_Lpeak.pdf', bbox_inches='tight')
# plt.tight_layout()
# #ax5.set_aspect('equal', adjustable='box')
#
ax6.legend(loc='best', bbox_to_anchor=(1.3, -0.2),
          fancybox=True, ncol=3, fontsize =10)
ax6.plot(np.arange(0, 0.8, 0.1),np.arange(0, 0.8, 0.1),'--',color='black',linewidth=1)

ax6.yaxis.set_minor_locator(AutoMinorLocator(5))
ax6.xaxis.set_minor_locator(AutoMinorLocator(5))
# x.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax6.set_ylabel(r'Prentice $\rm M_{Ni} \ (M_\odot$)')
ax6.set_xlabel(r'Arnett $\rm M_{Ni} \ (M_\odot$)')
ax6.set_xlim(0, 0.6)
ax6.set_ylim(0, 0.6)
ax6.xaxis.set_ticks(np.arange(0, 0.7, 0.1))
ax6.yaxis.set_ticks(np.arange(0, 0.7, 0.1))
ax6.set_aspect('equal', adjustable='box')
f6.savefig('./Plots/Ni_Prentice_Arnett1.pdf', bbox_inches='tight')
#
#
ax7.legend( loc='best', bbox_to_anchor=(1.3, -0.2),
          fancybox=True, ncol=3, fontsize =10)
#ax6.plot(np.arange(0, 0.8, 0.1),np.arange(0, 0.8, 0.1),'--',color='black',linewidth=1)
ax7.yaxis.set_minor_locator(AutoMinorLocator(5))
ax7.xaxis.set_minor_locator(AutoMinorLocator(5))
# x.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax7.set_ylabel(r'$\beta$ required')
ax7.set_xlabel(r'Tail $\rm M_{Ni} \ (M_\odot$)')
ax7.xaxis.set_ticks(np.arange(0, 0.4, 0.1))
ax7.set_xlim(0, 0.3)
ax7.set_ylim(0, 1.5)

#ax7.yaxis.set_ticks(np.arange(0, 0.7, 0.1))
#ax7.set_aspect('equal', adjustable='box')
plt.tight_layout()
f7.savefig('./Plots/Ni_Tail_BetaReq1.pdf', bbox_inches='tight')

#ax5.set_aspect('equal', adjustable='box')

ax8.legend(loc='best', bbox_to_anchor=(1.3, -0.2),
          fancybox=True, ncol=3, fontsize =10)
ax8.plot(np.arange(0, 0.8, 0.1),np.arange(0, 0.8, 0.1),'--',color='black',linewidth=1)

ax8.yaxis.set_minor_locator(AutoMinorLocator(5))
ax8.xaxis.set_minor_locator(AutoMinorLocator(5))
# x.xaxis.set_tick_params(width=1.5)
# ax.yaxis.set_tick_params(width=1.5)
ax8.set_ylabel(r'$\rm M_{Ni} \ (M_\odot$) from suggested $\beta$')
ax8.set_xlabel(r'Tail $\rm M_{Ni} \ (M_\odot$)')
ax8.set_xlim(0, 0.6)
ax8.set_ylim(0, 0.6)
ax8.xaxis.set_ticks(np.arange(0, 0.7, 0.1))
ax8.yaxis.set_ticks(np.arange(0, 0.7, 0.1))
ax8.set_aspect('equal', adjustable='box')
plt.tight_layout()
f8.savefig('./Plots/Ni_Tail_fromSugBeta1.pdf', bbox_inches='tight')
plt.show()