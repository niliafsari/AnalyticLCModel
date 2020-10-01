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

print np.median(lpeaks),np.median(tpeaks),np.median(lfac),np.median(lfac)