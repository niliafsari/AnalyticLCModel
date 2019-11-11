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
import random

# with open('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SN1993J_results.csv', mode='r') as csv_file:
#     csv_reader = csv.DictReader(csv_file)
#     line_count = 0
#     header=csv_reader.fieldnames

first=0
b=glob.glob("Data/*_results_wygoda.csv")
b= sorted(b)

b.remove('Data/iPTF13bvn_results_wygoda.csv')
b.insert(21,'Data/iPTF13bvn_results_wygoda.csv')

with open('Data/SNdata_wygoda.csv','w') as f:
    for file in b:
        #print file
        with open(file, mode='r') as csv_file:
            if first==0:
                csv_reader = csv.DictReader(csv_file)
                w = csv.DictWriter(f, fieldnames=csv_reader.fieldnames)
                w.writeheader()
                first=1
            line_count = 0
            for row in csv_file:
                if (line_count == 1) | ((line_count==0) & (file==b[0])):
                    temp=row.replace(", ",";").\
                        replace("[","").replace("]","").replace('"','').replace("(","").replace(")","")\
                        .replace("'","").replace(';UGC 1238',',UGC 1238')
                    f.write(temp)
                line_count = +1
#
types=[]
beta_req=[]
sn_names=[]
sn_type=[]
sn_host=[]
dist=[]
dist_e=[]
ebv_gal=[]
ebv_gal_e=[]
ebv_host=[]
ebv_host_e=[]
t0=[]
t0_e=[]
arnett_ni=[]
ni=[]
line=0
with open("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda.csv", "r") as f_input:
    csv_reader = csv.reader(f_input, delimiter=",")
    rows = list(csv_reader)
    for row in rows:
        if line==0:
            index_name=row.index('name')
            index_host = row.index('host')
            index_sn_type = row.index('sn_type_full')
            index_ebv_gal=row.index('ebv_gal')
            index_ebv_host=row.index('ebv_host')
            index_distance=row.index('distance')
            index_t0=row.index('t0')
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

            line=line+1
        else:
            line=line+1
            sn_names.append(row[index_name])
            sn_host.append(row[index_host])
            sn_type.append(row[index_sn_type])
            dist.append(np.round(float(row[index_distance].split(";")[0])/1e6,2))
            dist_e.append(np.round(float(row[index_distance].split(";")[1])/1e6,2))
            ebv_gal.append(float(row[index_ebv_gal].split(";")[0]))
            ebv_gal_e.append(float(row[index_ebv_gal].split(";")[1]))
            ebv_host.append(float(row[index_ebv_host].split(";")[0]))
            ebv_host_e.append(float(row[index_ebv_host].split(";")[1]))
            t0.append(float(row[index_t0].split(";")[0]))
            t0_e.append(float(row[index_t0].split(";")[1]))
            arnett_ni.append(float(row[index_arnett_ni_mass].split(";")[0]))
            ni.append(float(row[index_tail_ni_mass].split(";")[0]))

sn_names=np.array(sn_names)
ab = np.zeros((28,), dtype=[('SN name', 'U10'), ('Host','U16'),('Type','U5'),('Distance', np.float64),
                            ('Distance_e', np.float64),('E(B-V)_gal', np.float64),('E(B-V)_gal_e', np.float64),
                            ('E(B-V)_host', np.float64),('E(B-V)_host_e', np.float64), ('t0', np.float64),('t0_e', np.float64)])
ab['SN name'] = sn_names
ab['Host'] = np.array(sn_host)
ab['Type'] = np.array(sn_type)
ab['Distance'] = np.array(dist)
ab['Distance_e'] = np.array(dist_e)
ab['E(B-V)_gal'] = np.round(np.array(ebv_gal),4)
ab['E(B-V)_gal_e'] = np.round(np.array(ebv_gal_e),4)
ab['E(B-V)_host'] = np.round(np.array(ebv_host),2)
ab['E(B-V)_host_e'] = np.round(np.array(ebv_host_e),2)
ab['t0'] = np.round(np.array(t0),1)
ab['t0_e'] = np.round(np.array(t0_e),1)
np.savetxt("Data/table1.csv", ab, fmt="%10s & %12s & %6s & %10.1f $\pm$ %10.1f & %10.4f $\pm$ %10.4f & %10.2f $\pm$ %10.2f & %10.1f $\pm$ %10.1f", newline=' \\\\\n')
