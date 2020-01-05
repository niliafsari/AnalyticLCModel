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

with open('Data/SNdata_wygoda_26dec.csv','w') as f:
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
T0=[]
T0_e=[]
tpeak=[]
tpeak_e=[]
lpeak=[]
lpeak_e=[]
beta_req_e=[]
ni=[]
ni_e=[]
arnett_ni_e=[]
bc=[]
line=0
with open("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda_26dec.csv", "r") as f_input:
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
            index_bc=row.index('band_bc')
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
            bc.append(row[index_bc].replace(";","-"))
            beta_req.append(float(row[index_beta_req].split(";")[0]))
            beta_req_e.append(float(row[index_beta_req].split(";")[1]))
            arnett_ni.append(float(row[index_arnett_ni_mass].split(";")[0]))
            arnett_ni_e.append(float(row[index_arnett_ni_mass].split(";")[1]))
            ni.append(float(row[index_tail_ni_mass].split(";")[0]))
            ni_e.append(float(row[index_tail_ni_mass].split(";")[1]))
            lpeak.append(float(row[index_lpeak].split(";")[0]))
            lpeak_e.append(float(row[index_lpeak].split(";")[1]))
            tpeak.append(float(row[index_tpeak].split(";")[0]))
            tpeak_e.append(float(row[index_tpeak].split(";")[1]))
            T0.append(float(row[index_tail_meje].split(";")[0])*32)
            T0_e.append(float(row[index_tail_meje].split(";")[1])*32)


sn_names=np.array(sn_names)
ab = np.zeros((28,), dtype=[('SN name', 'U10'), ('Host','U16'),('Type','U5'),('Distance', np.float64),
                            ('Distance_e', np.float64),('E(B-V)_gal', np.float64),('E(B-V)_gal_e', np.float64),
                            ('E(B-V)_host', np.float64),('E(B-V)_host_e', np.float64), ('t0', np.float64),('t0_e', np.float64)])
res = np.zeros((28,), dtype=[('SN name', 'U10'), ('BC bands','U16'),('Lpeak',np.float64),('Lpeak_e',np.float64),('Tpeak', np.float64),('Tpeak_e', np.float64),
                            ('mni_tail', np.float64),('mni_tail_e', np.float64),('T0', np.float64),('T0_e', np.float64),('mni_arnett', np.float64),('mni_arnett_e', np.float64),
                            ('beta', np.float64),('beta_e', np.float64)])
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
res['SN name'] = sn_names
res['BC bands']=np.array(bc)
res['Lpeak']=np.round(np.array(lpeak),2)
res['Lpeak_e']=np.round(np.array(lpeak_e),2)
res['Tpeak']=np.round(np.array(tpeak),1)
res['Tpeak_e']=np.round(np.array(tpeak_e),1)
res['mni_tail']=np.round(np.array(ni),3)
res['mni_tail_e']=np.round(np.array(ni_e),3)
res['mni_arnett']=np.round(np.array(arnett_ni),2)
res['mni_arnett_e']=np.round(np.array(arnett_ni_e),2)
res['T0']=np.round(np.array(T0),1)
res['T0_e']=np.round(np.array(T0_e),1)
res['beta']=np.round(np.array(beta_req),2)
res['beta_e']=np.round(np.array(beta_req_e),2)

np.savetxt("Data/table1_26dec.csv", ab, fmt="%10s & %12s & %6s & %10.1f (%10.1f) & %10.4f (%10.4f) & %10.2f (%10.2f) & %10.1f (%10.1f)", newline=' \\\\\n')
np.savetxt("Data/table2_26dec.csv", res, fmt="%10s & $%6s$ & %10.2f (%10.2f) & %10.1f (%10.1f) & %10.3f (%10.3f) & %10.1f (%10.1f) & %10.2f (%10.2f) & %10.2f (%10.2f)", newline=' \\\\\n')
