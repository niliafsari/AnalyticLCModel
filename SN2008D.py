import urllib
import os
import glob
import subprocess
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.time import Time
from scipy.interpolate import UnivariateSpline
from scipy.integrate import trapz
import csv
import json
from pprint import pprint
import os.path
import sys



def nickel_mass(t_peak,L_peak,beta):
    return (L_peak*(beta**2)*(t_peak/8.8)**2)/(2*3.9e10*((0.83*(1-beta*t_peak/8.8)*np.exp(-beta*t_peak/8.8))+(26.56*(1-((1+beta*t_peak/111.3)*np.exp(-beta*t_peak/111.3))))))


sys.path.insert(0, '/home/afsari/')
from SNAP5.Analysis import *



# with open("logs/sn_names_color.txt") as f:
#     file_names = f.readlines()

#
#s&d 2011
coef = {'B': 3.626, 'V': 2.742, 'I': 1.505, 'i': 1.698}

ebv= "0.65 Mazzali, 0.62 (Soderburg), 0.8 (Malesani)"
z = "0.006494 (Mazzali)"
t0= "January 9 13:32:40 UT (Soderburg)"
d= "32 (Mazzali), 27(Soderburg)"
JD="2454475.0643"
DM="32.46 \pm 0.15 (modjaz) -> d=21 \pm 2"


Lbol=np.loadtxt("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/2008D/Lbol_2008D.csv",delimiter=",")
print Lbol
print Lbol.shape
index=np.argsort(Lbol[:,0])
Lbol=Lbol[index,:]
#z = np.polyfit(Lbol[:,0], Lbol[:,1],11)
#fit = np.poly1d(z)
fit= UnivariateSpline(Lbol[:,0], Lbol[:,1], s=0.01)
t_Lbol=np.arange(0,80,0.3)
L_Lbol=fit(t_Lbol)
#
# ax=plt.subplot(111)
# plt.plot(t_Lbol,L_Lbol,'-.',color='green')
# plt.plot(Lbol[:,0],Lbol[:,1],color='red')
# plt.show()

te=np.loadtxt("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/2008D/Te_2008D.csv",delimiter=",")
index=np.argsort(te[:,0])
te=te[index,:]

print te
#z = np.polyfit(Lbol[:,0], Lbol[:,1],11)
#fit = np.poly1d(z)
fit1= UnivariateSpline(te[:,0], te[:,1], s=0.01)
t_te=np.arange(0,10,0.3)
T_te=fit1(t_te)
#

ax=plt.subplot(111)
plt.plot(t_te,np.power(10,T_te),'-.',color='green')
plt.plot(te[:,0],np.power(10,te[:,1]),color='red')
plt.show()


vph=np.loadtxt("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/2008D/vph_2008D.csv",delimiter=",")
index=np.argsort(vph[:,0])
vph=vph[index,:]

print vph
#z = np.polyfit(Lbol[:,0], Lbol[:,1],11)
#fit = np.poly1d(z)
fit2= UnivariateSpline(vph[:,0], vph[:,1], s=0.01)
t_vph=np.arange(1,30,0.3)
V_vph=fit2(t_vph)
#

ax=plt.subplot(111)
plt.plot(t_vph,V_vph,'--',color='green')
plt.scatter(vph[:,0],vph[:,1],s=10,color='red')
plt.show()
t=4
tdelta=np.zeros(shape=(5,1))
tmin=4.3*np.sqrt(np.power(10,fit(t))/10**42)*np.power(np.power(10,fit1(t))/10**4,-2)*np.power(fit2(t)*10**3/10**4,-1)
tdelta[0]=t-tmin
#print tmin
t=5

tmin=4.3*np.sqrt(np.power(10,fit(t))/10**42)*np.power(np.power(10,fit1(t))/10**4,-2)*np.power(fit2(t)*10**3/10**4,-1)
tdelta[1]= t-tmin
#print tmin
t=7
tmin=4.3*np.sqrt(np.power(10,fit(t))/10**42)*np.power(np.power(10,fit1(t))/10**4,-2)*np.power(fit2(t)*10**3/10**4,-1)
print t-tmin
tdelta[2]=t-tmin
#print tmin
t=10
tmin=4.3*np.sqrt(np.power(10,fit(t))/10**42)*np.power(np.power(10,fit1(t))/10**4,-2)*np.power(fit2(t)*10**3/10**4,-1)
tdelta[3]= t-tmin
print t-tmin
t=1.7
tmin=4.3*np.sqrt(np.power(10,fit(t))/10**42)*np.power(np.power(10,fit1(t))/10**4,-2)*np.power(22*10**3/10**4,-1)
tdelta[4]=t-tmin
tdiff=np.mean(tdelta)

# magR=np.loadtxt("/home/afsari/PycharmProjects/typeIbcAnalysis/Data/2008D/Rmag_2008D.csv",delimiter=",")
# index=np.argsort(magR[:,0])
# magR=magR[index,:]
# magR[:,1]=22-(magR[:,1]-22)
# fit5= UnivariateSpline(magR[magR[:,0]>-0.1,0], magR[magR[:,0]>-0.1,1], s=0.05)
# t_mr=np.arange(0,40,0.1)
# mr=fit5(t_mr)
# ax=plt.subplot(111)
# #plt.plot(t_vph,V_vph,'--',color='green')
# plt.scatter(magR[:,0],magR[:,1],s=20)
# i=np.argmin(mr)
# print mr[i]
# print t_mr[i]
# print fit5(t_mr[i]+15)
# d_mr_15=fit5(t_mr[i]+15)-mr[i]
# plt.plot(t_mr,mr)
# plt.gca().invert_yaxis()
# plt.show()
# print d_mr_15
# t_rise=57.08-71.17*d_mr_15+32.98*d_mr_15**2
# print t_rise
#tdiff=7
Lbol[:,0]=Lbol[:,0]-tdiff
fit= UnivariateSpline(Lbol[:,0], Lbol[:,1], s=0.01)
t_Lbol=np.arange(0,40,0.1)
L_Lbol=fit(t_Lbol)


index=np.argmax(L_Lbol)
L_peak=L_Lbol[index]
t_peak=t_Lbol[index]
ax=plt.subplot(111)
plt.plot(t_Lbol,L_Lbol,'-.',color='green')
plt.plot(Lbol[:,0],Lbol[:,1],color='red')
plt.show()
print L_peak,t_peak
print nickel_mass(t_peak,10**L_peak,9/8)/2e33