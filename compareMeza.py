import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.ticker import AutoMinorLocator

df=pd.read_csv('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda_26dec.csv')

# print df[df['name']=='SN2008ax'].loc[['peakL']
df_meza=pd.read_csv('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/Meza2020.csv')

Nilou_new=df[df['name'].isin(df_meza['SN'].tolist())].reset_index()
Meza_new=df_meza[df_meza['SN'].isin(Nilou_new['name'].tolist())].reset_index()
#print [float(integer) for integer in Meza_new['L_p'].tolist()]
f=plt.figure(1,figsize=(4,4))
ax=plt.subplot(111)

ll=[float(integer.split(';')[0]) for integer in Nilou_new['peakL'].tolist()]
l=[np.log10(float(integer)*10**41) for integer in Meza_new['L_p'].tolist()]
plt.scatter(np.array(ll),np.array(l))
ax.plot(np.arange(41.5, 43.4, 0.2),np.arange(41.5, 43.4, 0.2),'--',color='black',linewidth=1)
# ax.set_xlim(41.5, 43.4)
# ax.set_ylim(41.5, 43.4)
# ax.set_xscale('log')
# ax.set_yscale('log')
ax.xaxis.set_ticks(np.arange(41.5, 43.4, 0.5))
ax.yaxis.set_ticks(np.arange(41.5, 43.4, 0.5))
plt.xlabel('L_p Nilou')
plt.ylabel('L_p Meza')
ax.set_aspect('equal', adjustable='box')
#ax.ticklabel_format(axis='both',style='plain')
#plt.ticklabel_format(style='plain', axis='both')
print np.median(np.array([float(integer.split(';')[0]) for integer in Nilou_new['arnett_ni_mass'].tolist()]))
print np.median(np.array([float(integer) for integer in Meza_new['Arnett Nickel'].tolist()]))

print np.median([10**float(integer.split(';')[0]) for integer in Nilou_new['peakL'].tolist()])
print np.median([float(integer)*10**41 for integer in Meza_new['L_p'].tolist()] )

f=plt.figure(2,figsize=(4,4))
ax=plt.subplot(111)
plt.scatter([float(integer.split(';')[0]) for integer in Nilou_new['peakt'].tolist()],[float(integer) for integer in Meza_new['t_p'].tolist()] )
ax.plot(np.arange(10, 25, 0.2),np.arange(10, 25, 0.2),'--',color='black',linewidth=1)
ax.set_xlim(10, 25)
ax.set_ylim(10, 25)
plt.xlabel('t_p Nilou')
plt.ylabel('t_p Meza')
ax.set_aspect('equal', adjustable='box')

f=plt.figure(3,figsize=(4,4))
ax=plt.subplot(111)
plt.scatter([float(integer.split(';')[0]) for integer in Nilou_new['arnett_ni_mass'].tolist()],[float(integer) for integer in Meza_new['Arnett Nickel'].tolist()] )
ax.plot(np.arange(0,0.7, 0.02),np.arange(0,0.7, 0.02),'--',color='black',linewidth=1)
ax.set_xlim(0,0.7)
ax.set_ylim(0, 0.7)
plt.xlabel('Arnett Nilou')
plt.ylabel('Arnett Meza')
ax.set_aspect('equal', adjustable='box')

# f=plt.figure(4,figsize=(4,4))
# ax=plt.subplot(111)
print Meza_new['SN']
print zip(Nilou_new['tail_ni_mass'].tolist(), Meza_new['Tail Nickel'].tolist())
# plt.scatter([float(integer.split(';')[0]) for integer in Nilou_new['tail_ni_mass'].tolist()],[float(integer) for integer in Meza_new['Tail Nickel'].tolist()] )
# ax.plot(np.arange(0,0.7, 0.02),np.arange(0,0.7, 0.02),'--',color='black',linewidth=1)
# ax.set_xlim(0,0.7)
# ax.set_ylim(0, 0.7)
# plt.xlabel('Tail Nilou')
# plt.ylabel('Tail Meza')
# ax.set_aspect('equal', adjustable='box')

plt.show()


# print df_meza.describe()