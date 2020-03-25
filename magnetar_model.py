import numpy as np
from scipy.optimize import fsolve
import math



def MagnetarLightCurve(t, Ep, tp):
    # Generate a light curve from the model of magnetar spin-down
    # Lsd(t) = (Ep/tp) * (1 + t/tp)^(-2)
    # Does not pass the diffusion process, nor leakage
    # t = time (after explosion in rest frame) corresponding to the calculated spin-down luminosity (day)
    # Ep = initial spin-down energy (erg)
    # tp = initial spin-down timescale (day)
    # output = total spin-down luminosity at a given time (erg/s)
    return (Ep / (tp * 24 * 60 * 60.0)) * (1.0 + t / tp) ** (-2)


def DiffusionMagnetarLightCurve(B, P, kappa, Esn, Mej,tmax):
    # Generate a light curve from the diffusion approximation with magnetar spin-down engine from 1 - tmax days
    # This calls function MagnetarLightCurve
    # By default, the function assumes small initial radius
    # Time step for integration is 1 day
    # Edge effect: set zero at 0 and tmax+1 days
    # Ep = initial spin-down energy (erg)
    # tp = spin-down timescale (day)
    # tmax = maximum time (after explosion in rest frame) corresponding to the calculated output luminosity (day)
    # Esn = Explosion energy (erg)
    # tLC = effective light curve timescale (i.e., = sqrt(2 * tDiff * tExp)) (day)
    # tDiff = diffusion timescale = tLC (default) (day)
    # r0 = 0 (default for assuming small initial radius) (cm)
    # v = ejecta velocity (km/s)
    # return(time grid in days, spin-down input luminosity in erg/s, spin-down diffused output luminosity in erg/s, fireball luminosity in erg/s)
    Msolar=2e33
    c=3e10
    Mns=1.4
    Rns=10
    tgrid = np.arange(0,tmax)  # generate time grids (day)
    Ep=2.5e52*P**(-2)*(Mns/1.4)**1.5
    v=((Ep+Esn)/(2*Mej*Msolar))**0.5
    tp=0.5*B**(-2)*P**2*(Mns/1.4)**1.5*(Rns/10)**(-6)
    tLC=((3/(4*np.pi))*Mej*Msolar*kappa/(v*c))**0.5/(3600.0*24)
    # values will be often used in the calculation
    a = (tgrid / tLC) ** 2
    b = tgrid / tLC
    c = np.exp(-a)
    d = np.exp(+a)
    ####

    LoutFireball = (Esn / (tLC * 24.0 * 60 * 60)) * c  # fireball term

    # magnetar term
    Linp = MagnetarLightCurve(tgrid, Ep, tp)
    integrand = np.multiply(np.multiply(d, b), Linp)
    integral = np.zeros(np.shape(tgrid))
    LoutMagnetar = np.zeros(np.shape(tgrid))
    for i in tgrid[1:]:
        integral[i] = 0.5 * (integrand[i-1] + integrand[i]) * 1.0
        LoutMagnetar[i] = LoutMagnetar[i-1] + integral[i]
    LoutMagnetar = (2.0 / tLC) * np.multiply(c, LoutMagnetar)
    return  tgrid[1:], Linp[1:], LoutMagnetar[1:], LoutFireball[1:]
from scipy.optimize import least_squares

def magnetar_eq(p, *mag_properties):
    Lpeak, tpeak, Mej=mag_properties
    P, B = p
    Esn=1e51
    kappa=0.1
    Msolar=2e33
    c=3e10
    Mns=1.4
    Rns=10
    Em=2.5e52*P**(-2)*(Mns/1.4)**1.5
    tm=0.5*B**(-2)*P**2*(Mns/1.4)**1.5*(Rns/10)**(-6)*24*3600
    v = ((Em + Esn) / (2 * Mej * Msolar)) ** 0.5
    tsn=((3/(4*np.pi))*Mej*Msolar*kappa/(v*c))**0.5
    #print Lpeak, tpeak,(tm*(np.sqrt(Em/(Lpeak*tm))-1)), 1.5*((Em*tm)/tsn**2)*(np.log(1+(tsn/tm))-(tsn/(tsn+tm)))
    # print 1.5*((Em*tm)/tsn**2)*(np.log(1+(tsn/tm))-(tsn/(tsn+tm)))
    return (Lpeak-1.5*((Em*tm)/tsn**2)*(np.log(1+(tsn/tm))-(tsn/(tsn+tm))), tpeak-(tm*(np.sqrt(Em/(Lpeak*tm))-1)))

import pandas as pd
df=pd.read_csv('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SNdata_wygoda_26dec.csv')

# print df[df['name']=='SN2008ax'].loc[['peakL']
df_lfrac=pd.read_csv('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/lfrac.dat')
Mej=2
df_save = pd.DataFrame(columns=['name', 'p', 'b','f'])
rows_list = []
lpeak_add_list=[]
for index, row in df.iterrows():
    if df_lfrac[df_lfrac['SN_name']==row['name']].iloc[0,1]>0:
        mag_properties= ((10**float(row['peakL'].split(';')[0]))*df_lfrac[df_lfrac['SN_name']==row['name']].iloc[0,1],float(row['peakt'].split(';')[0])*(3600.0*24),Mej)
        lpeak_add_list.append((10 ** float(row['peakL'].split(';')[0]))* (df_lfrac[df_lfrac['SN_name'] == row['name']].iloc[0, 1]))
        #print mag_properties, (10**float(row['peakL'].split(';')[0]))*df_lfrac[df_lfrac['SN_name']==row['name']].iloc[0,1]
        x, y =  fsolve(magnetar_eq, (10, 10),args=mag_properties)
        #res = least_squares(magnetar_eq, (40, 20), bounds=((1, 1), (100, 100)),args=mag_properties,max_nfev=100000)
        #p,b=res.x
        #print row['name'],x,y, (10**float(row['peakL'].split(';')[0])), df_lfrac[df_lfrac['SN_name']==row['name']].iloc[0,1],(10**float(row['peakL'].split(';')[0]))*df_lfrac[df_lfrac['SN_name']==row['name']].iloc[0,1]
        df_save.loc[index]=[row['name'],x,y, df_lfrac[df_lfrac['SN_name']==row['name']].iloc[0,1]]
        print row['name'], df_lfrac[df_lfrac['SN_name']==row['name']].iloc[0,1], x, y
df_save.to_csv('Data/magnetar_2.csv',index=False)


mag_5=pd.read_csv('Data/magnetar_5.csv')
mag_2=pd.read_csv('Data/magnetar_2.csv')

import matplotlib.pyplot as plt
import matplotlib
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
matplotlib.rcParams.update({'font.size': 18})
plt.rc('text', usetex=True)

#
name='SN2008ax'


import matplotlib.gridspec as gridspec
f = plt.figure(1,figsize=(6,9))
gs = gridspec.GridSpec(ncols=2, nrows=4)
gs.update(wspace=0.15, hspace=0.5)
ax=f.add_subplot(gs[0, 0])
plt.hist(mag_2['p'].tolist(),bins=15, color='orange', ec='none',label=r'$M_{\rm ej}=2 M_{\odot}$')
ax.set_xlabel(r'$P \ (\rm ms)$',fontsize=15)
ax.set_ylabel('Count',fontsize=15)
leg=ax.legend(frameon=False, fontsize=11)
for item in leg.legendHandles:
    item.set_visible(False)
ax=f.add_subplot(gs[0, 1])
plt.hist(mag_5['p'].tolist(),bins=15, color='orange', ec='none',label=r'$M_{\rm ej}=5 M_{\odot}$')
ax.set_xlabel(r'$P \ (\rm ms)$',fontsize=15)
# plt.gca().set_xlim([0,10])
leg=ax.legend(frameon=False, fontsize=11)
for item in leg.legendHandles:
    item.set_visible(False)
ax=f.add_subplot(gs[1, 0])
plt.hist(mag_2['b'].tolist(),bins=15, color='salmon', ec='none',label=r'$M_{\rm ej}=2 M_{\odot}$')
ax.set_xlabel(r'$B_{14}$ ',fontsize=15)
ax.set_ylabel('Count',fontsize=15)
# ax.set_xlim(0,30)
#plt.gca().set_xlim([0,70])
leg=ax.legend(frameon=False, fontsize=11)
for item in leg.legendHandles:
    item.set_visible(False)
ax=f.add_subplot(gs[1,1])
plt.hist(mag_5['b'].tolist(),bins=15, color='salmon', ec='none',label=r'$M_{\rm ej}=5 M_{\odot}$')
ax.set_xlabel(r'$B_{14}$ ',fontsize=15)
leg=ax.legend(frameon=False, fontsize=11)
# plt.gca().set_xlim([0,45])
for item in leg.legendHandles:
    item.set_visible(False)

ax=f.add_subplot(gs[2,:])
plt.hist(df_lfrac['lfac'].tolist(),bins=15, color='goldenrod', ec='none')
ax.set_xlabel(r'$f$ ',fontsize=15)
ax.set_ylabel('Count',fontsize=15)

ax=f.add_subplot(gs[3,:])
plt.hist(np.log10(lpeak_add_list),bins=15, color='gray', ec='none')
ax.set_xlabel(r'Log $f L_{\rm p} \ (\rm erg \ s^{-1})$ ',fontsize=15)
ax.set_ylabel('Count',fontsize=15)
# ax.set_ylim(0,6)
plt.tight_layout()
plt.savefig('./Plots/magnetars.pdf')


f=plt.figure(2,figsize=(6,4))
tail_data=np.load("Data/" +  name + "_tail_LC.npy")
peak_data=np.load("Data/" + name + "_peak_LC.npy")
print tail_data.shape
ax=plt.subplot(111)
plt.plot(tail_data[((tail_data[:,0] >= 60) & (tail_data[:,0] < 120)),0],tail_data[((tail_data[:,0] >= 60) & (tail_data[:,0] < 120)),1],'k',lw=2,label=r'SN $L_{\rm bol}$')
plt.plot(peak_data[:,0],peak_data[:,1],'k',lw=2)
t,le,lm,lf=DiffusionMagnetarLightCurve(mag_5[mag_5['name']==name].iloc[0,2],mag_5[mag_5['name']==name].iloc[0,1],0.1,1e51,5,100)
plt.plot(t,1.35*lm,'b',label=r'Magnetar $M_{\rm ej}=5 M_{\odot}, B_{14}=23, P=9~ms$',lw=2)
t,le,lm,lf=DiffusionMagnetarLightCurve(mag_2[mag_2['name']==name].iloc[0,2],mag_2[mag_2['name']==name].iloc[0,1],0.1,1e51,2,100)
plt.plot(t,1.35*lm,'orange',label=r'Magnetar $M_{\rm ej}=2 M_{\odot}, B_{14}=17, P=62~ms$',lw=2)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_yscale("log")
ax.tick_params(direction='in', which='both')
# plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
print np.max(lm)
ax.axhline(y=5.9e41, color='r',ls='--',lw=2)
# plt.annotate(r'$Additional peak power n$',
#              xy=(6.0 / 140, 0.4),
#              xycoords='axes fraction',
#              textcoords='offset points', fontsize=18)
plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
plt.legend(frameon=False, fontsize=11)
plt.gca().set_xlim([-5, 140])
plt.gca().set_ylim([5e40, 3e42])
plt.xlabel(r'Time (days)',fontsize=20)
plt.ylabel(r'$ L_{\rm bol} \ (\rm \ erg \ s^{-1})$',fontsize=20)
plt.tight_layout()
plt.savefig('./Plots/'+name+'_magnetar.pdf')
plt.show()

