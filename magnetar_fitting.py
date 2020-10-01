import numpy as np
from scipy.optimize import fsolve
import math
from valenti_ni56 import *
from scipy.optimize import curve_fit

def MagnetarLightCurve(t, Ep, tp):
    # Generate a light curve from the model of magnetar spin-down
    # Lsd(t) = (Ep/tp) * (1 + t/tp)^(-2)
    # Does not pass the diffusion process, nor leakage
    # t = time (after explosion in rest frame) corresponding to the calculated spin-down luminosity (day)
    # Ep = initial spin-down energy (erg)
    # tp = initial spin-down timescale (day)
    # output = total spin-down luminosity at a given time (erg/s)
    return (Ep / (tp * 24 * 60 * 60.0)) * (1.0 + t / tp) ** (-2)


name='SN2008ax'
tail_data=np.load("Data/" +  name + "_tail_LC.npy")
peak_data=np.load("Data/" + name + "_peak_LC.npy")

#print tail_data,tail_data[:,1]
print np.argwhere(np.isinf(tail_data))
tail_data=np.float64(tail_data)
print np.argwhere(np.isnan(tail_data))
x = tail_data[~np.isnan(tail_data[:,1]),0]
y = tail_data[~np.isnan(tail_data[:,1]),1]

popt, pcov = curve_fit(MagnetarLightCurve, tail_data[:,0],tail_data[:,1],p0=[10**51.0,10.0])
print popt

ax=plt.subplot(111)
plt.plot(tail_data[((tail_data[:,0] >= 60) & (tail_data[:,0] < 120)),0],tail_data[((tail_data[:,0] >= 60) & (tail_data[:,0] < 120)),1],'k',lw=2,label=r'SN $L_{\rm bol}$')
plt.plot(peak_data[:,0],peak_data[:,1],'k',lw=2)
plt.plot(np.arange(60,200),MagnetarLightCurve(np.arange(60,200),popt[0],popt[1]),'r',lw=1,label='magnetar',alpha=0.6)
plt.plot(np.arange(60,200),wygoda_bol(np.arange(60,200),0.06, 99.2/32),'y',lw=1,label='nickel',alpha=0.6)
ax.set_yscale("log")
plt.legend(frameon=False, fontsize=10)
# plt.gca().set_xlim([-5, 150])
# plt.gca().set_ylim([5e40, 3e42])
plt.xlabel(r'Time (days)',fontsize=20)
plt.ylabel(r'$ L_{\rm bol} \ (\rm \ erg \ s^{-1})$',fontsize=20)
plt.savefig('./Plots/'+name+'magnetar_tail.pdf')
plt.show()

