import numpy as np

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
    print Ep
    v=((Ep+Esn)/(2*Mej*Msolar))**0.5
    print v
    tp=0.5*B**(-2)*P**2*(Mns/1.4)**1.5*(Rns/10)**(-6)
    tLC=((3/(4*np.pi))*Mej*Msolar*kappa/(v*c))**0.5/(3600.0*24)
    print tLC
    # values will be often used in the calculation
    a = (tgrid / tLC) ** 2
    print a
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


name='SN2008ax'
tail_data=np.load("Data/" +  name + "_tail_LC.npy")
peak_data=np.load("Data/" + name + "_peak_LC.npy")
print tail_data.shape
ax=plt.subplot(111)
plt.plot(tail_data[((tail_data[:,0] >= 60) & (tail_data[:,0] < 120)),0],tail_data[((tail_data[:,0] >= 60) & (tail_data[:,0] < 120)),1],'k',lw=2,label=r'SN $L_{\rm bol}$')
plt.plot(peak_data[:,0],peak_data[:,1],'k',lw=2)
t,le,lm,lf=DiffusionMagnetarLightCurve(21,6,0.1,1e51,5,100)
plt.plot(t,lm,'b',label=r'Magentar $M_{\rm ej}=5 M_{\odot}, B_{14}=21, P=6~ms$')
t,le,lm,lf=DiffusionMagnetarLightCurve(24,40,0.1,1e51,2,100)
plt.plot(t,lm,'orange',label=r'Magentar $M_{\rm ej}=2 M_{\odot}, B_{14}=24, P=40~ms$')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.set_yscale("log")
ax.tick_params(direction='in', which='both')
# plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
ax.axhline(y=np.max(peak_data)*0.28, color='r',ls='--',lw=2)
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

