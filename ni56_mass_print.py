import numpy as np
from astropy.constants import M_sun
import matplotlib.pyplot as plt
import matplotlib
def nickel_mass(t_peak,L_peak,beta):
    return (L_peak*(beta**2)*(t_peak/8.8)**2)/(2*3.9e10*((0.83*(1-beta*t_peak/8.8)*np.exp(-beta*t_peak/8.8))+(26.56*(1-((1+beta*t_peak/111.3)*np.exp(-beta*t_peak/111.3))))))/M_sun.to("g")


#sn1987a
L_peak=41.85
t_peak=80
beta=0.82
print "M_ni(1987A)=",nickel_mass(t_peak,10**L_peak,beta).value

#sn2008d
L_peak=42.2
t_peak=18.5
beta=9/8
print "M_ni(2008D)=",nickel_mass(t_peak,10**L_peak,beta).value

#sn2016gkg
L_peak=42.3
t_peak=20
beta=0.82
print "M_ni(2016gkg)=",nickel_mass(t_peak,10**L_peak,beta).value


ax=plt.subplot(311)
beta=np.arange(0,5,0.05)
m56=np.zeros(shape=beta.shape)
for i,b in enumerate(beta):
    m56[i]=nickel_mass(t_peak, 10 ** L_peak, b).value
plt.plot(beta, m56, '--', color='blue',label="log(Lp)=42.2, tp=18.5")
plt.xlabel("Beta")
plt.ylabel("M_ni (M_sun)")
plt.legend()

ax1=plt.subplot(312)
#sensitivity to t_peak
t_peak=np.arange(t_peak-4,t_peak+4,0.5)
m56=np.zeros(shape=t_peak.shape)
for i,t in enumerate(t_peak):
    m56[i]=nickel_mass(t, 10 ** L_peak, 0.7).value

plt.plot(t_peak, m56, '--', color='blue',label="beta=0.7, log(Lp)=42.2")
plt.xlabel("t_peak (day)")
plt.ylabel("M_ni (M_sun)")
plt.legend()
ax2=plt.subplot(313)
#sensitivity to L_peak
t_peak=18.5
L_peak=np.arange(L_peak-0.5,L_peak+0.5,0.01)
m56=np.zeros(shape=L_peak.shape)
for i,L in enumerate(L_peak):
    m56[i]=nickel_mass(t_peak, 10 ** L, 0.7).value

plt.plot(L_peak, m56, '--', color='blue',label="beta=0.7, tp=18.5")
plt.xlabel("log(L_peak (erg/s))")
plt.ylabel("M_ni (M_sun)")
plt.legend()
plt.tight_layout()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(5.1, 8, forward=True)
fig.savefig('./Plots/sn2008D.pdf')
plt.show()

