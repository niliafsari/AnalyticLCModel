import urllib
import os
import glob
import subprocess
import commands
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.time import Time
import csv
import json
from pprint import pprint
import os.path
import sys
from scipy.interpolate import UnivariateSpline

def valenti_ni56(dt,Lbol,E51=0,Mej=0,Le=0):
    M_sun=2e33
    c=3e10
    tau_Ni=8.8*86400. # decay time of Ni56 in sec
    tau_Co=9.822e6 #decay time of Co56 in sec
    e_Ni=3.90e10 # erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    if (E51==0.0 or Mej==0.0):
        M_Ni = np.divide(Lbol,(e_Ni*np.exp(-dt*86400./tau_Ni)+ e_Co*(np.exp(-dt*86400./tau_Co) - np.exp(-dt*86400./tau_Ni))))
        M_Ni_e = np.divide(Le, (e_Ni * np.exp(-dt * 86400. / tau_Ni) + e_Co * (
                    np.exp(-dt * 86400. / tau_Co) - np.exp(-dt * 86400. / tau_Ni))))

    else:
        M_ej = Mej * M_sun
        E_51 = E51
        E_K = E51 * 1e51
        F = 32*(M_ej/M_sun)/np.sqrt(E_51)
        G = 16.1*F
        Epsilon =  e_Co*(np.exp(-dt*86400./tau_Co) - np.exp(-dt*86400./tau_Ni))
        S_Ni = e_Ni * np.exp(-dt * 86400. / tau_Ni)
        S_Co1 = 0.81*Epsilon*(1-np.exp(-(F/dt)**2.))
        S_Co2 = 0.164*Epsilon*(1-np.exp(-(F/dt)**2.))*(1-np.exp(-(G/dt)**2.))
        S_Co3 = 0.036*Epsilon*(1-np.exp(-(G/dt)**2.))
        M_Ni = np.divide(Lbol,S_Ni+S_Co1+S_Co2+S_Co3)
        M_Ni_e=np.divide(Le,S_Ni+S_Co1+S_Co2+S_Co3)
    return M_Ni/M_sun, M_Ni_e/M_sun

# ;free parameters
# M_Ni = double(MNi[0])*M_sun



#equations;;;;;;;;;;;;;;;;;

#Nickel source term:

#Cobalt Source terms:
#Epsilon = M_Ni[0]*

# F = 32*(M_ej[0]/M_sun)/sqrt(E_51[0])
# G = 16.1*F
#
# if keyword_set(self) then begin
#    F=Fs[0]
#    G=16.1*F
# endif
#
# S_Co1 = 0.81*Epsilon*(1-exp(-(F/dt)^2.))
# S_Co2 = 0.164*Epsilon*(1-exp(-(F/dt)^2.))*(1-exp(-(G/dt)^2.))
# S_Co3 = 0.036*Epsilon*(1-exp(-(G/dt)^2.))
#
# if keyword_set(fulltrap) then begin
#    S_Co1 = 0.81*Epsilon
#    S_Co2 = 0.164*Epsilon
#    S_Co3 = 0.036*Epsilon
# endif
#
# if keyword_set(KE_only) then begin
#    S_Co1 = 0.00
#    S_Co2 = 0.00
#    S_Co3 = 0.036*Epsilon
# endif