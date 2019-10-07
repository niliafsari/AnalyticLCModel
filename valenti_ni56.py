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
def valenti_bol(dt,M_Ni,MejE):
    M_sun=2e33
    M_Ni=M_Ni*M_sun
    c=3e10
    tau_Ni=8.8*86400. # decay time of Ni56 in sec
    tau_Co=9.822e6 #decay time of Co56 in sec
    e_Ni=3.90e10 # erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    if (MejE==0.0):
        Lbol =M_Ni*(e_Ni*np.exp(-dt*86400./tau_Ni)+ e_Co*(np.exp(-dt*86400./tau_Co) - np.exp(-dt*86400./tau_Ni)))
    else:
        F = 32*MejE
        G = 16.1*F
        Epsilon =  e_Co*(np.exp(-dt*86400./tau_Co) - np.exp(-dt*86400./tau_Ni))
        S_Ni = e_Ni * np.exp(-dt * 86400. / tau_Ni)
        S_Co1 = 0.81*Epsilon*(1-np.exp(-(F/dt)**2.))
        S_Co2 = 0.164*Epsilon*(1-np.exp(-(F/dt)**2.))*(1-np.exp(-(G/dt)**2.))
        S_Co3 = 0.036*Epsilon*(1-np.exp(-(G/dt)**2.))
        Lbol = M_Ni*(S_Ni+S_Co1+S_Co2+S_Co3)
    return Lbol

def valentiErr(p, t, L, L_err):
    err = (valenti_bol(t, p[0],p[1]) - L)/L_err
    return np.abs(err)

# ;free parameters
# M_Ni = double(MNi[0])*M_sun
def wheeler_ni56(dt,Lbol,E51=0,Mej=0,Le=0):
    M_sun=2e33
    c=3e10
    tau_Ni=8.8*86400. # decay time of Ni56 in sec
    tau_Co=111.3*86400 #decay time of Co56 in sec
    e_Ni=3.90e10 # erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    M_ej = Mej * M_sun
    E_51 = E51
    E_K = E51 * 1e51
    F = 32*(M_ej/M_sun)/np.sqrt(E_51)
    S=((e_Ni-e_Co)*np.exp(-dt*86400./tau_Ni) +e_Co*np.exp(-dt*86400./tau_Co) )
    S_gamma = (1-np.exp(-(F/dt)**2.))
    M_Ni = np.divide(Lbol,np.multiply(S,S_gamma))
    M_Ni_e=np.divide(Le,np.multiply(S,S_gamma))
    return M_Ni/M_sun, M_Ni_e/M_sun

def wheeler_bol(dt,M_Ni,MejE):
    M_sun=2e33
    M_Ni=M_Ni*M_sun
    c=3e10
    tau_Ni=8.8*86400. # decay time of Ni56 in sec
    tau_Co=111.3*86400  #decay time of Co56 in sec
    e_Ni=3.90e10 # erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    F = 32*MejE
    S=((e_Ni-e_Co)*np.exp(-dt*86400./tau_Ni) +e_Co*np.exp(-dt*86400./tau_Co))
    S_gamma = (1-np.exp(-(F/dt)**2.))
    Lbol=M_Ni*np.multiply(S,S_gamma)
    return Lbol

def wheelerErr(p, t, L, L_err):
    err = (wheeler_bol(t, p[0],p[1]) - L)/L_err
    return np.abs(err)

def nickel_mass_khatami_perturb(t_peak,L_peak,L_peak_e,beta):
    n=1000
    L= np.multiply(np.tile(L_peak_e,n), np.random.randn(n)+np.tile(L_peak,n))
    t_peak_t=np.tile(t_peak, n)
    MNi=np.divide(np.multiply(L_peak*(beta**2),(t_peak_t/8.8)**2),(2*3.9e10*((0.83*np.multiply((1-beta*t_peak_t/8.8),np.exp(-beta*t_peak_t/8.8)))+(26.56*(1-(np.multiply((1+beta*t_peak_t/111.3),np.exp(-beta*t_peak_t/111.3))))))))/M_sun.to("g").value
    return np.mean(MNi,axis=0), np.std(MNi,axis=0)

def nickel_mass_khatami(t_peak,L_peak,L_peak_err,beta):
    e_Ni=3.90e10 # erg/s/g energy produced by 1 gram of Ni
    e_Co=6.78e9 #erg/s/g energy produced by 1 gram of Co
    tau_Ni=8.8*86400. # decay time of Ni56 in sec
    tau_Co = 111.3 * 86400  #decay time of Co56 in sec
    M_sun = 2e33
    MNi=np.divide(L_peak*(beta**2)*(t_peak/8.8)**2,(2*e_Ni*(((1-(e_Co/e_Ni))*(1-np.multiply((1+beta*t_peak/8.8),np.exp(-beta*t_peak/8.8))))
                                                            +((e_Co*tau_Co**2/(e_Ni*tau_Ni**2))*(1-(np.multiply((1+beta*t_peak/111.3),np.exp(-beta*t_peak/111.3))))))))/M_sun
    MNi_err=np.divide(np.multiply(L_peak_err*(beta**2),(t_peak/8.8)**2),(2*e_Ni*(((1-(e_Co/e_Ni))*(1-np.multiply((1+beta*t_peak/8.8),np.exp(-beta*t_peak/8.8))))+((e_Co*tau_Co**2/(e_Ni*tau_Ni**2))*(1-(np.multiply((1+beta*t_peak/111.3),np.exp(-beta*t_peak/111.3))))))))/M_sun
    return MNi,MNi_err

def nickel_mass_khatami_err(beta,t_peak,L_peak,L_peak_err,Mni,Mni_err):
    err = (nickel_mass_khatami(t_peak,L_peak,L_peak_err,beta)[0] - Mni)/ np.sqrt(Mni_err**2+nickel_mass_khatami(t_peak,L_peak,L_peak_err,beta)[1]**2)
    return np.abs(err)

def khatami_err( beta,x,L_peak,L_peak_err):
    err = (nickel_mass_khatami(x[0], L_peak, L_peak_err, beta)[0] - x[1]) / nickel_mass_khatami(x[0], L_peak, L_peak_err, beta)[1]
    return np.abs(err)


def power_law(t, t0, a, n):
    return a * np.sign(t - t0) * (np.abs(t - t0)) ** n
def powerLawErr(p, t, f, f_err):
    err = (power_law(t, p[0],p[1],p[2]) - f)/1
    return np.abs(err)

def benzine(t, tfall,trise, A,B):
    t0=0
    return A*(np.exp(-(t-t0)/tfall)/(np.exp((t-t0)/trise)+1))+B

def benzineErr1(p, t, L, L_err):
    err = (benzine(t, p[0],p[1],p[2], p[3]) - L)/L_err
    return np.abs(err)
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