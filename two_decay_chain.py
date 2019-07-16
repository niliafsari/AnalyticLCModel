import numpy as np
from astropy.constants import M_sun

def nickel_mass(t_peak,L_peak,beta):
    return (L_peak*(beta**2)*(t_peak/8.8)**2)/(2*3.9e10*((0.83*(1-beta*t_peak/8.8)*np.exp(-beta*t_peak/8.8))+(26.56*(1-((1+beta*t_peak/111.3)*np.exp(-beta*t_peak/111.3))))))/M_sun.to("g")
