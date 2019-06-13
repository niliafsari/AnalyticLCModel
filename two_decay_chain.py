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
    return L_peak*beta**2*(t_peak/8.8)**2/(2*3.9e10*(0.83*(1-beta)*np.exp(-beta*t_peak/8.8))/(26.56*(1-(1+beta*t_peak/111.3))*np.exp(-beta*t_peak/111.3)))
