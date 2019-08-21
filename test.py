import csv
import requests
import numpy as np
import urllib
import os
import glob
import subprocess
import commands
import matplotlib.pyplot as plt
import matplotlib
from astropy.time import Time
import csv
import json
from pprint import pprint
import os.path
import sys
from scipy.interpolate import UnivariateSpline
from valenti_ni56 import *
from error_interpolation import *
from astropy.constants import M_sun
#from ni56_mass_print import *
#I = r - 1.2444*(r - i) - 0.3820 lupton 2005
sys.path.insert(0, '/home/afsari/')
from SNAP5.Analysis import *
from lyman_BC import  *

print valenti_ni56(13.461,10**42.299)
