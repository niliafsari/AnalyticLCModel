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


# with open('/home/afsari/PycharmProjects/typeIbcAnalysis/Data/SN1993J_results.csv', mode='r') as csv_file:
#     csv_reader = csv.DictReader(csv_file)
#     line_count = 0
#     header=csv_reader.fieldnames

first=0
with open('Data/SNdata.csv','w') as f:
    for file in glob.glob("Data/*_results.csv"):
        print file
        with open(file, mode='r') as csv_file:
            if first==0:
                csv_reader = csv.DictReader(csv_file)
                w = csv.DictWriter(f, fieldnames=csv_reader.fieldnames)
                w.writeheader()
                first=1
            line_count = 0
            for row in csv_file:
                if (line_count == 1) | ((line_count==0) & (file=='Data/SN2008ax_results.csv')):
                    temp=row.replace(", ",";").\
                        replace("[","").replace("]","").replace('"','').replace("(","").replace(")","")\
                        .replace("'","").replace(';UGC 1238',',UGC 1238')
                    f.write(temp)
                line_count = +1