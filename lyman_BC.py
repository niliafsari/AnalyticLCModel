import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.time import Time
import csv
import json
from pprint import pprint
import os.path

def lyman_BC(mag_x, mag_y, x,y, mag_x_e=0, mag_y_e=0):
    n=10000
    if (x=='B') & (y=='V'):
        if (mag_x_e==0) & (mag_y_e==0):
            mag_x_sam = mag_x
            mag_y_sam = mag_y
            rms = 0.109
            rms_sam = np.multiply(np.tile(rms, (n, np.shape(mag_y)[1])), np.random.randn(n, np.shape(mag_y)[1]))
        else:
            mag_x_sam=np.multiply(np.tile(mag_x_e,n),np.random.randn(n,np.shape(mag_x)[1])) + np.tile(mag_x,n)
            mag_x_sam = np.multiply(np.tile(mag_y_e,n), np.random.randn(n,np.shape(mag_y)[1])+np.tile(mag_x,n))
            rms=0.109
            rms_sam=np.multiply(np.tile(rms, (n,np.shape(mag_y)[1])), np.random.randn((n,np.shape(mag_y)[1])))
        BC = np.mean(rms_sam+-0.083 - 0.139 * (mag_x_sam -mag_y_sam) - 0.691 * (mag_x_sam -mag_y_sam) ** 2, axis=0)
        BC_e = np.std(rms_sam+-0.083 - 0.139 * (mag_x_sam -mag_y_sam) - 0.691 * (mag_x_sam -mag_y_sam) ** 2, axis=0)
    elif (x=='V') & (y=='I'):
        if (mag_x_e==0) & (mag_y_e==0):
            mag_x_sam = mag_x
            mag_y_sam = mag_y
            rms = 0.090
            rms_sam = np.multiply(np.tile(rms, (n, np.shape(mag_y)[1])), np.random.randn(n, np.shape(mag_y)[1]))
        else:
            mag_x_sam=np.multiply(np.tile(mag_x_e,n),np.random.randn(n,np.shape(mag_x)[1])) + np.tile(mag_x,n)
            mag_x_sam = np.multiply(np.tile(mag_y_e,n), np.random.randn(n,np.shape(mag_y)[1])+np.tile(mag_x,n))
            rms=0.090
            rms_sam=np.multiply(np.tile(rms, (n,np.shape(mag_y)[1])), np.random.randn((n,np.shape(mag_y)[1])))
        BC = np.mean(rms_sam+0.213 - 0.203  * (mag_x_sam -mag_y_sam) - 0.079* (mag_x_sam -mag_y_sam) ** 2, axis=0)
        BC_e = np.std(rms_sam+0.213 -  0.203  * (mag_x_sam -mag_y_sam) - 0.079* (mag_x_sam -mag_y_sam) ** 2, axis=0)
    elif (x == 'V') & (y == 'R'):
        if (mag_x_e==0) & (mag_y_e==0):
            mag_x_sam = mag_x
            mag_y_sam = mag_y
            rms =0.101
            rms_sam = np.multiply(np.tile(rms, (n, np.shape(mag_y)[1])), np.random.randn(n, np.shape(mag_y)[1]))
        else:
            mag_x_sam=np.multiply(np.tile(mag_x_e,n),np.random.randn(n,np.shape(mag_x)[1])) + np.tile(mag_x,n)
            mag_x_sam = np.multiply(np.tile(mag_y_e,n), np.random.randn(n,np.shape(mag_y)[1])+np.tile(mag_x,n))
            rms=0.101
            rms_sam=np.multiply(np.tile(rms, (n,np.shape(mag_y)[1])), np.random.randn((n,np.shape(mag_y)[1])))
        print rms_sam.shape
        BC = np.mean(rms_sam+0.197 -  0.183   * (mag_x_sam -mag_y_sam) - 0.419* (mag_x_sam -mag_y_sam) ** 2, axis=0)
        BC_e = np.std(rms_sam+0.197  -  0.183   * (mag_x_sam -mag_y_sam) -0.419* (mag_x_sam -mag_y_sam) ** 2, axis=0)
    else:
        raise ('error, not a valid band')
    Mbol = BC+ mag_x
    Mbol_e=np.sqrt(BC_e**2+mag_x_e**2)
    Msun = 4.74
    Lsun = 3.84e33
    lbol = Lsun * np.power(10, ((Msun - Mbol) / 2.5))
    lbol_e=lbol*np.log(10)*Mbol_e/2.5
    return lbol,lbol_e

def convertRtor(MV,MV_e,MR,MR_e):
    Mr=np.zeros(shape=np.shape(MV))
    Mr_e = np.zeros(shape=np.shape(MV_e))
    for i,color in enumerate(MV-MR):
        if color<=0.93:
            Mr[i]=MR[i]+0.267*color+0.088
            Mr_e[i]=np.sqrt(0.005**2*color**2+0.267**2*(MV_e[i]**2+MR_e[i]**2)+0.003**2+MR_e[i]**2)
        else:
            Mr[i]=MR[i]+0.77*color-0.37
            Mr_e[i]=np.sqrt(0.04**2*color**2+0.77**2*(MV_e[i]**2+MR_e[i]**2)+0.04**2+MR_e[i]**2)
    return Mr,Mr_e
def convertVtog(MV,MV_e,MR,MR_e,Mr,Mr_e):
    Mg=np.zeros(shape=np.shape(MV))
    Mg_e = np.zeros(shape=np.shape(MV_e))
    Mg=Mr+1.646*(MV-MR)-0.139
    Mg_e=np.sqrt(0.008**2*(MV-MR)**2+1.646**2*(MV_e**2+MR_e**2)+0.004**2+Mr_e**2)
    return Mg,Mg_e
