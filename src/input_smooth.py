#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 14:17:22 2021

@author: sunyuexi
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import UnivariateSpline
from scipy.signal import savgol_filter

PT_data=[]

with open('MR.txt') as f:
    for line in f:

        PT_data.append([elt.strip() for elt in line.split(',')])
PT_data = np.asfarray(PT_data,float)
PT=PT_data[:,1]

new_length=250

t_step=np.linspace(0,1,np.shape(PT)[0])
t_new=np.linspace(0,1,new_length)

f_PT=interpolate.interp1d(t_step,PT)
PT=f_PT(t_new)

PT_smooth = savgol_filter(PT, 79, 7)  #79 is the resize array number, must set in odd number and less or equal to current array size

plt.figure(3)
plt.plot(PT,label='original')
plt.plot(PT_smooth,label='smooth')
plt.title("compare smoothing graph")
plt.legend()

plt.savefig("output.png")
# print(np.shape(PT_smooth))

np.savetxt("smooth.txt",PT_smooth,delimiter=",")