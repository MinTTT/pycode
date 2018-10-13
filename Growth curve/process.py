# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 20:40:33 2017

@author: FuLab
"""

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import logistic as lg

# from scipy import interpolate
rawdata = pd.read_csv('20180810.csv')    # input data


# set a time serial
locNum = rawdata.iloc[:, 0].size
time_interval = 10
time = np.arange(9, locNum*time_interval, time_interval)    #set a time interval.
cl2 = rawdata[["F2","F3","F4","F5","F6"]]  # experimental group
cr = rawdata[["F7","F8","F9","F10","F11"]]      # blank

cr_cl2 = np.subtract(cl2, cr)
cr_cl2_mean = np.mean(cr_cl2, axis=1)
cr_cl2_std = cr_cl2.std(axis=1)

flg = lg.logisticleastsq(cr_cl2_mean, time)
print("logistic fitting: ", flg.x)
flg_y = lg.func(time, flg.x)
flg_diff_y = lg.diff(flg_y, flg.x)
flg_maxin = flg_diff_y.tolist().index(max(flg_diff_y))



import doubleing_time as dt

doubling_time, cr_cl2_logtime, fity=dt.doutime(time, cr_cl2_mean, flg_maxin)
annotation = ''.join(["Doubling time =", str(round(doubling_time, 3)), "min."])
plt.figure(figsize=(12, 6))



axe1 = plt.subplot(121)
axe2 = plt.subplot(122)
plt.sca(axe1)
plt.errorbar(time, cr_cl2_mean,cr_cl2_std, label="CL4")
plt.plot(cr_cl2_logtime, fity, label="Fitting", linewidth=4,
         color=(253/255, 242/255, 99/255))

plt.plot(time, flg_y, label="Logistic",
         color=(199/255, 65/255, 146/255))

#plt.scatter(time,cr_cl2_mean,alpha=0.5,
#            label="CL3",
#            color=(69/255,188/255,238/255) )

plt.annotate(annotation,(300,0.8))
plt.ylabel("$ABS_{600}$")
plt.xlabel("min")
plt.legend()
plt.sca(axe2)
plt.plot(cr_cl2_logtime, fity, label="Fitting", linewidth=6,
         color=(0/255, 0/255, 0/255))

plt.plot(time,flg_y,label="Logistic",
        color=(199/255, 65/255,146/255))
plt.errorbar(time, cr_cl2_mean, cr_cl2_std, label="CL4")

#plt.scatter(time,cr_cl2_mean,alpha=0.8,
#            label="CL3",
#            color=(69/255,188/255,238/255) )

plt.annotate(annotation,(300,0.8))
plt.ylabel("$ABS_{600}$")
plt.xlabel("min")
plt.legend()
plt.yscale("log")
plt.show()