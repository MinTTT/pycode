# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 14:33:38 2018
use histogram 2d
@author: FuLab
"""
import os
import numpy as np
import api
import matplotlib.pyplot as plt

data_path = r'E:\Pan\毕设\cytometry\20171124\Exp_20171124_1-2'
file_lis = [f.name[:-4] for f in os.scandir(path=data_path) if f.is_file()]       # capture file name
flowdatadic = {x: api.FCSParser(path=data_path+'\\'+str(x)+'.fcs').dataframe
               for x in file_lis}  # open all file and put into a dictionary

da1 = flowdatadic['cl2-g1']
da1 = np.log(da1)
da1 = da1.dropna(axis=0, how='any')
x = da1['FSC-A']
y = da1['FITC-A']
xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()
plt.figure(figsize=(8, 8))
c = plt.hexbin(x, y, gridsize=500, cmap='jet', bins='log')
plt.axis([xmin, xmax, ymin, ymax])
plt.colorbar(c)
plt.show()
