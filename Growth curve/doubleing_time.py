# -*- coding: utf-8 -*-
'''
This module need three parameters
x: time sequence, y: ABS or OD
max: logmarthtic midpoint
'''
import numpy as np
from scipy import optimize



def doutime(x,y,max):
   logtime=x[np.arange(max-38,max-10,1).tolist()]
   logod=y.iloc[np.arange(max-38,max-10,1).tolist()]
   
   def residuals(p):
    k,b=p
    return np.log(logod) - (k*logtime+b)

   r=optimize.leastsq(residuals,[0,-10])
   k,b=r[0]
   k=np.log(2)/k
   b=np.e**b
   print("Doubling time =",k,"min.")
   fity=b*(2**(logtime/k))
   return k,logtime,fity
    

