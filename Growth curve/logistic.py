# -*- coding: utf-8 -*-
import numpy as np
from scipy import optimize
def func(x,p):
    """
    (k*np.e**(r*t)/(np.e**(r*t)+k*c))
    """
    k, r, c =p
    return (k*np.e**(r*(x)))/(np.e**(r*(x)+k*c))

def func_error(p,y,x):
    return np.sum((y - func(x, p))**2)

def logisticleastsq(y,x):
    plsq=optimize.basinhopping(func_error,(0.4,0.001,200),
                               niter=50,
                               minimizer_kwargs={"method":"L-BFGS-B",
                                                 "args":(y,x)})
    
    return plsq

def diff(y,p):
    k=p[0]
    r=p[1]
    return r*y*(1-y/k)
