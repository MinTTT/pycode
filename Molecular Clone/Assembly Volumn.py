# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 23:16:02 2018

@author: FuLab

"""
import numpy as np

def VolmnCompute(con, length, molar, total_v, ratio):
    """
    Arguments:
    con: Fragments concentration. ng/ul
    length: the length of each fragments. bp
    molar: final molar of reaction system. fmol
    total_v: volum of reaction of reaction system. ul
    ratio: molar tatio of each fragments. number
    
    """
    Volum_ea = (molar * ratio * (length * 617.96 + 36.04) * 10**-6) / con

    return Volum_ea


print(VolmnCompute.__doc__)
con, length, molar, total_v, ratio =eval( input("Please input con, length, molar, total_v, ratio:"))
Volum = [VolmnCompute(i, length[con.index(i)], molar, total_v, ratio) for i in con ]
if np.sum(Volum) <= total_v :
    Volum[0] = Volum[0]/ratio 
    print(Volum)
else:
    print('concentrations of fragments are too low !')
        

    
