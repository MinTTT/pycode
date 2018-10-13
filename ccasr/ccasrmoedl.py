#-*- coding: utf-8 -*-
"""
Created on Thu Feb 15 21:44:04 2018

@author: pan_c
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def model_solv(y, kr):
    def expression_kindtic(w, t, kp1, kp2, Stot, K, n, leak, kdilm, kdil, y, ktrans, kr, beta):
        R, Rp, mRNA, GFP = w.tolist()
        dR = kr + kp2*Stot*(1-y)*Rp - (kp1*Stot*y+kdil)*R
        dRp = kp1*Stot*y*R - kp2*Stot*(1-y)*Rp - kdil*Rp
        dmRNA = beta*(Rp**n/(K**n+Rp**n))-kdilm*mRNA + leak
        dGFP = ktrans * mRNA - kdil * GFP
        return dR, dRp, dmRNA, dGFP

    args1 = 0.134366, 0.036881, 0.125480, 0.05827, 2, 0.000764, 0.030053, 0.000305, y, 0.004921, kr, 0.015811
    init_state = 0, 0, 0, 0
    t = np.arange(0, 40000, 0.01)
    result = odeint(expression_kindtic, init_state, t, args1)
    return result, t

# CcaS and CcaR expression level modeling
def CcaSR_Modeling(y, kr, Stot):
    def expression_kindtic(w, t, kp1, kp2, Stot, K, n, leak, kdilm, kdil, y, ktrans, kr, beta):
        R, Rp, mRNA, GFP = w.tolist()
        dR = kr + kp2*Stot*(1-y)*Rp - (kp1*Stot*y+kdil)*R
        dRp = kp1*Stot*y*R - kp2*Stot*(1-y)*Rp - kdil*Rp
        dmRNA = beta*(Rp**n/(K**n+Rp**n))-kdilm*mRNA + leak
        dGFP = ktrans * mRNA - kdil * GFP
        return dR, dRp, dmRNA, dGFP
    args1 = 0.134366, 0.036881, Stot, 0.05827, 2, 0.000764, 0.030053, 0.000305, y, 0.004921, kr, 0.015811
    init_state = 0, 0, 0, 0
    t = np.arange(0, 40000, 0.05)
    result = odeint(expression_kindtic, init_state, t, args1)
    return result, t

def full_modeling(y, kr, Stot, dRi, dRpi, dmrnai, dGFPi):
    def expression_kindtic(w, t, kp1, kp2, Stot, K, n, leak, kdilm, kdil, y, ktrans, kr, beta):
        R, Rp, mRNA, GFP = w.tolist()
        dR = kr + kp2*Stot*(1-y)*Rp - (kp1*Stot*y+kdil)*R
        dRp = kp1*Stot*y*R - kp2*Stot*(1-y)*Rp - kdil*Rp
        dmRNA = beta*(Rp**n/(K**n+Rp**n))-kdilm*mRNA + leak
        dGFP = ktrans * mRNA - kdil * GFP
        return dR, dRp, dmRNA, dGFP
    args1 = 0.134366, 0.036881, Stot, 0.05827, 2, 0.000764, 0.030053, 0.000305, y, 0.004921, kr, 0.015811
    init_state = dRi, dRpi, dmrnai, dGFPi
    t = np.arange(0, 40000, 0.05)
    result = odeint(expression_kindtic, init_state, t, args1)
    return result, t



if __name__ == '__main__':
#    # ratio of Sa/Sg
#    y = 1
#    kr = 0.0004  # ccar expression tate

#    #args kp1, kp2, Stot, K, n, leak, kdilm, kdil, y, ktrans, kr, beta
#    args1 = 0.036881, 0.134366, 0.125480, 0.05827, 2, 0.000764, 0.030053, 0.000305, y,      0.004921, kr, 0.015811
#    init_state = 0, 0, 0, 0
#    t = np.arange(0, 20000, 0.005)
#    result1 = odeint(expression_kindtic, init_state, t, args1)
    result1, t1 = model_solv(1, 0.0204)
    plt.figure(1,figsize=(16, 8))
    plt.plot(t1, result1[:, 0], label='R', linewidth=3)
    plt.plot(t1, result1[:, 1], label='Rp', linewidth=3)
    plt.plot(t1, result1[:, 2], label='mRNA', linewidth=3)
    plt.plot(t1, result1[:, 3], label='GFP', linewidth=3)
    plt.legend()
    plt.savefig("CcaRmodel.png")

    result2, t2 = CcaSR_Modeling(1, 0.0204, 0.500)
    plt.figure(2,figsize=(16,8))
    plt.plot(t2, result2[:, 0], label='R', linewidth=3)
    plt.plot(t2, result2[:, 1], label='Rp', linewidth=3)
    plt.plot(t2, result2[:, 2], label='mRNA', linewidth=3)
    plt.plot(t2, result2[:, 3], label='GFP', linewidth=3)
    plt.legend()
    plt.savefig("CcaSRmodel.png")

