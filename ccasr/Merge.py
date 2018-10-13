import ccasrmoedl as cc
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import responsetime as RG
import shiftfromgreen as GR
import CcaSR_Modeling as cSRM
from matplotlib import cm


def mergework(list):
    acc =40
    krseq = np.logspace(0, 3, acc) / 500000
    Stotseq = np.logspace(0, 3, acc) / 1800
    parmax = np.array([[[Stot, krseq] for krseq in krseq] for Stot in Stotseq])
    KSCom_Kr = parmax[:, :, 1].reshape(acc * acc)
    KSCom_Stot = parmax[:, :, 0].reshape(acc * acc)
    res_fc, res_red, res_green = cSRM.ccasrFC_Com(KSCom_Stot[list], KSCom_Kr[list])
    res_RtoG = RG.respon_time(KSCom_Kr[list], KSCom_Stot[list])
    res_GtoR = GR.shiftfromgreen(KSCom_Kr[list], KSCom_Stot[list])
    return  res_fc, res_RtoG, res_GtoR

acc =40
inputs = list(range(acc * acc))
pro_pool = mp.Pool(processes=14)
pro_pool_output = pro_pool.map(mergework, inputs)
pro_pool.close()
pro_pool.join()
pro_pool_output= np.array(pro_pool_output) #

## merge contour fig
krseq = np.logspace(0, 3, acc) / 500000
Stotseq = np.logspace(0, 3, acc) / 1800
parmax = np.array([[[Stot, krseq] for krseq in krseq] for Stot in Stotseq])
plt.figure(1,figsize=(8,8))
fig_fc = plt.contour(parmax[:,:,1], parmax[:,:,0], pro_pool_output[:, 0].reshape(acc,acc), 6.00, linewidths=4, linestyles='solid', colors='r')
fig_rg = plt.contour(parmax[:,:,1], parmax[:,:,0], pro_pool_output[:, 1].reshape(acc,acc), 0.2710, linewidths=4, linestyles='dashed', colors='y')
fig_gr = plt.contour(parmax[:,:,1], parmax[:,:,0], pro_pool_output[:, 2].reshape(acc,acc), 0.2778, linewidths=4, linestyles='dotted', colors='b')
background = plt.contourf(parmax[:,:,1], parmax[:,:,0], pro_pool_output[:, 0].reshape(acc,acc), acc, cmap=cm.PRGn)
plt.colorbar(background)
plt.xscale('log')
plt.yscale('log')
plt.savefig('Merge.png')
