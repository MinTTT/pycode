import ccasrmoedl as cc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import multiprocessing as mp

# definite a function that process the response time
def respon_time(kr, Stot):
    r1, t1 = cc.CcaSR_Modeling(1.00, kr, Stot)
    gfplevel = r1[-1, 3]
    halflevel = gfplevel / 2.0 # half expression level
    decatlevel = (np.array(r1[:, 3]) - halflevel)**2
    Index_of_response = np.argmin(decatlevel)
    respontime = t1[Index_of_response]
    respontime = 1/ np.log10(respontime)
    return respontime

# using multiprocessing
def multpor_resptime(list):
    acc = 100
    krseq = np.logspace(0, 3, acc) / 500000
    Stotseq = np.logspace(0, 3, acc) / 1800
    parmax = np.array([[[Stot, krseq] for krseq in krseq] for Stot in Stotseq])
    KSCom_Kr = parmax[:, :, 1].reshape(acc * acc)
    KSCom_Stot = parmax[:, :, 0].reshape(acc * acc)
    res_time = respon_time(KSCom_Kr[list], KSCom_Stot[list])
    return res_time

if __name__ == '__main__':
    acc = 100
    # multiprocessing
    inputs = list(range(acc*acc))
    #print(inputs)
    pro_pool = mp.Pool(processes=14)
    pro_pool_output = pro_pool.map(multpor_resptime, inputs)
    pro_pool.close()
    pro_pool.join()
    pro_pool_output = np.array(pro_pool_output)
    FC_mat = pro_pool_output.reshape(acc, acc)
    #print(FC_mat)
    # draw
    krseq = np.logspace(0, 3, acc) / 500000
    Stotseq = np.logspace(0, 3, acc) / 1800
    parmax = np.array([[[Stot, krseq] for krseq in krseq] for Stot in Stotseq])
    plt.figure(1,figsize=(8,8))
    aex = plt.gca()
    fig1 = aex.contourf(parmax[:,:,1], parmax[:,:,0], FC_mat, acc, cmap='jet' )
    plt.colorbar(fig1)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("ccaSR-responsetime.png")






