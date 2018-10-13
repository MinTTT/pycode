import multiprocessing as mp
import CcaSR_Modeling as cSRM
import numpy as np
import matplotlib.pyplot as plt


def mp_pro_SR(list):
    acc = 50
    krseq = np.logspace(0, 3, acc) / 500000
    Stotseq = np.logspace(0, 3, acc) / 1800
    parmax = np.array([[[Stot, krseq] for krseq in krseq] for Stot in Stotseq])
    KSCom_Kr = parmax[:, :, 1].reshape(acc * acc)
    KSCom_Stot = parmax[:, :, 0].reshape(acc * acc)
    #name = mp.current_process().name
    #print(name)
    res = cSRM.ccasrFC_Com(KSCom_Stot[list], KSCom_Kr[list])
    #print("finish: ", name)
    return res



# generate parameters sequences
acc = 50
# multiprocessing
inputs = list(range(acc*acc))
#print(inputs)
pro_pool = mp.Pool(processes=14)
pro_pool_output = pro_pool.map(mp_pro_SR, inputs)
pro_pool.close()
pro_pool.join()
pro_pool_output = np.array(pro_pool_output)
FC_mat = pro_pool_output[:, 0].reshape(acc, acc)
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
plt.savefig("ccaSR-contourf-mp.png")