import ccasrmoedl as cc
import numpy as np
import matplotlib.pyplot as plt

def ccasrFC_Com(Stot, Kr):
    y1 = 0.04
    y2 = 1
    r1, t1 = cc.CcaSR_Modeling(y1, Kr, Stot)
    r2, t2 = cc.CcaSR_Modeling(y2, Kr, Stot)
    r1 = (r1+0.1)*400
    r2 = (r2+0.1)*400
    fold = r2[-1,3] / r1[-1,3]
    # out put : fold_change, gfp in dark, gfp in green, time line
    return fold, r1[-1,3], r2[-1,3]

CcasrFC_Com_ufunc = np.frompyfunc(ccasrFC_Com, 2, 3)




if __name__ == '__main__':
    acc = 2
    krseq = np.logspace(0, 3, acc)/800000
    Stotseq = np.logspace(0, 3, acc)/1800
    parmax = np.array([[[Stot, krseq] for krseq in krseq]for Stot in Stotseq])
    KSCom_Kr = parmax[:, :, 1].reshape(acc*acc)
    KSCom_Stot = parmax[:, :, 0].reshape(acc*acc)
    results_of_foldC, re_of_dar_f, re_of_dar_g = CcasrFC_Com_ufunc(KSCom_Stot, KSCom_Kr)
    # fig. 1
    plt.figure(1, figsize=(8,8))
    c = plt.hexbin(KSCom_Kr, KSCom_Stot, C=results_of_foldC, cmap='inferno', xscale='log', yscale='log', gridsize=(acc-1, acc-1))
    cbar1 = plt.colorbar(c)
    plt.savefig("ccaSR-bexbin.png")
    #fig.2
    plt.figure(2,figsize=(8,8))
    aex = plt.gca()
    X = results_of_foldC.reshape(acc,acc)
    X = X.astype(np.float)
    heat = aex.imshow(X)
    aex.set_xticks(np.arange(len(krseq)))
    aex.set_yticks(np.arange(len(Stotseq)))
    aex.set_xticklabels(krseq)
    aex.set_yticklabels(Stotseq)
    plt.savefig("ccasR-hist2d.png")
    # fig.  1
    plt.figure(3,figsize=(8,8))
    fig1 = plt.contourf(KSCom_Kr.reshape(acc,acc), KSCom_Stot.reshape(acc,acc), X, 20, cmap='jet')
    plt.colorbar(fig1)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("ccaSR-contourf.png")
