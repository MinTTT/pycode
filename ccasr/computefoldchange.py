
import ccasrmoedl as cc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style

def foldchagne_compute(kr):
    y1 = 0.004
    y2 = 1
    r1,t1 = cc.model_solv(y1,kr)
    r2,t2 = cc.model_solv(y2,kr)
    r1 = (r1+0.1)*400
    r2 = (r2+0.1)*400
    fold = r2[-1,3] / r1[-1,3]
    return [fold, r1[-1,3], r2[-1,3]]
if __name__ == '__main__':
    krseq = np.logspace(0, 3, 200)/800000
    foldseq = [foldchagne_compute(kr) for kr in krseq]
    foldseq = np.array(foldseq)
    style.use('ggplot')
    plt.figure(figsize=(32, 16))
    axe1 = plt.subplot(121)  # Scatter of all samples.
    axe2 = plt.subplot(222)  # Scatter of red light inhibit expression level and fold change
    axe3 = plt.subplot(224)  # Scatter of green light induced expression level and fold change.
    plt.sca(axe1)
    c = plt.hexbin(foldseq[:, 1], foldseq[:, 2], C=foldseq[:, 0], cmap='jet', xscale='log', yscale='log')
    cbar1 = plt.colorbar(c, ax=axe1)
    cbar1.ax.tick_params(labelsize='25')
    cbar1.ax.set_ylabel('Fold-change', fontsize='25')
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    axe1.set_xlabel('R.I.T. (a.u.)', fontsize=25)
    axe1.set_ylabel('G.I.T. (a.u.)', fontsize=25)
    plt.grid(True)
    plt.sca(axe2)
    plt.scatter(foldseq[:, 1], foldseq[:, 0])
    plt.xscale('log')
    axe2.set_xlabel('R.I.T. (a.u.)', fontsize=25)
    axe2.set_ylabel('Fold-change', fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.sca(axe3)
    plt.scatter(foldseq[:, 2], foldseq[:, 0])
    plt.xscale('log')
    axe3.set_xlabel('G.I.T. (a.u.)', fontsize=25)
    axe3.set_ylabel('Fold-change', fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.savefig("simulation.png")
    plt.figure(2, figsize=(32,32))
    plt.scatter(krseq,foldseq[:,0])
    plt.xlabel('CcaR Expression Rate', fontsize=25)
    plt.ylabel('Fold-change', fontsize=25)
    plt.xlim(-0.0001,krseq[-1]+0.0001 )
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.savefig("ccasRfoldchange.png")