# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import style


class LibSample:
    def _init_(self, name, greenindu, redindu):
        self.name = name
        self.greenindu = greenindu
        self.redindu = redindu

    def fold_change(self):
        return self.greenindu / self.redindu


if __name__ == '__main__':
    data = pd.read_csv("data1.csv")
    data = data.iloc[:84, :4]
    data.iloc[40:84, 1:3] = data.iloc[40:84, 1:3]/3
    data.iloc[:, 3] = (data.iloc[:, 1])/(data.iloc[:, 2])  # deduct background?

    data.columns = ['tube_name', 'Green', 'Red', 'Fold_Change']
    style.use('ggplot')

    selected = data.loc[data.Fold_Change>8.9, "tube_name"]
    print(selected)
    plt.figure(figsize=(16, 8))
    axe1 = plt.subplot(121)  # Scatter of all samples.
    axe2 = plt.subplot(222)  # Scatter of red light inhibit expression level and fold change
    axe3 = plt.subplot(224)  # Scatter of green light induced expression level and fold change.
    plt.sca(axe1)
    c = plt.hexbin(data.iloc[:, 2], data.iloc[:, 1], C=data.iloc[:, 3], cmap='jet', xscale='log', yscale='log')
    cbar1 = plt.colorbar(c, ax=axe1)
    cbar1.ax.tick_params(labelsize='25')
    cbar1.ax.set_ylabel('Fold-change',fontsize='25')
    plt.xlim(80, data.iloc[:, 2].max()+100)
    plt.ylim(50, data.iloc[:, 1].max()+1000)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    axe1.set_xlabel('R.I.T.', fontsize=25)
    axe1.set_ylabel('G.I.T', fontsize=25)
    axe1.annotate("Origin \n RBS", fontsize='20',
                  xy=(data.iloc[83,2],data.iloc[83,1]), xycoords='data',
                  xytext=(data.iloc[83,2]-1000,data.iloc[83,1]-2500),textcoords="data",
                  arrowprops=dict(arrowstyle="fancy",
                                  color="0.5",
                                  shrinkB=5,
                                  connectionstyle="arc3,rad=0.3",),
                  )

    plt.grid(True)
    plt.sca(axe2)
    plt.scatter(data.iloc[:, 2], data.iloc[:, 3])
    plt.xscale('log')
    axe2.set_xlabel('R.I.T.', fontsize=25)
    axe2.set_ylabel('Fold-change', fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.sca(axe3)
    plt.scatter(data.iloc[:, 1], data.iloc[:, 3])
    plt.xscale('log')
    axe3.set_xlabel('G.I.T.', fontsize=25)
    axe3.set_ylabel('Fold-change', fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.figure(2)
    plt.hist(data.iloc[:, 3], bins=10, histtype='stepfilled')
    plt.xlabel('Fold-change', fontsize=25)
    plt.ylabel('Number', fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.show()
