import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import statsmodels.api as sm

qpdata_df = pd.read_csv("admin_2018-08-08 02-53-23_BR005256 -  Quantification Summary_0.csv")

def Sample_Statistic(Sample_Name):
    """
    This function is for reading the data in csv file which was export from Bio-Rad CFX Manager
    :param Sample_Name: input the sample name
    :return: a list contain cq-mean, cq-std, SQ-mean, cq-mean
    """
    qpdata_sd_df = qpdata_df.iloc[qpdata_df.Content.values == Sample_Name]
    sample_cq = qpdata_sd_df.Cq
    sample_SQ = qpdata_sd_df.SQ
    sample_cq_Mean = np.mean(sample_cq)
    sample_cq_Std = np.std(sample_cq)
    sample_SQ_Mean = np.mean(sample_SQ)
    sample_SQ_Std = np.std(sample_SQ)
    return np.array([sample_cq_Mean, sample_SQ_Mean, sample_cq_Std, sample_SQ_Std])

Stand_Curve_Sample = ['Std-01', 'Std-02', 'Std-03', 'Std-04', 'Std-05', 'Std-06', \
                      'Std-07', 'Std-08']
Stand_Curve_Values = np.array([Sample_Statistic(i) for i in Stand_Curve_Sample])

Untreat_Gro = [  'Std-04', 'Std-05', 'Std-07']
Treated_Gro = [  'Unkn-05', 'Unkn-07', 'Unkn-11']
Treated_Plus_Gro = [ 'Unkn-06', 'Unkn-08', 'Unkn-12']
ValuofUntre = np.array([Sample_Statistic(i) for i in Untreat_Gro])
ValuofTreatedGro = np.array([Sample_Statistic(i) for i in Treated_Gro])
ValueofTreatedGroPlus = np.array([Sample_Statistic(i) for i in Treated_Plus_Gro])
index = np.arange(len(Untreat_Gro))
width_of_bar = 0.24





x = np.log10(Stand_Curve_Values[:,1])
x = x[~np.isnan(x)]
x = sm.add_constant(x)
y = Stand_Curve_Values[:,0]
y = y[~np.isnan(y)]
model = sm.OLS(y, x)
results = model.fit()
r_squared = results.rsquared
cont = results.params[0]
slope = results.params[1]
y_predict = results.fittedvalues



plt.figure(figsize=(10,4),dpi=300)
axe1 = plt.subplot(1,2,1)
axe1.spines['top'].set_linewidth(3)
axe1.spines['bottom'].set_linewidth(3)
axe1.spines['right'].set_linewidth(3)
axe1.spines['left'].set_linewidth(3)
plt.plot(x[:, 1], y_predict, 'k--', linewidth=3, alpha=0.6, label='Fitting')
plt.plot(np.log10(Stand_Curve_Values[:,1]), Stand_Curve_Values[:,0], 'mo', ms=12, alpha= 0.5)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlabel('Log10([Rlative Phage Concentration]) (a.u.)',fontsize=11)
plt.ylabel('Ct',fontsize=11)
plt.legend(fontsize=9)
plt.text(4,25, r'$r^{2} = $ '+ '%.*f' % (4, r_squared), fontsize=12)

axe2 = plt.subplot(122)
axe2.spines['top'].set_linewidth(3)
axe2.spines['bottom'].set_linewidth(3)
axe2.spines['right'].set_linewidth(3)
axe2.spines['left'].set_linewidth(3)
plt.bar(index, ValuofUntre[:,1], width_of_bar, yerr = ValuofUntre[:,3], log=True,bottom=10**-1 ,\
        label='$DNase\:I^{-}$')
plt.bar(index+width_of_bar, ValuofTreatedGro[:,1], width_of_bar, yerr=ValuofTreatedGro[:,3], \
        log=True,bottom=10**-1,label='$Dnase\:I$')
plt.bar(index+width_of_bar*2, ValueofTreatedGroPlus[:,1], width_of_bar, \
        yerr=ValueofTreatedGroPlus[:,3], log=True,bottom=10**-1, label='$Dnase\:I^{+}$')
plt.xticks(index+width_of_bar, [ '$10^{9}$', '$10^{8}$', '$10^{6}$'])
plt.xlabel('Phage Concentration ($pfu\cdot mL^{-1}$)')
plt.ylabel('Relative Phage Concentration (a.u.)')
plt.legend()

plt.show()

