import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
import statsmodels.api as sm
from matplotlib import style
import matplotlib as mpl





def framean(x=np.array([1, 2, 3]), length=2):
    lengthofrm = x.size - length + 1
    Index = np.arange(lengthofrm)
    frameanlist = [np.sum(x[i:(i+length-1)])/length for i in Index]
    return frameanlist


def dataprocess(time, od600):
    """
    :param time: a list, time
    :param od600: dataframe, OD600
    :return: OD Mean, OD std, DO fitting,
    Fitting Parameters[r_squared, cont, slope, Doubling_Time,],
    interpolate time, interpolate OD
    """
    OD_Mean = np.mean(od600, axis=1)
    OD_Std = od600.std(axis=1)
    IndexofRC = np.where((OD_Mean > 0.02) & (OD_Mean < 0.4))
    x = time[IndexofRC]
    y = OD_Mean[IndexofRC[0]]
    #yer = OD_Std[IndexofRC[0]]
    xs = np.linspace(x[0], x[-1], 100)
    curve = interpolate.pchip(x, y.values)
    ys = curve(xs)
    logys = np.log10(ys)
    dlogys = np.gradient(logys, xs)
    frame_size = 10  #Frame Size!
    frameandlogys = framean(dlogys, frame_size)
    frameanxs = xs[0:(xs.size - frame_size + 1)]
    bottom_index = np.where(frameandlogys > frameandlogys[0] * 0.96)
    bottom_time = xs[(bottom_index[0][-1] + 100 / x.size * 0.3 * frame_size).astype(int)]
    exp_time = x[np.where(x < bottom_time)]
    exp_growth = y.iloc[np.where(x < bottom_time)[0]]

    if np.max(frameandlogys > frameandlogys[0] * 1.03):
        print("Diauxic growth!")

    # def residuals(p):
    #     k, b = p
    #     return np.log(exp_growth) - (k * exp_time + b)
    #
    # r = optimize.leastsq(residuals, [0, -6])
    # k, b = r[0]
    # dou_time = np.log(2) / k
    # print("Doubling time = %.*f min." % (2, dou_time))
    # fity = np.e ** (k * exp_time + b)

    # Start to calculate least square
    X = exp_time
    X = sm.add_constant(X)
    Y = np.log10(exp_growth)
    Model_Fit = sm.OLS(Y, X)
    Fit_Results = Model_Fit.fit()
    r_squared = Fit_Results.rsquared
    cont = Fit_Results.params[0]
    slope = Fit_Results.params[1]
    y_predict = Fit_Results.fittedvalues
    OD_predict = 10.0**y_predict
    Doubling_Time = np.log(2)/slope
    print(Fit_Results.summary())
    print("Doubling time = %.*f min." % (2, Doubling_Time))
    Fitting_Para = [r_squared, cont, slope, Doubling_Time]
    return OD_Mean, OD_Std, OD_predict, exp_time, Fitting_Para, xs, ys, frameanxs, frameandlogys

def draw_text(ax, str):
    from matplotlib.offsetbox import AnchoredText
    at = AnchoredText(str,
                      loc=2
                      )
    #at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)

if __name__ == '__main__':
    # from scipy import interpolate
    rawdata = pd.read_csv('20180810.csv')  # input data
    # set a time serial
    locNum = rawdata.iloc[:, 0].size
    time_interval = 10
    time = np.arange(9, locNum*time_interval, time_interval)    #set a time interval.
    cl2 = rawdata[["G2", "G3", "G4", "G5", "G6"]]  # experimental group
    cr = rawdata[["G7", "G8", "G9", "G10", "G11"]]      # blank
    SAMPLE_NAME = 'NONE'
    cr_cl2 = np.subtract(cl2, cr)

    cr_cl2_mean, cr_cl2_std, fity, exp_time, Parameters, xs, ys, frameanxs, frameandlogys = dataprocess(time, cr_cl2)
    style.use('bmh')
    mpl.rcParams['font.size'] = 16
    mpl.rcParams['lines.linewidth'] = 3.5
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['figure.subplot.wspace'] = 0.4

    plt.figure(figsize=(16,16))
    axe1 = plt.subplot(221)
    plt.errorbar(time, cr_cl2_mean, cr_cl2_std, label=SAMPLE_NAME, ms=8,\
                 elinewidth=1)
    plt.legend( )
    plt.yscale('log')
    plt.xlabel('Time(min)')
    plt.ylabel("ABS"+"$_{600}$")
    plt.tick_params(which="major", length=7, width=2)
    plt.tick_params(which="minor", length=3.5, width=2)
    plt.xticks()
    plt.yticks()
    axe2 = plt.subplot(222)
    plt.errorbar(time, cr_cl2_mean, cr_cl2_std, fmt='ko', alpha=0.2, ms=8)
    plt.plot(exp_time, fity, 'r--', label="Fitting")
    plt.xlim(xmax=exp_time[-1]+40, xmin=exp_time[0]-40)
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Time(min)',)
    plt.ylabel("ABS"+"$_{600}$",)
    draw_text(axe2,
              "Doubling Time: %.3f min \n $r^{2}$ : %.3f " % (Parameters[-1], Parameters[0]))

    plt.tick_params(which="major", length=7, width=2)
    plt.tick_params(which="minor", length=3.5, width=2)
    plt.xticks()
    plt.yticks()
    axe3 = plt.subplot(223)
    plt.plot(time, cr_cl2_mean, "ro")
    plt.plot(xs, ys, '--')
    plt.yscale('log')
    plt.xlabel('Time(min)')
    plt.ylabel("ABS"+"$_{600}$")
    plt.tick_params(which="major", length=7, width=2)
    plt.tick_params(which="minor", length=3.5, width=2)
    plt.xlim(xmax=xs[-1], xmin=xs[0])
    plt.ylim(ymin=ys[0], ymax=ys[-1])

    axe4 = plt.subplot(224)
    plt.plot(frameanxs, frameandlogys, 'g-', label="Growth Rate")
    plt.hlines(frameandlogys[0]*1.03, frameanxs[0], frameanxs[-1])
    plt.hlines(frameandlogys[0]*0.97, frameanxs[0], frameanxs[-1])
    plt.xlim(xmax=exp_time[-1]+40)
    plt.ylim(ymax=np.max(frameandlogys)*1.4, ymin=frameandlogys[-1]*0.8)
    plt.legend()
    plt.xlabel('Time(min)')
    plt.ylabel("ABS"+"$_{600}\cdot min^{-1}$")
    plt.tick_params(which="major", length=7, width=2)
    plt.tick_params(which="minor", length=3.5, width=2)
    plt.legend()
    plt.savefig('00')
