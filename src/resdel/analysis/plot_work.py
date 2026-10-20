# This file is a modified part of the PMX package. I have taken the original code from the PMX package 
# and modified it to fit my needs. The original code can be found in the analysis.py file of the PMX package. 
# Check the end section of the website https://degrootlab.github.io/pmx/examples/analysis.html for details. 
# This modified code saves the forward and reverse work distributions in a file, whose path is given by the user as input,
# instead of saving it as "Wdist.png" in the current working directory. The user can specify the path and filename for the output file.

from typing import Optional
from matplotlib import pyplot as plt
import numpy as np

def gauss_func(A, mean, dev, x):
    '''Given the parameters of a Gaussian and a range of the x-values, returns
    the y-values of the Gaussian function'''
    x = np.array(x)
    y = A*np.exp(-(((x-mean)**2.)/(2.0*(dev**2.))))
    return y

def data2gauss(data):
    '''Takes a one dimensional array and fits a Gaussian.

    Parameters
    ----------
    data : list
        1D array of values

    Returns
    -------
    float
        mean of the distribution.
    float
        standard deviation of the distribution.
    float
        height of the curve's peak.
    '''
    m = np.average(data)
    dev = np.std(data)
    A = 1./(dev*np.sqrt(2*np.pi))
    return m, dev, A


def plot_work_dist(wf=[], wr=[], fname: Optional[str] = 'Wdist.png', nbins=20, dG=None, dGerr=None,
                   units='kJ/mol', dpi=300, statesProvided='AB'):
    '''Plots forward and reverse work distributions. Optionally, it adds the
    estimate of the free energy change and its uncertainty on the plot.

    Parameters
    ----------
    wf : list
        list of forward work values.
    wr : list
        list of reverse work values.
    fname : str, optional
        filename of the saved image. Default is 'Wdist.png'.
    nbins : int, optional
        number of bins to use for the histogram. Default is 20.
    dG : float, optional
        free energy estimate.
    dGerr : float, optional
        uncertainty of the free energy estimate.
    units : str, optional
        the units of dG and dGerr. Default is 'kJ/mol'.
    dpi : int
        resolution of the saved image file.
    statesProvided: str
        work values for two states or only one

    Returns
    -------
    None

    '''

    def smooth(x, window_len=11, window='hanning'):

        if x.ndim != 1:
            raise ValueError("smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
            raise ValueError("Input vector needs to be bigger than "
                             "window size.")
        if window_len < 3:
            return x
        if window not in ['flat', 'hanning', 'hamming',
                          'bartlett', 'blackman']:
            raise ValueError("Window is on of 'flat', 'hanning', 'hamming', "
                             "'bartlett', 'blackman'")
        s = np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
        # moving average
        if window == 'flat':
            w = np.ones(window_len, 'd')
        else:
            w = eval('np.' + window + '(window_len)')
        y = np.convolve(w/w.sum(), s, mode='same')
        return y[window_len-1:-window_len+1]

    plt.figure(figsize=(8, 6))
    x1 = list(range(len(wf)))
    x2 = list(range(len(wr)))
    if 'A' in statesProvided:
        x1 = list(range(len(wf)))
    if 'B' in statesProvided:
        x2 = list(range(len(wr)))
    if len(x1) > len(x2):
        x = x1
    else:
        x = x2
    if 'A' in statesProvided:
        mf, devf, Af = data2gauss(wf)
    if 'B' in statesProvided:
        mb, devb, Ab = data2gauss(wr)

    if 'AB' in statesProvided:
        maxi = max(wf+wr)
        mini = min(wf+wr)
        plt.subplot(1, 2, 1)
        plt.plot(x1, wf, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
        plt.plot(x2, wr, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
        ### smoothing ###
        try:
            sm1 = smooth(np.array(wf))
            plt.plot(x1, sm1, 'g-', linewidth=3)
        except:
            print("Plotting: no smoothing for Wf")
        try:
            sm2 = smooth(np.array(wr))
            plt.plot(x2, sm2, 'b-', linewidth=3)
        except:
            print("Plotting: no smoothing for Wr")
    elif 'A' in statesProvided:
        maxi = max(wf)
        mini = min(wf)
        plt.subplot(1, 2, 1)
        plt.plot(x1, wf, 'g-', linewidth=2, label="Forward (0->1)", alpha=.3)
        ### smoothing ###
        try:
            sm1 = smooth(np.array(wf))
            plt.plot(x1, sm1, 'g-', linewidth=3)
        except:
            print("Plotting: no smoothing for Wf")
    elif 'B' in statesProvided:
        maxi = max(wr)
        mini = min(wr)
        plt.subplot(1, 2, 1)
        plt.plot(x2, wr, 'b-', linewidth=2, label="Backward (1->0)", alpha=.3)
        ### smoothing ###
        try:
            sm2 = smooth(np.array(wr))
            plt.plot(x2, sm2, 'b-', linewidth=3)
        except:
            print("Plotting: no smoothing for Wr")

    plt.legend(shadow=True, fancybox=True, loc='upper center',
               prop={'size': 12})
    plt.ylabel(r'W [kJ/mol]', fontsize=20)
    plt.xlabel(r'# Snapshot', fontsize=20)
    plt.grid(lw=2)
    plt.xlim(0, x[-1]+1)
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplot(1, 2, 2)
    if 'A' in statesProvided:
        plt.hist(wf, bins=nbins, orientation='horizontal', facecolor='green',
             alpha=.75, density=True)
    if 'B' in statesProvided:
        plt.hist(wr, bins=nbins, orientation='horizontal', facecolor='blue',
             alpha=.75, density=True)

    x = np.arange(mini, maxi, .5)

    if 'AB' in statesProvided:
        y1 = gauss_func(Af, mf, devf, x)
        y2 = gauss_func(Ab, mb, devb, x)
        plt.plot(y1, x, 'g--', linewidth=2)
        plt.plot(y2, x, 'b--', linewidth=2)
        size = max([max(y1), max(y2)])
    elif 'A' in statesProvided:
        y1 = gauss_func(Af, mf, devf, x)
        plt.plot(y1, x, 'g--', linewidth=2)
        size = max(y1)
    elif 'B' in statesProvided:
        y2 = gauss_func(Ab, mb, devb, x)
        plt.plot(y2, x, 'b--', linewidth=2)
        size = max(y2)

    res_x = [dG, dG]
    res_y = [0, size*1.2]
    if dG is not None and dGerr is not None:
        plt.plot(res_y, res_x, 'k--', linewidth=2,
                 label=r'$\Delta$G = %.2f $\pm$ %.2f %s' % (dG, dGerr, units))
        plt.legend(shadow=True, fancybox=True, loc='upper center',
                   prop={'size': 12})
    elif dG is not None and dGerr is None:
        plt.plot(res_y, res_x, 'k--', linewidth=2,
                 label=r'$\Delta$G = %.2f %s' % (dG, units))
        plt.legend(shadow=True, fancybox=True, loc='upper center',
                   prop={'size': 12})
    else:
        plt.plot(res_y, res_x, 'k--', linewidth=2)

    plt.xticks([])
    plt.yticks([])
    xl = plt.gca()
    for val in xl.spines.values():
        val.set_lw(2)
    plt.subplots_adjust(wspace=0.0, hspace=0.1)
    plt.savefig(fname, dpi=dpi)

