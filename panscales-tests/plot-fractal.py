import sys
import os
sys.path.insert(0, os.path.dirname(__file__) + "/../../panscales-main/2020-eeshower/submodules/AnalysisTools/python")
from hfile import *

import matplotlib
from matplotlib import pylab, mlab, pyplot
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
from matplotlib import gridspec
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib import rc
plt = pyplot
from pylab import *
from numpy import *

#from IPython.core.pylabtools import figsize, getfigs
from matplotlib import rc


font = {'color':  'black',
    'weight': 'normal',
    'size': 18
    }
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Palatino']
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.numpoints'] = 1
plt.rcParams['legend.handlelength'] = 1

#########################################
with PdfPages('fractal-lund.pdf') as pdf:
    data    = get_array_plus_comments('a', 'diff_hist:eps_area')
    fig,(ax,axr) = plt.subplots(figsize=(5,3.8),nrows=2,height_ratios=[2,1],sharex=True)   

    ax.grid(True, lw=0.5, ls=':', zorder=0)
    axr.grid(True, lw=0.5, ls=':', zorder=0)
    ax.set_title(r'$e^+e^-\to q\bar q$, $\sqrt{s}=2$ TeV',fontsize=16)
    ax.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )
    axr.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )
    ax.set_ylabel(r'$\ln\mathcal{A}$',fontsize=16)
    axr.set_xlabel(r'$\ln\varepsilon=\ln \frac{k_t}{Q}$',fontsize=16)
    # we take the log of both variables
    ax.plot(np.log(data[:,0]), np.log(data[:,3]),color='crimson')
    ax.set_xlim(-6,0)
    ax.minorticks_on()
    ax.legend(loc='best', edgecolor='white',borderpad=0.1, fontsize=14)
    ratio = -np.log(data[:,3])/np.log(data[:,0])
    axr.set_ylabel(r'$d$',fontsize=16)
    axr.set_ylim(1.99,2.05)
    axr.set_xlim(-6,0)
    # the dimension is just the log between the 2 things
    axr.plot(np.log(data[:,0]), ratio,color='navy')
    pdf.savefig(bbox_inches='tight')
    plt.close()
