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

def get_data(filename, histname):
    data = get_2darray(filename, "2d_hist_compact:"+histname)
    return data

eta_slices = [r'$0 < \eta < 1$',r'$1 < \eta < 2$', r'$2 < \eta < 3$', r'$3 < \eta < 4$']
colours = ['navy', 'crimson', 'gray', 'green']
files=['panglobal','panglobal-cffe','panlocal','pythia8']
showers=[r'PG$_\beta=0$',r'PG$^{\rm{CFFE}}_{\beta=0}$', r'PL$_{\beta=0.5}$',r'PY8']
#########################################
with PdfPages('exclusive-lund.pdf') as pdf:
    for f,s in zip(files,showers):
        ratio_xs = get_xsection('2024-10-29-results/'+f,'eta_lnkt_w_coll').xsc/get_xsection('2024-10-29-results/'+f,'eta_lnkt_wo_coll').xsc
        dataw    = get_data('2024-10-29-results/'+f, 'eta_lnkt_w_coll')
        datawo   = get_data('2024-10-29-results/'+f, 'eta_lnkt_wo_coll')
        fig = plt.figure(figsize=(5,3.8))
        ax = plt.gca()

        ax.grid(True, lw=0.5, ls=':', zorder=0)
        plt.title(s+r', $e^+e^-\to q\bar q$, $\sqrt{s}=2$ TeV',fontsize=16)
        ax.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True )
        ax.set_ylabel(r'veto / no veto ($5<\eta<6$)',fontsize=16)
        ax.set_xlabel(r'$\ln k_t$ [GeV]',fontsize=16)

        for i, l, c in zip([0,1,2,3], eta_slices,colours):
            ax.plot(dataw[:,:,1][i], (datawo[:,:,2][i]*ratio_xs/dataw[:,:,2][i]), label=l, color=c)

        ax.set_ylim(0.8,1.6)
        ax.set_xlim(1.5,6)
        ax.minorticks_on()
        ax.hlines(1,xmin=0.5,xmax=6,color='black',ls='--')
        ax.legend(loc='best', edgecolor='white',borderpad=0.1, fontsize=14)

        pdf.savefig(bbox_inches='tight')
        plt.close()
