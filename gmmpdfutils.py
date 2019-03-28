#!/usr/bin/env python

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
# Utilities for reading, manipulating, plotting ground motion residual data

def setupLogTPlot(Tlim=[0.04, 5.5], ylim=[-2,2], ylabel='ylabel'):
# use y instead of residuals...bc works for LLH plot too!
# NOTE: I assume PGA is denoted by T=-1, 
#       so ax_pga limits are hard-coded from -2 to 0
    '''Set up a plot with gridspec for two subplots: narrow strip for PGA, 
       and square-ish for PSA
       Overlay log x axis 
       Return figure and a dict of pga, psa, and logpsa axes'''
    fig = plt.figure() 

    gs_plots = gridspec.GridSpec(1,2, width_ratios=[0.05,1])
    gs_plots.update(wspace=0.0)
    gs_cb = gridspec.GridSpec(1,1)

    ax_pga = fig.add_subplot(gs_plots[0,0])
    ax_psa = fig.add_subplot(gs_plots[0,1])
    
    ax_psa.yaxis.set_ticks_position('both')
    ax_psa.tick_params(axis='y', labelleft='off')

    ax_pga.xaxis.set_visible(False)
    ax_pga.text(0.5, -0.02, 'PGA', ha='center', va='top', 
                transform=ax_pga.transAxes)

# Overlay log x-axis
    ax_psa_log = ax_psa.twiny()
    ax_psa_log.set_xlabel('Period (s)', fontsize='large')
    ax_psa_log.semilogx()
    logTlim=np.log10(np.array(Tlim))
    print(Tlim)
    ax_psa_log.set_xlim(Tlim[0], Tlim[1])
    ax_psa_log.xaxis.set_ticks_position('bottom')
    ax_psa_log.xaxis.set_label_position('bottom')
    formatter = ticker.FormatStrFormatter('%.1f')
    ax_psa_log.get_xaxis().set_major_formatter(formatter) 

# Limits and labels
# important to set these AFTER twinning!!
    ax_psa.set_xlim(logTlim[0], logTlim[1])
    ax_psa.set_ylim(ylim[0], ylim[1])
    ax_psa.xaxis.set_visible(False)

    ax_pga.set_ylabel(ylabel)
    ax_pga.set_xlim(-2,0)
    ax_pga.set_ylim(ylim[0], ylim[1])

    axdict = {'pga':ax_pga, 'psa':ax_psa, 'psalog':ax_psa_log}
    return fig, axdict 

###############################################################
def circleplot(modelname, T, logres, pdf, xbins, ybins, dx, dy):
    '''Quick and dirty figure just to check on raw data values:
       Plot values as circles so we can check hist2d results later'''
# Set up subplot
    fig, ax = plt.subplots(1)

# Plot data as colored circles, 
# note that we do not attempt to match vmax used in hist2d
    ax.scatter(T, logres, c=pdf, cmap='bone_r')

    
#####
# Plot bins to make sure we got them right
# xbins at y=1
# ybins at x=0
    ax.plot(xbins, np.ones(shape=len(xbins)), marker='|', mfc='k', mec='k', ms=15, linestyle='None')
    ax.plot(np.zeros(shape=len(ybins)), ybins, marker='_', mfc='k', mec='k', ms=15, linestyle='None')

# Plot dashed line at y=zero
    ax.plot([xbins.min(), xbins.max()],[0,0],linestyle='--', color='grey')

# Set axis limits and labels
    ax.set_ylim(logres.min(),logres.max())
    ax.set_ylabel('Residual')

    ax.set_xlim(xbins.min()-dy[0],xbins.max()+dy[-1])
    ax.set_xlabel('log10(period) (s)')
    
# Save plot
    fig.savefig(modelname+'circles.png')
    print('quick scatter plot saved to {0}circles.png'.format(modelname))
###############################################################

def xbinedges(X_in):
    '''Find bin edges for an array X:
       len(xbins) = len(X) + 1 
       first bin: pad on left/bottom by difference between 1st and 2nd pt
       last bin: pad on right/top by difference from 2nd to last pt'''

    X = X_in.copy()
    X.sort()
    dx = 0.5*np.diff(X)
    dx0 = dx[0]
    dx = np.append(dx, dx[-1])
    xbins = X + dx
    xbins = np.insert(xbins, 0, X[0]-dx0)
    return xbins

###############################################################
class Psa:
    def __init__(self, T, logres, pdf):
# naming convention for unique, sorted lists of vals, versus gridded vals needed for hist2d?
        self.T = T
        self.logT = np.log10(T) 
        self.logres = logres
        self.pdf = pdf

        self.T_unique = np.unique(self.T)
        self.logT_unique = np.unique(self.logT)
        self.logres_unique = np.unique(self.logres)

        self.Tbins = xbinedges(self.T_unique)
        self.logresbins = xbinedges(self.logres_unique)

    def xbinedges(self, X):
        '''Find bin edges for an array X:
           len(xbins) = len(X) + 1 
           first bin: pad on left/bottom by difference between 1st and 2nd pt
           last bin: pad on right/top by difference from 2nd to last pt'''

        dx = 0.5*np.diff(X)
        dx0 = dx[0]
        dx = np.append(dx, dx[-1])
        xbins = X + dx
        xbins = np.insert(xbins, 0, X[0]-dx0)
        return xbins


#class Pga: def __init__(self, ): # self.T = -1 # no period defined (no x axis on res plot) # basically a 1d hist on y axis so just need values of log residuals and pdf self.logres = self.pdf = class Pgv: def __init__(self, ): # like pga, no defined period, just need log residuals and pdf    
#       self.logres = 
#       self.pdf = 
#
#
#
#
#def readPDF(modelname):
#    
##    if 
#
## need to decide:
## is it better to split PSA, PGA, PGV early on?
##    return T_psa, logT_psa, logres_psa, pdf_psa, res_pga, pdf_pga, res_pgv, pdf_pgv
#
## or maybe we can design plotting stuff so that it auto identifies T>0?
# #   return T, logT, logres, pdf  # will need to take log of only T>0 anyway
#
## fewer variable names if split later, but more potential for confusion esp w/log of "negative" pds for pga+pgv if kept all together
#
## maybe define psa, pga, pgv CLASSES and return these??
#
#
#def readsigma(modelname):
#    return T_sigma, logT_sigma, mean_logres, std_logres, std_pred 
#
#
#
