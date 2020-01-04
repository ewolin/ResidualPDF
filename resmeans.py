#!/usr/bin/env python
#--------------------
# Make lots of figures for Hawaii ground motion model assessment paper
# EW 2018
#--------------------
# Expects to find in current working dir:
# *_sigma.out (one for each model)
#     period (s), mean residual, std of residuals, std of model predictions
#
# *_[volcanic,tectonic]_shakemap.txt (one pair of v,t for each model)
#     period (s), mean residual, std of residuals
#
# *_LLH.out (one for each model selected...in this case 7 out of starting 15)
#     period (s), LLH score
#
# Outputs 5 different plots, all with semilog x axes for [y] vs period
# PGA is shown in a separate little subplot in the left
# 1. Heatmap of mean residuals: 2 cols (PGA and PSAs) x n rows (one for each model)
# 2. Lines: mean residuals from A10 database, w/transparent shaded region showing +/- standard deviation
# 3. smt: same as #2, but for ShakeMap tectonic data
# 4. smv: same as #2, but for ShakeMap volcanic data
# 5. llh: log-likelihood weights for selected ground motion models
#--------------------
# to do:
# i don't love how much repeated code there is for reading and plotting input files, but hey it works for now.
#--------------------

import os
import sys
import numpy as np
from glob import glob

import itertools
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker

from gmmpdfutils import xbinedges, setupLogTPlot

# to use classic style:
#import matplotlib.style
#import matplotlib as mpl
#mpl.style.use('classic')

############################
# Input files
meanfiles = glob('*_sigma.out')
nfiles = len(meanfiles)

############################
# Set up figures using setupLogTPlot from gmmpdfutils.py
# good colormaps: bwr, seismic, RdBu_r
# other: PuOr_r

# Set up heatmap figure
# can't use my setupLogTPlot function, since we've got lots of rows...
fig_heat = plt.figure(figsize=(10,6))
gs_plots_heat = gridspec.GridSpec(nfiles,2, width_ratios=[0.05,1])
gs_cb_heat = gridspec.GridSpec(1,1)
gs_plots_heat.update(wspace=0.01)
 
# Set up lines figure
meanres_label = r'Mean residual log$_{10}$(obs)-log$_{10}$(pred)'
ylim_lines = [-2,2]
fig_lines, axdict_lines = setupLogTPlot(ylim=ylim_lines, 
                                        ylabel=meanres_label)
ax_lines_psa = axdict_lines['psa']
ax_lines_pga = axdict_lines['pga']

# Set up LLH figure
llh_label = 'Log-likelihood score'
ylim_llh = [0.5,2]
fig_llh, axdict_llh = setupLogTPlot(ylim=ylim_llh, ylabel=llh_label)
ax_llh_psa = axdict_llh['psa']
ax_llh_pga = axdict_llh['pga']

# Set up ShakeMap figures
# smt = shakemap tectonic
# smv = shakemap volcanic
ylim_sm = [-3,3]
fig_smt, axdict_smt = setupLogTPlot(ylim=ylim_sm, ylabel=meanres_label)
ax_smt_psa = axdict_smt['psa']
ax_smt_pga = axdict_smt['pga']

fig_smv, axdict_smv = setupLogTPlot(ylim=ylim_sm, ylabel=meanres_label)
ax_smv_psa = axdict_smv['psa']
ax_smv_pga = axdict_smv['pga']

##################################
# Color maps and symbol cycle
#colormap = plt.cm.nipy_spectral # for rainbow colors
#colormap = plt.cm.PuOr_r
colormap = plt.cm.coolwarm_r
cmap = [colormap(i) for i in np.linspace(0,1,nfiles)]
axdict_lines['psa'].set_prop_cycle('color', cmap)
axdict_lines['pga'].set_prop_cycle('color', cmap)
ax_llh_pga.set_prop_cycle('color', cmap)
ax_llh_psa.set_prop_cycle('color', cmap)
ax_smt_pga.set_prop_cycle('color', cmap)
ax_smt_psa.set_prop_cycle('color', cmap)
ax_smv_pga.set_prop_cycle('color', cmap)
ax_smv_psa.set_prop_cycle('color', cmap)
markers = itertools.cycle(['o', 'v', 'd', 's', '^'])

############################
# Calculate mean residuals across all periods and get sorted list of models 
# *_sigma.out files are list of:
# period, mean residual, std of residuals, std of model predictions
modelmeans = {}
for i,mf in enumerate(meanfiles):
    modelname = mf.split('_')[0]
    print(modelname)
    sigmafile = modelname+'_sigma.out'
    T_std, meanres, stdres, stdpred = np.loadtxt(sigmafile, unpack=True)
    i_psa, = np.where(T_std > 1e-5)
#    i_psa, = np.where(T_std>0.041)
#    i_psa, = np.where((T_std > 1e-5) & (T_std != 0.075) & (T_std != 4.0) )
    print(i_psa)
    meanres_psa = meanres[i_psa]
    print(meanres[i_psa])
    meanavg = np.mean(meanres_psa)
# try rms instead of mean residual:
#    meanavg = np.sqrt(np.mean(meanres_psa*meanres_psa))
    #print(modelname, meanavg)
    modelmeans[modelname] = meanavg
print('------')

# Get list of model names and average residual, then sort 
k = list(modelmeans.keys())
v = list(modelmeans.values())

i_sort = np.argsort(v)
v.sort()
newmods = []
for i in i_sort:
    newmods.append(k[i])

for i in range(len(i_sort)):
    print(newmods[i], v[i])


############################
# Plot models in order of decreasing mean residual 
# (from positive to negative)
n_llh = len(glob('*LLH.out'))
i_llh = 0

for i,modelname in enumerate(newmods[::-1]):
    marker = next(markers)
    # Read input file
    # should probably move reading input PSA vs PGA into a function
    sigmafile = modelname+'_sigma.out'
    T, meanres, stdres, stdpred = np.loadtxt(sigmafile, unpack=True)
    i_psa, = np.where(T > 1e-5)
    #i_psa, = np.where(T>0.041)
    #i_psa, = np.where((T > 1e-5) & (T != 0.075) & (T != 4.0 ))
    i_pga, = np.where(T==-1)
    T_psa = np.log10(T[i_psa])
    T_pga = T[i_pga]
    meanres_psa = meanres[i_psa]
    meanres_pga = meanres[i_pga]
    stdres_psa = stdres[i_psa]
    stdres_pga = stdres[i_pga]
    xbins = xbinedges(T_psa)

    # ShakeMap tectonic...kludgey AF
    smtfile = modelname+'_tectonic_shakemap.txt'
    if os.path.exists(smtfile):
        T_smt, meanres_smt, stdres_smt = np.loadtxt(smtfile, unpack=True)
        i_psa_smt, = np.where(T_smt!=-1)
        i_pga_smt, = np.where(T_smt==-1)
        T_psa_smt = np.log10(T_smt[i_psa_smt])
        T_pga_smt = T_smt[i_pga]
        meanres_psa_smt = meanres_smt[i_psa_smt]
        meanres_pga_smt = meanres_smt[i_pga_smt]
        stdres_psa_smt = stdres_smt[i_psa_smt]
        stdres_pga_smt = stdres_smt[i_pga_smt]
    else:
        print(smtfile, 'not found, will not plot ShakeMap tect residuals')
 
    # ShakeMap volcanic...kludgey AF
    smvfile = modelname+'_volcanic_shakemap.txt'
    if os.path.exists(smvfile):
        T_smv, meanres_smv, stdres_smv = np.loadtxt(smvfile, unpack=True)
        i_psa_smv, = np.where(T_smv!=-1)
        i_pga_smv, = np.where(T_smv==-1)
        T_psa_smv = np.log10(T_smv[i_psa_smv])
        T_pga_smv = T_smv[i_pga]
        meanres_psa_smv = meanres_smv[i_psa_smv]
        meanres_pga_smv = meanres_smv[i_pga_smv]
        stdres_psa_smv = stdres_smv[i_psa_smv]
        stdres_pga_smv = stdres_smv[i_pga_smv]
    else:
        print(smvfile, 'not found, will not plot ShakeMap volc residuals')

    # Plot LLH values: these only exist for 7 'selected' GMMs
    try:
        llhfile = modelname+'_LLH.out'
        T_llh, llh = np.loadtxt(llhfile, unpack=True)
#        i_psa_llh, = np.where(T_llh>0.041)
        i_pga_llh, = np.where(T_llh==-1)
        T_psa_llh = np.log10(T[i_psa_llh])
        T_pga_llh = T[i_pga_llh]
        llh_psa = llh[i_psa_llh]
        llh_pga = llh[i_pga_llh]
        plotllh = True
        i_llh += 1
        # LLH plot
        ax_llh_psa.plot(T_psa_llh, llh_psa, label=modelname, 
                        marker=marker, color=cmap[i])
        ax_llh_pga.plot(i_llh/n_llh-1.5, llh_pga, label=modelname, 
                        marker=marker, color=cmap[i])
    except:
        plotllh = False

    # Lines w/std
    ax_lines_psa.plot(T_psa, meanres_psa, label=modelname, marker=marker)
    ax_lines_psa.fill_between(T_psa, meanres_psa+stdres_psa, 
                              meanres_psa-stdres_psa, color=cmap[i], 
                              alpha=0.1)
    ax_lines_pga.plot(T_pga, meanres_pga, label=modelname, marker=marker)
    ax_lines_pga.fill_between([-2,0], np.repeat(meanres_pga+stdres_pga, 2),
                              np.repeat(meanres_pga-stdres_pga, 2), 
                              color=cmap[i], alpha=0.1)

    # Lines w/std for ShakeMap tectonic
    if os.path.exists(smvfile):
        psa_std_plus = meanres_psa_smt + stdres_psa_smt
        psa_std_minus = meanres_psa_smt - stdres_psa_smt
        pga_std_plus = np.repeat(meanres_pga_smt + stdres_pga_smt, 2)
        pga_std_minus = np.repeat(meanres_pga_smt - stdres_pga_smt, 2)
        ax_smt_psa.plot(T_psa_smt, meanres_psa_smt, label=modelname, 
                        marker=marker)
        ax_smt_psa.fill_between(T_psa_smt, psa_std_plus, psa_std_minus, 
                                color=cmap[i], alpha=0.1)
        ax_smt_pga.plot(T_pga_smt, meanres_pga_smt, label=modelname, 
                        marker=marker)
        ax_smt_pga.fill_between([-2,0], pga_std_plus, pga_std_minus,
                                color=cmap[i], alpha=0.1)

    # Lines w/std for ShakeMap volcanic
    if os.path.exists(smvfile):
        psa_std_plus = meanres_psa_smv + stdres_psa_smv
        psa_std_minus = meanres_psa_smv - stdres_psa_smv
        pga_std_plus = np.repeat(meanres_pga_smv + stdres_pga_smv, 2)
        pga_std_minus = np.repeat(meanres_pga_smv - stdres_pga_smv, 2)
        ax_smv_psa.plot(T_psa_smv, meanres_psa_smv, label=modelname, 
                        marker=marker)
        ax_smv_psa.fill_between(T_psa_smv, psa_std_plus, psa_std_minus, 
                                color=cmap[i], alpha=0.1)
        ax_smv_pga.plot(T_pga_smv, meanres_pga_smv, label=modelname, 
                        marker=marker)
        ax_smv_pga.fill_between([-2,0], pga_std_plus, pga_std_minus, 
                                color=cmap[i], alpha=0.1)

    # Heat map
    axes_heat = fig_heat.add_subplot(gs_plots_heat[i,1])
    axes_heat_pga = fig_heat.add_subplot(gs_plots_heat[i,0])
    counts, xedges, yedges, image = axes_heat.hist2d(T_psa, np.zeros(len(T_psa)), weights=meanres_psa, bins=[xbins, [-1,1]], cmap='bwr', vmin=-2, vmax=2)
    counts, xedges, yedges, image2 = axes_heat_pga.hist2d(T_pga, np.zeros(len(T_pga)), weights=meanres_pga, bins=[[-2,0], [-1,1]], cmap='bwr', vmin=-2, vmax=2)

    # Set heat map limits and overlay log x axis on each row
    axes_heat.set_xlim(xbins[0], xbins[-1])
    axes_heat_pga.xaxis.set_visible(False)
    axes_heat_pga.set_yticklabels([])
    axes_heat_pga.tick_params(axis='y', labelleft='off', size=0)
    axes_heat.yaxis.set_major_locator(plt.NullLocator())
    axes_heat.xaxis.set_major_locator(plt.NullLocator())
    axes_heat_pga.set_ylabel(modelname, rotation=0, ha='right', va='center')
    if i != nfiles-1:
        ax_log = axes_heat.twiny()
        ax_log.set_xlim(10**xbins[0], 10**xbins[-1])
        ax_log.semilogx()
        axes_heat.xaxis.set_visible(False) 
        ax_log.xaxis.set_ticks_position('bottom')
        plt.setp(ax_log.get_xticklabels(), visible=False)
    else:
        axtwin = axes_heat
        axlastpga = axes_heat_pga

##################################
# Add zero lines for reference
zeroline = {'ls':'--', 'color':'black', 'zorder':0, 'lw':3}

# Lines
ax_lines_psa.plot([xbins[0], xbins[-1]], [0,0], **zeroline) 
ax_lines_pga.plot([-1,1], [0,0], **zeroline) 

# ShakeMap tectonic
ax_smt_psa.plot([xbins[0], xbins[-1]], [0,0], **zeroline) 
ax_smt_pga.plot([-1,1], [0,0], **zeroline) 

# ShakeMap volcanic
ax_smv_psa.plot([xbins[0], xbins[-1]], [0,0], **zeroline) 
ax_smv_pga.plot([-1,1], [0,0], **zeroline) 

##################################
# Set limits and add legends 
# Lines
axdict_lines['psa'].set_xlim(xbins[0], xbins[-1])
axdict_lines['psalog'].set_xlim(10**xbins[0], 10**xbins[-1])
ax_lines_psa.legend(ncol=4, numpoints=1, fontsize='9', framealpha=0.5)

# LLH
ax_llh_psa.set_xlim(xbins[0], xbins[-1])
axdict_llh['psalog'].set_xlim(10**xbins[0], 10**xbins[-1])
ax_llh_psa.legend(ncol=4, numpoints=1, fontsize='9', framealpha=0.5)

# Shakemap tectonic
ax_smt_psa.set_xlim(xbins[0], xbins[-1])
axdict_smt['psalog'].set_xlim(10**xbins[0], 10**xbins[-1])
ax_smt_psa.legend(title='ShakeMap tectonic events', ncol=5, numpoints=1, 
                  fontsize='9', framealpha=0.5)

# Shakemap volcanic
ax_smv_psa.set_xlim(xbins[0], xbins[-1])
axdict_smv['psalog'].set_xlim(10**xbins[0], 10**xbins[-1])
ax_smv_psa.legend(title='ShakeMap volcanic events', ncol=5, numpoints=1, 
                  fontsize='9', framealpha=0.5)

##################################
# Finishing touches on heatmap axes 
# (similar to what's done in setupLogTPlot)
# Heat: add 'PGA' label
axlastpga.text(0.5, -0.2, 'PGA', ha='center', va='top', 
               transform=axlastpga.transAxes)

# Heat: labels, ticks, limits on bottom row
ax3 = axtwin.twiny()
ax3.set_xlim(10**xbins[0], 10**xbins[-1])
ax3.set_xlabel('Period (s)', fontsize='large')
ax3.semilogx()
ax3.xaxis.set_ticks_position('bottom')
ax3.xaxis.set_label_position('bottom')
formatter = ticker.FormatStrFormatter('%.1f')
ax3.get_xaxis().set_major_formatter(formatter) 

# Heat: colorbar
cblabel = r'Mean residual log$_{10}$(obs)-log$_{10}$(pred)'
cb = fig_heat.colorbar(image, shrink=0.7, ax=fig_heat.axes, 
                       orientation='vertical', extend='both', 
                       label=cblabel)
loc = ticker.MultipleLocator(base=0.5)
cb.locator = loc
cb.update_ticks()

##################################
# Save figures
print('saving figures...')
fig_heat.savefig('resmeans_heat.png')
#fig_heat.savefig('resmeans_heat.eps')
fig_lines.savefig('resmeans_lines.png')
#fig_lines.savefig('resmeans_lines.eps')
fig_smt.savefig('resmeans_smt.png')
#fig_smt.savefig('resmeans_smt.eps')
fig_smv.savefig('resmeans_smv.png')
#fig_smv.savefig('resmeans_smv.eps')
fig_llh.savefig('resmeans_llh.png')
#fig_llh.savefig('resmeans_llh.eps')
print('done')
#plt.show()
