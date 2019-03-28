#!/usr/bin/env python
#--------------------
# Plot residual PDFs for Hawaii ground motion model assessment paper
#--------------------
# Plot PDF of ground motion residuals: 
# x-axis = period
# NOTE: we convert to log(period) after reading for better behavior in hist2d, 
#       then overlay a semilog x axis
# y-axis = residual (log(obs)-log(pred))
# color = probability 
#
# input file *_PDF.out:
#    f, T, residual, probability
# input file *_sigma.out:
#    T, mean residual, std of residuals, std of model predictions
#--------------------
# to do: 
# -command line arg to turn axis label(s) and colorbar off (for multifig plots)
# EW 2018
#--------------------

############################
# Imports
import sys
import numpy as np

from matplotlib import ticker
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

#import gmmpdfutils

############################
# Usage and verbosity
if len(sys.argv) < 2:
#    print('Usage: respdf.py inputfile')
    print('Usage: {0} MODEL'.format(sys.argv[0]))
    print('Requires 2 input files in working directory:')
    print('MODEL_PDF.out must be a list of: f, T, logresidual, pdf')
    print('MODEL_sigma.out must be a list of: T, mean logresidual, sigma_obs, sigma_pred')

    print('Plots will be named Modelpdf.[png,eps]')
    print('Will also write Model_mean.txt for use with resmean.py')
    print('use -v for debugging output')
    sys.exit()

if '-v' in sys.argv:
    verbose = True
else:
    verbose = False

modelname = sys.argv[1]
pdffile = modelname+'_PDF.out'
sigmafile = modelname+'_sigma.out'
print('Plotting residuals vs. period for {0}'.format(modelname))
if verbose:
    print('Expecting to find {0} and {1} in current directory'.format(pdffile, sigmafile))

############################
# Read input files
# Identify PGA: f, T = -1
# Separate PGA from PSA
# and for PSA, convert T to log10(T) 
# (use log(pd) bc hist2d does NOT look good if log-scaled after plotting)
f, T_orig, logres_orig, pdf_orig = np.loadtxt(pdffile, unpack=True)

# PGA is at -1 so let's separate those from psa before we try to do log(T)
i_pga, = np.where(T_orig==-1)

# Define T, logres, and PDF for pga
T_pga = T_orig[i_pga]
logres_pga = logres_orig[i_pga]
pdf_pga = pdf_orig[i_pga]

i_psa, = np.where(T_orig!=-1) 
# actually, let's pick off unwanted short periods (0.0333 and 0.04 s)
#i_psa, = np.where(T_orig>0.041) 

# Define T, logres, and PDF for psa
T = T_orig[i_psa]
logres = logres_orig[i_psa]
pdf = pdf_orig[i_psa]
T = np.log10(T)

# Find unique values of period and log-residuals for histograms
Ts = np.unique(T)
lrs = np.unique(logres)
# find min and max values here...

############################
# Read mean, std of residuals, std of predictions from MODEL_sigma.out
T_std, meanres, stdres, stdpred = np.loadtxt(sigmafile, unpack=True)
#i_psa, = np.where(T_std>0.041) 
i_psa, = np.where(T_std!=-1) 
T_std[0] = 0.000001
T_std = np.log10(T_std)
T_std[0] = -1

print(i_psa, T_std[i_psa])
T_std_psa = T_std[i_psa]
meanres_psa = meanres[i_psa]
stdres_psa = stdres[i_psa]
stdpred_psa = stdpred[i_psa]

############################
# Find bin edges for 2d histogram
# len(xbins) = len(Ts) + 1 
# len(ybins) = len(lrs) + 1 
# first bin: pad on left/bottom by difference between 1st and 2nd pt
# last bin: pad on right/top by difference from 2nd to last pt 
dx = 0.5*np.diff(Ts)
dx = np.append(dx, dx[-1])
dx0 = dx[0]
xbins = Ts + dx
xbins = np.insert(xbins, 0, Ts[0]-dx0)
if verbose:
    print('xbins, len(xbins), len(Ts):')
    print(xbins, len(xbins), len(Ts))

dy = 0.5*np.diff(lrs)
dy = np.append(dy, dy[-1])
dy0 = dy[0]
ybins = lrs + dy
ybins = np.insert(ybins, 0, lrs[0]-dy0)

############################
# Quick and dirty figure
# Scatterplot just to check on raw data values
#gmmpdfutils.circleplot(modelname, T, logres, pdf, xbins, ybins, dx, dy)

############################
# From here on: make hist2d figure
# Set up hist2d figure
# Use 2 different grid specs: 
# -one with 2 subplots for pga + psa plot
# -one with 1 subplot for colorbar 
# use the gridspecs' update() function to:
# -get the pga and period plots to touch (wspace=0)
# -keep the colorbar offset by adjusting the right and left limits
fig = plt.figure()
gs_plots = gridspec.GridSpec(1,2, width_ratios=[0.05,1])
gs_cb = gridspec.GridSpec(1,1)

ax_pga = fig.add_subplot(gs_plots[0,0])
ax_logT = fig.add_subplot(gs_plots[0,1])
ax_cb = fig.add_subplot(gs_cb[0,0])

gs_plots.update(wspace=0.0)
gs_plots.update(right=0.82)
gs_cb.update(left=0.85)

#####
# Plot 2D histogram: 
# Since we already have the PDF values (output by Dan's matlab code),
# we just use them as weights
# Need to set 'image' variable for drawing colorbar later
# xedges and yedges should be identical to input xbins and ybins
# Dan requested set max probability to 0.35

# for PSA:
counts, xedges, yedges, image = ax_logT.hist2d(T,logres, weights=pdf, bins=[xbins, ybins],cmap='bone_r', vmax=0.35, zorder=1)

# for PGA:
c2, x2, y2, im2 = ax_pga.hist2d(T_pga, logres_pga, weights=pdf_pga, bins=[[-2,0], ybins], cmap='bone_r', vmax=0.35)

#####
# Plot line at zero to highlight positive vs negative residuals
# new! use stroke (instead of plotting 2 lines)
# not sure if this is any easier than just plotting 2 lines though :P
stroke=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()]
ax_logT.plot(ax_logT.get_xlim(),[0,0], lw=1, color='red',
         path_effects=stroke, linestyle='--')
ax_pga.plot(ax_pga.get_xlim(),[0,0], lw=1, color='red',
         path_effects=stroke, linestyle='--')

#####
# Add PGA label to PGA subplot
# transAxes uses axis coordinates instead of plot coordinates
#ax_pga.text(-1, -2.1, 'PGA', ha='center', va='top')
ax_pga.text(0.5, -0.02, 'PGA', ha='center', va='top',
            transform=ax_pga.transAxes)

#####
# Calculate average and standard deviation from PDF, for each period
#avgs = np.empty(len(Ts))
#stds = np.empty(len(Ts))
#for i,Ti in enumerate(Ts):
#    iwhere, = np.where(T==Ti)
#    if verbose:
#        print('Ti:', Ti)
#        print('iwhere:', iwhere)
#    avgres = np.average(logres[iwhere], weights=pdf[iwhere])
#    var = np.average((logres[iwhere]-avgres)**2, weights=pdf[iwhere])
#    std = np.sqrt(var)
#    stds[i] = std
#    avgs[i] = avgres
#    if verbose:
#        print('Ti, avgres, std:')
#        print(Ti,avgres, std)

#####
# Plot average and standard deviation
# remember:
# marker = symbol to use (circle, triangle, etc)
# ms = markersize (in pts)
# mec = markeredgecolor (fill color, -G in GMT)
# mfc = markerfacecolor (stroke color, half of -W in GMT)
# mew = markeredgewidth (weight of edge stroke, other half of -W in GMT)
# Set up line and marker styles
stroke2=[pe.Stroke(linewidth=4, foreground='w'), pe.Normal()]
lineeffects = {'lw':2, 'color':'black', 'path_effects':stroke2, 'alpha':0.7}
lineeffects2 = {'lw':2, 'color':'r', 'path_effects':stroke, 'alpha':0.7, 'linestyle':'-'}
avgmarker = {'marker':'o', 'ms':'10', 'mec':'k', 'mew':2, 'linestyle':'None'}

# Plot mean values for PSA and PGA
ax_logT.plot(T_std_psa, meanres_psa, **avgmarker, mfc='white', alpha=0.3)
ax_logT.plot(T_std_psa, meanres_psa, **avgmarker, mfc='None', label='mean')

ax_pga.plot(T_std[0],meanres[0], **avgmarker, mfc='white', alpha=0.3)
ax_pga.plot(T_std[0],meanres[0], **avgmarker, mfc='None', label='mean')

# Plot std devs of residuals across the whole bin, instead of only at average value
# Because len(xbin) = len(avg|std) + 1,
# need to repeat each xbin 3x, but chop off 1st and last value
# and repeat avgs and stds 3x, and duplicate last value
# from dan's new output file instead of my calculation from the pdf

# need to make sure bins and sigmas are in the same order...arrrrghhh
bin_repeat = np.repeat(xbins,3)[1:-1]
avg_repeat = np.append(np.repeat(meanres_psa,3), meanres_psa[-1])
std_repeat = np.append(np.repeat(stdres_psa,3), stdres_psa[-1])
dum = np.append(np.repeat(10**T_std_psa,3), 10**T_std_psa[-1])
std_plus = np.zeros(shape=len(avg_repeat))
std_minus = np.zeros(shape=len(avg_repeat))
for i in range(len(avg_repeat)):
    print(dum[i], avg_repeat[i], std_repeat[i], 
          avg_repeat[i]+std_repeat[i], avg_repeat[i]-std_repeat[i])
    std_plus[i] = avg_repeat[i]+std_repeat[i]
    std_minus[i] = avg_repeat[i]-std_repeat[i]
#    ax_logT.plot(bin_repeat[i], std_plus[i], 'bx')
#    ax_logT.plot(bin_repeat[i], std_minus[i], 'bx')
if verbose:
    print('len of xbins, meanres_psa, stdres_psa:')
    print(len(xbins), len(meanres_psa), len(stdres_psa))
    print('len of repeated xbins, repeated means, repeated stds:')
    print(len(bin_repeat), len(avg_repeat), len(std_repeat))
    print(T_std)


#ax_logT.plot(bin_repeat[::-1], avg_repeat-std_repeat, **lineeffects, label='residual std', zorder=8)
#ax_logT.plot(bin_repeat[::-1], avg_repeat+std_repeat, **lineeffects, zorder=8)
ax_logT.plot(bin_repeat, std_plus, **lineeffects, label='residual std', zorder=8)
ax_logT.plot(bin_repeat, std_minus, **lineeffects, zorder=8)

# Now plot stds for model predictions...again repeat vals so we plot across the whole bin
bin_repeat = np.repeat(xbins,3)[1:-1]
std_repeat = np.append(np.repeat(stdpred_psa,3), stdpred_psa[-1])

ax_logT.plot(bin_repeat[::-1], std_repeat, **lineeffects2, label='predicted std', zorder=3) 
ax_logT.plot(bin_repeat[::-1], -1*std_repeat, **lineeffects2, zorder=3)

# Plot stds of residuals and model predictions for PGA
ax_pga.plot([-2,0], np.repeat(meanres[0]+stdres[0], 2), **lineeffects) 
ax_pga.plot([-2,0], np.repeat(meanres[0]-stdres[0], 2), **lineeffects)
ax_pga.plot([-2,0], np.repeat(stdpred[0], 2), **lineeffects2) 
ax_pga.plot([-2,0], np.repeat(-1*stdpred[0], 2), **lineeffects2)

# Legend
# Turn off by commenting out line below
#ax_logT.legend(numpoints=1, ncol=2, loc='lower center', fontsize='small')

#####
# Labels and such
# make labels big and make ticks larger
# hide linear log10(period) axis (we'll plot a log axis later)
ax_logT.xaxis.set_visible(False)
ax_logT.tick_params(axis='y', labelleft=False)
ylabel = r'Residual log$_{10}$(obs)-log$_{10}$(pred)'
ax_pga.set_ylabel(ylabel, fontsize='large')
ax_pga.tick_params(axis='y', labelsize='large', direction='in', which='both')
ax_pga.xaxis.set_visible(False)
ax_logT.set_title(modelname, fontsize='xx-large', x=0.5, y=0.9)

#####
# Add axes with logarithmic period x-axis: ax_Tlog
ax_Tlog = ax_logT.twiny()
ax_Tlog.set_xlabel('Period (s)', fontsize='large')
ax_Tlog.semilogx()
ax_Tlog.tick_params(labelsize='large', direction='in', which='both')
ax_Tlog.tick_params(axis='x', which='major', length=15)
ax_Tlog.tick_params(axis='x', which='minor', length=8)
ax_Tlog.tick_params(labelbottom=True, labeltop=False)
ax_Tlog.xaxis.set_ticks_position('both')
ax_Tlog.xaxis.set_label_position('bottom')

# Use decimal values (0.1, 1.0) instead of powers (10^-1, 10^0)
formatter = ticker.FormatStrFormatter('%.1f')
ax_Tlog.get_xaxis().set_major_formatter(formatter) 

# Set axis limits:
ax_Tlog.set_xlim(10**(np.min(xbins)), 10**np.max(xbins))
ax_logT.set_xlim(xbins[0], xbins[-1])

ax_logT.set_ylim(-2,2)
ax_pga.set_ylim(-2,2)

#####
# Draw colorbar 
# need to pass both lin and log axes to the 'ax' keyword arg to get proper colorbar placement
#fig.colorbar(image, ax=[ax_logT,ax_Tlog], label='Probability')

# Turn off colorbar by commenting out line below and uncommenting set_visible command
#fig.colorbar(image, cax=ax_cb, label='Probability')
ax_cb.set_visible(False)

#####
# Finish up and save
# tight layout doesn't work with the color bar (boooooo)
# if you want to trim off the whitespace,
# use this bash command after processing:
# convert -trim pdf.png meh.png
fig.savefig(modelname+'pdf.eps')
print('residual histogram saved to {0}pdf.eps'.format(modelname))
fig.savefig(modelname+'pdf.png')
fig.savefig(modelname+'pdf.pdf')
print('residual histogram saved to {0}pdf.png'.format(modelname))
#plt.show()


