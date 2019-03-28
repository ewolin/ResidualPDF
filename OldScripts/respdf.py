#!/usr/bin/env python
#####
# Plot PDF of ground motion residuals: 
# x-axis = period (actually hist2d uses log(period) but we plot a logarithmic axis on top...MAKE SURE YOUR X AXES HAVE THE SAME RANGE)
# y-axis = residual (ln(obs)-ln(pred))
# color = probability (calculated separately)
# input file:
# f, T, residual, probability

# to do: 
# -finish cleaning up organization/imports/etc
# x write name of output files to screen
# x implement verbose switch
# -switch to turn axis label(s) and colorbar off (for multifig plots)
# -add in PGA and PGV: autodetect negative periods, plot separately?
# E. Wolin July 2018
#####

#####
# Imports
import sys
import numpy as np

from matplotlib import ticker
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

#####
# Usage and verbosity
if len(sys.argv) < 2:
    print('Usage: respdf.py inputfile')
    print('Input file must be a list of: f, T, logresidual, pdf')
    print('Naming convention for input files: Model_whatever.etc')
    print('Plots will be named Modelpdf.[png,eps]')
    print('Will also write Model_mean.txt for use with resmean.py')
    print('use -v for debugging output')
    sys.exit()

if '-v' in sys.argv:
    verbose = True
else:
    verbose = False

modelname = sys.argv[1].split('_')[0]
print('Plotting residuals vs. period for {0}'.format(modelname))

#####
# Read input file and convert to log(pd)
# (use log(pd) bc hist2d does NOT look good if log-scaled after plotting)
f, T, logres_orig, pdf_orig = np.loadtxt(sys.argv[1], unpack=True)

# NEW! PGA is at -1 so let's peel that off before we do log10
i_pga, = np.where(T==-1)
#i_other, = np.where(T!=-1) 
# pick off unwanted short periods; leave 2nd lowest so we can stick PGA in
i_other, = np.where(T>0.035) 
logres_pga = logres_orig[i_pga]
pdf_pga = pdf_orig[i_pga]


T = T[i_other]
logres = logres_orig[i_other]
pdf = pdf_orig[i_other]
T = np.log10(T)

# for now let's just stick PGA at the 2nd lowest period? or lowest?
n_pga = len(i_pga)
#logres[n_pga:2*n_pga] = np.zeros(n_pga)#logres_pga
pdf[len(T)-n_pga:] = pdf_pga # shortest period
#pdf[len(T)-2*n_pga:len(T)-n_pga] = pdf_pga # 2nd shortest period
#print(len(T))
#print(len(T)-n_pga)
#print(T[len(T)-n_pga:])

# Find unique values of period and log-residuals
Ts = np.unique(T)
lrs = np.unique(logres)
# find min and max values here...

#####
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

###############################################################
# Quick and dirty figure just to check on raw data values
#####
# Plot values as circles so we can check hist2d results later
fig, ax = plt.subplots(1)
ax.scatter(T, logres, c=pdf, cmap='bone_r')
ax.set_ylim(logres.min(),logres.max())
ax.set_ylabel('Residual')
ax.set_xlabel('log10(period) (s)')

#####
# Plot bins to make sure we got them right
# xbins at y=1
# ybins at x=0
ax.plot(xbins, np.ones(shape=len(xbins)), marker='|', mfc='k', mec='k', ms=15, linestyle='None')
ax.plot(np.zeros(shape=len(ybins)), ybins, marker='_', mfc='k', mec='k', ms=15, linestyle='None')
# Plot dashed line at y=zero
ax.plot([xbins.min(), xbins.max()],[0,0],linestyle='--', color='grey')

ax.set_xlim(xbins.min()-dy[0],xbins.max()+dy[-1])
fig.savefig(modelname+'circles.png')
print('quick scatter plot saved to {0}circles.png'.format(modelname))
###############################################################

###############################################################
# From here on: make hist2d figure
#####
# Set up hist2d figure
fig2, ax2 = plt.subplots(1, figsize=(10,8))

# Old: Rough imshow plot but this doesn't work well for unevenly sampled x axis
#plt.imshow(pdf.reshape(21,80).T, interpolation='none', cmap='magma_r')

#####
# Plot 2D histogram: 
# Since we already have the PDF values (output by Dan's matlab code),
# we just use them as weights
# Need to set 'image' variable for drawing colorbar later
counts, xedges, yedges, image = ax2.hist2d(T,logres, weights=pdf, bins=[xbins, ybins],cmap='bone_r', vmax=0.35)

#####
# Plot line at zero to highlight positive vs negative residuals
# new! use stroke (instead of plotting 2 lines)
# not sure if this is any easier than just plotting 2 though :P
stroke=[pe.Stroke(linewidth=3, foreground='w'), pe.Normal()]
ax2.plot(ax2.get_xlim(),[0,0], lw=1, color='k',
         path_effects=stroke, linestyle='--')
# plot standard deviation 
# (dummy value for now, read input file after Dan makes it)
ax2.plot(ax2.get_xlim(),[-0.7,-0.7], lw=1, color='k',
         path_effects=stroke, linestyle='--')
ax2.plot(ax2.get_xlim(),[0.7,0.7], lw=1, color='k',
         path_effects=stroke, linestyle='--')


#####
# Plot line to divide PGA (and PGV in future?) residuals 
# from period residuals

#ax2.axvline(x=Ts[0])
ax2.axvline(x=xbins[1], color='black', linewidth=3, zorder=20)
ax2.text(Ts[0], -2.1, 'PGA', ha='center', va='top')

#####
# Calculate average and standard deviation for each period
avgs = np.empty(len(Ts))
stds = np.empty(len(Ts))
for i,Ti in enumerate(Ts):
    iwhere, = np.where(T==Ti)
    if verbose:
        print('Ti:', Ti)
        print('iwhere:', iwhere)
    avgres = np.average(logres[iwhere], weights=pdf[iwhere])
    var = np.average((logres[iwhere]-avgres)**2, weights=pdf[iwhere])
    std = np.sqrt(var)
    stds[i] = std
    avgs[i] = avgres
    if verbose:
        print('Ti, avgres, std:')
        print(Ti,avgres, std)

# Write average residual to output file for use with resmeans.py
avg_outfile = open(modelname+'_mean.txt', 'w')
avg_out_header = '#log10(period) period log10(leftbinedge) log10(rightbinedge) binaverage binstd \n'
avg_outfile.write(avg_out_header)
for i in range(len(avgs)):
    avg_outfile.write('{0} {1} {2} {3} {4} {5}\n'.format(Ts[i], 10**Ts[i], xbins[i], xbins[i+1], avgs[i], stds[i]))
avg_outfile.close()
print('mean residuals written to {0}_mean.txt'.format(modelname))

#####
# Plot average and standard deviation
# remember:
# marker = symbol to use (circle, triangle, etc)
# ms = markersize (in pts)
# mec = markeredgecolor (fill color, -G in GMT)
# mfc = markerfacecolor (stroke color, half of -W in GMT)
# mew = markeredgewidth (weight of edge stroke, other half of -W in GMT)
# Set up line and marker styles
stroke2=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()]
lineeffects = {'lw':3, 'color':'grey', 'path_effects':stroke2, 'alpha':0.7}
#avgmarker = {'marker':'o', 'ms':'10', 'mec':'k', 'mew':2, 'linestyle':'None'}
#avgmarker = {'marker':'o', 'ms':'10', 'mec':'k', 'mew':2, 'color':'k'}
avgmarker = {'marker':'o', 'ms':'10', 'mec':'k', 'mew':2}

# Plot average values 
ax2.plot(Ts[0],avgs[0], **avgmarker, mfc='white', alpha=0.3, linestyle='None')
ax2.plot(Ts[0],avgs[0], **avgmarker, mfc='None', label='mean', linestyle='None')
ax2.plot(Ts[1:],avgs[1:], **avgmarker, mfc='white', alpha=0.3, linestyle='None')
ax2.plot(Ts[1:],avgs[1:], **avgmarker, mfc='None', label='mean', linestyle='None')
#ax2.plot(Ts[1:],avgs[1:], **avgmarker, mfc='None', label='mean', color='k', lw='2')

# Plot std devs across the whole bin, instead of only at average value
# Because len(xbin) = len(avg|std) + 1,
# need to repeat each xbin 3x, but chop off 1st and last value
# and repeat avgs and stds 3x, and duplicate last value
bin_repeat = np.repeat(xbins,3)[1:-1]
avg_repeat = np.append(np.repeat(avgs,3), avgs[-1])
std_repeat = np.append(np.repeat(stds,3), stds[-1])
ax2.plot(bin_repeat, avg_repeat+std_repeat, **lineeffects, label=r'± 1$\sigma$') 
ax2.plot(bin_repeat, avg_repeat-std_repeat, **lineeffects)

#ax2.plot(bin_repeat[2:], avg_repeat[2:], 'k-', lw=2)

# Legend
#ax2.legend(numpoints=1)

#####
# Labels and such
# make labels big and make ticks larger
# hide linear log10(period) axis (we'll plot a log axis later)
#ax2.set_xlabel('log10 period (s)', fontsize='large')
#ax2.tick_params(labelsize='large', direction='inout', which='both')
ax2.xaxis.set_visible(False)
#ax2.set_ylabel('Residual ln(obs)-ln(pred)', fontsize='large')
ax2.set_ylabel('Residual log10(obs)-log10(pred)', fontsize='large')
ax2.tick_params(axis='y', labelsize='large', direction='inout', which='both')
ax2.tick_params(axis='y', which='major', length=10)
ax2.set_title(modelname, fontsize='xx-large')

#####
# Draw logarithmic period axis
ax3 = ax2.twiny()
ax3.set_xlabel('Period (s)', fontsize='large')
ax3.semilogx()
ax3.tick_params(labelsize='large', direction='inout', which='both')
ax3.tick_params(axis='x', which='major', length=15)
ax3.tick_params(axis='x', which='minor', length=8)
ax3.tick_params(labelbottom=True, labeltop=False)
ax3.xaxis.set_ticks_position('both')
ax3.xaxis.set_label_position('bottom')

# Use decimal values (0.1, 1.0) instead of powers (10^-1, 10^0)
formatter = ticker.FormatStrFormatter('%.1f')
ax3.get_xaxis().set_major_formatter(formatter) 

# Set axis limits:
ax3.set_xlim(10**(np.min(xbins)), 10**np.max(xbins))
ax2.set_xlim(xbins[0], xbins[-1])

# Skip first 2 pds and only plot residuals from +/-3 (Dan requested) 
#ax3.set_xlim(10**xbins[2], 10**np.max(xbins))
#ax2.set_xlim(xbins[2], xbins[-1])

#ax2.set_ylim(-3,3)
ax2.set_ylim(-2,2)

#####
# Draw colorbar 
# need to pass both lin and log axes to the 'ax' keyword arg to get proper colorbar placement
fig2.colorbar(image, ax=[ax2,ax3], label='Probability')

#####
# Finish up and save
#fig2.tight_layout()
# tight layout doesn't work with the color bar (boooooo)
# if you want to trim off the whitespace,
# use this bash command after processing:
# convert -trim pdf.png meh.png
fig2.savefig(modelname+'pdf.eps')
print('residual histogram saved to {0}pdf.eps'.format(modelname))
fig2.savefig(modelname+'pdf.png')
fig2.savefig(modelname+'pdf.pdf')
print('residual histogram saved to {0}pdf.png'.format(modelname))
#plt.show()


