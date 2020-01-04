#!/usr/bin/env python
# Plot to demonstrate meaning of residual PDF plots
# Plot histogram of residuals, 
# and 2 Gaussians to show model and residual distribution

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib import colors

############################
# Check for correct # of args
if len(sys.argv) != 3:
    print('Usage: resdemo.py MODEL period')

# Read input files
try: 
    modelname = sys.argv[1]
    pdffile='{0}_PDF.out'.format(modelname)
    sigmafile='{0}_sigma.out'.format(modelname)
    f, T_orig, logres_orig, pdf_orig = np.loadtxt(pdffile, unpack=True)
    T_std, meanres1, stdres1, stdpred1 = np.loadtxt(sigmafile, unpack=True)
except Exception as e:
    print('Error: ',e)
    print('MODEL should be the name of a ground motion model')
    print('with MODEL_PDF.out and MODEL_sigma.out files in cwd')
    sys.exit()
   
# Set specified period T_ex
# if not given, print a list of periods to choose
try:
    if len(sys.argv) != 3:
        print('choose from periods found in {0}_PDF.out:'.format(modelname))
        print(np.unique(T_orig))
        sys.exit()
    T_ex = float(sys.argv[2])
except Exception as e:
    print(e)
    sys.exit()

# Pull values at period T_ex out of *_PDF.out file 
iwhere, = np.where(T_orig == T_ex)
T = T_orig[iwhere]
logres = logres_orig[iwhere]
pdf = pdf_orig[iwhere]
print('sum of all bins at T={0} (should be 1.0):'.format(T_ex), sum(pdf))
print(sum(pdf*0.25))

# Pull mean and sigmas at period T_ex out of *_sigma.out file
iwhere, = np.where(T_std == T_ex)
meanres = meanres1[iwhere]
stdres = stdres1[iwhere]
stdpred = stdpred1[iwhere]

meanres=-0.2989
stdres=0.2950
print('mean residual: {0}'.format(meanres))
print('std of observed residuals: {0}'.format(stdres))
print('std of model predictions:  {0}'.format(stdpred))

############################
# Set some plot parameters
stdmarker = {'marker':'|', 'ms':'20', 'mew':2, 'linestyle':'None'}
avgmarker = {'marker':'o', 'ms':'10', 'mec':'k', 'mew':2, 'linestyle':'None'}

############################
# Plot Gaussian showing model predictions and observed residuals
# for model predictions, mean is always 0; sigma varies
binwidth=0.25 # so PDF has same area as histogram (bins sum to 1, need to multiply by width)

x = np.linspace(-2,2,100)

w_model = norm(loc=0, scale=stdpred)
plt.plot(x, binwidth*w_model.pdf(x), label=modelname, color='r', lw=1)
plt.plot(-1*stdpred,binwidth*w_model.pdf(stdpred), **stdmarker, 
         mec='r')#, label=r'Predicted 1$\sigma$')
plt.plot(stdpred,binwidth*w_model.pdf(stdpred), **stdmarker, mec='r')

w_res = norm(loc=meanres, scale=stdres)

plt.plot(x,binwidth*w_res.pdf(x), label='data', color='grey', lw=1)
plt.plot(meanres-stdres,binwidth*w_res.pdf(meanres-stdres), 
         **stdmarker, mec='grey')#, label=r'Residual 1$\sigma$')
plt.plot(meanres+stdres,binwidth*w_res.pdf(meanres-stdres), 
         **stdmarker, mec='grey')
plt.plot(meanres,binwidth*w_res.pdf(meanres), 
         **avgmarker, mfc='None', label='Mean residual')
plt.plot(meanres,binwidth*w_res.pdf(meanres), **avgmarker, mfc='w', alpha=0.3)

############################
# Plot histogram of observed residuals
N, bins, patches = plt.hist(logres, len(logres), weights=pdf, zorder=0)
fracs = N #/ N.max()

# uncomment for debugging
#print('---pdf---')
#for i in range(len(pdf)):
#    print(logres[i], pdf[i])
#print('---pdf---')

# Get fancy: fill each histogram bin with the same color used 
# in residual PDF colorbar (which is normalized btw 0 and 0.35)
# we need to normalize the data to 0..1 for the full range of the colormap
norm = colors.Normalize(0, 0.35) 
for thisfrac, thispatch in zip(fracs, patches):
    color = plt.cm.bone_r(norm(thisfrac))
    thispatch.set_facecolor(color)
    thispatch.set_zorder(0)

############################
# Add title, labels, and legend, set limits, save to png and eps
plt.xlim(-2,2)
plt.xlabel('Residual log10(obs)-log10(pred)', fontsize='large')
#plt.ylim(0,.5)
plt.ylim(0,0.5)

plt.plot([0,0],[0,1], 'r--', lw=3)
plt.title('{0}, T={1} s'.format(modelname, T_ex))
plt.legend(numpoints=1)
plt.savefig('LLH_{0}_{1}.png'.format(modelname,T_ex))
plt.savefig('LLH_{0}_{1}.eps'.format(modelname,T_ex))
print('saved', 'LLH_{0}_{1}.png'.format(modelname,T_ex))
