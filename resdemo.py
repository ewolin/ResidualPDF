#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from scipy.stats import norm
import sys

def mygauss(x,mean,sigma):
#    scale = 1./(sigma*np.sqrt(2*np.pi))
#    g = scale*np.exp(-0.5*((x-mean)/sigma)**2)
    scale = 1./(np.sqrt(2*np.pi*sigma**2))
    g = scale*np.exp(-1*(x-mean)**2/(2*sigma**2))
    return g
    

modelname = sys.argv[1]
T_ex = float(sys.argv[2])

x = np.linspace(-2,2,100)

stdmarker = {'marker':'|', 'ms':'20', 'mew':2, 'linestyle':'None'}
avgmarker = {'marker':'o', 'ms':'10', 'mec':'k', 'mew':2, 'linestyle':'None'}

############################
# Read input file
pdffile='{0}_PDF.out'.format(modelname)
sigmafile='{0}_sigma.out'.format(modelname)
#T_ex = 0.1 
#T_ex = 1.0 
# Find PDF 
f, T_orig, logres_orig, pdf_orig = np.loadtxt(pdffile, unpack=True)
iwhere, = np.where(T_orig == T_ex)
T = T_orig[iwhere]
print(T)
logres = logres_orig[iwhere]
pdf = pdf_orig[iwhere]
#pdf = 0.25*norm(loc=0.076073, scale=1.04046/2).pdf(logres)
print(sum(pdf))

# Find mean and sigma
T_std, meanres1, stdres1, stdpred1 = np.loadtxt(sigmafile, unpack=True)
iwhere, = np.where(T_std == T_ex)
meanres = meanres1[iwhere]
stdres = stdres1[iwhere]
stdpred = stdpred1[iwhere]
print(meanres, stdres, stdpred)
############################
scale=0.25
#scale=1
#scale=1./0.25

# Plot model
# mean is always 0, sigma varies
#stdpred = 0.4
#w_model = norm(loc=0, scale=stdpred)
#w_model = mygauss(x,0,stdpred)

#plt.plot(x,w_model, label='model', color='r', lw=3)
#plt.plot(-1*stdpred,w_model(-1*stdpred), **stdmarker, mec='r')
#plt.plot(stdpred,w_model(stdpred), **stdmarker, mec='r')

w_model = norm(loc=0, scale=stdpred/2)
plt.plot(x,scale*w_model.pdf(x), label=modelname, color='r', lw=1)
plt.plot(-1*stdpred,scale*w_model.pdf(-1*stdpred), **stdmarker, mec='r', label=r'Predicted 1$\sigma$')
plt.plot(stdpred,scale*w_model.pdf(stdpred), **stdmarker, mec='r')


# Plot data
w_res = norm(loc=meanres, scale=stdres/2)

plt.plot(x,scale*w_res.pdf(x), label='data', color='grey', lw=1)
plt.plot(meanres-stdres,scale*w_res.pdf(meanres-stdres), **stdmarker, mec='grey', label=r'Residual 1$\sigma$')
plt.plot(meanres+stdres,scale*w_res.pdf(meanres-stdres), **stdmarker, mec='grey')
plt.plot(meanres,scale*w_res.pdf(meanres), **avgmarker, mfc='None', label='Mean residual')
plt.plot(meanres,scale*w_res.pdf(meanres), **avgmarker, mfc='w', alpha=0.3)


N, bins, patches = plt.hist(logres, len(logres), weights=pdf, zorder=0)
fracs = N #/ N.max()

print('---pdf---')
print(pdf)
print('---pdf---')
#print(N)
sumN = 0.
for i in range(len(N)):
    sumN += N[i]
    print(i, bins[i], N[i], sumN)

# we need to normalize the data to 0..1 for the full range of the colormap
norm = colors.Normalize(0, 0.35) #fracs.min(), fracs.max())
#norm = colors.Normalize(0, 1)

for thisfrac, thispatch in zip(fracs, patches):
#    color = plt.cm.bone_r(set_clim)
    color = plt.cm.bone_r(norm(thisfrac))
    thispatch.set_facecolor(color)
    thispatch.set_zorder(0)




# Finish up
plt.xlim(-2,2)
plt.xlabel('Residual log10(obs)-log10(pred)', fontsize='large')
plt.ylim(0,.5)
plt.plot([0,0],[0,1], 'r--', lw=3)
plt.title('T={0} s'.format(T_ex))
plt.legend(numpoints=1)
#plt.savefig('meh.png')
plt.savefig('LLH_{0}_{1}.png'.format(modelname,T_ex))
plt.savefig('LLH_{0}_{1}.eps'.format(modelname,T_ex))
print('saved', 'LLH_{0}_{1}.png'.format(modelname,T_ex))
