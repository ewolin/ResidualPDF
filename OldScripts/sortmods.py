#!/usr/bin/env python

import numpy as np
from glob import glob

meanfiles = glob('*_sigma.out')
modelmeans = {}
for i,mf in enumerate(meanfiles):
    modelname = mf.split('_')[0]
    sigmafile = modelname+'_sigma.out'
    T_std, meanres, stdres, stdpred = np.loadtxt(sigmafile, unpack=True)
    i_psa, = np.where(T_std>0.041)
#    T_std_psa = np.log10(T_std[i_psa])
    meanres_psa = meanres[i_psa]
#    stdres_psa = stdres[i_psa]
#    stdpred_psa = stdpred[i_psa]
    meanavg = np.mean(meanres_psa)
    print(modelname, meanavg)
    modelmeans[modelname] = meanavg
print('------')

k = list(modelmeans.keys())
v = list(modelmeans.values())

i_sort = np.argsort(v)

v.sort()

newmods = []

for i in i_sort:
    newmods.append(k[i])

for i in range(len(i_sort)):
    print(newmods[i], v[i])
