#!/usr/bin/env python

# Convert raw output from Matlab as list of:
#gmm, PGAmean, PGAstd, PSA03m, PSA03s, PSA1m, PSA1s, PSA3m, PSA3s
# to input format more similar to *_sigma.out
# pd, mean, std 
# Expects Matlab raw output to be in text file called
# MODEL_[volcanic, tectonic]_shakemap.txt

import sys

if len(sys.argv) != 2:
    print('Usage: shakemapraw2pdf.py MODEL [volcanic | tectonic]')
    print('Expects raw Matlab output for volcanic/tectonic events')
    print('to be saved in Shakemap_volcanic.txt or Shakemap_tectonic.txt')
    sys.exit

# specify tectonic or volcanic as command line arg
evtype = sys.argv[1]

infile = open('Shakemap_{0}.txt'.format(evtype), 'r')

header = infile.readline()
lines = infile.readlines()

for line in lines:

    l = line.strip().split()
    
    m = [ k.strip(',') for k in l ]
    
    modelname = m[0]
    print(modelname)
    
    # periods from Shakemap: PGA (-1), 0.3, 1.0, 3.0
    # WARNING hard-coded
    pds = [-1, 0.3, 1.0, 3.0]
    
    filename = '{0}_{1}_shakemap.txt'.format(modelname, evtype)
    print(filename)
    outfile = open(filename, 'w')
    
    for i in range(len(pds)):
        j = 2*i+1
        print(pds[i], float(m[j]), float(m[j+1]))
        outfile.write('{0} {1} {2}\n'.format(pds[i], float(m[j]), float(m[j+1])))
    outfile.close
    
