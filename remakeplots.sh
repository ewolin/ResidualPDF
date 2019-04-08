#!/bin/bash

for f in $(ls *_PDF.out)
do
    modelname=$(echo $f | cut -d _ -f 1)
#    python $HOME/code/ResidualPDF/respdf.py $modelname
    python $HOME/code/ResidualPDF/respdf.py $modelname --tmin 0.04
done
