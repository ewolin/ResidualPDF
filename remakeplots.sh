#!/bin/bash

for f in $(ls *_PDF.out)
do
    modelname=$(echo $f | cut -d _ -f 1)
    python respdf.py $modelname
done
