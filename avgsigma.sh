#!/bin/bash

modelname=$1
sigmafile=${1}_sigma.out

awk 'BEGIN {n=0; sum=0} {sum+=($3-$4); n+=1} END {print sum/n}' $sigmafile

awk 'BEGIN {n=0; sum=0} {sum+=sqrt(($3-$4)*($3-$4)); n+=1} END {print sum/n}' $sigmafile
