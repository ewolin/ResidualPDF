#!/bin/bash

convert resmeans_smv.png -gravity West -chop 30x0 -gravity East -chop 30x0 trim_west.png
convert resmeans_smt.png -gravity East -chop 30x0 -gravity West -chop 30x0 trim_east.png
convert +append trim_east.png trim_west.png shakemap.png
open shakemap.png

convert resmeans_smt.png -gravity North -chop 0x40 trim_north2.png
convert resmeans_smv.png -gravity North -chop 0x40 trim_north.png

convert -append trim_north2.png trim_north.png +repage shakemap-vstack-whitespace.png
convert shakemap-vstack-whitespace.png -gravity West -chop 30x0 -gravity East -chop 30x0 shakemap-vstack.png
open shakemap-vstack.png
