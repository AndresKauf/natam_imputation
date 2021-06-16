#!/bin/bash

#Script para simular datos en las proporciones presentes en 1KGP

seed=113

stdpopsim HomSap -s $seed -g HapMapII_GRCh37 -c chr9 -o ../sim.data.ts -d AmericanAdmixture_4B11 1322 1006 6000 1294

tskit vcf --ploidy 2 ../sim.data.ts > ../sim.data.vcf

exit 0 
