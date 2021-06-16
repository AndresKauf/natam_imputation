#!/bin/bash

ind=$1
impute=/data/users/ajimenez/bin/impute_v2.3.2_x86_64_static/impute2
f=/data/users/ajimenez/projects/sim_imp/imputation_run/chr_index/chr9_index.txt
refhap=/data/users/ajimenez/projects/sim_imp2000/imputation_run/ref_panels/ref.148nats.haps
refleg=/data/users/ajimenez/projects/sim_imp2000/imputation_run/ref_panels/ref.148nats.legend
outdir=/data/users/ajimenez/projects/sim_imp2000/imputation_run/results/ref148/ind${ind}
targ=/data/users/ajimenez/projects/sim_imp2000/imputation_run/targets/admix${ind}.target.haps
map=/data/reference_panels/genetic_map_b37/genetic_map_chr9_combined_b37.txt

for i in {1..57}
do
    l=`awk -v ind=$i 'NR==ind' $f | cut -f1`
    ri=`awk -v ind=$i 'NR==ind' $f | cut -f2`
    out=${outdir}/ind${ind}_chunk_${i}
    $impute \
      -use_prephased_g \
      -known_haps_g $targ \
      -m $map \
      -h $refhap \
      -l $refleg \
      -Ne 20000 \
      -k_hap 2796 \
      -int $l $ri \
      -o $out
done

exit 0