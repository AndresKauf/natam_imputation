#!/bin/bash

ind=$1

array_pos=/data/users/ajimenez/projects/sim_imp2000/sim_array/simulated_pos_updateprob.txt
array_nopos=/data/users/ajimenez/projects/sim_imp2000/sim_array/simulated_nopos.txt
include=/data/users/ajimenez/projects/sim_imp2000/imputation_run/ind_index/admix${ind}.txt
tarhaps=/data/users/ajimenez/projects/sim_imp2000/imputation_run/targets/admix${ind}.target.haps
tarsamp=/data/users/ajimenez/projects/sim_imp2000/imputation_run/targets/admix${ind}.target.sample
valhaps=/data/users/ajimenez/projects/sim_imp2000/imputation_run/validation/admix${ind}.val.haps
valsamp=/data/users/ajimenez/projects/sim_imp2000/imputation_run/validation/admix${ind}.val.sample

#make target
shapeit -convert --input-haps /data/users/ajimenez/projects/sim_imp2000/sim.data.all --include-ind $include --include-snp $array_pos --output-haps $tarhaps $tarsamp

#make validation
shapeit -convert --input-haps /data/users/ajimenez/projects/sim_imp2000/sim.data.all --include-ind $include --include-snp $array_nopos --output-haps $valhaps $valsamp

exit 0 
