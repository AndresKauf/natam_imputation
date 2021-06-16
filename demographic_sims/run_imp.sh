#!/bin/bash

refs=(100 134)
for r in "${refs[@]}"
do
  ref_script=./imp${r}.sh
  parallel -j20 $ref_script :::: numbers.txt
done

exit 0
