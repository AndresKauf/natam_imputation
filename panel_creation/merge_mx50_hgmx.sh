#!/bin/bash

chr=$1
dir=/data/users/ajimenez/data/genomic_data/natives_panel
merge_file=${dir}/mergefile.chr${chr}.txt
out_raw=${dir}/raw/natpan.raw.chr${chr}.vcf.gz
#keep_pos=${dir}/keep_chr${chr}.txt
#out_freeze=${dir}/hgdp_sgmx.freeze.chr${chr}.vcf.gz

bcftools merge -l $merge_file --missing-to-ref | bgzip -c > $out_raw

tabix -p vcf $out_raw

#vcftools --gzvcf $out_raw --positions $keep_pos --recode --recode-INFO-all --stdout | bgzip -c > $out_freeze

#tabix -p vcf $out_freeze

exit 0
