#!/bin/bash

chr=$1
hgdp_dir=/data/users/ajimenez/data/genomic_data/HGDP_americans/lift_hg19
vcf=${hgdp_dir}/rejected.chr${chr}.hg19.vcf.gz
s_dir=/data/users/ajimenez/data/genomic_data/HGDP_americans/lift_hg19
out_biallelic=${s_dir}/chr${chr}.bisnps
out_indels=${s_dir}/chr${chr}.indels

###############################
### Get only biallelic SNPs ###
###############################

#Save indels and others to a list
vcftools --gzvcf $vcf --remove-indels --min-alleles 2 --max-alleles 2 --removed-sites --out $out_indels

vcftools --gzvcf $vcf --remove-indels --min-alleles 2 --max-alleles 2 --kept-sites --out $out_biallelic

exit 0
