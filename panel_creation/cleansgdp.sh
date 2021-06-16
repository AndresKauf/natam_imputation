#!/bin/bash


name=$1
simons_dir=/data/users/ajimenez/data/genomic_data/simons_vcfs
vcf=${simons_dir}/${name}.annotated.nh2.variants.vcf.gz
out_dir=/data/users/ajimenez/data/genomic_data/simons_americans
out_biallelic=${out_dir}/${name}.biallelic.vcf.gz
out_indels=${out_dir}/${name}.indels
out_missing=${out_dir}/${name}.missing
out_nomiss=${out_dir}/${name}.nomissing.biallelic.vcf.gz
out_kept=${out_dir}/${name}

###############################
### Get only biallelic SNPs ###
###############################

#Save indels and others to a list
vcftools --gzvcf $vcf --remove-indels --min-alleles 2 --max-alleles 2 --removed-sites --out $out_indels

#Remove said indels and create index
vcftools --gzvcf $vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout | bgzip -c > $out_biallelic

tabix -p vcf $out_biallelic

################################
### Get no-missing positions ###
################################

#Save missing positions to a list

vcftools --gzvcf $out_biallelic --max-missing 1.0 --removed-sites --out $out_missing

#Remove missing positions

vcftools --gzvcf $out_biallelic --max-missing 1.0 --recode --recode-INFO-all --stdout | bgzip -c > $out_nomiss

tabix -p vcf $out_nomiss
#Get the kept sites

vcftools --gzvcf $out_nomiss --kept-sites --out $out_kept

exit 0
