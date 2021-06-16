#!/bin/bash


chr=$1
mex_dir=/data/users/ajimenez/data/genomic_data/mx_bio50g/raw
vcf=${mex_dir}/mx.50g.raw.chr${chr}.vcf.gz
out_dir=/data/users/ajimenez/data/genomic_data/mx_bio50g/biallelic
s_dir=/data/users/ajimenez/data/genomic_data/mx_bio50g/kept_n_removed
out_biallelic=${out_dir}/mx.50g.biallelic.chr${chr}vcf.gz
out_indels=${s_dir}/chr${chr}.indels
out_missing=${s_dir}/chr${chr}.missing
out_nomiss=${out_dir}/mx.50g.nomissing.biallelic.chr${chr}.vcf.gz
out_kept=${s_dir}/chr${chr}

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

#Rm biallelic file

rm $out_biallelic*

exit 0
