#!/bin/bash


chr=$1
hgdp_dir=/data/users/ajimenez/data/genomic_data/HGDP_americans/raw
vcf=${hgdp_dir}/hgdp_amrs.full.chr${chr}.hg19.vcf.gz
out_dir=/data/users/ajimenez/data/genomic_data/HGDP_americans/biallelic
s_dir=/data/users/ajimenez/data/genomic_data/HGDP_americans/kept_n_removed
out_biallelic=${out_dir}/hgdp_amrs.biallelic.chr${chr}vcf.gz
out_indels=${s_dir}/chr${chr}.indels
out_missing=${s_dir}/chr${chr}.missing
out_nomiss=${out_dir}/hgdp_amrs.nomissing.biallelic.chr${chr}.vcf.gz
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

vcftools --gzvcf $vcf --max-missing 1.0 --removed-sites --out $out_missing

#Remove missing positions

vcftools --gzvcf $vcf --max-missing 1.0 --recode --recode-INFO-all --stdout | bgzip -c > $out_nomiss

tabix -p vcf $out_nomiss
#Get the kept sites

vcftools --gzvcf $out_nomiss --kept-sites --out $out_kept

#Rm biallelic file

#rm $out_biallelic

exit 0
