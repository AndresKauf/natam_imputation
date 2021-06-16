#!/bin/bash

num=$1
hap=../imputation_run/ref_panels/ref.${num}nats.haps
leg=../imputation_run/ref_panels/ref.${num}nats.legend
samp=../imputation_run/ref_panels/ref.${num}nats.sample
include=../imputation_run/ref_panels/panel.${num}.nats.txt

#hap=../imputation_run/ref_panels/ref.0nats.downs.haps
#leg=../imputation_run/ref_panels/ref.0nats.downs.legend
#samp=../imputation_run/ref_panels/ref.0nats.downs.sample
#include=../imputation_run/ref_panels/panel.0.nats.downsample.txt
shapeit -convert --input-haps ../sim.data.all --include-ind $include --output-ref $hap $leg $samp
#temp1=../imputation_run/ref_panels/temp1.${num}.vcf
#temp2=../imputation_run/ref_panels/temp2.${num}.vcf
#temp3=../imputation_run/ref_panels/temp3.${num}vcf
#temp4=../imputation_run/ref_panels/temp4.${num}

#vcftools --vcf ../sim.data.newhead.vcf --keep $include --recode --recode-INFO-all --stdout > $temp1

#bcftools +fill-tags $temp1 -- -t AF,AN,AC > $temp2

#vcftools --vcf $temp2 --non-ref-ac 1 --recode --recode-INFO-all --stdout > $temp3

#bcftools convert --hapsample $temp4 $temp3

#gunzip ${temp4}.hap.gz

#mv ${temp4}.hap ${temp4}.haps

#shapeit -convert --input-haps $temp4 --output-ref $hap $leg $samp
 
#rm $temp1*
#rm $temp2*
#rm $temp3*
#rm $temp4*

exit 0

