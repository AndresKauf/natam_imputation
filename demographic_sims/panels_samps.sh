#!/bin/bash

cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.100.nats.txt

for i in {1..100}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.100.nats.txt
done

cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.134.nats.txt

for i in {1..134}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.134.nats.txt
done

cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.200.nats.txt

for i in {1..200}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.200.nats.txt
done

cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.400.nats.txt

for i in {1..400}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.400.nats.txt
done

cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.600.nats.txt

for i in {1..600}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.600.nats.txt
done

cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.800.nats.txt

for i in {1..800}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.800.nats.txt
done

cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.1000.nats.txt

for i in {1..1000}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.1000.nats.txt
done


cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.1500.nats.txt

for i in {1..1500}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.1500.nats.txt
done


cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.2000.nats.txt

for i in {1..2000}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.2000.nats.txt
done


cp ../imputation_run/ref_panels/panel.0.nats.txt ../imputation_run/ref_panels/panel.3000.nats.txt

for i in {1..3000}
do
    echo NAT_${i} >> ../imputation_run/ref_panels/panel.3000.nats.txt
done

exit 0

