#!/bin/bash
set -euo pipefail

# Simple example
echo "Running simple example"
medicc2 simple_example/simple_example.tsv output
echo ""

# Testing example
echo "Running testing example"
medicc2 testing_example/testing_example.tsv output
echo ""

# WGD example
echo "Running WGD example"
medicc2 WGD_example/WGD_example.tsv output
echo ""

# OV03-04
echo "Running OV03-04 example"
medicc2 OV03-04/OV03-04_descr.txt output -i fasta --normal-name "OV03-04_diploid" --prefix "OV03-04"
echo ""

# Gundem 2015 patients
for file in $(ls gundem_et_al_2015 | grep PTX)
do
    patient=$(echo $file | egrep "PTX0.." -o)
    echo "Running Gundem $patient"
    medicc2 gundem_et_al_2015/$file "output_gundem_et_al_2015" --prefix $patient
    echo ""
done
