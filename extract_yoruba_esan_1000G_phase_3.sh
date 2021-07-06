#!/bin/bash

default_path_1000G_phase3="/mnt/lachance_projects/1000_Genomes_Project/Phase3_vcf/"
path_1000G_phase3=${1:-$default_path_1000G_phase3}
echo "Expect VCF files in directory ${path_1000G_phase3}"

# filter out children and uncertain samples
awk -F '\t' '{if ($8 !~ /child/ && $NF == 0) print $2}' data/pedigree_YRI.ped > data/yoruba_samples_1000G_phase3_pedigree_filtered.txt
awk -F '\t' '{if ($8 !~ /child/ && $NF == 0) print $2}' data/pedigree_ESN.ped > data/esan_samples_1000G_phase3_pedigree_filtered.txt
# iterate over all chromosomes
for file in ${path_1000G_phase3}*.vcf.gz; do
  # create output name
  out_file_yoruba=$(echo $file | awk -F '/' '{print $NF}' | sed 's/ALL/YRI/')
  out_file_esan=$(echo $file | awk -F '/' '{print $NF}' | sed 's/ALL/ESN/')
  echo $out_file_yoruba
  echo $out_file_esan
  # extract samples
  bcftools view -Oz -S data/yoruba_samples_1000G_phase3_pedigree_filtered.txt --force-samples -o data/${out_file_yoruba} --threads 16 $file
  bcftools view -Oz -S data/esan_samples_1000G_phase3_pedigree_filtered.txt --force-samples -o data/${out_file_esan} --threads 16 $file

done
