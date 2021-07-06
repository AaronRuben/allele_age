#!/bin/bash
if [[ -f data/heterozygous_counts_plink2.tab ]]; then

  rm data/heterozygous_counts_plink2.tab
fi

for file in data/YRI*vcf.gz; do
  if [[ "$file" != *chrX* ]] && [[ "$file" != *chrY* ]] && [[ "$file" != *chrMT* ]] && [[ "$file" != *wgs* ]]; then
    filebase=$(echo $file | sed 's/.vcf.gz/_heterozygosity_plink2/')
    # get het count using plink
    plink2 --vcf $file --sample-counts cols=hom,het --out $filebase
    # extract counts
    if [[ -f data/heterozygous_counts_plink2.tab.tmp ]]; then
      cut -f3 ${filebase}.scount | paste data/heterozygous_counts_plink2.tab.tmp - |
        sponge data/heterozygous_counts_plink2.tab.tmp
    else
      cut -f3 ${filebase}.scount > data/heterozygous_counts_plink2.tab.tmp
    fi
  fi
done
# sum counts
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; print sum}' data/heterozygous_counts_plink2.tab.tmp |
  sponge data/heterozygous_counts_plink2.tab.tmp
cut -f1 ${filebase}.scount | paste - data/heterozygous_counts_plink2.tab.tmp |
  sponge data/heterozygous_counts_plink2.tab.tmp
awk -F '\t' 'BEGIN{OFS="\t"; print "#IID", "HET_CT", "H"}{ if (NR cols!= 1) print $0,$2/(3200000000)}' data/heterozygous_counts_plink2.tab.tmp |
  sponge data/heterozygous_counts_plink2.tab.tmp
awk -F '\t' 'BEGIN{OFS="\t"; print "#IID", "HET_CT", "H", "theta"}{if (NR != 1) print $0,-$3/($3 - 1)}' data/heterozygous_counts_plink2.tab.tmp |
  sponge data/heterozygous_counts_plink2.tab.tmp
awk -F '\t' 'BEGIN{OFS="\t"; print "#IID", "HET_CT", "H", "theta", "Ne"}{if (NR != 1) print $0, $4 / (4 * 1.2 * 1e-8)}' data/heterozygous_counts_plink2.tab.tmp |
  sponge data/heterozygous_counts_plink2.tab.tmp

for file in data/ESN*vcf.gz; do
  if [[ "$file" != *chrX* ]] && [[ "$file" != *chrY* ]] && [[ "$file" != *chrMT* ]] && [[ "$file" != *wgs* ]]; then
    filebase=$(echo $file | sed 's/.vcf.gz/_heterozygosity_plink2/')
    # get het count using plink
    plink2 --vcf $file --sample-counts cols=hom,het --out $filebase
    # extract counts
    if [[ -f data/heterozygous_counts_plink2.tab ]]; then
      cut -f3 ${filebase}.scount | paste data/heterozygous_counts_plink2.tab - | sponge data/heterozygous_counts_plink2.tab
    else
      cut -f3 ${filebase}.scount > data/heterozygous_counts_plink2.tab
    fi
  fi
done
# sum counts
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; print sum}' data/heterozygous_counts_plink2.tab |
  sponge data/heterozygous_counts_plink2.tab
cut -f1 ${filebase}.scount | paste - data/heterozygous_counts_plink2.tab |
  sponge data/heterozygous_counts_plink2.tab
awk -F '\t' 'BEGIN{OFS="\t"; print "#IID", "HET_CT", "H"}{ if (NR cols!= 1) print $0,$2/(3100000000)}' data/heterozygous_counts_plink2.tab |
  sponge data/heterozygous_counts_plink2.tab
awk -F '\t' 'BEGIN{OFS="\t"; print "#IID", "HET_CT", "H", "theta"}{if (NR != 1) print $0,-$3/($3 - 1)}' data/heterozygous_counts_plink2.tab |
  sponge data/heterozygous_counts_plink2.tab
awk -F '\t' 'BEGIN{OFS="\t"; print "#IID", "HET_CT", "H", "theta", "Ne"}{if (NR != 1) print $0, $4 / (4 * 1.2 * 1e-8)}' data/heterozygous_counts_plink2.tab |
  sponge data/heterozygous_counts_plink2.tab

# merge datasets
tail -n+2 data/heterozygous_counts_plink2.tab.tmp >> data/heterozygous_counts_plink2.tab
rm data/heterozygous_counts_plink2.tab.tmp
