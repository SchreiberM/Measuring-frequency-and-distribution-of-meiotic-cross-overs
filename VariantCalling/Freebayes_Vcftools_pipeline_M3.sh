#!/bin/bash -e

# Tools required: freebayes, vcftools
# to use freebayes with --variant-input we went back to version v1.0.2 (see GitHub Issues)

## User specified parameters
prefix=$1
bam_file=$2
reference_sequence=$3
numThreads=$4
Background=$5

#=================================
# genotype calling
#=================================

freebayes \
--fasta-reference $reference_sequence \
--min-alternate-count 1 \
--min-coverage 1 \
--min-alternate-fraction 0 \
--min-mapping-quality 30 \
--use-best-n-alleles 4 \
--no-population-priors \
--legacy-gls \
--variant-input $Background \
--only-use-input-alleles \
--vcf $prefix.vcf \
--bam $bam_file 


# For M3 plants and 2x coverage we used these parameters:
vcftools --vcf $prefix.vcf --minDP 2 --recode --stdout 2>/dev/null | grep -v "#" | awk -v OFS="\t" '{print $1,$2,$10}' | cut -d":" -f1 | sed 's#1/1#0#' | sed 's#0/1#0.5#' | sed 's#0/0#0#' | awk '{if ($3==0.5) print $0"\t""blue"; else print $0"\t""grey"}' >  $prefix.filtered_2.tsv

# For M3 plants and 4x coverage we used these parameters:
vcftools --vcf $prefix.vcf --minDP 4 --recode --stdout 2>/dev/null | grep -v "#" | awk -v OFS="\t" '{print $1,$2,$10}' | cut -d":" -f1 | sed 's#1/1#0#' | sed 's#0/1#0.5#' | sed 's#0/0#0#' | awk '{if ($3==0.5) print $0"\t""blue"; else print $0"\t""grey"}' >  $prefix.filtered_3.tsv


