#!/bin/bash

## User specified parameters
min_alternate_count=$1
min_coverage=$2
min_alt_frac=$3
prefix=$4
bamFile=$5
backgroundBamFile=$6

## Parameters untouched
reference_sequence=/mnt/shared/projects/barley/201906_MorexGenomeV2/genome/splitChromos/111018_Morex_pseudomolecules_parts.fasta

#=================================
# genotype calling
#=================================

freebayes \
--fasta-reference $reference_sequence \
--min-alternate-count $min_alternate_count \
--min-coverage $min_coverage \
--min-alternate-fraction $min_alt_frac \
--min-mapping-quality 30 \
--no-indels \
--no-mnps \
--no-complex \
--use-best-n-alleles 4 \
--no-population-priors \
--legacy-gls \
--vcf $prefix.vcf \
-b $backgroundBamFile \
-b $bamFile 

bcftools stats $prefix.vcf > $prefix.stats

#=================================
# Convert to bcf and index
#=================================

bcftools convert --output-type b --output $prefix.bcf $prefix.vcf
bcftools index $prefix.bcf
# Remove vcf file
rm $prefix.vcf
